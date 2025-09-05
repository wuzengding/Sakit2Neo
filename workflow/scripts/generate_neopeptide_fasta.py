#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generates a FASTA file of neoantigen peptides based on a manually reviewed
variant report and a PHASED VCF file.

This script reads the "Somatic & Noise Event Report" sheet, selects variants
where 'Manual_Select' is 'yes', and retrieves full annotation and PHASE
information from the corresponding VEP-annotated and phased VCF file.

It uses an external UniProt FASTA file to obtain full-length protein sequences,
ensuring that the generated peptide windows are based on the correct biological
context. It also intelligently reconstructs local haplotype sequences for somatic
variants by incorporating phased 'cis' germline variants.
"""

import re
import pandas as pd
import pysam
import argparse
import logging
import sys
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# --- Configuration for Logging ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)


def load_uniprot_sequences(fasta_file: str) -> Dict[str, str]:
    """
    Parses a UniProt FASTA file and returns a dictionary mapping
    UniProt ID (without isoform/version) to its canonical protein sequence.
    """
    logging.info(f"Loading UniProt protein sequences from {fasta_file}...")
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = "".join(current_seq)
                
                # Parse new ID: >sp|Q15111|... or >tr|A0A...|...
                match = re.search(r'>[a-z]{2}\|([A-Z0-9]+)', line)
                if match:
                    #uniprot_id_full = match.group(1)
                    #current_id = uniprot_id_full.split('-')[0] # Use base ID
                    current_id = line.split("|")[1]
                    current_seq = []
                else:
                    current_id = None
            elif current_id:
                current_seq.append(line.strip())
    
    if current_id:
        sequences[current_id] = "".join(current_seq)
        
    logging.info(f"Loaded {len(sequences)} canonical protein sequences.")
    return sequences


def get_csq_header(vcf_reader: pysam.VariantFile) -> List[str]:
    """Parses VEP CSQ header to get the field order."""
    if "CSQ" not in vcf_reader.header.info:
        raise ValueError("CSQ field not found in VCF header.")
    csq_header_str = vcf_reader.header.info["CSQ"].description
    match = re.search(r"Format: ([\w|]+)", csq_header_str)
    if match: return match.group(1).split("|")
    raise ValueError("Could not parse CSQ header format string.")


class SequenceProcessor:
    def __init__(self, window_size=12):
        self.window_size = window_size
        self.amino_acid_3to1 = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': '*'
        }

    def convert_to_short_notation(self, hgvsp: str) -> str:
        """Converts p.Arg273His to p.R273H."""
        if not hgvsp or 'p.' not in hgvsp or hgvsp == '/': return "/"
        mutation = hgvsp.split('p.')[-1]
        match = re.match(r'([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)', mutation)
        if match:
            ref, pos, alt = match.groups()
            # The line below is the corrected one
            return f"p.{self.amino_acid_3to1.get(ref, '?')}{pos}{'*' if alt in ['Ter', '*'] else self.amino_acid_3to1.get(alt, '?')}"
        return f"p.{mutation}"

    def parse_protein_change(self, hgvsp: str) -> Optional[Tuple[int, str, str]]:
        """Parses HGVSp string to get (position, ref_aa, alt_aa)."""
        if not hgvsp or hgvsp == '/': return None
        match = re.search(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|\*)', hgvsp)
        if match:
            ref, pos, alt = match.groups()
            return (int(pos), self.amino_acid_3to1.get(ref, '?'), '*' if alt in ['Ter', '*'] else self.amino_acid_3to1.get(alt, '?'))
        match = re.search(r'p\.([A-Z])(\d+)([A-Z]|\*)', hgvsp)
        if match:
            ref, pos, alt = match.groups()
            return (int(pos), ref, alt)
        return None

    def get_sequence_window(self, protein_seq: str, pos: int) -> Tuple[str, int]:
        """Gets a sequence window around a mutation site."""
        pos_0based = pos - 1
        if not (0 <= pos_0based < len(protein_seq)):
            logging.warning(f"Position {pos} is out of range for seq length {len(protein_seq)}")
            return "", -1
        start = max(0, pos_0based - self.window_size)
        end = min(len(protein_seq), pos_0based + self.window_size + 1)
        return protein_seq[start:end], pos_0based - start

    def get_full_variant_info(self, vcf_record, target_transcript_id: str, target_alt: str, csq_header: List[str], uniprot_db: Dict) -> Optional[Dict]:
        """
        Extracts all details for a variant, including phase and protein sequence from UniProt DB.
        """
        all_csq = [dict(zip(csq_header, e.split("|"))) for e in vcf_record.info.get("CSQ", [])]
        target_csq = next((c for c in all_csq if c.get("Feature") == target_transcript_id and c.get("Allele") == target_alt), None)

        if not target_csq:
            logging.warning(f"Could not find annotation for transcript {target_transcript_id} in VCF record at {vcf_record.chrom}:{vcf_record.pos}")
            return None

        # --- CORRECTED AND ROBUST UniProt ID EXTRACTION ---
        uniprot_id = None
        # 1. Prioritize the SWISSPROT field if it exists
        if 'SWISSPROT' in target_csq and target_csq['SWISSPROT']:
            uniprot_id = target_csq['SWISSPROT']
        # 2. Fallback to TREMBL if SWISSPROT is not available
        elif 'TREMBL' in target_csq and target_csq['TREMBL']:
            uniprot_id = target_csq['TREMBL']
        # 3. Last resort: regex search on other fields that VEP might use
        else:
            uniprot_regex = r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})'
            # Check the specific UniProt accession field VEP sometimes uses
            if 'Uniprot_ACC' in target_csq and target_csq['Uniprot_ACC']:
                 match = re.search(uniprot_regex, target_csq['Uniprot_ACC'])
                 if match: uniprot_id = match.group(1)

        # Clean the ID by removing any version numbers (.123)
        if uniprot_id:
            uniprot_id = uniprot_id.split('.')[0]

        protein_sequence = uniprot_db.get(uniprot_id, "")
        if uniprot_id and not protein_sequence:
            logging.warning(f"UniProt ID '{uniprot_id}' found for transcript {target_transcript_id}, but no matching sequence in the provided FASTA.")

        hgvsp_full = target_csq.get("HGVSp", "/")
        protein_change = self.parse_protein_change(hgvsp_full)
        if not protein_change: return None
        pos, ref_aa, alt_aa = protein_change

        phase_set, genotype_phase = None, None
        tumor_sample_name = next((s for s in vcf_record.samples if 'tumor' in s.lower()), None)
        if tumor_sample_name:
            sample_data = vcf_record.samples[tumor_sample_name]
            if 'PS' in sample_data and sample_data['PS'] is not None:
                phase_set = sample_data['PS']
            if sample_data.get('GT') and '|' in sample_data['GT']:
                try: genotype_phase = sample_data['GT'].split('|').index('1')
                except (ValueError, IndexError): genotype_phase = None
        
        return {
            "gene": target_csq.get("SYMBOL", "/"), "transcript_id": target_transcript_id,
            "hgvsp_full": hgvsp_full.split(':')[-1] if ':' in hgvsp_full else hgvsp_full,
            "p_short": self.convert_to_short_notation(hgvsp_full), "aa_pos": pos,
            "aa_ref": ref_aa, "aa_alt": alt_aa, "uniprot_id": uniprot_id or "",
            "protein_sequence": protein_sequence, "phase_set": phase_set, "haplotype": genotype_phase
        }

    def build_local_haplotype_sequence(self, somatic_variant: Dict, all_phased_variants: Dict) -> str:
        """Reconstructs the local protein sequence by incorporating phased cis germline variants."""
        base_sequence = somatic_variant["protein_sequence"]
        somatic_ps, somatic_haplotype = somatic_variant["phase_set"], somatic_variant["haplotype"]

        if not base_sequence or somatic_ps is None or somatic_haplotype is None:
            return base_sequence
        
        haplotype_seq_list = list(base_sequence)
        for germline_variant in all_phased_variants.get("Nearby Germline", []):
            if germline_variant["phase_set"] == somatic_ps and germline_variant["haplotype"] == somatic_haplotype:
                pos_0based = germline_variant["aa_pos"] - 1
                if 0 <= pos_0based < len(haplotype_seq_list):
                    original_aa, new_aa = haplotype_seq_list[pos_0based], germline_variant["aa_alt"]
                    haplotype_seq_list[pos_0based] = new_aa
                    logging.info(f"Applied cis germline {germline_variant['p_short']}. Changed AA at {pos_0based+1} from {original_aa} to {new_aa}.")
        return "".join(haplotype_seq_list)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report_xlsx', required=True, help='Input Excel report file.')
    parser.add_argument('--vcf_file', required=True, help='The PHASED VEP-annotated VCF file.')
    parser.add_argument('--uniprot_fasta', required=True, help='UniProt canonical protein sequence FASTA file.')
    parser.add_argument('--output_fasta', required=True, help='Output FASTA file.')
    parser.add_argument('--window_size', type=int, default=12, help='Amino acid window size.')
    args = parser.parse_args()

    uniprot_db = load_uniprot_sequences(args.uniprot_fasta)
    processor = SequenceProcessor(window_size=args.window_size)

    logging.info(f"Reading variants from {args.report_xlsx}...")
    try:
        df = pd.read_excel(args.report_xlsx, sheet_name="Somatic & Noise Event Report")
        selected_df = df[df['Manual_Select'].astype(str).str.lower() == 'yes'].copy()
    except Exception as e:
        logging.error(f"Failed to read Excel file: {e}"); sys.exit(1)

    if selected_df.empty:
        logging.warning("No variants marked 'yes'. Output FASTA will be empty.")
        open(args.output_fasta, 'w').close(); sys.exit(0)
    
    selected_lookup = defaultdict(list)
    for _, row in selected_df.iterrows():
        selected_lookup[f"{row['CHROM']}:{row['POS']}"].append((row['Transcript'], row['REF'], row['ALT'], row['Somatic_Status']))

    logging.info(f"Found {len(selected_df)} variants selected for FASTA generation.")
    logging.info(f"Extracting full annotations from {args.vcf_file}...")
    
    all_phased_variants = defaultdict(list)
    with pysam.VariantFile(args.vcf_file) as vcf:
        csq_header = get_csq_header(vcf)
        for record in vcf:
            key = f"{record.chrom}:{record.pos}"
            if key in selected_lookup:
                for transcript, ref, alt, status in selected_lookup[key]:
                    if record.ref == ref and alt in record.alts:
                        full_info = processor.get_full_variant_info(record, transcript, alt, csq_header, uniprot_db)
                        if full_info:
                            all_phased_variants[status].append(full_info)

    somatic_variants = all_phased_variants.get("Somatic", [])
    logging.info(f"Writing {len(somatic_variants)} somatic peptide sequences to {args.output_fasta}...")
    
    with open(args.output_fasta, 'w') as f_out:
        for somatic_var in somatic_variants:
            try:
                haplotype_seq = processor.build_local_haplotype_sequence(somatic_var, all_phased_variants)
                ref_window, rel_pos = processor.get_sequence_window(haplotype_seq, somatic_var["aa_pos"])

                if not ref_window or rel_pos < 0:
                    logging.warning(f"Could not generate window for {somatic_var['p_short']}. Using artificial X-padded sequence.")
                    padding = 'X' * processor.window_size
                    ref_window = padding + somatic_var["aa_ref"] + padding
                    rel_pos = processor.window_size
                
                alt_aa = somatic_var["aa_alt"]
                # Check for stop codon ('*')
                if alt_aa == '*':
                    var_window = ref_window[:rel_pos]
                else:
                    var_window = ref_window[:rel_pos] + alt_aa + ref_window[rel_pos + 1:]
                f_out.write(f">Ref|{somatic_var['transcript_id']}\n{ref_window}\n")
                uniprot_id_str = somatic_var.get('uniprot_id', '')
                header = (f">Var|{somatic_var['transcript_id']}|{somatic_var['gene']}|"
                          f"{somatic_var['hgvsp_full']}|{somatic_var['p_short']}|"
                          f"{uniprot_id_str}\n")
                f_out.write(header)
                f_out.write(f"{var_window}\n")

            except Exception as e:
                logging.error(f"Failed to write FASTA for {somatic_var.get('p_short')}: {e}")

    logging.info("FASTA generation complete.")


if __name__ == "__main__":
    main()