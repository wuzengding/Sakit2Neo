#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
filter_vcf_for_ase.py: A robust VCF filter for GATK ASEReadCounter compatibility.
"""

import pysam
import argparse
import sys
from collections import OrderedDict
import logging

# --- Configure logging ---
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    stream=sys.stderr
)

def get_samples_and_tumor(vcf_header):
    """Robustly identifies all sample names and the tumor sample name."""
    all_samples = list(vcf_header.samples)
    if all_samples:
        tumor_name = all_samples[-1]
        logging.info(f"Identified samples via pysam: {all_samples}, Tumor: '{tumor_name}'.")
        return all_samples, tumor_name
    logging.warning("pysam header.samples is empty. Falling back to manual #CHROM parsing.")
    header_str = str(vcf_header)
    chrom_line = next((line for line in header_str.strip().split('\n') if line.startswith('#CHROM')), None)
    if not chrom_line: raise ValueError("VCF header does not contain a '#CHROM' line.")
    header_fields = chrom_line.strip().split('\t')
    if len(header_fields) < 10: raise ValueError("VCF '#CHROM' line has no sample columns.")
    all_samples = header_fields[9:]
    tumor_name = all_samples[-1]
    logging.info(f"Identified samples via fallback: {all_samples}, Tumor: '{tumor_name}'.")
    return all_samples, tumor_name

def add_to_buffer(buffer, record, tumor_name):
    """Adds a record to the buffer, resolving conflicts based on tumor AD."""
    key = (record.chrom, record.start)
    try:
        ad = record.samples[tumor_name].get('AD')
        new_vd = ad[1] if ad and len(ad) > 1 and ad[1] is not None else 0
    except (KeyError, IndexError): new_vd = 0
    
    current = buffer.get(key)
    if not current or new_vd > current.get('vd', 0):
        buffer[key] = {'record': record, 'vd': new_vd}

def flush_buffer(buffer, vcf_out, current_chrom, current_pos_0based, max_dist=1000):
    """Writes records from the buffer to the output VCF."""
    keys_to_flush = [
        k for k in list(buffer.keys())
        if k[0] != current_chrom or (current_pos_0based - k[1]) > max_dist
    ]
    for key in sorted(keys_to_flush):
        if key in buffer:
            vcf_out.write(buffer[key]['record'])
            del buffer[key]

def main(args):
    """Main function to execute the VCF filtering."""
    logging.info(f"Starting VCF processing.")
    logging.info(f"Input VCF: {args.input_vcf}")
    logging.info(f"Output VCF: {args.output_vcf}")

    try:
        vcf_in = pysam.VariantFile(args.input_vcf)
        vcf_out = pysam.VariantFile(args.output_vcf, 'w', header=vcf_in.header)
        all_samples, tumor_name = get_samples_and_tumor(vcf_in.header)
    except Exception as e:
        logging.error(f"Initialization failed: {e}")
        sys.exit(1)

    conflict_buffer = OrderedDict()
    last_chrom = None
    processed_records = 0

    for record in vcf_in:
        processed_records += 1
        if processed_records % 20000 == 0: logging.info(f"Processed {processed_records} records...")
        
        if last_chrom and record.chrom != last_chrom:
            flush_buffer(conflict_buffer, vcf_out, " flushing_all", float('inf'))
        flush_buffer(conflict_buffer, vcf_out, record.chrom, record.start, max_dist=args.buffer_size)
        last_chrom = record.chrom

        records_to_process = []
        
        # --- NEW SIMPLIFIED AND ROBUST LOGIC ---
        try:
            # Condition 1: Is this a true MNP that needs decomposition?
            is_true_mnp = ('clustered_events' in record.info and len(record.ref) > 1 and 
                           len(record.alts) == 1 and len(record.ref) == len(record.alts[0]))
            
            if is_true_mnp:
                logging.info(f"Decomposing MNP at {record.chrom}:{record.pos} ({record.ref} -> {record.alts[0]})")
                for i, (ref_base, alt_base) in enumerate(zip(record.ref, record.alts[0])):
                    if ref_base != alt_base:
                        # Create a copy of the original record to modify
                        snp_rec = record.copy()
                        snp_rec.start += i
                        snp_rec.ref = ref_base
                        snp_rec.alts = (alt_base,)
                        # Remove the tag that is no longer relevant
                        if 'clustered_events' in snp_rec.info:
                            del snp_rec.info['clustered_events']
                        records_to_process.append(snp_rec)
            
            # Condition 2: Is this a multi-allelic site?
            elif len(record.alts) > 1:
                tumor_ad = record.samples[tumor_name]['AD']
                if tumor_ad is None or len(tumor_ad) < 2: 
                    raise ValueError("Tumor AD field is missing or malformed.")

                alt_ads = tumor_ad[1:]
                max_ad, best_alt_idx = -1, -1
                for i, d in enumerate(alt_ads):
                    if d is not None and d > max_ad:
                        max_ad, best_alt_idx = d, i
                if best_alt_idx == -1: 
                    raise ValueError("No valid alternate allele found in tumor AD.")
                
                # Create a copy and modify its attributes
                new_rec = record.copy()
                new_rec.alts = (record.alts[best_alt_idx],)

                # Fix INFO fields where Number='A' or 'R'
                for k in list(new_rec.info):
                    try:
                        meta = vcf_in.header.info[k]
                        if meta.number in ('A', 'R') and isinstance(new_rec.info[k], (list, tuple)):
                            new_rec.info[k] = (new_rec.info[k][best_alt_idx],)
                    except (KeyError, IndexError):
                        # If something goes wrong with a field, just remove it to be safe
                        del new_rec.info[k]

                # Fix FORMAT fields for ALL samples
                for sample in all_samples:
                    # GT: rebuild based on presence of the chosen ALT
                    original_gt = record.samples[sample].get('GT', (None, None))
                    new_rec.samples[sample]['GT'] = (0, 1) if original_gt and (best_alt_idx + 1) in original_gt else (0, 0)
                    
                    # AD: rebuild to be biallelic
                    original_ad = record.samples[sample].get('AD')
                    if original_ad and len(original_ad) > best_alt_idx + 1:
                        new_rec.samples[sample]['AD'] = (original_ad[0], original_ad[best_alt_idx + 1])
                
                # Remove the multiallelic tag
                if 'multiallelic' in new_rec.info:
                    del new_rec.info['multiallelic']

                records_to_process.append(new_rec)
            
            # Condition 3: A simple, biallelic record
            else:
                records_to_process.append(record)

            # --- Universal Buffering ---
            for rec in records_to_process:
                add_to_buffer(conflict_buffer, rec, tumor_name)

        except Exception as e:
            logging.warning(f"Skipping record at {record.chrom}:{record.pos} due to processing error: {e}")
            # Add original record to buffer to avoid losing it completely if processing fails
            add_to_buffer(conflict_buffer, record, tumor_name)

    flush_buffer(conflict_buffer, vcf_out, "flushing_all", float('inf'))
    vcf_in.close()
    vcf_out.close()
    logging.info(f"Processing complete. Output written to {args.output_vcf}.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A robust VCF filter for GATK ASEReadCounter compatibility.")
    parser.add_argument("-i", "--input-vcf", required=True)
    parser.add_argument("-o", "--output-vcf", required=True)
    parser.add_argument("-b", "--buffer-size", type=int, default=1000)
    args = parser.parse_args()
    main(args)