#!/usr/bin/env python3

"""
Generate a comprehensive, annotated Excel report for somatic variants.

This script performs a deep analysis of a single VEP-annotated VCF from Mutect2,
classifying each variant as Somatic, Germline, or Noise based on a detailed,
hierarchical rule-set. It calculates an evidence score and provides detailed
justification for each classification.

Key Features:
- Classifies variants from a single input VCF.
- Handles multi-allelic sites by creating a row for each alternate allele.
- Calculates a "Coding Ratio" for transcripts affected by each variant.
- Identifies germline variants within a 54 exonic bp radius of Somatic and Noise hits.
- Produces a polished, formatted Excel report ready for manual review.
"""

import pandas as pd
import pysam
import logging
import argparse
import sys
import os
import re
import yaml
from collections import defaultdict
from openpyxl import Workbook
from openpyxl.styles import PatternFill, Font
from openpyxl.utils.dataframe import dataframe_to_rows

try:
    import pyensembl
except ImportError:
    pyensembl = None


def get_args():
    """
    Determines if running under Snakemake or standalone and returns
    an appropriate configuration object.
    """
    if 'snakemake' in globals():
        logging.info("Running in Snakemake mode.")
        if not pyensembl: raise ImportError("pyensembl is required but not installed.")
        return snakemake

    logging.info("Running in Standalone mode.")
    if not pyensembl: sys.exit("ERROR: pyensembl is required. Please `pip install pyensembl`.")
        
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    # The --germline-vcf argument is removed from help text as it's no longer used.
    parser.add_argument("--somatic-vcf", required=True, help="The primary Mutect2 VEP-annotated VCF.")
    parser.add_argument("--varscan-vcf", required=False, help="[OPTIONAL] VarScan2 VCF for validation.")
    parser.add_argument("--rna-counts-tumor", required=True, help="Tumor RNA ASE counts TSV.")
    parser.add_argument("--cancer-genes", required=True, help="List of known cancer genes.")
    parser.add_argument("--gtf", required=True, help="Reference GTF for pyensembl.")
    parser.add_argument("--fasta", required=True, help="Reference FASTA.")
    parser.add_argument("--output-xlsx", required=True, help="Output Excel report path.")
    parser.add_argument("--sample-id", required=True, help="Sample identifier (e.g., 'PATIENT1').")
    parser.add_argument("--tiering-params-yaml", required=True, help="YAML with somatic tiering parameters.")
    parser.add_argument("--log", help="Log file path. Logs to console if not provided.")
    args = parser.parse_args()

    class MockSnakemake:
        def __init__(self):
            self.input = {
                'somatic_vcf': args.somatic_vcf,
                'varscan_vcf': args.varscan_vcf if args.varscan_vcf else [],
                'rna_counts_tumor': args.rna_counts_tumor, 'cancer_genes': args.cancer_genes,
                'gtf': args.gtf, 'fasta': args.fasta
            }
            self.output = {'xlsx_report': args.output_xlsx}
            self.params = {'sample_id': args.sample_id}
            self.log = [args.log] if args.log else None
            with open(args.tiering_params_yaml, 'r') as f:
                config_yaml = yaml.safe_load(f)
                self.params['somatic_tiering_params'] = config_yaml.get("somatic_tiering")
                if not self.params['somatic_tiering_params']:
                    sys.exit(f"ERROR: 'somatic_tiering' not found in {args.tiering_params_yaml}")
    return MockSnakemake()


def get_csq_header(vcf_reader):
    if "CSQ" not in vcf_reader.header.info:
        raise ValueError("CSQ field not found in VCF header.")
    match = re.search(r"Format: ([\w|]+)", vcf_reader.header.info["CSQ"].description)
    if match: return match.group(1).split("|")
    raise ValueError("Could not parse CSQ header format string.")


def get_safe_info_field(info, key, default):
    if key in info:
        val = info.get(key)
        return val[0] if isinstance(val, (list, tuple)) else val
    return default


def get_safe_format_field(sample_format, key, default):
    val = sample_format.get(key)
    if val is None: return default
    return val[0] if isinstance(val, (list, tuple)) else val


def calculate_coding_ratio(csq_entries):
    """Calculates the ratio of protein-coding transcripts."""
    if not csq_entries: return "0//0"
    total_transcripts = len(csq_entries)
    coding_transcripts = sum(1 for e in csq_entries if e.get("BIOTYPE") == "protein_coding")
    return f"{coding_transcripts}//{total_transcripts}"


def parse_vcf_record(record, alt_allele_index, csq_header, sample_id):
    """Parses a single allele from a VCF record."""
    alt_allele = record.alts[alt_allele_index]
    chrom, pos, ref = record.chrom, record.pos, record.ref
    variant_key = f"{chrom}-{pos}-{ref}-{alt_allele}"
    info = record.info

    if f"{sample_id}_normal" in record.samples:
        n_name, t_name = f"{sample_id}_normal", f"{sample_id}_tumor"
    else:
        n_name, t_name = "NORMAL", "TUMOR"

    n_fmt, t_fmt = record.samples.get(n_name, {}), record.samples.get(t_name, {})
    n_dp, t_dp = get_safe_format_field(n_fmt, "DP", 0), get_safe_format_field(t_fmt, "DP", 0)
    n_ad, t_ad = n_fmt.get("AD", (0, 0)), t_fmt.get("AD", (0, 0))

    if isinstance(n_ad, (list, tuple)) and len(n_ad) > alt_allele_index + 1:
        n_alt_reads = n_ad[alt_allele_index + 1]
    else: n_alt_reads = 0
    
    if isinstance(t_ad, (list, tuple)) and len(t_ad) > alt_allele_index + 1:
        t_alt_reads = t_ad[alt_allele_index + 1]
    else: t_alt_reads = 0
    
    n_ref_reads = n_ad[0] if isinstance(n_ad, (list, tuple)) else 0
    t_ref_reads = t_ad[0] if isinstance(t_ad, (list, tuple)) else 0
    
    n_vaf = n_alt_reads / n_dp if n_dp > 0 else 0
    t_vaf = t_alt_reads / t_dp if t_dp > 0 else 0
    
    all_csq = [dict(zip(csq_header, e.split("|"))) for e in info.get("CSQ", [])]
    
    allele_csq = [e for e in all_csq if e.get("Allele") == alt_allele]

    canonical = next((e for e in allele_csq if e.get("CANONICAL") == "YES"), None)
    if not canonical:
        canonical = next((e for e in allele_csq if e.get("BIOTYPE") == "protein_coding"), allele_csq[0] if allele_csq else {})
    if not canonical: canonical = {}

    hgvsc = canonical.get("HGVSc", "/")
    hgvsp = canonical.get("HGVSp", "/")
    
    return {
        "VariantKey": variant_key, "CHROM": chrom, "POS": pos, "ID": record.id, "REF": ref, "ALT": alt_allele,
        "Normal_DP": n_dp, "Normal_AD": f"{n_ref_reads},{n_alt_reads}", "Normal_VAF": f"{n_vaf:.4f}",
        "Tumor_DP": t_dp, "Tumor_AD": f"{t_ref_reads},{t_alt_reads}", "Tumor_VAF": f"{t_vaf:.4f}",
        "TLOD": get_safe_info_field(info, "TLOD", 0.0), "gnomAD_AF": get_safe_info_field(info, "gnomAD_AF", 0.0),
        "COSMIC_CNT": get_safe_info_field(info, "COSMIC_CNT", 0), "Gene": canonical.get("SYMBOL", "/"),
        "Transcript": canonical.get("Feature", "/"), "HGVSc": hgvsc.split(':')[-1] if ':' in hgvsc else hgvsc,
        "HGVSp": hgvsp.split(':')[-1] if ':' in hgvsp else hgvsp, "Consequence": canonical.get("Consequence", "/"),
        "Impact": canonical.get("Impact", "/"), "coding_ratio": calculate_coding_ratio(allele_csq)
    }


def classify_variant_(variant, rna_support, varscan_support, params):
    """
    Classifies a variant as Somatic, Germline, or Noise using a hierarchical,
    depth-aware rule-set and calculates an evidence score.
    """
    p = {
        'MIN_TUMOR_DP': 3, 
        'MIN_NORMAL_DP': 8,
        'MIN_NORMAL_VD': 3,
        'MIN_NORMAL_VAF_GERMLINE': 0.20,
        'MAX_GNOMAD_AF_GERMLINE': 0.01,
        'NORMAL_VAF_FILTERS': [
            {'max_depth': 10, 'max_vd': 1},
            {'max_depth': 50, 'max_vd': 2},
            {'max_depth': 100, 'max_vd': 3},
            {'max_depth': 500, 'max_vd': 5},
            {'max_depth': 1000, 'max_vd': 8}
        ],
        'MIN_TUMOR_VAF_SOMATIC': 0.01,
        'MIN_TLOD_SOMATIC': 6.0,
        'MAX_GNOMAD_AF_SOMATIC': 0.001,
        'MIN_RNA_ALT_READS_SUPPORT': 3,
        'SCORE_BASE_SOMATIC': 5,
        'SCORE_BONUS_COSMIC': 3,
        'SCORE_BONUS_RNA': 2,
        'SCORE_BONUS_VARSCAN': 1,
        'SCORE_PENALTY_GERMLINE_GNOMAD': -2
    }
    n_vaf, t_vaf = float(variant["Normal_VAF"]), float(variant["Tumor_VAF"])
    n_dp, t_dp = variant["Normal_DP"], variant["Tumor_DP"]
    n_ad, t_ad = int(variant["Normal_AD"].split(",")[1]), int(variant["Tumor_AD"].split(",")[1])
    gnomad_af = float(variant["gnomAD_AF"]) if variant["gnomAD_AF"] != "/" else 0.0
    tlod = float(variant["TLOD"])
    cosmic_cnt = int(variant["COSMIC_CNT"])
    rna_alt_reads = rna_support.get('alt_reads', 0)
    details, score = [], 0

    if n_vaf >= p['MIN_NORMAL_VAF_GERMLINE'] and n_dp >= p['MIN_NORMAL_DP']:
        details.append(f"High_Normal_VAF({n_vaf:.2f})")
        if gnomad_af > p['MAX_GNOMAD_AF_GERMLINE']:
            details.append(f"Common_gnomAD_AF({gnomad_af:.4f})")
            score += p['SCORE_PENALTY_GERMLINE_GNOMAD']
        return "Germline", score, "; ".join(details)

    is_somatic_candidate = True
    reasons_for_failure = []
    
    pass_normal_vaf_check = False
    for rule in p['NORMAL_VAF_FILTERS']:
        if n_dp <= rule['max_depth']:
            if n_vd <= rule['max_vd']:
                pass_normal_vaf_check = True
            details.append(f"N_Filter(DP<={rule['max_depth']},VD<={rule['max_vd']})")
            break
    if not pass_normal_vaf_check:
        is_somatic_candidate = False
        reasons_for_failure.append(f"Fail_N_VAF({n_vaf:.2f})")


    if n_vd == 0:
        if t_vaf < p['MIN_TUMOR_VAF_SOMATIC']:
            if t_vd < p['MIN_TUMOR_DP']*2:
                is_somatic_candidate = False
                reasons_for_failure.append("filter1")
        elif t_vaf < p['MIN_TUMOR_VAF_SOMATIC']*2: ##这里需要两倍VAF
            if t_vd < p['MIN_TUMOR_DP']:
                is_somatic_candidate = False
                reasons_for_failure.append("filter2")
    elif n_vd <= p['MIN_NORMAL_VD']: ## 0<n_vd<3
        if t_vaf < p['MIN_TUMOR_VAF_SOMATIC']:
            if t_vd < p['MIN_TUMOR_DP']*2 or tlod < p['MIN_TLOD_SOMATIC']:
                is_somatic_candidate = False
                reasons_for_failure.append("filter3")
        elif t_vaf < p['MIN_TUMOR_VAF_SOMATIC']*2: ##这里需要两倍VAF
            if t_vd < p['MIN_TUMOR_DP'] or tlod < p['MIN_TLOD_SOMATIC']:
                is_somatic_candidate = False 
                reasons_for_failure.append("filter4")
    elif n_vd >= p['MIN_NORMAL_VD']: ## 0<n_vd<3
        if t_vaf < p['MIN_TUMOR_VAF_SOMATIC']*2:
            if t_vd < p['MIN_TUMOR_DP']*3 or tlod < p['MIN_TLOD_SOMATIC']*10:
                is_somatic_candidate = False
                reasons_for_failure.append("filter5")
        elif t_vaf < p['MIN_TUMOR_VAF_SOMATIC']*3: ##这里需要两倍VAF
            if t_vd < p['MIN_TUMOR_DP']*3 or tlod < p['MIN_TLOD_SOMATIC']*20:
                is_somatic_candidate = False
                reasons_for_failure.append("filter6")      
            
    if is_somatic_candidate:
        status = "Somatic"
        score += p['SCORE_BASE_SOMATIC']
        if not details: details.append("Pass_Core_Filters")
        if cosmic_cnt > 0:
            details.append(f"COSMIC(x{cosmic_cnt})")
            score += p['SCORE_BONUS_COSMIC']
        if rna_alt_reads >= p['MIN_RNA_ALT_READS_SUPPORT']:
            details.append(f"RNA_Support({rna_alt_reads})")
            score += p['SCORE_BONUS_RNA']
        if varscan_support:
            details.append("VarScan_Support")
            score += p['SCORE_BONUS_VARSCAN']
        return status, score, "; ".join(details)
    else:
        return "Noise", -1, "; ".join(reasons_for_failure)

def classify_variant(variant, rna_support, varscan_support, params):
    """
    Classifies a variant as Somatic, Germline, or Noise using a hierarchical,
    depth-aware rule-set based on normal sample variant depth (VD).
    """
    # --- Classification Parameters (Easy to edit) ---
    p = {
        'MIN_TUMOR_DP': 3, #
        'MIN_NORMAL_DP_OVERALL': 8,
        'MIN_NORMAL_VD_LOW_THRESHOLD': 3, # 您代码中的 n_vd <= 3 的阈值
        'MIN_NORMAL_VAF_GERMLINE': 0.20,
        'MAX_GNOMAD_AF_GERMLINE': 0.01,
        'NORMAL_AD_FILTERS': [ # 使用VD（即AD[1]）而不是VAF
            {'max_depth': 10, 'max_vd': 1},
            {'max_depth': 50, 'max_vd': 2},
            {'max_depth': 100, 'max_vd': 3},
            {'max_depth': 500, 'max_vd': 5},
            {'max_depth': 9999, 'max_vd': 8}
        ],
        'MIN_TUMOR_VAF_SOMATIC': 0.01, # 您代码中的基准VAF
        'MIN_TLOD_SOMATIC': 6.0,
        'MAX_GNOMAD_AF_SOMATIC': 0.001,
        'MIN_RNA_ALT_READS_SUPPORT': 3,
        'SCORE_BASE_SOMATIC': 5,
        'SCORE_BONUS_COSMIC': 3,
        'SCORE_BONUS_RNA': 2,
        'SCORE_BONUS_VARSCAN': 1,
        'SCORE_PENALTY_GERMLINE_GNOMAD': -2
    }
    
    # --- Data Extraction & Renaming for Clarity ---
    n_vaf, t_vaf = float(variant["Normal_VAF"]), float(variant["Tumor_VAF"])
    n_dp, t_dp = variant["Normal_DP"], variant["Tumor_DP"]
    # Use VD (Variant Depth) as the primary variable name, matching your logic
    n_vd = int(variant["Normal_AD"].split(",")[1])
    t_vd = int(variant["Tumor_AD"].split(",")[1])
    gnomad_af = float(variant["gnomAD_AF"]) if variant["gnomAD_AF"] != "/" else 0.0
    tlod = float(variant["TLOD"])
    cosmic_cnt = int(variant["COSMIC_CNT"])
    rna_alt_reads = rna_support.get('alt_reads', 0)
    
    details, score = [], 0

    # =================================================================
    # HIERARCHICAL CLASSIFICATION LOGIC
    # =================================================================

    # --- Step 1: Check for clear signs of being a GERMLINE variant ---
    if n_vaf >= p['MIN_NORMAL_VAF_GERMLINE'] and n_dp >= p['MIN_NORMAL_DP_OVERALL']:
        details.append(f"High_Normal_VAF({n_vaf:.2f})")
        if gnomad_af > p['MAX_GNOMAD_AF_GERMLINE']:
            details.append(f"Common_gnomAD_AF({gnomad_af:.4f})")
            score += p['SCORE_PENALTY_GERMLINE_GNOMAD']
        return "Germline", score, "; ".join(details)

    # --- Step 2: Apply your detailed, n_vd-based filtering logic ---
    is_somatic_candidate = True
    reasons_for_failure = []
    
    # First, apply the depth-aware VD filter for the normal sample
    pass_normal_vd_check = False
    for rule in p['NORMAL_AD_FILTERS']:
        if n_dp <= rule['max_depth']:
            if n_vd <= rule['max_vd']:
                pass_normal_vd_check = True
            details.append(f"N_Filter(DP<={rule['max_depth']},VD<={rule['max_vd']})")
            break
    if not pass_normal_vd_check:
        is_somatic_candidate = False
        reasons_for_failure.append(f"Fail_Normal_Contamination(VD={n_vd})")
    
    # Then, apply the tiered tumor evidence checks
    base_vaf = p['MIN_TUMOR_VAF_SOMATIC']
    base_dp = p['MIN_TUMOR_DP'] # Using MIN_TUMOR_DP from your parameters
    base_tlod = p['MIN_TLOD_SOMATIC']

    if n_vd == 0:
        if t_vaf < base_vaf and t_vd < base_dp * 2:
            is_somatic_candidate = False
            reasons_for_failure.append("Low_Tumor_Evidence_Vs_Clean_Normal")
        elif t_vaf < base_vaf * 2 and t_vd < base_dp:
            is_somatic_candidate = False
            reasons_for_failure.append("Insufficient_Tumor_Depth_Vs_Clean_Normal")
    
    elif n_vd <= p['MIN_NORMAL_VD_LOW_THRESHOLD']: # 0 < n_vd <= 3
        if t_vaf < base_vaf and (t_vd < base_dp * 2 or tlod < base_tlod):
            is_somatic_candidate = False
            reasons_for_failure.append("Low_Tumor_Evidence_Vs_Minor_Normal_Noise")
        elif t_vaf < base_vaf * 2 and (t_vd < base_dp or tlod < base_tlod):
            is_somatic_candidate = False
            reasons_for_failure.append("Insufficient_Tumor_Signal_Vs_Minor_Normal_Noise")
            
    elif n_vd > p['MIN_NORMAL_VD_LOW_THRESHOLD']: # n_vd > 3
        if t_vaf < base_vaf * 2 and (t_vd < base_dp * 3 or tlod < base_tlod * 10):
            is_somatic_candidate = False
            reasons_for_failure.append("Insufficient_Tumor_Evidence_Vs_High_Normal_Noise")
        elif t_vaf < base_vaf * 3 and (t_vd < base_dp * 3 or tlod < base_tlod * 20):
            is_somatic_candidate = False
            reasons_for_failure.append("Indistinguishable_From_High_Normal_Noise")
    
    # Final check on gnomAD AF
    if gnomad_af > p['MAX_GNOMAD_AF_SOMATIC']:
        is_somatic_candidate = False
        reasons_for_failure.append(f"Fail_gnomAD_AF({gnomad_af:.4f})")
            
    # --- Step 3: Final Classification ---
    if is_somatic_candidate:
        status = "Somatic"
        score += p['SCORE_BASE_SOMATIC']
        if not details: details.append("Pass_Core_Filters")
        if cosmic_cnt > 0:
            details.append(f"COSMIC(x{cosmic_cnt})")
            score += p['SCORE_BONUS_COSMIC']
        if rna_alt_reads >= p['MIN_RNA_ALT_READS_SUPPORT']:
            details.append(f"RNA_Support({rna_alt_reads})")
            score += p['SCORE_BONUS_RNA']
        if varscan_support:
            details.append("VarScan_Support")
            score += p['SCORE_BONUS_VARSCAN']
        return status, score, "; ".join(details)
    else:
        # If it failed any check, it's Noise. The reasons are already collected.
        return "Noise", -1, "; ".join(reasons_for_failure)

def find_nearby_germline_exon_aware(somatic_row, germline_df, ensembl_data):
    """Finds germline variants within 54 exonic bp of a somatic variant."""
    chrom, pos, tx_id = somatic_row["CHROM"], somatic_row["POS"], somatic_row["Transcript"]
    if not tx_id or not tx_id.startswith("ENST") or tx_id == "/": return []
    try:
        transcript = ensembl_data.transcript_by_id(tx_id.split('.')[0])
        somatic_offset = transcript.spliced_offset(pos)
    except Exception: return []
    nearby_germline = []
    candidate_germlines = germline_df[(germline_df["CHROM"] == chrom) & (germline_df["POS"].between(pos - 10000, pos + 10000))]
    for _, germline_row in candidate_germlines.iterrows():
        try:
            if transcript.contains(chrom, germline_row["POS"]):
                germline_offset = transcript.spliced_offset(germline_row["POS"])
                if abs(somatic_offset - germline_offset) <= 54:
                    nearby_germline.append(germline_row.to_dict())
        except Exception: continue
    return nearby_germline


def main():
    config = get_args()
    log_file = config.log[0] if config.log else None
    log_params = {"level": logging.INFO, "format": "%(asctime)s - %(levelname)s - %(message)s", "force": True}
    if log_file: log_params["filename"] = log_file
    else: log_params["stream"] = sys.stdout
    logging.basicConfig(**log_params)

    somatic_vcf, varscan_vcf = config.input['somatic_vcf'], config.input['varscan_vcf']
    rna_counts, cancer_genes_file = config.input['rna_counts_tumor'], config.input['cancer_genes']
    gtf, fasta = config.input['gtf'], config.input['fasta']
    output_xlsx, sample_id, tiering_params = config.output['xlsx_report'], config.params['sample_id'], config.params['somatic_tiering_params']
    
    logging.info("Initializing pyensembl...")
    pyensembl_cache = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "..", "resources", "pyensembl_cache")
    logging.info(f"Setting pyensembl cache directory to: {pyensembl_cache}")
    os.makedirs(pyensembl_cache, exist_ok=True)
    os.environ['PYENSEBL_CACHE_DIR'] = pyensembl_cache
    try:
        ensembl_data = pyensembl.Genome(reference_name='GRCh38_local', annotation_name='gencode_local', gtf_path_or_url=gtf)
        ensembl_data.index(); logging.info("pyensembl initialized successfully.")
    except Exception as e:
        logging.error(f"Failed to initialize pyensembl. Exon-aware search will be skipped. Error: {e}")
        ensembl_data = None

    logging.info("Loading auxiliary data...")
    cancer_genes = set(line.strip() for line in open(cancer_genes_file))
    varscan_variants = set()
    if varscan_vcf:
        with pysam.VariantFile(varscan_vcf) as vcf:
            for r in vcf: varscan_variants.add(f"{r.chrom}-{r.pos}-{r.ref}-{r.alts[0]}")
    rna_support = {f"{r['contig']}-{r['position']}-{r['refAllele']}-{r['altAllele']}": {'alt_reads': r['altCount']} for _, r in pd.read_csv(rna_counts, sep='\t').iterrows()}

    logging.info("Parsing and classifying variants from primary VCF...")
    all_variants = []
    with pysam.VariantFile(somatic_vcf) as vcf:
        csq_header = get_csq_header(vcf)
        for record in vcf:
            for i in range(len(record.alts)):
                try:
                    parsed = parse_vcf_record(record, i, csq_header, sample_id)
                    all_variants.append(parsed)
                except Exception as e: logging.warning(f"Could not parse record {record.chrom}:{record.pos}: {e}")
    if not all_variants:
        logging.warning("No variants found. Creating empty report."); Workbook().save(output_xlsx); return
        
    df = pd.DataFrame(all_variants)
    df["VarScan_Support"] = df["VariantKey"].isin(varscan_variants) if varscan_vcf else False
    classification_results = df.apply(lambda r: classify_variant(r, rna_support.get(r["VariantKey"], {}), r["VarScan_Support"], tiering_params), axis=1)
    df[["Somatic_Status", "Evidence_Score", "Evidence_Details"]] = pd.DataFrame(classification_results.tolist(), index=df.index)

    # --- MODIFIED: New report generation logic starts here ---
    logging.info("Constructing final report based on new logic...")
    somatic_df = df[df["Somatic_Status"] == "Somatic"].copy()
    noise_df = df[df["Somatic_Status"] == "Noise"].copy()
    germline_df = df[df["Somatic_Status"] == "Germline"].copy()
    
    report_rows = []

    # Process Somatic events
    for _, somatic_row in somatic_df.iterrows():
        somatic_event = somatic_row.to_dict()
        # The "one-vote veto" logic for Manual_Select
        is_protein_affecting = somatic_event.get("HGVSp", "/") not in ["/", ""]
        somatic_event["Manual_Select"] = "yes" if is_protein_affecting else "no"
        report_rows.append(somatic_event)
        
        if ensembl_data:
            nearby_germlines = find_nearby_germline_exon_aware(somatic_row, germline_df, ensembl_data)
            for germline_variant in nearby_germlines:
                germline_variant["Somatic_Status"] = "Nearby Germline"
                # Inherit the 'yes'/'no' from the parent somatic event
                germline_variant["Manual_Select"] = somatic_event["Manual_Select"]
                report_rows.append(germline_variant)

    # Process Noise events
    for _, noise_row in noise_df.iterrows():
        noise_event = noise_row.to_dict()
        noise_event["Manual_Select"] = "no"
        report_rows.append(noise_event)
        
        if ensembl_data:
            nearby_germlines = find_nearby_germline_exon_aware(noise_row, germline_df, ensembl_data)
            for germline_variant in nearby_germlines:
                germline_variant["Somatic_Status"] = "Nearby Germline"
                germline_variant["Manual_Select"] = "no"
                report_rows.append(germline_variant)

    final_df_sheet1 = pd.DataFrame(report_rows)
    if not final_df_sheet1.empty:
        final_df_sheet1.sort_values(by=["CHROM", "POS"], inplace=True)
    
    logging.info("Generating Excel report...")
    wb = Workbook()
    ws1 = wb.active; ws1.title = "Somatic & Noise Event Report"
    
    final_cols = ["CHROM", "POS", "REF", "ALT", "Gene", "Transcript", "HGVSc", "HGVSp", "Consequence", "Impact", "coding_ratio", "Tumor_DP", "Tumor_AD", "Tumor_VAF", "Normal_DP", "Normal_AD", "Normal_VAF", "gnomAD_AF", "COSMIC_CNT", "Somatic_Status", "Manual_Select", "Evidence_Score", "Evidence_Details"]
    ws1.append(final_cols)
    for cell in ws1[1]: cell.font = Font(bold=True)
    
    somatic_fill = PatternFill("solid", fgColor="DCDCDC") # Light Grey

    if not final_df_sheet1.empty:
        for _, row in final_df_sheet1.reindex(columns=final_cols).fillna('/').iterrows():
            ws1.append(list(row))
            if row["Somatic_Status"] == "Somatic":
                for cell in ws1[ws1.max_row]: cell.fill = somatic_fill
            
    ws2 = wb.create_sheet(title="Pathogenic Germline Variants")
    pathogenic_germline_df = germline_df[germline_df["Impact"] == "HIGH"].copy()
    ws2.append(final_cols)
    for cell in ws2[1]: cell.font = Font(bold=True)
    if not pathogenic_germline_df.empty:
        for r in dataframe_to_rows(pathogenic_germline_df.reindex(columns=final_cols).fillna('/'), index=False, header=False):
            ws2.append(r)
        
    wb.save(output_xlsx)
    logging.info(f"Report generation complete. Saved to {output_xlsx}")

if __name__ == "__main__":
    main()