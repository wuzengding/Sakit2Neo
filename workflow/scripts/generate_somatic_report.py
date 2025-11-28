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


# VEP Consequence Severity Order (from most to least severe)
# Based on Ensembl's recommendations and general impact on protein function.
CONSEQUENCE_SEVERITY_ORDER = [
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'transcript_amplification',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
    'splice_donor_5th_base_variant',
    'splice_region_variant',
    'splice_donor_region_variant',
    'splice_polypyrimidine_tract_variant',
    'incomplete_terminal_codon_variant',
    'start_retained_variant',
    'stop_retained_variant',
    'synonymous_variant',
    'coding_sequence_variant',
    'mature_miRNA_variant',
    '5_prime_UTR_variant',
    '3_prime_UTR_variant',
    'non_coding_transcript_exon_variant',
    'intron_variant',
    'NMD_transcript_variant',
    'non_coding_transcript_variant',
    'coding_transcript_variant', # Added for completeness
    'upstream_gene_variant',
    'downstream_gene_variant',
    'TFBS_ablation',
    'TFBS_amplification',
    'TF_binding_site_variant',
    'regulatory_region_ablation',
    'regulatory_region_amplification',
    'regulatory_region_variant',
    'feature_elongation',
    'feature_truncation',
    'intergenic_variant',
    'sequence_variant' # Should be last
]

# Create a mapping for quick lookups. The lower the number, the more severe.
CONSEQUENCE_SEVERITY_MAP = {
    consequence: i for i, consequence in enumerate(CONSEQUENCE_SEVERITY_ORDER)
}

# Define a set of consequences that DEFINITELY alter the protein sequence.
# This excludes synonymous variants and other non-coding/benign changes.
PROTEIN_ALTERING_CONSEQUENCES = {
    'transcript_ablation',
    'splice_acceptor_variant',
    'splice_donor_variant',
    'stop_gained',
    'frameshift_variant',
    'stop_lost',
    'start_lost',
    'inframe_insertion',
    'inframe_deletion',
    'missense_variant',
    'protein_altering_variant',
}


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

def calculate_AAchaged_ratio(csq_entries):
    if not csq_entries: return "0//0"
    total_transcripts = len(csq_entries)
    protein_altering_count = 0
    for e in csq_entries:
        consequences = e.get("Consequence", "").split('&')
        if any(c in PROTEIN_ALTERING_CONSEQUENCES for c in consequences):
            protein_altering_count += 1
    return f"{protein_altering_count}//{total_transcripts}"

def get_most_severe_consequence_rank(csq_entry):
    """
    Finds the most severe consequence in a CSQ entry (which can have multiple, e.g., "A&B").
    Returns its rank from our map. Lower is more severe.
    """
    consequences = csq_entry.get("Consequence", "").split('&')
    # Assign a very high (not severe) rank as default
    min_rank = len(CONSEQUENCE_SEVERITY_ORDER) 
    for c in consequences:
        rank = CONSEQUENCE_SEVERITY_MAP.get(c, min_rank)
        if rank < min_rank:
            min_rank = rank
    return min_rank

def get_amino_acid_sub_status(consequences):
    """
    根据VEP的consequence列表，判断氨基酸的变化子类型。

    Args:
        consequences (list): 从VEP CSQ字段解析出的consequence字符串列表。

    Returns:
        str: 'No_AA_Altered', 'Stop_Gained', 'Long_AA_Altered', or 'Single_AA_Altered'.
    """
    is_protein_altering = any(c in PROTEIN_ALTERING_CONSEQUENCES for c in consequences)
    
    if not is_protein_altering:
        return "No_AA_Altered"
    
    # 如果是蛋白质改变类型，再进行细分
    if 'stop_gained' in consequences:
        return "Stop_Gained"
    elif any(c in ['frameshift_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_lost'] for c in consequences):
        return "Long_AA_Altered"
    else:
        # 包括 missense_variant, inframe_insertion, inframe_deletion 等
        return "Single_AA_Altered"

def parse_vcf_record(record, alt_allele_index, csq_header, sample_id):
    """Parses a single allele from a VCF record."""
    alt_allele = record.alts[alt_allele_index]
    chrom, pos, ref = record.chrom, record.pos, record.ref
    #print("chrom", "pos", chrom, pos, record.chrom, record.pos)
    variant_key = f"{chrom}-{pos}-{ref}-{alt_allele}"
    info = record.info

    if f"{sample_id}_normal" in record.samples:
        n_name, t_name = f"{sample_id}_normal", f"{sample_id}_tumor"
    else:
        n_name, t_name = "NORMAL", "TUMOR"

    #print(n_name, t_name)
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

    # 一个record只有一个"CSQ"字符串，但是一个CSQ，可能存在多个Allele，"
    # "比如 “chr1	1312198	.	TG	T,TGG,TGGG,TGGGG”这样的复杂突变"
    all_csq = [dict(zip(csq_header, e.split("|"))) for e in info.get("CSQ", [])]
    # 所以用Allele来选择对应突变的CSQ内容
    allele_csq = [e for e in all_csq if e.get("Allele") == alt_allele]
    
    # --- 新代码 ---
    if not allele_csq:
        canonical = {}
    else:
        # Sort all annotations for this allele based on severity, then by CANONICAL tag as a tie-breaker
        sorted_csq = sorted(
            allele_csq,
            key=lambda csq: (
                get_most_severe_consequence_rank(csq),
                csq.get("CANONICAL") != "YES" # False (is CANONICAL) comes before True
            ) ## key 返回的是一个元组（Tuple），例如 (3, False) 或 (10, True)
        )
        # The best annotation is the first one after sorting
        canonical = sorted_csq[0]

    hgvsc = canonical.get("HGVSc", "/")
    hgvsp = canonical.get("HGVSp", "/")
    
    return {
        "VariantKey": variant_key, "CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt_allele,
        "Normal_DP": n_dp, "Normal_AD": f"{n_ref_reads},{n_alt_reads}", "Normal_VAF": f"{n_vaf:.4f}",
        "Tumor_DP": t_dp, "Tumor_AD": f"{t_ref_reads},{t_alt_reads}", "Tumor_VAF": f"{t_vaf:.4f}",
        "TLOD": get_safe_info_field(info, "TLOD", 0.0), "gnomAD_AF": get_safe_info_field(info, "gnomAD_AF", 0.0),
        "COSMIC_CNT": get_safe_info_field(info, "COSMIC_CNT", 0), "Gene": canonical.get("SYMBOL", "/"),
        "Transcript": canonical.get("Feature", "/"), "HGVSc": hgvsc.split(':')[-1] if ':' in hgvsc else hgvsc,
        "HGVSp": hgvsp.split(':')[-1] if ':' in hgvsp else hgvsp, "Consequence": canonical.get("Consequence", "/"),
        "Impact": canonical.get("IMPACT", "/"), 
        "Coding_ratio": calculate_coding_ratio(allele_csq), "AA_chaged_ratio": calculate_AAchaged_ratio(allele_csq),
        "all_csq": allele_csq
    }

def classify_variant_(variant, median_T_dp, median_N_dp, params):
    """
    Classifies a variant as Somatic, Germline, or Noise using a hierarchical,
    depth-aware rule-set based on normal sample variant depth (VD).
    """
    # --- Classification Parameters (Easy to edit) ---
    p = {
        'median_dp_ratio':0.05,
        'median_dp_threshold': 25,
        
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
        'MIN_RNA_READS_SUPPORT': 1,
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
    rna_support_reads = int(variant["RNA_Support_Reads"])
    varscan_support = variant["VarScan_Support"]
    
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
        if rna_support_reads >= p['MIN_RNA_READS_SUPPORT']:
            details.append(f"RNA_Support({rna_support_reads})")
            score += p['SCORE_BONUS_RNA']
        if varscan_support:
            details.append("VarScan_Support")
            score += p['SCORE_BONUS_VARSCAN']
        return status, score, "; ".join(details)
    else:
        # If it failed any check, it's Noise. The reasons are already collected.
        return "Noise", -1, "; ".join(reasons_for_failure)

def classify_variant_new(variant, params):
    """
    根据思维导图中的规则对变异进行分类 (增加了RNA捞回逻辑并重构了代码)。
    
    返回:
        tuple: (primary_status, sub_status, evidence_details)
    """
    # --- 参数提取 ---
    p = params

    # --- 数据提取 ---
    n_vaf, t_vaf = float(variant["Normal_VAF"]), float(variant["Tumor_VAF"])
    n_dp, t_dp = variant["Normal_DP"], variant["Tumor_DP"]
    n_vd = int(variant["Normal_AD"].split(",")[1])
    t_vd = int(variant["Tumor_AD"].split(",")[1])
    gnomad_af = float(variant["gnomAD_AF"]) if variant["gnomAD_AF"] != "/" else 0.0
    tlod = float(variant["TLOD"])
    consequences = variant.get("Consequence", "").split('&')
    rna_support_reads = int(variant.get("RNA_Support_Reads",0))
    
    details = []

    # --- 分类逻辑 (最终版) ---

    # === 规则 1: 优先判断明确的 Germline/LOH 信号 ===
    if n_vaf >= p.get('min_normal_vaf_germline', 0.18) and n_dp >= p.get('min_normal_dp_germline', 10):
        details.append(f"Germline_Signal(N_VAF={n_vaf:.2f}, N_DP={n_dp})")
        
        if (0.3 <= n_vaf <= 0.7) and (t_vaf > 0.9 or t_vaf < 0.1):
            primary_status = "LOH"
            details.append(f"LOH_Shift(T_VAF={t_vaf:.2f})")
        else:
            primary_status = "Germline"
            
        # 调用辅助函数获取子状态
        sub_status = get_amino_acid_sub_status(consequences)
        return primary_status, sub_status, "; ".join(details)

    # === 规则 2: 噪音 (Noise) 判定与RNA捞回 ===
    is_noise = False
    noise_reasons = []
    
    if t_vd < p.get('min_tumor_vd_somatic', 3):
        is_noise = True; noise_reasons.append(f"Low_T_VD({t_vd})")
    if t_vaf < p.get('min_tumor_vaf_somatic', 0.02):
        is_noise = True; noise_reasons.append(f"Low_T_VAF({t_vaf:.2f})")
    if n_vd > p.get('max_normal_vd_somatic', 3):
        is_noise = True; noise_reasons.append(f"High_N_VD_Contamination({n_vd})")
    if tlod < p.get('min_tlod_somatic', 6.0):
        is_noise = True; noise_reasons.append(f"Low_TLOD({tlod})")
    if gnomad_af > p.get('max_gnomad_af_somatic', 0.001):
        is_noise = True; noise_reasons.append(f"High_gnomAD_AF({gnomad_af:.4f})")
        
    # RNA捞回
    if is_noise:
        # 检查是否有足够的RNA证据来“捞回”这个变异
        if rna_support_reads >= p.get('min_rna_reads_for_rescue', 2):
            details.append(f"Rescued_from_Noise(RNA_Reads={rna_support_reads})")
            # 在捞回后，重新进行简化的 Germline vs Somatic 判断
            # 这里我们使用一个稍微宽松的标准来判断是否是Germline，因为它最初信号不强
            # 例如，只要n_vaf > 0.1 就可以认为是Germline的迹象
            if n_vaf >= p.get('min_rescue_n_vaf_germline', 0.1):
                primary_status = "Germline"
                details.append("Classified_as_Germline_Post_Rescue")
            else:
                primary_status = "Somatic"
                details.append("Classified_as_Somatic_Post_Rescue")
            
            sub_status = get_amino_acid_sub_status(consequences)
            return primary_status, sub_status, "; ".join(details)
            
        else:
            # 如果没有足够的RNA证据，则确定为Noise
            return "Noise", None, "; ".join(noise_reasons)

    # === 规则 3: 如果不是Germline也不是Noise，则认为是Somatic ===
    # 这个代码块现在会处理“原生Somatic”和“被RNA捞回的Somatic”
    primary_status = "Somatic"
    if not details: # 如果不是被捞回的，就加上默认的通过信息
        details.append("Pass_Somatic_Filters")
    
    # 调用辅助函数获取子状态
    sub_status = get_amino_acid_sub_status(consequences)
            
    return primary_status, sub_status, "; ".join(details)

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

def load_auxiliary_data(config):
    """加载所有辅助数据文件。"""
    logging.info("Loading auxiliary data...")
    varscan_vcf_path = config.input.get('varscan_vcf')
    rna_counts_path = config.input['rna_counts_tumor']

    varscan_variants = set()
    if varscan_vcf_path:
        with pysam.VariantFile(varscan_vcf_path) as vcf:
            for r in vcf:
                varscan_variants.add(f"{r.chrom}-{r.pos}-{r.ref}-{r.alts[0]}")
    
    rna_support_df = pd.read_csv(rna_counts_path, sep='\t')
    rna_support = {
        f"{r['contig']}-{r['position']}-{r['refAllele']}-{r['altAllele']}": r['altCount']
        for _, r in rna_support_df.iterrows()
    }
    
    return varscan_variants, rna_support


def parse_and_classify_variants(config, varscan_variants, rna_support):
    """解析VCF，添加辅助信息，并对变异进行分类。"""
    logging.info("Parsing and classifying variants from primary VCF...")
    somatic_vcf_path = config.input['somatic_vcf']
    sample_id = config.params['sample_id']
    tiering_params = config.params['somatic_tiering_params']

    all_variants = []
    with pysam.VariantFile(somatic_vcf_path) as vcf:
        csq_header = get_csq_header(vcf)
        for record in vcf:
            for i in range(len(record.alts)):
                try:
                    parsed = parse_vcf_record(record, i, csq_header, sample_id)
                    all_variants.append(parsed)
                except Exception as e:
                    logging.warning(f"Could not parse record {record.chrom}:{record.pos}: {e}")
    
    if not all_variants:
        logging.warning("No variants found in VCF.")
        return pd.DataFrame()

    df = pd.DataFrame(all_variants)
    df["VarScan_Support"] = df["VariantKey"].isin(varscan_variants)
    df['RNA_Support_Reads'] = df['VariantKey'].map(rna_support).fillna(0).astype(int)

    classification_results = df.apply(
        lambda r: classify_variant_new(r, tiering_params), 
        axis=1
    )
    df[["Primary_Status", "Sub_Status", "Evidence_Details"]] = pd.DataFrame(
        classification_results.tolist(), index=df.index
    )
    
    return df

def parse_hgvsp_for_aa_pos(hgvsp_str):
    """从HGVSp字符串中解析出氨基酸位置。例如 'p.Val600Glu' -> 600"""
    if not isinstance(hgvsp_str, str) or not hgvsp_str.startswith('p.'):
        return None
    match = re.search(r'p\.[A-Z][a-z]{2}(\d+)', hgvsp_str)
    if match:
        return int(match.group(1))
    return None

def find_companion_germline_map(df):
    """
    使用稳健的蛋白质感知方法，查找与Somatic变异相邻的Germline变异。
    逻辑更新：不再只取最近的一个，而是取18aa范围内所有的Somatic，
    以便在计算Manual_Select时，只要周围有一个高质量Somatic，该Germline就能被选中。
    """
    logging.info("Finding nearby germline variants using a robust protein-aware method...")
    
    protein_map_list = []
    # 筛选出包含有效 'all_csq' 列表的行
    valid_rows = df[df['all_csq'].apply(lambda x: isinstance(x, list) and len(x) > 0)]
    
    for _, row in valid_rows.iterrows():
        if pd.isna(row.get('Primary_Status')): continue

        for csq in row['all_csq']:
            hgvsp = csq.get('HGVSp')
            if hgvsp and ':' in hgvsp:
                protein_id = hgvsp.split(':')[0]
                aa_pos = parse_hgvsp_for_aa_pos(hgvsp.split(':')[1])
                if protein_id and aa_pos is not None:
                    protein_map_list.append({
                        'VariantKey': row['VariantKey'],
                        'Primary_Status': row['Primary_Status'],
                        'Sub_Status': row['Sub_Status'],
                        'protein_id': protein_id,
                        'aa_pos': aa_pos
                    })

    if not protein_map_list:
        return {}

    protein_map_df = pd.DataFrame(protein_map_list)
    
    # 1. 筛选 Somatic：必须改变氨基酸
    somatic_map = protein_map_df[
        (protein_map_df['Primary_Status'] == 'Somatic') & 
        (protein_map_df['Sub_Status'] != 'No_AA_Altered')
    ]

    # 2. 筛选 Germline：必须改变氨基酸
    germline_map = protein_map_df[
        (protein_map_df['Primary_Status'] == 'Germline') & 
        (protein_map_df['Sub_Status'] != 'No_AA_Altered')
    ]

    if somatic_map.empty or germline_map.empty:
        return {}
        
    merged = pd.merge(somatic_map, germline_map, on='protein_id', suffixes=('_somatic', '_germline'))
    merged['distance'] = abs(merged['aa_pos_somatic'] - merged['aa_pos_germline'])
    
    # === 核心修改点 ===
    # 保留所有距离 <= 18 的组合，而不是只保留最近的一个
    companions = merged[merged['distance'] <= 18]
    
    if companions.empty:
        return {}

    # 聚合：一个 Germline 可能对应多个 Somatic
    companion_dict = defaultdict(set)
    for _, row in companions.iterrows():
        g_key = row['VariantKey_germline']
        s_key = row['VariantKey_somatic']
        companion_dict[g_key].add(s_key)
    
    # 转换为列表返回
    result_map = {k: list(v) for k, v in companion_dict.items()}

    logging.info(f"Found {len(result_map)} germline variants paired with somatic events.")
    return result_map


def generate_reports(df, companion_germline_map, config):
    """根据分类结果和伴侣映射，生成两个Excel报告。"""
    logging.info("Generating final Excel reports...")
    output_filtered_xlsx = config.output['xlsx_report_filter']
    output_neopeptides_xlsx = config.output['xlsx_report_Somatic']

    # === 定义浅黄色样式 (Light Yellow) ===
    # FFF2CC 是 Excel 中标准的 "浅黄色 40%"，视觉效果比较柔和
    germline_highlight_fill = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")
    
    # === 1. 初始化新列 ===
    df['Manual_Select'] = 'no'
    df['Companion_To'] = '/'  # 新增列：用于标记Germline对应的Somatic Key

    # === 2. 处理 Somatic 的 Manual_Select ===
    # 规则：Somatic 且 RNA_Support_Reads > 0 -> yes
    somatic_yes_mask = (df['Primary_Status'] == 'Somatic') & (df['RNA_Support_Reads'] > 0)
    df.loc[somatic_yes_mask, 'Manual_Select'] = 'yes'
    
    # 创建一个快速查询 Somatic Manual_Select 状态的字典
    # Key: VariantKey, Value: 'yes'/'no'
    somatic_manual_status_map = df[df['Primary_Status'] == 'Somatic'].set_index('VariantKey')['Manual_Select'].to_dict()

    # === 3. 处理 Germline 的 Manual_Select 和 Companion 信息 ===
    # 遍历我们在 find_companion_germline_map 中生成的字典
    for germline_key, somatic_partners in companion_germline_map.items():
        # 找到对应的 Germline 行索引
        germline_idx = df[df['VariantKey'] == germline_key].index
        
        if not germline_idx.empty:
            # 3.1 填充 Companion_To 列 (把列表转换成字符串，如 "key1;key2")
            partners_str = ";".join(somatic_partners)
            df.loc[germline_idx, 'Companion_To'] = partners_str
            
            # 3.2 判断 Manual_Select
            # 规则：如果伴侣 Somatic 中有任意一个是 'yes'，则该 Germline 为 'yes'
            is_companion_yes = False
            for s_key in somatic_partners:
                if somatic_manual_status_map.get(s_key) == 'yes':
                    is_companion_yes = True
                    break
            
            if is_companion_yes:
                df.loc[germline_idx, 'Manual_Select'] = 'yes'

    # === 4. 定义最终输出的列顺序 (Requirements 2 & 3) ===
    # 加入了 'Manual_Select' 在 VarScan_Support 之前
    # 加入了 'Companion_To' (建议放在 Gene 附近或者最后，这里放在最后方便查看)
    final_cols = [
        "VariantKey", "CHROM", "POS", "REF", "ALT", 
        "Gene", "Transcript", "HGVSc", "HGVSp", 
        "Consequence", "Impact", 
        "Coding_ratio", "AA_chaged_ratio", "Manual_Select", # <-- 插入在这里
        "VarScan_Support", "RNA_Support_Reads", 
        "Normal_DP", "Normal_AD", "Normal_VAF", 
        "Tumor_DP", "Tumor_AD", "Tumor_VAF", 
        "TLOD", "gnomAD_AF", "COSMIC_CNT", 
        "Primary_Status", "Sub_Status", "Evidence_Details",
        "Companion_To" # <-- 新增列，方便溯源
    ]

    # 确保 DataFrame 只包含存在的列 (防止某些列未生成报错)
    available_cols = [c for c in final_cols if c in df.columns]

    # === 5. 分配 Sheet 逻辑 (保持之前的修正) ===
    somatic_neopeptides_sheets = defaultdict(list)
    filtered_variants_sheets = defaultdict(list)

    # 5.1 基础分类
    filtered_variants_sheets['Noise'] = df[df['Primary_Status'] == 'Noise']
    filtered_variants_sheets['LOH'] = df[df['Primary_Status'] == 'LOH']
    filtered_variants_sheets['Somatic_without_AA_altered'] = df[(df['Primary_Status'] == 'Somatic') & (df['Sub_Status'] == 'No_AA_Altered')]
    filtered_variants_sheets['Germline_without_AA_altered'] = df[(df['Primary_Status'] == 'Germline') & (df['Sub_Status'] == 'No_AA_Altered')]

    # 5.2 改变AA的 Somatic
    somatic_aa_altered = df[(df['Primary_Status'] == 'Somatic') & (df['Sub_Status'] != 'No_AA_Altered')]
    somatic_neopeptides_sheets['low scale variants'].append(somatic_aa_altered[somatic_aa_altered['Sub_Status'].isin(['Single_AA_Altered', 'Stop_Gained'])])
    somatic_neopeptides_sheets['Large scale variants'].append(somatic_aa_altered[somatic_aa_altered['Sub_Status'] == 'Long_AA_Altered'])

    # 5.3 改变AA的 Germline
    germline_aa_altered = df[(df['Primary_Status'] == 'Germline') & (df['Sub_Status'] != 'No_AA_Altered')]
    
    isolated_germlines = []
    for _, row in germline_aa_altered.iterrows():
        if row['VariantKey'] in companion_germline_map:
            # 有伴侣的 Germline
            # 根据伴侣的类型决定去哪个 Sheet (这里取第一个伴侣的类型简单判断，或者根据优先级)
            # 为了简单起见，我们查看其伴侣Somatic Sub_Status
            somatic_partners = companion_germline_map[row['VariantKey']]
            
            # 获取伴侣的 Sub_Status (从 df 中反查)
            # 这里取伴侣中最严重的类型决定 Sheet
            partner_statuses = df[df['VariantKey'].isin(somatic_partners)]['Sub_Status'].tolist()
            
            if any(s == 'Long_AA_Altered' for s in partner_statuses):
                target_sheet = 'Large scale variants'
            else:
                target_sheet = 'low scale variants'
                
            somatic_neopeptides_sheets[target_sheet].append(pd.DataFrame([row]))
        else:
            isolated_germlines.append(row)
    
    if isolated_germlines:
        filtered_variants_sheets['Germline_with_AA_altered'] = pd.DataFrame(isolated_germlines)
        
    # === 6. 写入 Somatic 报告并应用高亮样式 ===
    with pd.ExcelWriter(output_neopeptides_xlsx, engine='openpyxl') as writer:
        logging.info(f"Writing to {output_neopeptides_xlsx}...")
        for sheet_name, df_list in somatic_neopeptides_sheets.items():
            if df_list:
                # 准备数据：合并、去重、排序、重置索引
                # 重置索引非常重要，这样我们遍历数据行时才能对应到 Excel 的行号
                final_df = pd.concat(df_list).drop(columns=['all_csq'], errors='ignore')
                final_df = final_df.reindex(columns=available_cols).fillna('/').sort_values(by=["CHROM", "POS"]).reset_index(drop=True)
                
                # 写入 Excel
                final_df.to_excel(writer, sheet_name=sheet_name, index=False)
                
                # 获取 Worksheet 对象
                ws = writer.sheets[sheet_name]
                
                # --- 应用高亮逻辑 ---
                # 遍历刚才写入的 DataFrame
                for i, row in final_df.iterrows():
                    # 判断条件：是 Germline 且 Companion_To 不为 '/'
                    if row['Primary_Status'] == 'Germline' and row.get('Companion_To', '/') != '/':
                        # Excel 行号 = DataFrame索引(i) + 标题行(1) + 1 = i + 2
                        excel_row_idx = i + 2
                        # 遍历该行的所有列进行填色
                        for col_idx in range(1, len(available_cols) + 1):
                            cell = ws.cell(row=excel_row_idx, column=col_idx)
                            cell.fill = germline_highlight_fill

    # === 7. 写入 Filtered 报告 (保持原样，无需高亮) ===
    with pd.ExcelWriter(output_filtered_xlsx, engine='openpyxl') as writer:
        logging.info(f"Writing to {output_filtered_xlsx}...")
        for sheet_name, df_data in filtered_variants_sheets.items():
            if not df_data.empty:
                final_df = df_data.drop(columns=['all_csq'], errors='ignore')
                final_df = final_df.reindex(columns=available_cols).fillna('/').sort_values(by=["CHROM", "POS"])
                final_df.to_excel(writer, sheet_name=sheet_name, index=False)

def main():
    """主流程协调器。"""
    config = get_args()
    log_file = config.log[0] if config.log else None
    log_params = {"level": logging.INFO, "format": "%(asctime)s - %(levelname)s - %(message)s", "force": True}
    if log_file: 
        log_params["filename"] = log_file
    else: 
        log_params["stream"] = sys.stdout
    logging.basicConfig(**log_params)

    # 步骤 1: 加载辅助数据
    varscan_variants, rna_support = load_auxiliary_data(config)

    # 步骤 2: 解析VCF并分类变异
    classified_df = parse_and_classify_variants(config, varscan_variants, rna_support)
    if classified_df.empty:
        logging.info("No variants to process. Exiting.")
        return

    # 步骤 3: 查找伴侣Germline变异
    companion_map = find_companion_germline_map(classified_df)

    # 步骤 4: 生成报告
    generate_reports(classified_df, companion_map, config)

    logging.info("Script finished successfully.")

if __name__ == "__main__":
    main()