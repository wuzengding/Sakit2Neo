import pandas as pd
import json
import argparse
import sys
import os

# === 定义所有预期的列名，确保输出结构固定 ===
DNA_COLS = [
    'Sample', 'Type',
    'total_reads', 'total_bases', 'q30_bases', 'q30_rate', 'gc_content', 
    'total_bases_after', 'filtered_ratio', 'insert_size', 
    'read1_mean_length', 'read2_mean_length', 
    'PF_READS_ALIGNED', 'PCT_PF_READS_ALIGNED', 'STRAND_BALANCE', 'PCT_SOFTCLIP', 
    'PCT_USABLE_BASES_ON_TARGET', 'FOLD_ENRICHMENT', 
    'MEAN_TARGET_COVERAGE', 'MEDIAN_TARGET_COVERAGE', 'PCT_TARGET_BASES_20X',
    'PCT_TARGET_BASES_50X', 'PCT_TARGET_BASES_100X', 'PCT_TARGET_BASES_250X'
]

RNA_COLS = [
    'Sample', 'Type',
    'total_reads', 'total_bases', 'q30_bases', 'q30_rate', 'gc_content', 
    'total_bases_after', 'filtered_ratio', 'insert_size', 
    'read1_mean_length', 'read2_mean_length', 
    'PF_ALIGNED_BASES', 'RIBOSOMAL_BASES', 'CODING_BASES', 'UTR_BASES', 
    'INTRONIC_BASES', 'INTERGENIC_BASES', 'PCT_RIBOSOMAL_BASES', 
    'PCT_CODING_BASES', 'PCT_UTR_BASES', 'PCT_INTRONIC_BASES', 
    'PCT_INTERGENIC_BASES', 'PCT_MRNA_BASES'
]

def is_valid_file(path):
    """检查路径是否有效"""
    if not path: return False
    if path == "None": return False
    if not os.path.exists(path): return False
    # 检查文件是否为空
    if os.path.getsize(path) == 0: return False
    return True

def parse_fastp_json(json_path):
    if not is_valid_file(json_path): return {}
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        summary = data.get('summary', {})
        before = summary.get('before_filtering', {})
        after = summary.get('after_filtering', {})
        total_bases = before.get('total_bases', 0)
        total_bases_after = after.get('total_bases', 0)
        filtered_ratio = total_bases_after / total_bases if total_bases > 0 else 0

        return {
            'total_reads': before.get('total_reads'),
            'total_bases': total_bases,
            'q30_bases': before.get('q30_bases'),
            'q30_rate': before.get('q30_rate'),
            'gc_content': before.get('gc_content'),
            'total_bases_after': total_bases_after,
            'filtered_ratio': filtered_ratio,
            'read1_mean_length': before.get('read1_mean_length'),
            'read2_mean_length': before.get('read2_mean_length'),
            'insert_size_fastp': data.get('insert_size', {}).get('peak')
        }
    except Exception as e:
        print(f"Warning: Failed to parse fastp {json_path}: {e}")
        return {}

def parse_mosdepth_summary(summary_path):
    if not is_valid_file(summary_path): return {}
    try:
        df = pd.read_csv(summary_path, sep='\t')
        total_row = df[df['chrom'] == 'total_region']
        if not total_row.empty:
            return {'mean_depth': total_row.iloc[0]['mean']}
    except Exception:
        pass
    return {}

def read_picard_metrics(file_path):
    if not is_valid_file(file_path): return pd.DataFrame()
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
        start_line = 0
        for i, line in enumerate(lines):
            if not line.startswith('#') and line.strip() and ("CATEGORY" in line or "BAIT_SET" in line or "PF_BASES" in line):
                start_line = i
                break
        return pd.read_csv(file_path, sep='\t', skiprows=start_line)
    except Exception:
        return pd.DataFrame()

def parse_picard_alignment(file_path):
    df = read_picard_metrics(file_path)
    if df.empty: return {}
    
    # 尝试寻找 PAIR 或 UNPAIRED
    if 'CATEGORY' in df.columns:
        row = df[df['CATEGORY'] == 'PAIR']
        if row.empty: row = df[df['CATEGORY'] == 'UNPAIRED']
        if row.empty: row = df.iloc[0:1] # Fallback
    else:
        row = df.iloc[0:1]
        
    if row.empty: return {}
    r = row.iloc[0]
    return {
        'PF_READS_ALIGNED': r.get('PF_READS_ALIGNED'),
        'PCT_PF_READS_ALIGNED': r.get('PCT_PF_READS_ALIGNED'),
        'STRAND_BALANCE': r.get('STRAND_BALANCE'),
        'PCT_SOFTCLIP': r.get('PCT_SOFTCLIP'),
    }

def parse_picard_hs(file_path):
    df = read_picard_metrics(file_path)
    if df.empty: return {}
    r = df.iloc[0]
    return {
        'PCT_USABLE_BASES_ON_TARGET': r.get('PCT_USABLE_BASES_ON_TARGET'),
        'FOLD_ENRICHMENT': r.get('FOLD_ENRICHMENT'),
        'MEAN_TARGET_COVERAGE': r.get('MEAN_TARGET_COVERAGE'),
        'MEDIAN_TARGET_COVERAGE': r.get('MEDIAN_TARGET_COVERAGE'),
        'PCT_TARGET_BASES_20X': r.get('PCT_TARGET_BASES_20X'),
        'PCT_TARGET_BASES_50X': r.get('PCT_TARGET_BASES_50X'),
        'PCT_TARGET_BASES_100X': r.get('PCT_TARGET_BASES_100X'),
        'PCT_TARGET_BASES_250X': r.get('PCT_TARGET_BASES_250X')
    }

def parse_picard_rna(file_path):
    df = read_picard_metrics(file_path)
    if df.empty: return {}
    r = df.iloc[0]
    return {
        'PF_ALIGNED_BASES': r.get('PF_ALIGNED_BASES'),
        'RIBOSOMAL_BASES': r.get('RIBOSOMAL_BASES'),
        'CODING_BASES': r.get('CODING_BASES'),
        'UTR_BASES': r.get('UTR_BASES'),
        'INTRONIC_BASES': r.get('INTRONIC_BASES'),
        'INTERGENIC_BASES': r.get('INTERGENIC_BASES'),
        'PCT_RIBOSOMAL_BASES': r.get('PCT_RIBOSOMAL_BASES'),
        'PCT_CODING_BASES': r.get('PCT_CODING_BASES'),
        'PCT_UTR_BASES': r.get('PCT_UTR_BASES'),
        'PCT_INTRONIC_BASES': r.get('PCT_INTRONIC_BASES'),
        'PCT_INTERGENIC_BASES': r.get('PCT_INTERGENIC_BASES'),
        'PCT_MRNA_BASES': r.get('PCT_MRNA_BASES')
    }

def parse_picard_insert_metrics(file_path):
    df = read_picard_metrics(file_path)
    if df.empty: return {}
    return {'insert_size': df.iloc[0].get('MEDIAN_INSERT_SIZE')}

def main():
    parser = argparse.ArgumentParser()
    # 所有参数设为非必须，默认 None
    # DNA Inputs
    parser.add_argument("--dna-tumor-fastp", default=None)
    parser.add_argument("--dna-tumor-mosdepth", default=None)
    parser.add_argument("--dna-tumor-align", default=None)
    parser.add_argument("--dna-tumor-hs", default=None)
    #parser.add_argument("--dna-tumor-insert", default=None)
    
    parser.add_argument("--dna-normal-fastp", default=None)
    parser.add_argument("--dna-normal-mosdepth", default=None)
    parser.add_argument("--dna-normal-align", default=None)
    parser.add_argument("--dna-normal-hs", default=None)
    parser.add_argument("--dna-normal-insert", default=None)

    # RNA Inputs
    parser.add_argument("--rna-tumor-fastp", default=None)
    parser.add_argument("--rna-tumor-metrics", default=None)
    
    parser.add_argument("--rna-normal-fastp", default=None)
    parser.add_argument("--rna-normal-metrics", default=None)
    
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # === 1. 构建 DNA Sheet ===
    dna_rows = []
    
    def process_dna(sample_type, fastp_f, mosdepth_f, align_f, hs_f):
        # 即使文件全不存在，也生成一行包含 Sample 和 Type 的数据
        res = {'Sample': sample_type, 'Type': 'DNA'}
        res.update(parse_fastp_json(fastp_f))
        res.update(parse_mosdepth_summary(mosdepth_f))
        res.update(parse_picard_alignment(align_f))
        res.update(parse_picard_hs(hs_f))
        #ins = parse_picard_insert_metrics(insert_f)
        if 'insert_size_fastp' in res:
             res['insert_size'] = res['insert_size_fastp']
        return res

    dna_rows.append(process_dna('tumor', args.dna_tumor_fastp, args.dna_tumor_mosdepth, args.dna_tumor_align, args.dna_tumor_hs))
    dna_rows.append(process_dna('normal', args.dna_normal_fastp, args.dna_normal_mosdepth, args.dna_normal_align, args.dna_normal_hs))

    dna_df = pd.DataFrame(dna_rows)
    # 核心：使用 reindex 强制列对齐，缺失值填 "/"
    dna_df = dna_df.reindex(columns=DNA_COLS, fill_value='/')

    # === 2. 构建 RNA Sheet ===
    rna_rows = []
    
    def process_rna(sample_type, fastp_f, metrics_f):
        res = {'Sample': sample_type, 'Type': 'RNA'}
        res.update(parse_fastp_json(fastp_f))
        res.update(parse_picard_rna(metrics_f))
        if 'insert_size_fastp' in res:
             res['insert_size'] = res['insert_size_fastp']
        return res

    rna_rows.append(process_rna('tumor', args.rna_tumor_fastp, args.rna_tumor_metrics))
    rna_rows.append(process_rna('normal', args.rna_normal_fastp, args.rna_normal_metrics))

    rna_df = pd.DataFrame(rna_rows)
    # 核心：使用 reindex 强制列对齐，缺失值填 "/"
    rna_df = rna_df.reindex(columns=RNA_COLS, fill_value='/')

    # === 3. 写入 Excel ===
    # 确保输出目录存在
    out_dir = os.path.dirname(args.output)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    with pd.ExcelWriter(args.output, engine='openpyxl') as writer:
        dna_df.to_excel(writer, sheet_name='DNA', index=False)
        rna_df.to_excel(writer, sheet_name='RNA', index=False)
        # 即使为空也创建一个 ISO sheet
        pd.DataFrame().to_excel(writer, sheet_name='ISO', index=False)

if __name__ == "__main__":
    main()