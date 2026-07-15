#!/usr/bin/env python3
import pandas as pd
from pyfaidx import Fasta
import sys
import os
import argparse

# ================= 1. 基础生物学工具函数 =================
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq.upper()))

def translate(seq):
    seq = seq.upper()
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*': break
        protein += aa
    return protein

def get_first_mismatch(wt_prot, mut_prot):
    min_len = min(len(wt_prot), len(mut_prot))
    for i in range(min_len):
        if wt_prot[i] != mut_prot[i]:
            return i
    if len(mut_prot) > len(wt_prot):
        return min_len
    return -1

# ================= 2. 核心处理类 =================
class ASProcessor:
    def __init__(self, excel_path, gtf_path, fasta_path):
        self.excel_path = excel_path
        self.gtf_path = gtf_path
        self.fasta_path = fasta_path
        self.df = pd.DataFrame()
        self.target_transcripts = set()
        self.cds_dict = {}
        
    def load_and_filter_data(self):
        """读取指定 Sheet，并只保留 Manual_Select == yes 的行"""
        print(f"[*] 尝试加载 Excel: {self.excel_path}")
        try:
            raw_df = pd.read_excel(self.excel_path, sheet_name="Large scale variants")
        except ValueError:
            print("[警告] 未找到名为 'Large scale variants' 的 Sheet。将输出空文件。")
            return False
        
        if raw_df.empty:
            print("[警告] 'Large scale variants' Sheet 为空。将输出空文件。")
            return False
            
        if 'Manual_Select' not in raw_df.columns:
            print("[错误] 未找到 'Manual_Select' 列！请检查输入文件格式。")
            return False
            
        # 过滤: 将 Manual_Select 转换为小写并去除空格，筛选值为 'yes'
        self.df = raw_df[raw_df['Manual_Select'].astype(str).str.lower().str.strip() == 'yes'].copy()
        
        if self.df.empty:
            print("[INFO] 'Large scale variants' 中没有标记为 'yes' 的变异。将输出空文件。")
            return False
            
        self.target_transcripts = set(self.df['Transcript'].dropna().unique())
        print(f"[*] 过滤完毕，共有 {len(self.df)} 个变异需要处理，涉及 {len(self.target_transcripts)} 个独特转录本。")
        return True

    def _parse_gtf(self):
        print("[*] 正在解析 GTF 文件...")
        with open(self.gtf_path, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'CDS': continue
                
                info = parts[8]
                if 'transcript_id' not in info: continue
                # 获取带有版本号的 ID，例如 "ENST00000394420.9"
                full_enst = info.split('transcript_id "')[1].split('"')[0]
                # 关键修改：提取版本号前的部分，即 "ENST00000394420"
                base_enst = full_enst.split('.')[0]
                
                # 现在用 base_enst 与 Excel 中的 ID 进行匹配
                if base_enst in self.target_transcripts:
                    if base_enst not in self.cds_dict: self.cds_dict[base_enst] = []
                    self.cds_dict[base_enst].append({
                        'chrom': parts[0], 'strand': parts[6],
                        'start': int(parts[3]), 'end': int(parts[4])
                    })
        
        # 同时修改后面排序的字典 key
        for enst, blocks in self.cds_dict.items():
            strand = blocks[0]['strand']
            if strand == '+': blocks.sort(key=lambda x: x['start'])
            else: blocks.sort(key=lambda x: x['end'], reverse=True)
        print("[*] GTF 解析完成！")
    
    def _get_sequence_from_blocks(self, blocks):
        seq = ""
        for b in blocks:
            segment = str(self.fasta[b['chrom']][b['start']-1 : b['end']])
            if b['strand'] == '-': segment = reverse_complement(segment)
            seq += segment
        return seq

    def process_and_export(self, out_wt, out_mut, out_prot, out_peptide):
        print("[*] 正在加载参考基因组 (pyfaidx)...")
        self.fasta = Fasta(self.fasta_path)
        self._parse_gtf()
        
        with open(out_wt, "w") as f_wt_cds, \
             open(out_mut, "w") as f_mut_cds, \
             open(out_prot, "w") as f_mut_prot, \
             open(out_peptide, "w") as f_peptide:

            for index, row in self.df.iterrows():
                enst = row['Transcript']
                gene = row['Gene']
                event_type = row['EventAnnotation']
                trunc_or_ext = row['Truncation_or_Extension']
                if pd.isna(event_type) or pd.isna(trunc_or_ext): continue
                
                length = int(row['Length'])
                pos = int(row['POS'])
                
                if enst not in self.cds_dict:
                    print(f"[警告] 转录本 {enst} 未在 GTF 的 CDS 中找到，跳过。")
                    continue
                    
                wt_blocks = self.cds_dict[enst]
                strand = wt_blocks[0]['strand']
                
                # --- 1. 定位受影响的 Block ---
                target_block_idx = -1
                min_dist = float('inf')
                for i, b in enumerate(wt_blocks):
                    if event_type == "alt-3'":
                        dist = abs(b['start'] - pos) if strand == '+' else abs(b['end'] - pos)
                    elif event_type == "alt-5'":
                        dist = abs(b['end'] - pos) if strand == '+' else abs(b['start'] - pos)
                    else:
                        dist = min(abs(b['start'] - pos), abs(b['end'] - pos))
                    if dist < min_dist:
                        min_dist = dist
                        target_block_idx = i
                        
                if target_block_idx == -1: continue
                
                # --- 2. 修改 Mutant Blocks 坐标 ---
                mut_blocks = [dict(b) for b in wt_blocks]
                target_b = mut_blocks[target_block_idx]
                
                if event_type == "alt-3'":
                    if strand == '+': 
                        target_b['start'] = target_b['start'] - length if trunc_or_ext == 'Extension' else target_b['start'] + length
                    else:             
                        target_b['end'] = target_b['end'] + length if trunc_or_ext == 'Extension' else target_b['end'] - length
                elif event_type == "alt-5'":
                    if strand == '+': 
                        target_b['end'] = target_b['end'] + length if trunc_or_ext == 'Extension' else target_b['end'] - length
                    else:             
                        target_b['start'] = target_b['start'] - length if trunc_or_ext == 'Extension' else target_b['start'] + length

                # --- 3. 获取序列并翻译 ---
                wt_dna = self._get_sequence_from_blocks(wt_blocks)
                mut_dna = self._get_sequence_from_blocks(mut_blocks)
                wt_prot = translate(wt_dna)
                mut_prot = translate(mut_dna)
                
                idx = get_first_mismatch(wt_prot, mut_prot)
                if idx == -1: continue 

                # --- 4. 核心截取逻辑 (带 'A' 占位) ---
                var_peptide = ""
                ref_peptide = ""
                
                # 统一将保留长度从 8/12 调整为 12
                flank_len = 12
                
                if length % 3 != 0:
                    # ====== 发生移码 (Frameshift) ======
                    # Var: 左侧 12 个 WT aa + 新生的移码序列
                    var_peptide = mut_prot[max(0, idx - flank_len) : ]
                    # 移码产生的全新序列长度
                    mut_tail_len = len(mut_prot) - idx
                    # Ref: 对应左侧 12 个 WT aa + 用 A 占位变异部分
                    # 注意：分歧点之前的序列必须完全一致，所以前半部分取 wt_prot
                    ref_peptide = wt_prot[max(0, idx - flank_len) : idx] + 'A' * mut_tail_len
                    
                else:
                    # ====== 未发生移码 (In-frame) ======
                    aa_len = length // 3
                    if trunc_or_ext == 'Extension':
                        # Var: 左侧 12 WT aa + 插入的新 AA + 右侧 12 WT aa
                        var_peptide = wt_prot[max(0, idx - flank_len) : idx] + mut_prot[idx : idx + aa_len] + wt_prot[idx : idx + flank_len]
                        # Ref: 左侧 12 WT aa + 缺失 AA 数量的 A 占位 + 右侧 12 WT aa
                        ref_peptide = wt_prot[max(0, idx - flank_len) : idx] + 'A' * aa_len + wt_prot[idx : idx + flank_len]
                    else: # Truncation
                        # Var: 左侧 12 WT aa + 跳过缺失部分后右侧 12 WT aa
                        var_peptide = wt_prot[max(0, idx - flank_len) : idx] + wt_prot[idx + aa_len : idx + aa_len + flank_len]
                        # Ref: 左侧 12 WT aa + 缺失 AA 数量的 A 占位 + 右侧 12 WT aa
                        ref_peptide = wt_prot[max(0, idx - flank_len) : idx] + 'A' * aa_len + wt_prot[idx + aa_len : idx + aa_len + flank_len]
                        
                # --- 5. 写入 ---
                f_wt_cds.write(f">{enst}|WT_CDS|{gene}\n{wt_dna}\n")
                f_mut_cds.write(f">{enst}|MUT_CDS|{gene}|{event_type}|{trunc_or_ext}|Len_{length}\n{mut_dna}\n")
                f_mut_prot.write(f">{enst}|MUT_PROT|{gene}|{event_type}|{trunc_or_ext}\n{mut_prot}\n")
                f_peptide.write(f">Ref|{enst}\n{ref_peptide}\n")
                f_peptide.write(f">Var|{enst}|{gene}|{event_type}|-\n{var_peptide}\n")

        print("[*] 成功写入所有 4 个文件。")


def create_empty_outputs(outputs):
    """创建空文件防止报错"""
    for file_path in outputs:
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, 'w') as f:
            pass 

# ================= 3. 灵活执行接口 (Snakemake & Standalone) =================
def main():
    # 判定是否在 Snakemake 环境中运行
    if "snakemake" in globals():
        print("[INFO] 检测到 Snakemake 环境，使用自动注入的参数运行。")
        
        # 捕捉并重定向 Snakemake 日志
        if hasattr(snakemake, 'log') and snakemake.log:
            sys.stderr = open(snakemake.log[0], "w")
            sys.stdout = sys.stderr

        excel_file = snakemake.input.manual_check_file
        gtf_file = snakemake.params.gtf
        fasta_file = snakemake.params.fasta
        
        out_wt = snakemake.output.wt_cds
        out_mut = snakemake.output.mut_cds
        out_prot = snakemake.output.mut_prot
        out_peptide = snakemake.output.peptide

    else:
        print("[INFO] 独立命令行模式运行。")
        parser = argparse.ArgumentParser(description="Generate Neoantigen Peptides from AS events")
        parser.add_argument("-i", "--input_excel", required=True, help="Input manually checked Excel report")
        parser.add_argument("-g", "--gtf", required=True, help="Reference GTF annotation file")
        parser.add_argument("-f", "--fasta", required=True, help="Reference Genome FASTA file")
        parser.add_argument("--out_wt_cds", required=True, help="Output path for wildtype CDS fasta")
        parser.add_argument("--out_mut_cds", required=True, help="Output path for mutant CDS fasta")
        parser.add_argument("--out_mut_prot", required=True, help="Output path for mutant protein faa")
        parser.add_argument("--out_peptide", required=True, help="Output path for neo-peptide faa")
        
        args = parser.parse_args()
        
        excel_file = args.input_excel
        gtf_file = args.gtf
        fasta_file = args.fasta
        out_wt = args.out_wt_cds
        out_mut = args.out_mut_cds
        out_prot = args.out_mut_prot
        out_peptide = args.out_peptide

    # ================= 公共执行逻辑 =================
    # 确保输出目录存在
    os.makedirs(os.path.dirname(out_wt), exist_ok=True)
    os.makedirs(os.path.dirname(out_mut), exist_ok=True)
    os.makedirs(os.path.dirname(out_prot), exist_ok=True)
    os.makedirs(os.path.dirname(out_peptide), exist_ok=True)

    processor = ASProcessor(excel_file, gtf_file, fasta_file)
    
    if processor.load_and_filter_data():
        processor.process_and_export(out_wt, out_mut, out_prot, out_peptide)
    else:
        print("[INFO] 未找到需处理的有效变异，生成空的占位文件结束任务。")
        create_empty_outputs([out_wt, out_mut, out_prot, out_peptide])

if __name__ == "__main__":
    main()