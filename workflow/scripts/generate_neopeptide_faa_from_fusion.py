#!/usr/bin/env python3
import csv
import argparse
import sys
import re

# ================= 1. 基础辅助功能 =================

def clean_sequence(seq):
    """移除终止密码子星号"""
    return seq.replace('*', '')

def clean_gene_name(gene_raw):
    """
    Arriba 有时输出带有位置或多个基因的格式，如 'CDHR3(31874),SYPL1(22199)'。
    我们只取第一个主基因，并去除括号内容。
    """
    gene = gene_raw.split(',')[0]
    gene = re.sub(r'\(.*?\)', '', gene)
    return gene

def score_gene_name(gene_name):
    """
    为基因名打分，以在去冗余时优先选择经典的蛋白质编码基因。
    """
    if not gene_name or gene_name == '.':
        return 0
    name_upper = gene_name.upper()
    # 假基因、克隆、长链非编码 RNA 给予低分
    if '.' in name_upper or name_upper.startswith('LINC') or name_upper.startswith('AC0') or name_upper.startswith('AL'):
        return 1
    # 经典的标准的 Gene Symbol 给予高分
    return 10

def generate_ref_and_var(raw_seq_part):
    """
    核心升级：字符级精准映射。
    大写代表野生型 (保留)，小写代表变异/插入/移码 (Var转大写，Ref转A)。
    """
    ref_seq = ""
    var_seq = ""
    for char in raw_seq_part:
        if char.islower():
            ref_seq += "A"
            var_seq += char.upper()
        else:
            ref_seq += char
            var_seq += char
    return ref_seq, var_seq

# ================= 2. 核心处理与去重逻辑 =================

def process_fusion_peptides(input_tsv, output_faa, confidence_levels, max_aa_len, window_size):
    target_conf = [c.lower() for c in confidence_levels]
    
    # 使用字典进行去重：key = var_seq (最终突变肽段), value = fusion_info_dict
    unique_peptides = {}
    
    try:
        with open(input_tsv, 'r', encoding='utf-8') as tsv_file:
            reader = csv.DictReader(tsv_file, delimiter='\t')

            # --- 第一阶段：读取、解析、去重 ---
            for row in reader:
                confidence = row.get('confidence', '').lower()
                if confidence not in target_conf: continue
                
                reading_frame = row.get('reading_frame', '').lower()
                if reading_frame == 'stop-codon' or reading_frame == '.': continue
                
                peptide_raw = row.get('peptide_sequence', '')
                if not peptide_raw or '|' not in peptide_raw: continue
                
                gene1 = clean_gene_name(row.get('#gene1', 'Unknown'))
                gene2 = clean_gene_name(row.get('gene2', 'Unknown'))
                
                # 分割左右序列并清除星号 (保留大小写特征)
                parts = peptide_raw.split('|')
                left_full = parts[0]
                right_full = clean_sequence(parts[1])
                
                # 截取窗口：左侧取最后 window_size 个字符
                left_flank = left_full[-window_size:] if len(left_full) >= window_size else left_full
                
                # 右侧截取逻辑：严格基于 Arriba 的 reading_frame
                if reading_frame == 'in-frame':
                    # In-frame: 右侧取 window_size 个字符 (通常包含少许小写连接符 + 大写野生型Gene2)
                    right_flank = right_full[:window_size]
                else:
                    # Out-of-frame: 右侧全取 (因为移码后全是小写的全新序列)
                    right_flank = right_full
                
                # 精准生成 Ref 和 Var
                ref_left, var_left = generate_ref_and_var(left_flank)
                ref_right, var_right = generate_ref_and_var(right_flank)
                
                final_ref_seq = ref_left + ref_right
                final_var_seq = var_left + var_right
                
                # 计算打分用于优选
                conf_score = 2 if confidence == 'high' else 1
                filter_score = 1 if row.get('filters', '.') == '.' else 0
                try:
                    read_evidence = int(row.get('split_reads1', 0)) + int(row.get('split_reads2', 0)) + int(row.get('discordant_mates', 0))
                except ValueError:
                    read_evidence = 0
                
                # 评估基因组合的“经典”程度
                gene_quality_score = score_gene_name(gene1) + score_gene_name(gene2)
                
                fusion_record = {
                    'gene1': gene1,
                    'gene2': gene2,
                    'ref_seq': final_ref_seq,
                    'var_seq': final_var_seq,
                    'aa_length': len(final_var_seq),
                    'conf_score': conf_score,
                    'filter_score': filter_score,
                    'read_evidence': read_evidence,
                    'gene_quality': gene_quality_score
                }
                
                # 去冗余逻辑：如果该肽段已存在，对比并保留质量更好的记录
                if final_var_seq in unique_peptides:
                    existing = unique_peptides[final_var_seq]
                    # 优先级：基因经典度 > Filter无警告 > 表达量
                    if (gene_quality_score, filter_score, read_evidence) > \
                       (existing['gene_quality'], existing['filter_score'], existing['read_evidence']):
                        unique_peptides[final_var_seq] = fusion_record
                    else:
                        # 基因名不变的情况下，如果有相同的肽段，我们把 Read evidence 取最大值以防低估
                        existing['read_evidence'] = max(existing['read_evidence'], read_evidence)
                else:
                    unique_peptides[final_var_seq] = fusion_record

        # --- 第二阶段：提取字典的 value 为列表，进行贪心排序 ---
        candidate_fusions = list(unique_peptides.values())
        
        # 按照 [置信度, 无质控警告, 支持Read数量] 进行降序排列
        candidate_fusions.sort(
            key=lambda x: (x['conf_score'], x['filter_score'], x['read_evidence']), 
            reverse=True
        )
        
        # --- 第三阶段：动态截断输出 ---
        entry_count = 0
        current_total_aa = 0
        
        with open(output_faa, 'w', encoding='utf-8') as faa_file:
            for fusion in candidate_fusions:
                # 检查加上当前序列后，总长度是否会超过最大限制
                if current_total_aa + fusion['aa_length'] > max_aa_len:
                    print(f"INFO: Reached max AA limit ({max_aa_len}). Capping output.")
                    break 
                
                current_total_aa += fusion['aa_length']
                entry_count += 1
                fusion_id = f"fusion{entry_count:03d}"
                
                faa_file.write(f">Ref|-\n{fusion['ref_seq']}\n")
                faa_file.write(f">Var|-|{fusion['gene1']}-{fusion['gene2']}|{fusion_id}|-\n{fusion['var_seq']}\n")
                
        print(f"Success: {entry_count} unique, high-quality fusions processed into '{output_faa}'.")
        print(f"Total novel amino acids exported: {current_total_aa} / {max_aa_len}")

    except FileNotFoundError:
        print(f"Error: Input file '{input_tsv}' not found.", file=sys.stderr)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)

# ================= 3. 命令行接口 =================

def main():
    parser = argparse.ArgumentParser(
        description="Convert Arriba fusion results to neo-peptide FASTA with deduplication, smart A-padding, and dynamic capping."
    )
    parser.add_argument("-i", "--input", required=True, help="Input Arriba fusion TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output FAA file path")
    parser.add_argument("-c", "--confidence", nargs='+', default=["high", "medium"],
                        help="Confidence levels to extract. Default: high medium")
    parser.add_argument("-m", "--max_aa_len", type=int, default=500,
                        help="Maximum total length of variant amino acids to output. Default: 500")
    parser.add_argument("-w", "--window_size", type=int, default=18,
                        help="Number of amino acids to extract on each side of the fusion point. Default: 18")

    args = parser.parse_args()
    process_fusion_peptides(args.input, args.output, args.confidence, args.max_aa_len, args.window_size)

if __name__ == "__main__":
    main()