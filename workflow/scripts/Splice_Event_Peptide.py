import pandas as pd
from Bio import SeqIO
import re
import os
import numpy as np
import collections
from typing import Tuple, List, Dict
from Bio.Seq import Seq
import argparse
from collections import Counter
import sys

# =============================================================================
# Helper Functions
# =============================================================================

def calculate_similarity(alignment_file):
    alignment = SeqIO.parse(alignment_file, "fasta")
    splice_seq = ""
    raw_seq = ""
    for record in alignment:
        if "Splice" in record.id:
            splice_seq = str(record.seq)
        elif "Database" in record.id:
            raw_seq = str(record.seq)
    if len(splice_seq) != len(raw_seq):
        return 1001
    identical_count = 0
    splice_length = 0
    for i in range(len(splice_seq)):
        splice_aa = splice_seq[i]
        raw_aa = raw_seq[i]
        if splice_aa == '-' and raw_aa == '-': continue
        if splice_aa != '-': splice_length += 1
        if splice_aa != '-' and raw_aa != '-' and splice_aa == raw_aa: identical_count += 1
    similarity = identical_count / splice_length if splice_length > 0 else 0
    return round(similarity,4)

def get_ungapped_maps(seq: str) -> Tuple[Dict[int, int], Dict[int, int], str]:
    ungap_to_align = {}
    align_to_ungap = {}
    ungapped_seq = []    
    ungapped_idx = 0
    for aligned_idx, aa in enumerate(seq):
        if aa != '-':
            ungap_to_align[ungapped_idx] = aligned_idx
            align_to_ungap[aligned_idx] = ungapped_idx
            ungapped_seq.append(aa)
            ungapped_idx += 1
    return ungap_to_align, align_to_ungap, ''.join(ungapped_seq)

def aligned_to_ungapped(aligned_pos: int, align_to_ungap: Dict[int, int]) -> int:
    return align_to_ungap.get(aligned_pos, -1)

def ungapped_to_aligned(ungapped_pos: int, ungap_to_align: Dict[int, int]) -> int:
    return ungap_to_align.get(ungapped_pos, -1)

def read_aln_file(file_path: str) -> Tuple[str, str, Dict[int, int], Dict[int, int], str, int]:
    records = list(SeqIO.parse(file_path, "fasta"))
    raw_aligned, splice_aligned = None, None
    for record in records:
        if 'Database' in record.id: raw_aligned = str(record.seq)
        elif 'Splice' in record.id: splice_aligned = str(record.seq)
    if not raw_aligned or not splice_aligned: raise ValueError(f"Missing seq in {file_path}")
    ungap_to_align, align_to_ungap, splice_ungapped = get_ungapped_maps(splice_aligned)
    return (raw_aligned, splice_aligned, ungap_to_align, align_to_ungap, splice_ungapped, len(splice_aligned))

def is_new_peptide_position(raw_aa: str, splice_aa: str) -> bool:
    return splice_aa != '-' and (raw_aa == '-' or splice_aa != raw_aa)

def calculate_new_ratio(raw_seq: str, splice_seq: str, start: int, end: int) -> float:
    if start >= end or start < 0 or end > len(raw_seq): return 0.0
    new_count = 0
    splice_valid_count = 0
    for i in range(start, end):
        raw_aa = raw_seq[i]
        splice_aa = splice_seq[i]
        if splice_aa != '-':
            splice_valid_count += 1
            if raw_aa == '-' or splice_aa != raw_aa: new_count += 1
    return new_count / splice_valid_count if splice_valid_count > 0 else 0.0

def merge_nearby_regions(regions, raw_aligned, splice_aligned, distance_threshold, new_ratio_threshold):
    if not regions: return []
    sorted_regions = sorted(regions, key=lambda x: x[3])
    merged = [list(sorted_regions[0])]
    for current in sorted_regions[1:]:
        curr_align_s, curr_align_e, curr_pep, curr_ungap_s, curr_ungap_e = current
        last_align_s, last_align_e, last_pep, last_ungap_s, last_ungap_e = merged[-1]
        distance = curr_ungap_s - last_ungap_e - 1
        if distance < distance_threshold:
            merged_align_s = last_align_s
            merged_align_e = curr_align_e
            new_ratio = calculate_new_ratio(raw_aligned, splice_aligned, merged_align_s-1, merged_align_e)
            if new_ratio >= new_ratio_threshold:
                merged_pep = splice_aligned[merged_align_s-1:merged_align_e].replace('-', '')
                merged[-1] = [merged_align_s, merged_align_e, merged_pep, last_ungap_s, curr_ungap_e]
                continue
        merged.append(list(current))
    return [(r[0], r[1], r[2]) for r in merged]

def expand_regions(regions, splice_ungapped, ungap_to_align, align_to_ungap, valid_length_threshold, max_total_expand=4, bidirectional_threshold=2):
    expanded_regions = []
    ungap_total = len(splice_ungapped)
    for start_align_1, end_align_1, peptide in regions:
        if len(peptide) >= valid_length_threshold:
            expanded_regions.append((start_align_1, end_align_1, peptide))
            continue
        start_align_0 = start_align_1 - 1
        end_align_0 = end_align_1 - 1
        ungap_positions = []
        for align_pos in range(start_align_0, end_align_0 + 1):
            ungap_pos = aligned_to_ungapped(align_pos, align_to_ungap)
            if ungap_pos != -1: ungap_positions.append(ungap_pos)
        if not ungap_positions:
            expanded_regions.append((start_align_1, end_align_1, peptide))
            continue
        start_ungap_0 = min(ungap_positions)
        end_ungap_0 = max(ungap_positions)
        upstream_space = start_ungap_0
        downstream_space = (ungap_total - 1) - end_ungap_0
        expand_up, expand_down = 0, 0
        if upstream_space >= bidirectional_threshold and downstream_space >= bidirectional_threshold:
            expand_up = bidirectional_threshold
            expand_down = bidirectional_threshold
        else:
            if upstream_space <= downstream_space:
                expand_up = min(upstream_space, max_total_expand)
                remaining = max_total_expand - expand_up
                expand_down = min(downstream_space, remaining)
            else:
                expand_down = min(downstream_space, max_total_expand)
                remaining = max_total_expand - expand_down
                expand_up = min(upstream_space, remaining)
        new_start_ungap = max(0, start_ungap_0 - expand_up)
        new_end_ungap = min(ungap_total - 1, end_ungap_0 + expand_down)
        new_end_ungap = max(new_start_ungap, new_end_ungap)
        expanded_peptide = splice_ungapped[new_start_ungap : new_end_ungap + 1]
        new_start_align_0 = ungapped_to_aligned(new_start_ungap, ungap_to_align)
        new_end_align_0 = ungapped_to_aligned(new_end_ungap, ungap_to_align)
        if new_start_align_0 == -1: new_start_align_0 = start_align_0
        if new_end_align_0 == -1: new_end_align_0 = end_align_0
        expanded_regions.append((new_start_align_0 + 1, new_end_align_0 + 1, expanded_peptide))
    return expanded_regions

def find_new_peptide_regions(raw_aligned, splice_aligned, ungap_to_align, align_to_ungap, splice_ungapped, total_length, distance_threshold, new_ratio_threshold, valid_length_threshold, max_total_expand=4, bidirectional_threshold=2):
    new_regions = []
    n = len(splice_aligned)
    in_new_region = False
    region_start_align_0 = 0
    for i in range(n):
        raw_aa = raw_aligned[i]
        splice_aa = splice_aligned[i]
        is_new = is_new_peptide_position(raw_aa, splice_aa)
        if is_new and not in_new_region:
            in_new_region = True
            region_start_align_0 = i
        elif not is_new and in_new_region:
            in_new_region = False
            region_end_align_0 = i - 1
            peptide = splice_aligned[region_start_align_0:region_end_align_0+1].replace('-', '')
            ungap_start = aligned_to_ungapped(region_start_align_0, align_to_ungap)
            ungap_end = aligned_to_ungapped(region_end_align_0, align_to_ungap)
            if ungap_start != -1 and ungap_end != -1 and peptide:
                new_regions.append((region_start_align_0+1, region_end_align_0+1, peptide, ungap_start, ungap_end))
    if in_new_region:
        region_end_align_0 = n - 1
        peptide = splice_aligned[region_start_align_0:region_end_align_0+1].replace('-', '')
        ungap_start = aligned_to_ungapped(region_start_align_0, align_to_ungap)
        ungap_end = aligned_to_ungapped(region_end_align_0, align_to_ungap)
        if ungap_start != -1 and ungap_end != -1 and peptide:
            new_regions.append((region_start_align_0+1, region_end_align_0+1, peptide, ungap_start, ungap_end))
    merged_regions = merge_nearby_regions(new_regions, raw_aligned, splice_aligned, distance_threshold, new_ratio_threshold)
    expanded_regions = expand_regions(merged_regions, splice_ungapped, ungap_to_align, align_to_ungap, valid_length_threshold, max_total_expand, bidirectional_threshold)
    return new_regions, merged_regions, expanded_regions

# =============================================================================
# Main Execution
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--Splice_File", help="Trans Output File (txt)", required=True)
    parser.add_argument("-o", "--outputdir", help="output dir", required=True)
    parser.add_argument("-n", "--sample_name", help="sample name", required=True)
    parser.add_argument("-g", "--gtf_bed_file_dir", help="gtf bed file dir", required=True)
    parser.add_argument("-b", "--bed_file_intersect", help="bed file intersect with Paper", required=True)
    args = parser.parse_args()

    Splice_File = args.Splice_File
    outputdir = args.outputdir
    sample_name = args.sample_name
    
    # 基础路径
    file_prefix = f"{outputdir}/{sample_name}_Alt_Splice_full_length"

    # 读取输入文件
    if not os.path.exists(Splice_File) or os.stat(Splice_File).st_size == 0:
        print("Empty or missing splice input.")
        # 创建空的输出文件以满足 Snakemake 要求
        pd.DataFrame().to_csv(f"{outputdir}/{sample_name}_Full_Length_Protein_Peptide_info.txt", sep="\t", index=False)
        open(f"{outputdir}/{sample_name}_Full_Length_Protein_Peptide_seq.faa", 'w').close()
        sys.exit(0)

    alt_splice_df = pd.read_csv(Splice_File, sep="\t", header=0)
    
    # 1. 初始化关键列 (防止 KeyError)
    # 这一步非常重要！确保后续逻辑总能找到这些列
    init_cols = [
        "Peptide_list_merge", "TransLateSoftware_Used", 
        "Peptide_Confidence_Highest", "Peptide_Confidence_List",
        "Peptide_Seq_Confidence4", "Peptide_Seq_Confidence3", "Peptide_Seq_Confidence2", "Peptide_Seq_Confidence1",
        "Peptide_Software_Confidence4", "Peptide_Software_Confidence3", "Peptide_Software_Confidence2", "Peptide_Software_Confidence1"
    ]
    for col in init_cols:
        alt_splice_df[col] = "-"
    alt_splice_df["Peptide_Confidence_Highest"] = 0

    # 2. 软件结果处理
    software = ["TRANSAID","gmst","TransDecoder","3Codon"]
    for sf in software:
        protein_file = f"{file_prefix}_Protein_{sf}.txt"
        
        # 初始化该软件特定的列
        alt_splice_df[f"{sf}_peptide_list_merge"] = "-"
        
        if not os.path.exists(protein_file): continue
            
        try:
            protein_data = pd.read_csv(protein_file, sep="\t", header=0)
            if "alt_splice_num" not in protein_data.columns: continue
            
            alt_data = protein_data[protein_data["alt_splice_num"]>=1]
            if alt_data.empty: continue
            
            alt_data = alt_data[["alt_splice_name"]].drop_duplicates()
            alt_splice_peptide_info = []
            
            inputdir = f"{outputdir}/{sample_name}/{sf}"
            
            for i in alt_data.index:
                alt_splice_name = alt_data.loc[i,"alt_splice_name"]
                try:
                    trans_id = re.split(r"_", alt_splice_name)[-2]
                except: continue
                    
                alt_splice_aln_file = f"{inputdir}/{alt_splice_name}.faa.aln"
                
                # 检查比对文件是否存在且不为空
                if os.path.exists(alt_splice_aln_file) and os.stat(alt_splice_aln_file).st_size > 0:
                    try:
                        similarity = calculate_similarity(alt_splice_aln_file)
                        (raw_aligned, splice_aligned, ungap_to_align, align_to_ungap, splice_ungapped, total_length) = read_aln_file(alt_splice_aln_file)
                        
                        new_peptide_region, _, new_peptide_region_merge_expand = find_new_peptide_regions(
                                raw_aligned, splice_aligned, ungap_to_align, align_to_ungap, splice_ungapped, total_length,
                                5, 0.5, 9, 4, 2) # 使用默认参数
                        
                        peptide_list_str = "-"
                        if len(new_peptide_region) > 0:
                            peptides = [f"{p[0]}_{p[1]}:{p[2]}" for p in new_peptide_region_merge_expand]
                            peptide_list_str = "|".join(peptides)
                        
                        alt_splice_peptide = pd.DataFrame({
                            "alt_splice_name": [alt_splice_name],
                            f"{sf}_trans_id": [trans_id],
                            f"{sf}_similarity": [similarity],
                            f"{sf}_peptide_list_merge": [peptide_list_str]
                        })
                        alt_splice_peptide_info.append(alt_splice_peptide)
                    except Exception as e:
                        # print(f"Error processing {alt_splice_name} in {sf}: {e}") # 可选：打印警告
                        pass
            
            if alt_splice_peptide_info:
                sf_df = pd.concat(alt_splice_peptide_info, axis=0)
                # Merge logic: if column exists, update it; otherwise add it
                # 这里使用 merge 后 update 的逻辑比较复杂，简单起见：
                # 先 merge 到临时 df，然后 update 原 df 的特定列
                temp_merged = pd.merge(alt_splice_df[["alt_splice_name"]], sf_df, on="alt_splice_name", how="left")
                
                # Update main dataframe columns for this software
                col_name = f"{sf}_peptide_list_merge"
                # 创建映射字典
                mapping = sf_df.set_index("alt_splice_name")[col_name].to_dict()
                alt_splice_df[col_name] = alt_splice_df["alt_splice_name"].map(mapping).fillna("-")

        except Exception as e:
            print(f"Error reading/processing protein file {protein_file}: {e}")
            continue

    # 3. 汇总统计 (Confidence Calculation)
    for i in alt_splice_df.index:
        all_peptides = []
        found_in_software = {} # peptide -> list of softwares
        
        for sf in software:
            col = f"{sf}_peptide_list_merge"
            if col in alt_splice_df.columns:
                val = str(alt_splice_df.loc[i, col])
                if val != "-":
                    # 提取肽段序列 (格式: start_end:SEQUENCE)
                    # 使用正则提取冒号后的序列部分
                    seqs = [x.split(':')[-1] for x in val.split('|') if ':' in x]
                    # 长度过滤
                    seqs = [s for s in seqs if len(s) >= 8]
                    
                    all_peptides.extend(seqs)
                    for s in seqs:
                        if s not in found_in_software: found_in_software[s] = []
                        found_in_software[s].append(sf)
        
        if all_peptides:
            counts = Counter(all_peptides)
            # 更新最高置信度
            max_conf = max(counts.values())
            alt_splice_df.loc[i, "Peptide_Confidence_Highest"] = max_conf
            
            # 更新 Confidence List (去重后的 count 列表字符串)
            conf_list = sorted(list(set(counts.values())), reverse=True)
            alt_splice_df.loc[i, "Peptide_Confidence_List"] = "|".join(map(str, conf_list))
            
            # 填充各级置信度详情
            # 这里的逻辑是：找出出现次数为 N 的所有肽段，并记录发现它们的软件
            for n in range(1, 5):
                peps_n = [p for p, c in counts.items() if c == n]
                if peps_n:
                    alt_splice_df.loc[i, f"Peptide_Seq_Confidence{n}"] = "|".join(peps_n)
                    
                    # 记录软件来源
                    softs_n = []
                    for p in peps_n:
                        softs_n.extend(found_in_software[p])
                    
                    softs_n = sorted(list(set(softs_n)))
                    if len(softs_n) == 4:
                        alt_splice_df.loc[i, f"Peptide_Software_Confidence{n}"] = "All_Software"
                    else:
                        alt_splice_df.loc[i, f"Peptide_Software_Confidence{n}"] = "|".join(softs_n)
            
            # 填充主 Peptide_list_merge (取所有软件结果的并集)
            # 或者根据某种策略（例如取置信度最高的肽段）
            # 这里简单取并集
            unique_peps = sorted(list(set(all_peptides)))
            alt_splice_df.loc[i, "Peptide_list_merge"] = "|".join(unique_peps)

    # 4. 排序与输出 Info 文件
    alt_splice_df = alt_splice_df.sort_values(by=["Peptide_Confidence_Highest", "Cancer_SampleRatio"], ascending=[False, False])
    outfile_info = f"{outputdir}/{sample_name}_Full_Length_Protein_Peptide_info.txt"
    alt_splice_df.to_csv(outfile_info, sep="\t", index=False)

    # 5. 输出 FASTA 文件
    outfile_seq = f"{outputdir}/{sample_name}_Full_Length_Protein_Peptide_seq.faa"
    
    # 筛选出有结果的行
    valid_df = alt_splice_df[alt_splice_df["Peptide_list_merge"] != "-"]
    
    with open(outfile_seq, "w") as fa:
        for i in valid_df.index:
            alt_name = valid_df.loc[i, "alt_splice_name"]
            peptides = valid_df.loc[i, "Peptide_list_merge"].split("|")
            
            for pep in peptides:
                if len(pep) < 8 or len(pep) > 300: continue
                # header: >AltSpliceName||PeptideSeq
                header = f">{alt_name}||{pep}"
                fa.write(f"{header}\n{pep}\n")