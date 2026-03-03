import pandas as pd
import pybedtools
import os
import sys
import numpy as np

# =============================================================================
# 配置
# =============================================================================
output_dir = os.path.dirname(snakemake.output.tsv_final)
if not os.path.exists(output_dir): os.makedirs(output_dir)
pybedtools.helpers.set_tempdir(output_dir)

# 新增：最小 Reads 支持数过滤，防止 Normal 只有 1 条 read 造成的假阳性
MIN_TOTAL_READS = 5  # (IJC + SJC) 必须大于这个数

# =============================================================================
# 辅助函数
# =============================================================================

def parse_inc_level(inc_level_str):
    try:
        vals = [float(x) for x in str(inc_level_str).split(',') if x != 'NA']
        return np.mean(vals) if vals else 0.0
    except:
        return 0.0

def safe_int(val):
    """将浮点或字符串转为整数，空值转为 -1"""
    try:
        if pd.isna(val) or val == "":
            return -1
        return int(float(val))
    except:
        return -1

def get_target_coords(row, as_type):
    """
    统一不同 AS 类型的核心变异区域坐标，用于显示和 IGV。
    返回 (Start, End)
    """
    s, e = -1, -1
    try:
        if as_type == 'SE':
            s, e = row['exonStart_0base'], row['exonEnd']
        elif as_type == 'MXE':
            # MXE 有两个外显子，取第一个做代表，或者取范围
            s, e = row['1stExonStart_0base'], row['1stExonEnd']
        elif as_type == 'RI':
            s, e = row['riExonStart_0base'], row['riExonEnd']
        elif as_type == 'A3SS':
            # A3SS 的核心差异区是 LongExonStart 到 ShortExonStart 之间
            s, e = row['longExonStart_0base'], row['shortES']
        elif as_type == 'A5SS':
            # A5SS 的核心差异区是 ShortExonEnd 到 LongExonEnd 之间
            s, e = row['shortEE'], row['longExonEnd']
    except:
        pass
    return safe_int(s), safe_int(e)

def get_junctions_for_query(row, as_type):
    """提取用于 GTEx 查询的 Junctions"""
    junctions = []
    chrom = row['chr']
    try:
        if as_type == 'SE':
            junctions.append((chrom, safe_int(row['upstreamEE']), safe_int(row['exonStart_0base']), 'J1'))
            junctions.append((chrom, safe_int(row['exonEnd']), safe_int(row['downstreamES']), 'J2'))
        elif as_type == 'MXE':
            junctions.append((chrom, safe_int(row['upstreamEE']), safe_int(row['1stExonStart_0base']), 'J1'))
            junctions.append((chrom, safe_int(row['1stExonEnd']), safe_int(row['downstreamES']), 'J2'))
        elif as_type == 'A3SS':
            junctions.append((chrom, safe_int(row['flankingEE']), safe_int(row['longExonStart_0base']), 'J_Long'))
            junctions.append((chrom, safe_int(row['flankingEE']), safe_int(row['shortES']), 'J_Short'))
        elif as_type == 'A5SS':
            junctions.append((chrom, safe_int(row['longExonEnd']), safe_int(row['flankingES']), 'J_Long'))
            junctions.append((chrom, safe_int(row['shortEE']), safe_int(row['flankingES']), 'J_Short'))
        elif as_type == 'RI':
            junctions.append((chrom, safe_int(row['upstreamEE']), safe_int(row['downstreamES']), 'J_Skip'))
    except:
        pass
    return junctions

def process_all_events(
    rmats_files, gtex_bed_file, 
    output_final, output_bed, output_report,
    min_dpsi, max_fdr, 
    gtex_low_thresh, gtex_high_thresh
):
    all_candidates = []
    
    # 1. 读取并标准化列
    for as_type, filepath in rmats_files.items():
        print(f"[*] Processing {as_type}: {filepath}")
        try:
            df = pd.read_csv(filepath, sep='\t')
            if df.empty: continue
        except: continue
            
        df['AS_Type'] = as_type
        
        # 统一坐标列 (解决空列和浮点数问题)
        # 我们创建统一的 Target_Start 和 Target_End 列
        coords = df.apply(lambda r: get_target_coords(r, as_type), axis=1)
        df['Target_Start'] = coords.apply(lambda x: x[0])
        df['Target_End'] = coords.apply(lambda x: x[1])
        
        # 格式化 IGV 坐标 string
        df['IGV_Coords'] = df.apply(
            lambda r: f"{r['chr']}:{r['Target_Start']}-{r['Target_End']}", axis=1
        )
        
        # 计算 PSI 和 Stats
        has_paired = 'IncLevelDifference' in df.columns
        df['TumorPSI'] = df['IncLevel1'].apply(parse_inc_level)
        df['NormalPSI'] = df['IncLevel2'].apply(parse_inc_level) if 'IncLevel2' in df.columns else 0.0
        
        # 增加总 Reads 数计算 (用于过滤低质量数据)
        # rMATS 的计数列是逗号分隔的字符串，如 "10,12"
        def sum_counts(s):
            try: return sum([int(x) for x in str(s).split(',')])
            except: return 0
            
        df['Tumor_Reads'] = df['IJC_SAMPLE_1'].apply(sum_counts) + df['SJC_SAMPLE_1'].apply(sum_counts)
        if has_paired:
            df['Normal_Reads'] = df['IJC_SAMPLE_2'].apply(sum_counts) + df['SJC_SAMPLE_2'].apply(sum_counts)
            df['dPSI'] = df['IncLevelDifference']
            df['FDR'] = df['FDR'].fillna(1.0)
        else:
            df['Normal_Reads'] = 0
            df['dPSI'] = df['TumorPSI']
            df['FDR'] = 0.0
            
        df['UID'] = df['AS_Type'] + "_" + df['GeneID'].astype(str) + "_" + df['ID'].astype(str)
        
        # 只保留核心列和通用列，防止 concat 时产生大量 NaN
        # (可选：保留所有列，但最后输出时清理)
        all_candidates.append(df)

    if not all_candidates:
        pd.DataFrame().to_csv(output_final, sep='\t'); open(output_bed, 'w').close(); return

    full_df = pd.concat(all_candidates, ignore_index=True)
    
    # 2. GTEx 查询准备 (使用物理文件)
    junc_list = []
    uid_to_juncs = {}
    tmp_query_file = os.path.join(output_dir, "query_candidates.tmp.bed")
    
    for idx, row in full_df.iterrows():
        juncs = get_junctions_for_query(row, row['AS_Type'])
        uid_to_juncs[row['UID']] = juncs
        for (chrom, s, e, label) in juncs:
            if e > s and s > 0:
                name = f"{row['UID']}|{label}"
                junc_list.append(f"{chrom}\t{s}\t{e}\t{name}\t0\t+")
    
    with open(tmp_query_file, 'w') as f:
        f.write('\n'.join(junc_list) + '\n')
    
    # 3. 执行 bedtools intersect
    gtex_hits = {}
    if os.path.exists(gtex_bed_file) and os.path.getsize(tmp_query_file) > 0:
        print("[*] Intersecting with GTEx DB...")
        try:
            # 加载 GTEx DB (假设已是大文件，bedtools 会处理)
            # 这里的 gtex_bed_file 是你的 sort 过的 bed
            # 使用 pybedtools 或直接 subprocess 调用
            # 注意：pybedtools 对于大文件内存消耗可能较大，但 intersect 是流式的
            
            a = pybedtools.BedTool(tmp_query_file)
            b = pybedtools.BedTool(gtex_bed_file)
            
            # -wa -wb -f 0.95 -r
            # 输出: query_chr, q_s, q_e, q_name, ..., db_chr, db_s, db_e, db_count
            # GTEx bed 是 4 列: chr, s, e, count
            # query bed 是 6 列
            # 结果列数 = 6 + 4 = 10
            # count 在第 10 列 (index 9)
            
            for feat in a.intersect(b, wa=True, wb=True, f=0.95, r=True):
                q_name = feat.name # query name
                try:
                    count = int(feat[9]) # 第10列
                except: count = 0
                gtex_hits[q_name] = max(gtex_hits.get(q_name, 0), count)
                
        except Exception as e:
            print(f"[!] Bedtools error: {e}")
    
    # 4. 判定与过滤
    results = []
    # 定义输出列顺序 (美化)
    out_cols = [
        'UID', 'geneSymbol', 'AS_Type', 'IGV_Coords', 'Tier', 'Filter_Status', 'Filter_Reason',
        'dPSI', 'FDR', 'TumorPSI', 'NormalPSI', 'Tumor_Reads', 'Normal_Reads', 'GTEx_Max_Count',
        'chr', 'strand', 'Target_Start', 'Target_End'
    ]
    
    for idx, row in full_df.iterrows():
        uid = row['UID']
        status = "PASS"
        reason = []
        tier = "Unknown"
        
        # A. 覆盖度过滤 (新增)
        # 如果是配对，要求 Tumor 和 Normal 总 Reads 数不能太少
        # 如果 Normal Reads 太少 (如 < 5)，计算出的 PSI 极不可靠，容易假阳性
        if has_paired:
            if row['Normal_Reads'] < MIN_TOTAL_READS:
                status = "FILTERED"
                reason.append(f"Low_Normal_Depth({row['Normal_Reads']})")
            if row['Tumor_Reads'] < MIN_TOTAL_READS:
                status = "FILTERED"
                reason.append(f"Low_Tumor_Depth({row['Tumor_Reads']})")
        else:
            if row['Tumor_Reads'] < MIN_TOTAL_READS:
                status = "FILTERED"
                reason.append(f"Low_Tumor_Depth({row['Tumor_Reads']})")

        # B. 统计学过滤
        if abs(row['dPSI']) < min_dpsi:
            status = "FILTERED"
            reason.append(f"Low_dPSI({row['dPSI']:.2f})")
        if has_paired and row['FDR'] > max_fdr:
            status = "FILTERED"
            reason.append(f"High_FDR")
            
        # C. 背景分级
        max_gtex_count = 0
        juncs = uid_to_juncs.get(uid, [])
        for (_, _, _, label) in juncs:
            key = f"{uid}|{label}"
            count = gtex_hits.get(key, 0)
            max_gtex_count = max(max_gtex_count, count)
            
        if max_gtex_count == 0:
            tier = "Tier1_HighlySpecific"
        elif max_gtex_count <= gtex_low_thresh:
            tier = "Tier2_RareBackground"
        elif max_gtex_count <= gtex_high_thresh:
            tier = "Tier3_CommonBackground"
            if status == "PASS": status = "FILTERED"; reason.append("GTEx_Common")
        else:
            tier = "Tier4_Ubiquitous"
            if status == "PASS": status = "FILTERED"; reason.append("GTEx_Ubiquitous")
            
        # 更新 row
        row['Filter_Status'] = status
        row['Filter_Reason'] = ";".join(reason) if reason else "Pass"
        row['GTEx_Max_Count'] = max_gtex_count
        row['Tier'] = tier
        results.append(row)
        
    final_df_all = pd.DataFrame(results)
    
    # 整理列顺序
    # 提取存在的列
    final_cols = [c for c in out_cols if c in final_df_all.columns]
    # 把其他 rMATS 原生列放到后面
    other_cols = [c for c in final_df_all.columns if c not in final_cols]
    final_df_all = final_df_all[final_cols + other_cols]

    # 5. 写出
    print(f"[*] Writing Report...")
    final_df_all.to_csv(output_report, sep='\t', index=False)
    
    pass_df = final_df_all[final_df_all['Filter_Status'] == "PASS"]
    pass_df.to_csv(output_final, sep='\t', index=False)
    
    # BED (用 Target Start/End)
    with open(output_bed, 'w') as f:
        for idx, row in pass_df.iterrows():
            name = f"{row['geneSymbol']}_{row['UID']}_{row['Tier']}"
            score = row['dPSI']
            try:
                s, e = int(row['Target_Start']), int(row['Target_End'])
                if e > s:
                    f.write(f"{row['chr']}\t{s}\t{e}\t{name}\t{score:.3f}\t{row['strand']}\n")
            except: pass

    if os.path.exists(tmp_query_file): os.remove(tmp_query_file)
    pybedtools.cleanup(remove_all=True)

if __name__ == "__main__":
    files = {
        'SE': snakemake.input.se,
        'MXE': snakemake.input.mxe,
        'RI': snakemake.input.ri,
        'A3SS': snakemake.input.a3ss,
        'A5SS': snakemake.input.a5ss
    }
    
    process_all_events(
        rmats_files     = files,
        gtex_bed_file   = snakemake.input.gtex_db,
        output_final    = snakemake.output.tsv_final,
        output_bed      = snakemake.output.bed_final,
        output_report   = snakemake.output.report_full,
        min_dpsi        = float(snakemake.params.min_dpsi),
        max_fdr         = float(snakemake.params.max_fdr),
        gtex_low_thresh = int(snakemake.params.gtex_low_freq),
        gtex_high_thresh= int(snakemake.params.gtex_high_freq)
    )