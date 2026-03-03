import pandas as pd
import pybedtools
import os
import numpy as np

# =============================================================================
# 工具函数
# =============================================================================
def parse_psi(psi_str):
    try:
        vals = [float(x) for x in str(psi_str).split(',') if x != 'NA']
        return np.mean(vals) if vals else 0.0
    except: return 0.0

def sum_counts(count_str):
    try:
        return sum([int(x) for x in str(count_str).split(',')])
    except: return 0

def get_tumor_target_junctions(row, as_type):
    """
    确定需要去 GTEx 查询的 Junctions。
    逻辑：找出肿瘤偏好的形态 (Tumor Form)。
    - 对于 SE/MXE/AltSS: 如果 GTEx 有 Tumor Form 的 Junction -> 过滤 (Bad)。
    - 对于 RI (IR): 如果 GTEx 有 Spliced Form 的 Junction -> 证明正常人是切除的 -> Tumor Retention 是特异的 (Good)。
    """
    juncs = [] # list of (chr, start, end, label)
    chrom = row['chr']
    dpsi = row['dPSI']
    
    # 判断肿瘤偏好 Inclusion (>0) 还是 Skipping/Splicing (<0)
    # 注意：单样本模式下 dPSI = TumorPSI，通常 > 0
    tumor_prefers_inclusion = dpsi > 0
    
    if as_type == 'SE':
        if tumor_prefers_inclusion:
            # Tumor 保留外显子 -> 查 Inclusion Juncs
            juncs.append((chrom, int(row['upstreamEE']), int(row['exonStart_0base']), 'J1_Inc'))
            juncs.append((chrom, int(row['exonEnd']), int(row['downstreamES']), 'J2_Inc'))
        else:
            # Tumor 跳过外显子 -> 查 Skipping Junc
            juncs.append((chrom, int(row['upstreamEE']), int(row['downstreamES']), 'J_Skip'))

    elif as_type == 'MXE':
        if tumor_prefers_inclusion: # Exon 1
            juncs.append((chrom, int(row['upstreamEE']), int(row['1stExonStart_0base']), 'J1_E1'))
            juncs.append((chrom, int(row['1stExonEnd']), int(row['downstreamES']), 'J2_E1'))
        else: # Exon 2
            juncs.append((chrom, int(row['upstreamEE']), int(row['2ndExonStart_0base']), 'J1_E2'))
            juncs.append((chrom, int(row['2ndExonEnd']), int(row['downstreamES']), 'J2_E2'))
            
    elif as_type in ['A3SS', 'A5SS']:
        if tumor_prefers_inclusion: # Long form
            if as_type == 'A3SS': juncs.append((chrom, int(row['flankingEE']), int(row['longExonStart_0base']), 'J_Long'))
            else: juncs.append((chrom, int(row['longExonEnd']), int(row['flankingES']), 'J_Long'))
        else: # Short form
            if as_type == 'A3SS': juncs.append((chrom, int(row['flankingEE']), int(row['shortExonStart_0base']), 'J_Short'))
            else: juncs.append((chrom, int(row['shortExonEnd']), int(row['flankingES']), 'J_Short'))

    elif as_type == 'RI':
        # RI 特殊处理
        if tumor_prefers_inclusion:
            # Tumor = Retention. 
            # 我们查 Spliced Junction。如果 GTEx 有 -> 说明正常是切除的 -> Tumor Retention 特异。
            # 标记为 Proxy，后续逻辑反转。
            juncs.append((chrom, int(row['upstreamEE']), int(row['downstreamES']), 'J_Spliced_Proxy'))
        else:
            # Tumor = Spliced. 通常不关注，但如果关注，逻辑同 SE Skipping。
            juncs.append((chrom, int(row['upstreamEE']), int(row['downstreamES']), 'J_Spliced'))

    return juncs

def run_filter():
    # 0. 准备
    params = snakemake.params
    output_dir = os.path.dirname(snakemake.output.tsv_strict)
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    pybedtools.helpers.set_tempdir(output_dir)

    # 1. 读取数据
    as_files = {'SE': snakemake.input.se, 'MXE': snakemake.input.mxe, 'RI': snakemake.input.ri, 
                'A3SS': snakemake.input.a3ss, 'A5SS': snakemake.input.a5ss}
    all_dfs = []
    for atype, path in as_files.items():
        if os.path.exists(path) and os.path.getsize(path) > 0:
            df = pd.read_csv(path, sep='\t')
            if not df.empty:
                df['AS_Type'] = atype
                all_dfs.append(df)
    
    if not all_dfs:
        for o in snakemake.output: open(o, 'w').close()
        return

    full_df = pd.concat(all_dfs, ignore_index=True)
    
    # 2. 统计量计算 (兼容单/双样本)
    has_normal = 'IncLevel2' in full_df.columns
    
    full_df['Tumor_PSI'] = full_df['IncLevel1'].apply(parse_psi)
    full_df['Tumor_Count'] = full_df['IJC_SAMPLE_1'].apply(sum_counts) + full_df['SJC_SAMPLE_1'].apply(sum_counts)
    
    if has_normal:
        full_df['Normal_PSI'] = full_df['IncLevel2'].apply(parse_psi)
        full_df['dPSI'] = full_df['IncLevelDifference']
    else:
        full_df['Normal_PSI'] = 0.0
        # 单样本模式下，假设 dPSI = Tumor_PSI (如果关注 Inclusion)
        # 简单起见，我们主要筛选 Tumor Inclusion 事件
        full_df['dPSI'] = full_df['Tumor_PSI'] 
        full_df['FDR'] = 0.0 

    full_df['UID'] = full_df['AS_Type'] + "_" + full_df['ID'].astype(str)

    # 3. 准备 GTEx 查询
    junc_bed_lines = []
    uid_to_target_juncs = {}
    
    for idx, row in full_df.iterrows():
        targets = get_tumor_target_junctions(row, row['AS_Type'])
        uid_to_target_juncs[row['UID']] = targets
        for (c, s, e, label) in targets:
            if e > s:
                junc_bed_lines.append(f"{c}\t{s}\t{e}\t{row['UID']}|{label}\t0\t+")

    # 4. 执行查询
    gtex_hits = {}
    if os.path.exists(snakemake.input.gtex_db) and junc_bed_lines:
        tmp_bed = os.path.join(output_dir, "query.bed")
        with open(tmp_bed, 'w') as f: f.write('\n'.join(junc_bed_lines) + '\n')
        
        a = pybedtools.BedTool(tmp_bed)
        b = pybedtools.BedTool(snakemake.input.gtex_db)
        
        # 精确匹配内含子坐标
        res = a.intersect(b, wa=True, wb=True, f=0.99, r=True)
        for feat in res:
            q_name = feat.name 
            try: count = int(feat.fields[-1]); gtex_hits[q_name] = max(gtex_hits.get(q_name, 0), count)
            except: pass
        if os.path.exists(tmp_bed): os.remove(tmp_bed)

    # 5. 综合判定
    results = []
    for idx, row in full_df.iterrows():
        reasons = []
        dpsi = row['dPSI']
        uid = row['UID']
        
        # 基础过滤
        # RI 使用 Retention count (IJC), 其他使用对应形态的 count
        if dpsi > 0: tumor_reads = sum_counts(row['IJC_SAMPLE_1'])
        else: tumor_reads = sum_counts(row['SJC_SAMPLE_1'])
            
        if tumor_reads < params.min_tumor_inc_reads: reasons.append("LowTumorSupport")
        if has_normal:
            if row.get('FDR', 1.0) > params.max_fdr: reasons.append("HighFDR")
            if abs(dpsi) < params.min_dpsi: reasons.append("LowDeltaPSI")

        # GTEx 背景判定
        target_juncs = uid_to_target_juncs.get(uid, [])
        
        # --- RI 特殊逻辑 ---
        if row['AS_Type'] == 'RI' and dpsi > 0:
            # 对于 RI Retention，我们查的是 Proxy (Spliced Junction)
            # 如果 GTEx 有 Splicing (Count > 10) -> 说明正常是切除的 -> Tumor Retention 是特异的
            proxy_count = 0
            for (_, _, _, label) in target_juncs:
                proxy_count = max(proxy_count, gtex_hits.get(f"{uid}|{label}", 0))
            
            if proxy_count > 10:
                tier = "Tier1_Specific_IR" # Good! GTEx supports splicing, Tumor has retention
            else:
                tier = "Tier3_Ambiguous_IR" # Bad. GTEx has no splicing evidence.
                reasons.append("NoSplicedEvidenceInGTEx")
        
        # --- 其他类型正常逻辑 (SE, MXE...) ---
        else:
            max_gtex = 0
            for (_, _, _, label) in target_juncs:
                max_gtex = max(max_gtex, gtex_hits.get(f"{uid}|{label}", 0))
            
            if max_gtex == 0: tier = "Tier1_HighlySpecific"
            elif max_gtex <= params.gtex_low_freq: tier = "Tier2_RareBackground"
            else: 
                tier = "Tier3_CommonBackground"
                reasons.append(f"InGTEx({max_gtex})")

        row['Filter_Status'] = "PASS" if not reasons else "FILTERED"
        row['Tier'] = tier
        results.append(row)

    # 6. 输出
    final_df = pd.DataFrame(results)
    pass_df = final_df[final_df['Filter_Status'] == "PASS"]
    pass_df.to_csv(snakemake.output.tsv_all, sep='\t', index=False)
    
    strict_df = pass_df[pass_df['Tier'].str.contains("Tier1")]
    strict_df.to_csv(snakemake.output.tsv_strict, sep='\t', index=False)
    
    # BED Output
    with open(snakemake.output.bed, 'w') as f:
        for _, r in strict_df.iterrows():
            # 简化：只输出第一个特征 Junction
            juncs = uid_to_target_juncs.get(r['UID'], [])
            if juncs:
                c, s, e, lbl = juncs[0]
                name = f"{r['geneSymbol']}_{r['UID']}_{lbl}"
                f.write(f"{c}\t{s}\t{e}\t{name}\t{abs(r['dPSI']):.2f}\t{r['strand']}\n")

if __name__ == "__main__":
    run_filter()