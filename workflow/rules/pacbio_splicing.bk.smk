# workflow/rules/rna_splicing.smk
import os
import pandas as pd

# =============================================================================
# 辅助函数：动态判断样本配对情况 & 获取输入
# =============================================================================

def has_paired_normal(sample_id):
    """
    查询 samples 表，检查该病人是否有 normal 的 RNA 数据。
    返回 True 表示有配对对照，False 表示只有 Tumor 数据。
    """
    # 假设 samples 全局变量已由 Snakemake 加载 configfile 后生成
    # 必须确保 samples.csv 中有 sample_id, sampletype, data_type 列
    try:
        subset = samples[(samples['sample_id'] == sample_id) & 
                         (samples['sampletype'] == 'normal') & 
                         (samples['data_type'] == 'rna')]
        return not subset.empty
    except Exception as e:
        # 如果查找失败，默认假设为无对照模式，避免报错中断
        print(f"[WARNING] Could not determine paired status for {sample_id}: {e}")
        return False

def get_sqanti_input_gtf(wildcards):
    """
    自适应输入 GTF：
    - 有配对：使用 StringTie2 Merge 后的 GTF (包含 T+N+Isoseq)
    - 无配对：使用 Tumor 单独的 Hybrid GTF
    """
    if has_paired_normal(wildcards.sample):
        # 对应 rna_assembly.smk 中 stringtie_merge 的输出
        return f"rna/assembly/{wildcards.sample}_merged.assembly.gtf.gz"
    else:
        # 对应 rna_assembly.smk 中 stringtie_assemble 的输出
        return f"rna/assembly/{wildcards.sample}_tumor.assembly.gtf.gz"

def get_rmats_inputs(wildcards):
    """
    获取 rMATS 所需的 BAM 文件：
    总是返回 Tumor BAM，如果有 Normal 则一并返回。
    """
    inputs = {
        "gtf": f"rna/splicing/{wildcards.sample}_sqanti3/{wildcards.sample}_corrected.gtf",
        "bam_tumor": f"rna/align/{wildcards.sample}_tumor/Aligned.sortedByCoord.out.bam"
    }
    if has_paired_normal(wildcards.sample):
        inputs["bam_normal"] = f"rna/align/{wildcards.sample}_normal/Aligned.sortedByCoord.out.bam"
    return inputs


# =============================================================================
# Rule 0: 准备 GTEx v8 (hg38) 数据库 (本地处理版)
# =============================================================================
rule splice_00_process_local_gtex_v8:
    """
    直接处理本地下载好的 GTEx v8 junctions.bgz 文件。
    提取 chr, start, end, count 用于后续过滤。
    """
    input:
        bgz = "/home/bioinfo/05.pipeline_dev/SakitNeo/resources/recount3_db/junctions.bgz"
    output:
        gtex_bed = config["splicing_params"]["filtering_db"]["gtex_junctions"]
    log: "logs/splicing/process_gtex_v8.log"
    conda: "../../envs/python_filter.yaml"
    shell:
        """
        echo "[*] Processing local GTEx v8 junctions..." > {log}
        
        # Snaptron/Recount3 GTEx v2 结构:
        # Col 2: chr
        # Col 3: start
        # Col 4: end
        # Col 13: sample_count
        
        # 逻辑：
        # 1. 过滤样本数 >= 10 (有效过滤背景)
        # 2. 规范化 chr 名称 (如果缺失 chr 前缀则补全)
        
        gzip -dc {input.bgz} | \
        awk -F "\\t" 'NR>1 && $13 >= 10 {{
            chrom = $2;
            if (index(chrom, "chr") != 1) {{ chrom = "chr"chrom }};
            print chrom"\\t"$3"\\t"$4"\\t"$13
        }}' | \
        sort -k1,1 -k2,2n > {output.gtex_bed} 2>> {log}
        
        # 检查输出文件是否成功生成且不为空
        if [ ! -s {output.gtex_bed} ]; then
            echo "[ERROR] Output file is empty! Check awk column selection." >> {log}
            exit 1
        fi
        
        echo "[*] Success. Generated GTEx DB lines:" >> {log}
        wc -l {output.gtex_bed} >> {log}
        """

# =============================================================================
# Rule 1: SQANTI3 质控 (极速版：跳过 ORF)
# =============================================================================
rule splice_01_sqanti3_qc:
    """
    对 StringTie2 Hybrid 组装的 GTF 进行结构纠正和分类。
    关键优化：启用 --skipORF，避免 TransDecoder 单线程运行数小时的瓶颈。
    """
    input:
        gtf = get_sqanti_input_gtf,
        ref_gtf = config["reference"]["gtf"],
        ref_fasta = config["reference"]["genome"]
    output:
        gtf = "rna/splicing/{sample}_sqanti3/{sample}_corrected.gtf",
        classification = "rna/splicing/{sample}_sqanti3/{sample}_classification.txt",
        # 即使跳过 ORF，SQANTI3 仍会生成修正后的 FASTA
        fasta = "rna/splicing/{sample}_sqanti3/{sample}_corrected.fasta"
    params:
        out_dir = "rna/splicing/{sample}_sqanti3",
        out_prefix = "{sample}",
        sqanti_script = os.path.join(config["tools_path"]["sqanti3_dir"], "sqanti3_qc.py")
    log: "logs/splicing/sqanti3_{sample}.log"
    threads: 24
    conda: "../../envs/sqanti3.yaml"
    shell:
        """
        # 1. 环境准备
        SQ_DIR="{config[tools_path][sqanti3_dir]}"
        export PYTHONPATH="$SQ_DIR:${{PYTHONPATH:-}}"
        
        # [Fix] 解决 libcurl.so.4 not found 问题
        # 自动获取当前 Conda 环境的 lib 目录并加入 LD_LIBRARY_PATH
        if [ -n "$CONDA_PREFIX" ]; then
            export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        fi
        
        # 2. 处理输入 GTF
        INPUT_GTF="{input.gtf}"
        WORK_GTF="{params.out_dir}/input_temp.gtf"
        
        mkdir -p {params.out_dir}
        
        if [[ "$INPUT_GTF" == *.gz ]]; then
            echo "[*] Decompressing input GTF..." > {log}
            gunzip -c "$INPUT_GTF" > $WORK_GTF
        else
            cp "$INPUT_GTF" $WORK_GTF
        fi

        # 3. 清理旧文件
        rm -rf {params.out_dir}/*corrected* {params.out_dir}/TD2*

        # 4. 执行 SQANTI3
        echo "[*] Running SQANTI3..." >> {log}
        
        python3 {params.sqanti_script} \
            --isoforms $WORK_GTF \
            --refGTF {input.ref_gtf} \
            --refFasta {input.ref_fasta} \
            -d {params.out_dir} \
            -o {params.out_prefix} \
            -t {threads} \
            --skipORF \
            --report skip >> {log} 2>&1
            
        # 5. 清理临时解压文件
        rm -f $WORK_GTF
        """

# =============================================================================
# Rule 2: rMATS 差异分析 (自适应：单样本 vs 配对)
# =============================================================================
rule splice_02_rmats_adaptive:
    """
    利用二代数据的高深度进行定量。
    - 有配对：Tumor vs Normal (输出 FDR, DeltaPSI)
    - 无配对：Tumor Only (输出 PSI)
    """
    input:
        unpack(get_rmats_inputs) # 自动解包 bam_tumor 和可选的 bam_normal
    output:
        se = "rna/splicing/{sample}_rmats/SE.MATS.JC.txt",
        mxe  = "rna/splicing/{sample}_rmats/MXE.MATS.JC.txt",
        ri   = "rna/splicing/{sample}_rmats/RI.MATS.JC.txt",
        a3ss = "rna/splicing/{sample}_rmats/A3SS.MATS.JC.txt",
        a5ss = "rna/splicing/{sample}_rmats/A5SS.MATS.JC.txt"
    params:
        out_dir = "rna/splicing/{sample}_rmats",
        tmp_dir = "rna/splicing/{sample}_rmats/tmp",
        read_len = config["splicing_params"]["rna_read_length"],
        # [Fix]: 在 lambda 中不能使用 params.out_dir，改用 wildcards (w) 重构路径
        # 这里仅作为一个非空标记使用，具体的路径拼接在 shell 中通过 params.out_dir 完成更安全
        rmats_b2_arg = lambda w, input: "yes" if hasattr(input, "bam_normal") else ""
    log: "logs/splicing/rmats_{sample}.log"
    threads: 40
    conda: "../../envs/rmats.yaml"
    shell:
        """
        mkdir -p {params.out_dir} {params.tmp_dir}
        
        # 1. 准备 BAM 列表 (Tumor)
        echo {input.bam_tumor} > {params.out_dir}/b1.txt
        
        # 2. 处理 B2 (Normal) 参数逻辑
        B2_CMD=""
        
        # 判断 rmats_b2_arg 是否非空 (即是否存在 bam_normal)
        if [ -n "{params.rmats_b2_arg}" ]; then
            # 只有在有 Normal 时才写入 b2.txt
            echo {input.bam_normal} > {params.out_dir}/b2.txt
            B2_CMD="--b2 {params.out_dir}/b2.txt"
            echo "[INFO] Running rMATS in Paired Mode (Tumor vs Normal)" >> {log}
        else
            echo "[INFO] Running rMATS in Single Mode (Tumor Only)" >> {log}
        fi
        
        # 3. 运行 rMATS
        rmats.py \
            --b1 {params.out_dir}/b1.txt \
            $B2_CMD \
            --gtf {input.gtf} \
            --od {params.out_dir} \
            --tmp {params.tmp_dir} \
            -t paired \
            --readLength {params.read_len} \
            --variable-read-length \
            --allow-clipping \
            --nthread {threads} \
            --novelSS > {log} 2>&1
        """
# =============================================================================
# Rule 3: 肿瘤特异性过滤 (Python 脚本 Filter & Tiering)
# ---------------------------------------------------------------------
rule splice_03_filter_tsas_strict:
    """
    参考 IRIS (PNAS 2023) 逻辑的肿瘤特异性剪接筛选。
    - 强制 RI 使用 JCEC 文件进行内含子保留分析
    - 肿瘤端 Presence: Reads >= 5
    - 正常端 Absence: Reads <= 1 & PSI <= 0.02
    - GTEx 背景: 严格匹配 (0 count)
    """
    input:
        se   = "rna/splicing/{sample}_rmats/SE.MATS.JC.txt",
        mxe  = "rna/splicing/{sample}_rmats/MXE.MATS.JC.txt",
        ri   = "rna/splicing/{sample}_rmats/RI.MATS.JCEC.txt",
        a3ss = "rna/splicing/{sample}_rmats/A3SS.MATS.JC.txt",
        a5ss = "rna/splicing/{sample}_rmats/A5SS.MATS.JC.txt",
        gtex_db = config["splicing_params"]["filtering_db"]["gtex_junctions"]
    output:
        tsv_strict = "rna/splicing/{sample}_final/Tumor_Specific_Strict_Tier1.tsv",
        tsv_all    = "rna/splicing/{sample}_final/Tumor_Specific_Filtered_All.tsv",
        report     = "rna/splicing/{sample}_final/Full_Filter_Report.tsv",
        bed        = "rna/splicing/{sample}_final/Tumor_Specific_Strict.bed"
    params:
        min_tumor_inc_reads = 5,   # 肿瘤中支持该形态的 Reads 必须 >= 5 (IRIS Presence)
        max_normal_inc_reads = 1,  # 正常组中支持该形态的 Reads 必须 <= 1 (IRIS Absence)
        max_normal_psi = 0.02,     # 正常组 PSI 背景阈值
        min_dpsi = 0.1,            # 最小差异
        max_fdr = 0.05,
        gtex_low_freq = 10         # GTEx 罕见背景阈值
    log: "logs/splicing/filter_tsas_strict_{sample}.log"
    conda: "../../envs/python_filter.yaml"
    script:
        "../scripts/filter_splice_events.py"