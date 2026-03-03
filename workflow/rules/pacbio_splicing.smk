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
    try:
        # 确保 samples 是全局 pandas DataFrame
        subset = samples[(samples['sample_id'] == sample_id) & 
                         (samples['sampletype'] == 'normal') & 
                         (samples['data_type'] == 'rna')]
        return not subset.empty
    except Exception as e:
        print(f"[WARNING] Could not determine paired status for {sample_id}: {e}")
        return False

def get_sqanti_input_gtf(wildcards):
    """
    自适应输入 GTF：
    - 有配对：使用 StringTie2 Merge 后的 GTF (T+N+Isoseq)
    - 无配对：使用 Tumor 单独的 Hybrid GTF
    """
    if has_paired_normal(wildcards.sample):
        return f"rna/assembly/{wildcards.sample}_merged.assembly.gtf.gz"
    else:
        return f"rna/assembly/{wildcards.sample}_tumor.assembly.gtf.gz"

def get_rmats_inputs(wildcards):
    """
    获取 rMATS 所需的 BAM 文件。
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
    处理本地下载的 Snaptron/GTEx v2 junctions.bgz。
    Snaptron 格式通常为: id, chr, start, end, ..., coverage_sum (Col 11)
    """
    input:
        bgz = "resources/recount3_db/junctions.bgz" # 请确保路径正确
    output:
        gtex_bed = config["splicing_params"]["filtering_db"]["gtex_junctions"]
    log: "logs/splicing/process_gtex_v8.log"
    conda: "../../envs/python_filter.yaml"
    shell:
        """
        echo "[*] Processing local GTEx v8 junctions..." > {log}
        
        # 使用 gzip -dc 处理
        # 提取 chr(2), start(3), end(4), count(11)
        # 过滤掉 count < 10 的极低频噪音
        
        gzip -dc {input.bgz} | \
        awk -F "\\t" 'NR>1 && $11 >= 10 {{
            chrom = $2;
            if (index(chrom, "chr") != 1) {{ chrom = "chr"chrom }};
            print chrom"\\t"$3"\\t"$4"\\t"$11
        }}' | \
        sort -k1,1 -k2,2n > {output.gtex_bed} 2>> {log}
        
        if [ ! -s {output.gtex_bed} ]; then
            echo "[ERROR] Output file is empty! Check input file format." >> {log}
            exit 1
        fi
        """

# =============================================================================
# Rule 1: SQANTI3 质控 (极速版：跳过 ORF)
# =============================================================================
rule splice_01_sqanti3_qc:
    """
    StringTie2 GTF 结构纠正。启用 --skipORF 提速。
    解决 libcurl 缺失问题。
    """
    input:
        gtf = get_sqanti_input_gtf,
        ref_gtf = config["reference"]["gtf"],
        ref_fasta = config["reference"]["genome"]
    output:
        gtf = "rna/splicing/{sample}_sqanti3/{sample}_corrected.gtf",
        classification = "rna/splicing/{sample}_sqanti3/{sample}_classification.txt",
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
        
        # [Fix] 解决 libcurl.so.4 not found
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

        # 4. 执行 SQANTI3 (--skipORF 核心提速)
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
            
        rm -f $WORK_GTF
        """

# =============================================================================
# Rule 2: rMATS (自适应模式)
# =============================================================================
rule splice_02_rmats_adaptive:
    input:
        unpack(get_rmats_inputs)
    output:
        se = "rna/splicing/{sample}_rmats/SE.MATS.JC.txt",
        mxe  = "rna/splicing/{sample}_rmats/MXE.MATS.JC.txt",
        ri   = "rna/splicing/{sample}_rmats/RI.MATS.JCEC.txt",
        a3ss = "rna/splicing/{sample}_rmats/A3SS.MATS.JC.txt",
        a5ss = "rna/splicing/{sample}_rmats/A5SS.MATS.JC.txt"
    params:
        out_dir = "rna/splicing/{sample}_rmats",
        tmp_dir = "rna/splicing/{sample}_rmats/tmp",
        read_len = config["splicing_params"]["rna_read_length"],
        # 仅作为标记传递给 shell，具体路径在 shell 中处理
        has_normal = lambda w, input: "yes" if hasattr(input, "bam_normal") else ""
    log: "logs/splicing/rmats_{sample}.log"
    threads: 40
    conda: "../../envs/rmats.yaml"
    shell:
        """
        mkdir -p {params.out_dir} {params.tmp_dir}
        echo {input.bam_tumor} > {params.out_dir}/b1.txt
        
        B2_CMD=""
        if [ -n "{params.has_normal}" ]; then
            echo {input.bam_normal} > {params.out_dir}/b2.txt
            B2_CMD="--b2 {params.out_dir}/b2.txt"
            echo "[INFO] Mode: Paired (Tumor vs Normal)" >> {log}
        else
            echo "[INFO] Mode: Single (Tumor Only)" >> {log}
        fi
        
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
            --novelSS >> {log} 2>&1
        """

# =============================================================================
# Rule 3: 肿瘤特异性过滤 (Python Script)
# =============================================================================
rule splice_03_filter_tsas_strict:
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
        min_tumor_inc_reads = 5,
        max_normal_inc_reads = 1,
        max_normal_psi = 0.05,
        min_dpsi = 0.1,
        max_fdr = 0.05,
        gtex_low_freq = 10
    log: "logs/splicing/filter_tsas_strict_{sample}.log"
    conda: "../../envs/python_filter.yaml"
    script:
        "../../scripts/filter_splice_events.py"