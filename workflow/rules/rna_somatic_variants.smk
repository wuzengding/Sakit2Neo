# rna_somatic_variants.smk
# ========================
# 辅助函数
# ========================

def get_rna_tumor_bam(wildcards):
    """获取处理好的 Tumor RNA BAM"""
    return f"rna/snv_validation/{wildcards.sample}_tumor_analysis_ready.bam"

def get_rna_normal_bam(wildcards):
    """
    获取处理好的 Normal RNA BAM。
    注意：通常 Somatic RNA 分析建议使用 Tumor RNA vs Normal DNA (WES/WGS)。
    如果是 Tumor RNA vs Normal RNA，请确保 Normal 样本也经过了相同的 SplitNCigarReads 处理。
    """
    return f"rna/snv_validation/{wildcards.sample}_normal_analysis_ready.bam"

def get_rna_tumor_bai(wildcards):
    return f"rna/snv_validation/{wildcards.sample}_tumor_analysis_ready.bai"

def get_rna_normal_bai(wildcards):
    return f"rna/snv_validation/{wildcards.sample}_normal_analysis_ready.bai"

# 检查PON文件是否存在 (复用 DNA 逻辑)
def check_pon_exists():
    pon_path = config["reference"].get("pon", "")
    if pon_path and os.path.exists(pon_path):
        return True
    return False

# ========================
# MuTect2 RNA 变异检测 (Scatter per Chromosome)
# ========================
# 注意：这里我们直接在 GATK 命令中使用 -L {wildcards.chromosome}，
# 而不是像 DNA 流程那样读取物理分割后的 BAM 文件。这样更节省空间且适合 RNA。

rule mutect2_rna_somatic:
    input:
        tumor_bam = get_rna_tumor_bam,
        normal_bam = get_rna_normal_bam,
        #tumor_bai = get_rna_tumor_bai,
        #normal_bai = get_rna_normal_bai,
        ref = config["reference"]["genome"],
        ref_idx = config["reference"]["genome"] + ".fai",
        intervals = config["reference"]["interval_list"], # 或者是外显子 bed
        gnomad = config["reference"]["gnomad"],
        pon = config["reference"]["pon"] if check_pon_exists() else []
    output:
        vcf = "rna/variants/mutect2/scatter/{sample}.{chromosome}.vcf.gz",
        stats = "rna/variants/mutect2/scatter/{sample}.{chromosome}.vcf.gz.stats",
        f1r2 = "rna/variants/mutect2/scatter/{sample}.{chromosome}.f1r2.tar.gz"
    params:
        tumor_name = "{sample}_tumor",
        normal_name = "{sample}_normal",
        # RNA 特定参数: 不使用 soft-clipped bases 以减少剪接位点假阳性
        extra_args = "--dont-use-soft-clipped-bases true", 
        pon_arg = lambda wildcards, input: f"--panel-of-normals {input.pon}" if check_pon_exists() else ""
    log:
        "logs/rna_mutect2/{sample}.{chromosome}.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 24000,
        time = lambda wildcards, attempt: attempt * 360,
        threads = 2
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_mb}m" Mutect2 \
            -R {input.ref} \
            -I {input.tumor_bam} -tumor {params.tumor_name} \
            -I {input.normal_bam} -normal {params.normal_name} \
            -L chr{wildcards.chromosome} \
            --germline-resource {input.gnomad} \
            {params.pon_arg} \
            {params.extra_args} \
            --genotype-germline-sites true \
            --f1r2-tar-gz {output.f1r2} \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} \
            2> {log}
        """

# ========================
# 获取污染信息 (Pileup)
# ========================
rule get_rna_pileup_summaries:
    input:
        bam = "rna/snv_validation/{sample}_{sampletype}_analysis_ready.bam",
        #bai = "rna/snv_validation/{sample}_{sampletype}_analysis_ready.bai",
        sites = config["reference"]["gnomad"]
    output:
        pileup = "rna/variants/mutect2/scatter/{sample}_{sampletype}.{chromosome}.pileups.table"
    log:
        "logs/rna_mutect2/{sample}_{sampletype}.{chromosome}.pileup.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.sites} \
            -L chr{wildcards.chromosome} \
            -O {output.pileup} \
            2> {log}
        """

# ========================
# 计算污染率 (Gather Pileups & Calculate)
# ========================
rule calculate_rna_contamination:
    input:
        tumor_pileups = expand("rna/variants/mutect2/scatter/{{sample}}_tumor.{chromosome}.pileups.table",
                        chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        normal_pileups = expand("rna/variants/mutect2/scatter/{{sample}}_normal.{chromosome}.pileups.table",
                        chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        contamination = "rna/variants/mutect2/{sample}.contamination.table",
        segments = "rna/variants/mutect2/{sample}.segments.table"
    log:
        "logs/rna_mutect2/{sample}.contamination.log"
    conda:
        "../../envs/gatk.yaml"
    params:
        tempdir = "rna/variants/mutect2"
    shell:
        """
        # 合并分散的 pileup tables
        TUMOR_MERGED=$(mktemp --tmpdir={params.tempdir} --suffix='.pileups.table')
        NORMAL_MERGED=$(mktemp --tmpdir={params.tempdir} --suffix='.pileups.table')

        (head -n 1 {input.tumor_pileups[0]} && tail -n +2 -q {input.tumor_pileups}) > $TUMOR_MERGED
        (head -n 1 {input.normal_pileups[0]} && tail -n +2 -q {input.normal_pileups}) > $NORMAL_MERGED

        gatk CalculateContamination \
            -I $TUMOR_MERGED \
            -matched $NORMAL_MERGED \
            -O {output.contamination} \
            --tumor-segmentation {output.segments} \
            2> {log}
        """

# ========================
# 学习读取方向模型 (Artifacts)
# ========================
rule learn_rna_read_orientation_model:
    input:
        f1r2s = expand("rna/variants/mutect2/scatter/{{sample}}.{chromosome}.f1r2.tar.gz", 
                chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        model = "rna/variants/mutect2/{sample}.read-orientation-model.tar.gz"
    log:
        "logs/rna_mutect2/{sample}.orientation_model.log"
    conda:
        "../../envs/gatk.yaml"
    shell:
        """
        INPUT_ARGS=$(for f1r2 in {input.f1r2s}; do echo -n "-I $f1r2 "; done)
        gatk LearnReadOrientationModel $INPUT_ARGS -O {output.model} 2> {log}
        """

# ========================
# FilterMutectCalls (RNA)
# ========================
rule filter_rna_mutect_calls:
    input:
        vcf = "rna/variants/mutect2/scatter/{sample}.{chromosome}.vcf.gz",
        vcf_tbi = "rna/variants/mutect2/scatter/{sample}.{chromosome}.vcf.gz.tbi",
        stats = "rna/variants/mutect2/scatter/{sample}.{chromosome}.vcf.gz.stats",
        contamination = "rna/variants/mutect2/{sample}.contamination.table",
        segments = "rna/variants/mutect2/{sample}.segments.table",
        orientation_model = "rna/variants/mutect2/{sample}.read-orientation-model.tar.gz",
        ref = config["reference"]["genome"]
    output:
        vcf = "rna/variants/mutect2/scatter/{sample}.{chromosome}.filtered.vcf.gz",
        filtering_stats = "rna/variants/mutect2/scatter/{sample}.{chromosome}.filtering.stats"
    log: "logs/rna_mutect2/{sample}.{chromosome}.filter.log"
    conda: "../../envs/gatk.yaml"
    shell:
        """
        gatk FilterMutectCalls \
            -R {input.ref} \
            -V {input.vcf} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --orientation-bias-artifact-priors {input.orientation_model} \
            --stats {input.stats} \
            --filtering-stats {output.filtering_stats} \
            -O {output.vcf} \
            2> {log}
        """

# ========================
# 聚合 RNA Mutect2 VCFs
# ========================
ruleorder: gather_rna_vcfs > index_vcf
rule gather_rna_vcfs:
    input:
        vcfs = expand("rna/variants/mutect2/scatter/{{sample}}.{chromosome}.filtered.vcf.gz",
               chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        vcfs_tbi = expand("rna/variants/mutect2/scatter/{{sample}}.{chromosome}.filtered.vcf.gz.tbi",
               chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        vcf = "rna/variants/mutect2/{sample}.mutect2.vcf.gz",
        vcf_idx = "rna/variants/mutect2/{sample}.mutect2.vcf.gz.tbi"
    log: "logs/rna_mutect2/{sample}.gather.log"
    conda: "../../envs/gatk.yaml"
    shell:
        """
        INPUT_ARGS=$(for vcf in {input.vcfs}; do echo -n "-I $vcf "; done)
        gatk GatherVcfs $INPUT_ARGS -O {output.vcf} 2> {log}
        gatk IndexFeatureFile -I {output.vcf} 2>> {log}
        """

# ========================
# VarScan2 RNA 变异检测
# ========================
rule varscan2_rna_somatic_scatter:
    input:
        tumor_bam = get_rna_tumor_bam,
        normal_bam = get_rna_normal_bam,
        ref = config["reference"]["genome"]
    output:
        snp_vcf = "rna/variants/varscan2/scatter/{sample}.{chromosome}.snp.vcf.gz",
        indel_vcf = "rna/variants/varscan2/scatter/{sample}.{chromosome}.indel.vcf.gz"
    params:
        prefix_uncompressed = "rna/variants/varscan2/scatter/{sample}.{chromosome}",
        min_coverage = config["varscan2"]["min_coverage"],
        min_var_freq = config["varscan2"]["min_var_freq"],
        p_value = config["varscan2"]["p_value"],
        somatic_p = config["varscan2"]["somatic_p_value"]
    log:
        "logs/rna_varscan2/{sample}.{chromosome}.log"
    conda:
        "../../envs/varscan.yaml"
    resources:
        mem_mb = 12000
    shell:
        """
        # 使用 samtools -r 限制区域，避免读取整个BAM
        # 注意这里加上 chr 前缀
        (samtools mpileup \
            -f {input.ref} \
            -r chr{wildcards.chromosome} \
            -q 1 -Q 13 \
            {input.normal_bam} {input.tumor_bam} \
        | varscan somatic \
            /dev/stdin \
            {params.prefix_uncompressed} \
            --mpileup 1 \
            --min-coverage {params.min_coverage} \
            --min-var-freq {params.min_var_freq} \
            --p-value {params.p_value} \
            --somatic-p-value {params.somatic_p} \
            --strand-filter 1 \
            --output-vcf 1 \
        ) > /dev/null 2> {log}

        bgzip {params.prefix_uncompressed}.snp.vcf
        tabix -p vcf {output.snp_vcf}

        bgzip {params.prefix_uncompressed}.indel.vcf
        tabix -p vcf {output.indel_vcf}
        """

# ========================
# VarScan2 聚合规则 (复用 DNA 逻辑，修改路径)
# ========================
rule varscan2_rna_gather:
    input:
        snp_vcfs = expand("rna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz",
                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        indel_vcfs = expand("rna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz",
                            chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        ref_fai = config["reference"]["genome"] + ".fai"
    output:
        snp_vcf = "rna/variants/varscan2/{sample}.snp.vcf",
        indel_vcf = "rna/variants/varscan2/{sample}.indel.vcf",
        snp_vcf_gz = "rna/variants/varscan2/{sample}.snp.vcf.gz",
        indel_vcf_gz = "rna/variants/varscan2/{sample}.indel.vcf.gz"
    log:
        "logs/rna_varscan2/{sample}.gather.log"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        # (代码逻辑同 DNA 脚本，只是文件路径变为 rna/variants/varscan2/...)
        
        SNP_LIST_FILE={wildcards.sample}.rna_snp_list.txt
        INDEL_LIST_FILE={wildcards.sample}.rna_indel_list.txt

        for f in {input.snp_vcfs}; do if [ -f "$f" ]; then echo "$f" >> $SNP_LIST_FILE; fi; done
        bcftools concat -a -f $SNP_LIST_FILE > {output.snp_vcf} 2> {log}
        bgzip -k {output.snp_vcf} && tabix -p vcf {output.snp_vcf_gz}

        for f in {input.indel_vcfs}; do if [ -f "$f" ]; then echo "$f" >> $INDEL_LIST_FILE; fi; done
        bcftools concat -a -f $INDEL_LIST_FILE > {output.indel_vcf} 2>> {log}
        bgzip -k {output.indel_vcf} && tabix -p vcf {output.indel_vcf_gz}

        rm -f $SNP_LIST_FILE $INDEL_LIST_FILE
        """

rule varscan2_rna_reformat:
    input:
        snp_vcf = "rna/variants/varscan2/{sample}.snp.vcf",
        indel_vcf = "rna/variants/varscan2/{sample}.indel.vcf"
    output:
        vcf = "rna/variants/varscan2/{sample}.varscan2.vcf.gz"
    params:
        prefix = "rna/variants/varscan2/{sample}",
        script_dir = workflow.basedir
    conda:
        "../../envs/python.yaml"
    shell:
        """
        # 1. 合并 (增加 -h 修复染色体名称前带文件名的问题)
        # 增加对 ## 头的提取，确保 Python 脚本能识别 VCF 格式
        (
            grep -h "^##" {input.snp_vcf} || true;
            grep -h "^#CHROM" {input.snp_vcf} || true;
            grep -h -v "^#" {input.snp_vcf} {input.indel_vcf} | sort -k1,1V -k2,2n 
        ) > {params.prefix}.combined.vcf

        # 2. 调用 Python 脚本 (传入 3 个参数)
        # 第 1 参数: 输入的合并文本
        # 第 2 参数: 脚本内部需要的中间文本路径 (脚本处理完后会自动删除它)
        # 第 3 参数: 最终生成的压缩 VCF (.gz) —— 直接指向 Snakemake 的 output
        python {params.script_dir}/scripts/reformat_varscan_vcf.py \
            {params.prefix}.combined.vcf \
            {params.prefix}.reformatted.vcf \
            {output.vcf}

        # 3. 清理合并产生的临时文本文件
        rm -f {params.prefix}.combined.vcf
        """
