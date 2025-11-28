# workflow/rules/somatic_variants.smk
# ========================
# 检查PON文件是否存在的函数
# ========================
def check_pon_exists():
    """检查PON文件是否存在"""
    pon_path = config["reference"].get("pon", "")
    if pon_path and os.path.exists(pon_path):
        return True
    return False

def get_pon_input():
    """根据PON是否存在返回相应的输入"""
    if check_pon_exists():
        return config["reference"]["pon"]
    else:
        return []

# ========================
# MuTect2 变异检测 (有PON情况)
# ========================
rule mutect2_somatic_with_pon:
    input:
        tumor_bam = get_tumor_bam_scattered,
        normal_bam = get_normal_bam_scattered,
        tumor_bai = get_tumor_bai_scattered,
        normal_bai = get_normal_bai_scattered,
        ref = config["reference"]["genome"],
        ref_idx = config["reference"]["genome"] + ".fai",
        intervals = config["reference"]["interval_list"],
        gnomad = config["reference"]["gnomad"],
        pon = config["reference"]["pon"]
    output:
        vcf = "dna/variants/mutect2/scatter/{sample}.{chromosome}.withpon.vcf.gz",
        stats = "dna/variants/mutect2/scatter/{sample}.{chromosome}.withpon.vcf.gz.stats"
    params:
        tumor_name = "{sample}_tumor",
        normal_name = "{sample}_normal",
        tumor_lod = config["mutect2"]["tumor_lod"],
        normal_lod = config["mutect2"]["normal_lod"],
        af_threshold = config["mutect2"]["af_of_alleles_not_in_resource"],
        max_reads = config["mutect2"]["max_reads_per_alignment_start"]
    log:
        "logs/mutect2/{sample}.{chromosome}.mutect2_with_pon.log"
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
            -L {input.intervals} \
            --germline-resource {input.gnomad} \
            --panel-of-normals {input.pon} \
            --genotype-germline-sites true \
            --tumor-lod-to-emit {params.tumor_lod} \
            --normal-lod {params.normal_lod} \
            --af-of-alleles-not-in-resource {params.af_threshold} \
            --max-reads-per-alignment-start {params.max_reads} \
            --dont-use-soft-clipped-bases \
            --native-pair-hmm-threads {threads} \
            -O {output.vcf} \
            2> {log}
        """

# ========================
# MuTect2 变异检测 (无PON情况 - 原始调用)
# ========================
rule mutect2_somatic_no_pon:
    input:
        tumor_bam = get_tumor_bam_scattered,
        normal_bam = get_normal_bam_scattered,
        tumor_bai = get_tumor_bai_scattered,
        normal_bai = get_normal_bai_scattered,
        ref = config["reference"]["genome"],
        ref_idx = config["reference"]["genome"] + ".fai",
        intervals = config["reference"]["interval_list"],
        gnomad = config["reference"]["gnomad"]
    output:
        vcf_raw = "dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz",
        stats_raw = "dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz.stats",
        f1r2 = "dna/variants/mutect2/scatter/{sample}.{chromosome}.f1r2.tar.gz"
    params:
        tumor_name = "{sample}_tumor",
        normal_name = "{sample}_normal",
        tumor_lod = config["mutect2"]["tumor_lod"],
        normal_lod = config["mutect2"]["normal_lod"],
        af_threshold = config["mutect2"]["af_of_alleles_not_in_resource"],
        max_reads = config["mutect2"]["max_reads_per_alignment_start"]
    log:
        "logs/mutect2/{sample}.{chromosome}.mutect2_no_pon.log"
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
            -L {input.intervals} \
            --germline-resource {input.gnomad} \
            --genotype-germline-sites true \
            --tumor-lod-to-emit {params.tumor_lod} \
            --normal-lod {params.normal_lod} \
            --af-of-alleles-not-in-resource {params.af_threshold} \
            --max-reads-per-alignment-start {params.max_reads} \
            --dont-use-soft-clipped-bases \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2} \
            -O {output.vcf_raw} \
            2> {log}
        """

# ========================
# 获取污染信息 (无PON情况需要)
# ========================
rule get_pileup_summaries:
    input:
        bam = "dna/aligned/temp/{sample}_{sampletype}_{chromosome}.recal.bam",
        bai = "dna/aligned/temp/{sample}_{sampletype}_{chromosome}.recal.bam.bai",
        intervals = config["reference"]["interval_list"],
        sites = config["reference"]["gnomad"]
    output:
        pileup = "dna/variants/mutect2/scatter/{sample}_{sampletype}.{chromosome}.pileups.table"
    log:
        "logs/mutect2/{sample}_{sampletype}.{chromosome}.pileup.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60,
        threads = 2
    shell:
        """
        gatk GetPileupSummaries \
            -I {input.bam} \
            -V {input.sites} \
            -L {input.intervals} \
            -O {output.pileup} \
            2> {log}
        """

# ========================
# 计算污染率
# ========================
rule calculate_contamination:
    input:
        tumor_pileups = expand("dna/variants/mutect2/scatter/{{sample}}_tumor.{chromosome}.pileups.table",
                        chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        normal_pileups = expand("dna/variants/mutect2/scatter/{{sample}}_normal.{chromosome}.pileups.table",
                        chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        contamination = "dna/variants/mutect2/{sample}.contamination.table",
        segments = "dna/variants/mutect2/{sample}.segments.table"
    log:
        "logs/mutect2/{sample}.contamination.log"
    conda:
        "../../envs/gatk.yaml"
    params:
        tempdir = "dna/variants/mutect2"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = lambda wildcards, attempt: attempt * 30,
        threads = 1
    shell:
        """
        # 创建临时的合并文件
        TUMOR_MERGED=$(mktemp --tmpdir={params.tempdir} --suffix='.pileups.table')
        NORMAL_MERGED=$(mktemp --tmpdir={params.tempdir} --suffix='.pileups.table')

        # 安全地合并肿瘤pileup文件
        (head -n 1 {input.tumor_pileups[0]} && tail -n +2 -q {input.tumor_pileups}) > $TUMOR_MERGED

        # 安全地合并正常样本pileup文件
        (head -n 1 {input.normal_pileups[0]} && tail -n +2 -q {input.normal_pileups}) > $NORMAL_MERGED

        # 使用合并后的临时文件运行GATK
        gatk CalculateContamination \
            -I $TUMOR_MERGED \
            -matched $NORMAL_MERGED \
            -O {output.contamination} \
            --tumor-segmentation {output.segments} \
            2> {log}
        """

# ========================
# 学习读取方向模型
# ========================
rule learn_read_orientation_model:
    input:
        f1r2s = expand("dna/variants/mutect2/scatter/{{sample}}.{chromosome}.f1r2.tar.gz", 
                chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        model = "dna/variants/mutect2/{sample}.read-orientation-model.tar.gz"
    log:
        "logs/mutect2/{sample}.orientation_model.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60,
        threads = 1
    shell:
        """
        INPUT_ARGS=$(for f1r2 in {input.f1r2s}; do echo -n "-I $f1r2 "; done)
        gatk LearnReadOrientationModel \
             $INPUT_ARGS \
            -O {output.model} \
            2> {log}
        """

# ========================
# FilterMutectCalls 过滤 (无PON情况)
# ========================
rule filter_mutect_calls:
    input:
        vcf = "dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz",
        vcf_tbi = "dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz.tbi",
        stats = "dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz.stats",
        contamination = "dna/variants/mutect2/{sample}.contamination.table",
        segments = "dna/variants/mutect2/{sample}.segments.table",
        orientation_model = "dna/variants/mutect2/{sample}.read-orientation-model.tar.gz",
        ref = config["reference"]["genome"]
    output:
        vcf = "dna/variants/mutect2/scatter/{sample}.{chromosome}.nopon.filtered.vcf.gz",
        filtering_stats = "dna/variants/mutect2/scatter/{sample}.{chromosome}.nopon.filtering.stats"
    log: "logs/mutect2/{sample}.{chromosome}.filter_mutect.log"
    conda: "../../envs/gatk.yaml"
    resources: mem_mb = 8000, time = 45
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
# 最终聚合VCF文件 - Gather
# ========================
rule gather_vcfs_with_pon:
    input:
        vcfs = expand("dna/variants/mutect2/scatter/{{sample}}.{chromosome}.withpon.vcf.gz",
               chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        vcfs_tbi = expand("dna/variants/mutect2/scatter/{{sample}}.{chromosome}.withpon.vcf.gz.tbi",
               chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.withpon.vcf.gz"
    log: "logs/gatk/{sample}.gather_withpon.log"
    conda: "../../envs/gatk.yaml"
    resources: mem_mb = 8000, time = 60
    shell:
        """
        INPUT_ARGS=$(for vcf in {input.vcfs}; do echo -n "-I $vcf "; done)
        gatk GatherVcfs $INPUT_ARGS -O {output.vcf} 2> {log}
        """

rule gather_vcfs_no_pon:
    input:
        vcfs = expand("dna/variants/mutect2/scatter/{{sample}}.{chromosome}.nopon.filtered.vcf.gz",
               chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        vcfs_tbi = expand("dna/variants/mutect2/scatter/{{sample}}.{chromosome}.nopon.filtered.vcf.gz.tbi",
               chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.nopon.vcf.gz",
    log: 
        "logs/gatk/{sample}.gather_nopon.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = 8000, time = 60
    shell:
        """
        INPUT_ARGS=$(for vcf in {input.vcfs}; do echo -n "-I $vcf "; done)
        gatk GatherVcfs $INPUT_ARGS -O {output.vcf} 2> {log}
        """

# ========================
# 条件规则选择器
# ========================
def get_mutect2_input(wildcards):
    """根据PON是否存在决定使用哪个规则的输出"""
    if check_pon_exists():
        return {
            "vcf": f"dna/variants/mutect2/{wildcards.sample}.mutect2.withpon.vcf.gz",
            "vcf_idx": f"dna/variants/mutect2/{wildcards.sample}.mutect2.withpon.vcf.gz.tbi", 
        }
    else:
        return {
            "vcf": f"dna/variants/mutect2/{wildcards.sample}.mutect2.nopon.vcf.gz",
            "vcf_idx": f"dna/variants/mutect2/{wildcards.sample}.mutect2.nopon.vcf.gz.tbi",
        }


# ========================
# 规则聚合器 - 确保选择正确的流程
# ========================
ruleorder: mutect2_caller > index_vcf
rule mutect2_caller:
    input:
        unpack(get_mutect2_input)
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.vcf.gz",
        vcf_idx = "dna/variants/mutect2/{sample}.mutect2.vcf.gz.tbi"
    log:
        "logs/mutect2/{sample}.finalize.log"
    run:
        # 获取输出文件所在的目录
        output_dir = os.path.dirname(output.vcf)
        
        # 计算输入文件相对于输出目录的相对路径
        target_vcf = os.path.relpath(input.vcf, output_dir)
        target_idx = os.path.relpath(input.vcf_idx, output_dir)

        # 在shell中执行创建链接的命令
        # 注意：这里需要重定向日志
        shell(f"ln -sf {target_vcf} {output.vcf} 2> {log}")
        shell(f"ln -sf {target_idx} {output.vcf}.tbi 2>> {log}")

# ========================
# 为 VCF.GZ 文件创建索引
# ========================
rule index_vcf:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.vcf.gz.tbi"
    log:
        "logs/gatk/index_{prefix}.vcf.log"
    conda:
        "../../envs/gatk.yaml"  # 使用与GATK相同的环境
    resources:
        mem_mb = 4000
    shell:
        # GATK的IndexFeatureFile会自动生成.tbi后缀的索引文件
        "gatk IndexFeatureFile -I {input} 2> {log}"

# ========================
# Strelka2 变异检测 (辅助验证)
# ========================
rule strelka2_somatic:
    input:
        tumor_bam = get_tumor_bam,
        normal_bam = get_normal_bam,
        tumor_bai = get_tumor_bai,
        normal_bai = get_normal_bai,
        ref = config["reference"]["genome"],
        ref_idx = config["reference"]["genome"] + ".fai",
        bed = config["reference"]["capture_kit_bedgz"]
    output:
        vcf_snv = "dna/variants/strelka2/{sample}/results/variants/somatic.snvs.vcf.gz",
        vcf_indel = "dna/variants/strelka2/{sample}/results/variants/somatic.indels.vcf.gz",
        vcf_combined = "dna/variants/strelka2/{sample}.strelka2.vcf.gz"
    params:
        rundir = "dna/variants/strelka2/{sample}",
        tumor_name = "{sample}_tumor",
        normal_name = "{sample}_normal"
    log:
        config = "logs/strelka2/{sample}.config.log",
        run = "logs/strelka2/{sample}.run.log",
        merge = "logs/strelka2/{sample}.merge.log"
    conda:
        "../../envs/strelka.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        time = lambda wildcards, attempt: attempt * 240,
        threads = 16
    shell:
        """
        # 配置Strelka2
        configureStrelkaSomaticWorkflow.py \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam} \
            --referenceFasta {input.ref} \
            --targeted --exome \
            --callRegions {input.bed} \
            --runDir {params.rundir} \
            2> {log.config}
        
        # 运行分析
        {params.rundir}/runWorkflow.py \
            -m local \
            -j {threads} \
            2> {log.run}
        
        # 合并SNV和INDEL结果
        bcftools concat \
            {output.vcf_snv} {output.vcf_indel} \
            -a -O z -o {output.vcf_combined} \
            2> {log.merge}
        
        # 创建索引
        tabix -p vcf {output.vcf_combined}
        """

#rule varscan2_call:
#    input:
#        tumor_bam = get_tumor_bam,
#        normal_bam = get_normal_bam,
#        tumor_bai = get_tumor_bai,
#        normal_bai = get_normal_bai,
#        ref = config["reference"]["genome"],
#        intervals = config["reference"]["capture_kit"]
#    output:
#        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf",
#        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf"
#    params:
#        prefix = "dna/variants/varscan2/{sample}",
#        min_coverage = config["varscan2"]["min_coverage"],
#        min_var_freq = config["varscan2"]["min_var_freq"],
#        p_value = config["varscan2"]["p_value"],
#        somatic_p = config["varscan2"]["somatic_p_value"],
#        tumor_name = "{sample}_tumor",
#        normal_name = "{sample}_normal"
#    log:
#        mpileup = "logs/varscan2/{sample}.mpileup.log",
#        varscan = "logs/varscan2/{sample}.varscan.log"
#    conda:
#        "../../envs/varscan.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * 16000,
#        time = lambda wildcards, attempt: attempt * 240,
#    shell:
#        """
#        # 1. 生成mpileup文件
#        samtools mpileup \
#            -f {input.ref} \
#            -l {input.intervals} \
#            -q 1 -Q 13 \
#            {input.normal_bam} {input.tumor_bam} \
#            > {params.prefix}.mpileup \
#            2> {log.mpileup}
#        
#        # 2. VarScan2体细胞调用
#        varscan somatic \
#            {params.prefix}.mpileup \
#            {params.prefix} \
#            --mpileup 1 \
#            --min-coverage {params.min_coverage} \
#            --min-var-freq {params.min_var_freq} \
#            --p-value {params.p_value} \
#            --somatic-p-value {params.somatic_p} \
#            --strand-filter 1 \
#            --output-vcf 1 \
#            2> {log.varscan}
#        
#        # 3. 清理中间文件
#        #rm -f {params.prefix}.mpileup
#        """

# ========================
# VarScan2 变异检测 (按染色体分散，并输出压缩VCF)
# ========================
ruleorder: varscan2_somatic_scatter > index_vcf
rule varscan2_somatic_scatter:
    input:
        tumor_bam = get_tumor_bam_scattered,
        normal_bam = get_normal_bam_scattered,
        tumor_bai = get_tumor_bai_scattered,
        normal_bai = get_normal_bai_scattered,
        ref = config["reference"]["genome"],
        intervals = config["reference"]["capture_kit"]
    output:
        # **修改输出为 .vcf.gz 和 .vcf.gz.tbi**
        snp_vcf = "dna/variants/varscan2/scatter/{sample}.{chromosome}.snp.vcf.gz",
        snp_tbi = "dna/variants/varscan2/scatter/{sample}.{chromosome}.snp.vcf.gz.tbi",
        indel_vcf = "dna/variants/varscan2/scatter/{sample}.{chromosome}.indel.vcf.gz",
        indel_tbi = "dna/variants/varscan2/scatter/{sample}.{chromosome}.indel.vcf.gz.tbi"
    params:
        prefix_uncompressed = "dna/variants/varscan2/scatter/{sample}.{chromosome}",
        # ... 其他params保持不变 ...
        min_coverage = config["varscan2"]["min_coverage"],
        min_var_freq = config["varscan2"]["min_var_freq"],
        p_value = config["varscan2"]["p_value"],
        somatic_p = config["varscan2"]["somatic_p_value"]
    log:
        "logs/varscan2/{sample}.{chromosome}.log"
    conda:
        # **重要：确保环境中包含 bgzip 和 tabix (通常samtools或htslib提供)**
        "../../envs/varscan.yaml" 
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 120
    shell:
        """
        # 1. 运行 VarScan2，生成临时的未压缩VCF
        (samtools mpileup \
            -f {input.ref} \
            -l {input.intervals} \
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

        # 2. 压缩和索引 SNP VCF
        bgzip {params.prefix_uncompressed}.snp.vcf
        tabix -p vcf {output.snp_vcf}

        # 3. 压缩和索引 INDEL VCF
        bgzip {params.prefix_uncompressed}.indel.vcf
        tabix -p vcf {output.indel_vcf}
        """

# ========================
# 最终聚合VarScan2的VCF文件 - Gather (处理压缩VCF)
# ========================
#rule varscan2_gather_vcfs:
#    input:
#        # **输入现在是 .vcf.gz 和 .vcf.gz.tbi**
#        snp_vcfs = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz", 
#                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        snp_tbis = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz.tbi", 
#                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        indel_vcfs = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz", 
#                            chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        indel_tbis = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz.tbi", 
#                            chromosome=[c.replace("chr", "") for c in get_chromosomes()])
#    output:
#        # **输出也应该是标准的压缩格式**
#        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf.gz",
#        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf.gz"
#    log:
#        snp_log = "logs/varscan2/{sample}.gather_snp.log",
#        indel_log = "logs/varscan2/{sample}.gather_indel.log"
#    conda:
#        "../../envs/gatk.yaml" # 使用GATK环境
#    resources:
#        mem_mb = 8000, 
#        time = 60
#    shell:
#        """
#        # 和您的 mutect2 gather 规则几乎完全一样
#        SNP_INPUT_ARGS=$(for vcf in {input.snp_vcfs}; do echo -n "-I $vcf "; done)
#        gatk --java-options "-Xmx{resources.mem_mb}m" GatherVcfs $SNP_INPUT_ARGS -O {output.snp_vcf} 2> {log.snp_log}
#        
#        INDEL_INPUT_ARGS=$(for vcf in {input.indel_vcfs}; do echo -n "-I $vcf "; done)
#        gatk --java-options "-Xmx{resources.mem_mb}m" GatherVcfs $INDEL_INPUT_ARGS -O {output.indel_vcf} 2> {log.indel_log}
#        """

## ========================
## 最终聚合VarScan2的VCF文件 - Gather (合并、修复、压缩)
## ========================
#rule varscan2_gather_vcfs:
#    input:
#        snp_vcfs = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz", 
#                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        snp_idx = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz.tbi", 
#                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        indel_vcfs = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz", 
#                            chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        indel_idx = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz.tbi", 
#                            chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
#        ref_fai = config["reference"]["genome"] + ".fai"
#    output:
#        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf.gz",
#        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf.gz"
#    log:
#        "logs/varscan2/{sample}.gather.log"
#    conda:
#        "../../envs/bcftools.yaml" # 确保有 bcftools, bgzip, tabix
#    resources:
#        mem_mb = 8000, 
#        time = 60
#    shell:
#        """
#        # 定义一些临时文件名
#        COMBINED_SNP_VCF=dna/variants/varscan2/{wildcards.sample}.combined.snp.vcf
#        COMBINED_INDEL_VCF=dna/variants/varscan2/{wildcards.sample}.combined.indel.vcf
#        HEADER_FILE=dna/variants/varscan2/{wildcards.sample}.header.txt
#
#        # 步骤 1: 使用 bcftools 合并所有 SNP VCF 文件
#        # -a 允许不同样本集, -f 允许合并列表中的第一个文件不存在(以防某个染色体没变异)
#        bcftools concat -a -f {input.snp_vcfs} > $COMBINED_SNP_VCF 2> {log}
#
#        # 步骤 2: 使用 bcftools 合并所有 INDEL VCF 文件
#        bcftools concat -a -f {input.indel_vcfs} > $COMBINED_INDEL_VCF 2>> {log}
#
#        # 步骤 3: 从参考基因组 .fai 文件生成 contig 头部信息
#        awk 'BEGIN{{FS="\\t"}}; {{print "##contig=<ID="$1",length="$2">"}}' {input.ref_fai} > $HEADER_FILE
#
#        # 步骤 4: 修复头部、压缩并索引最终的 SNP VCF
#        bcftools reheader -h $HEADER_FILE $COMBINED_SNP_VCF | bgzip > {output.snp_vcf}
#        tabix -p vcf {output.snp_vcf}
#
#        # 步骤 5: 修复头部、压缩并索引最终的 INDEL VCF
#        bcftools reheader -h $HEADER_FILE $COMBINED_INDEL_VCF | bgzip > {output.indel_vcf}
#        tabix -p vcf {output.indel_vcf}
#
#        # 步骤 6: 清理临时合并文件和头部文件
#        rm -f $COMBINED_SNP_VCF $COMBINED_INDEL_VCF $HEADER_FILE
#        """

rule varscan2_gather_vcfs:
    input:
        snp_vcfs = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz",
                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        snp_idx = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.snp.vcf.gz.tbi",
                          chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        indel_vcfs = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz",
                            chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        indel_idx = expand("dna/variants/varscan2/scatter/{{sample}}.{chromosome}.indel.vcf.gz.tbi",
                            chromosome=[c.replace("chr", "") for c in get_chromosomes()]),
        ref_fai = config["reference"]["genome"] + ".fai"
    output:
        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf",
        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf",
        snp_vcf_gz = "dna/variants/varscan2/{sample}.snp.vcf.gz",
        indel_vcf_gz = "dna/variants/varscan2/{sample}.indel.vcf.gz"
    log:
        "logs/varscan2/{sample}.gather.log"
    conda:
        "../../envs/bcftools.yaml"
    resources:
        mem_mb = 8000,
        time = 60
    shell:
        """
        # 定义临时文件
        SNP_LIST_FILE={wildcards.sample}.snp_list_actual.txt
        INDEL_LIST_FILE={wildcards.sample}.indel_list_actual.txt
        CONTIG_HEADER_FRAGMENT=dna/variants/varscan2/{wildcards.sample}.contig_header_frag.txt
        FULL_HEADER_TEMP=dna/variants/varscan2/{wildcards.sample}.full_header.txt

        # 从参考基因组 .fai 文件生成 contig 头部信息片段
        awk 'BEGIN{{FS="\\t"}}; {{print "##contig=<ID="$1",length="$2">"}}' {input.ref_fai} > $CONTIG_HEADER_FRAGMENT

        # --- 处理 SNP VCFs ---
        # 过滤掉不存在的 SNP VCF 文件，只将存在的写入列表文件
        for f in {input.snp_vcfs}; do
            if [ -f "$f" ]; then
                echo "$f" >> $SNP_LIST_FILE
            fi
        done

        # 确保至少有一个 SNP VCF 文件存在以提取头部
        if [ ! -s $SNP_LIST_FILE ]; then
            echo "Error: No SNP VCF files found for {wildcards.sample}. This should not happen." >&2
            exit 1
        fi

        # 1. 使用 bcftools concat 合并所有存在的 SNP VCF 文件
        bcftools concat -a -f $SNP_LIST_FILE > {output.snp_vcf} 2> {log}
        bgzip -k {output.snp_vcf}
        tabix -p vcf {output.snp_vcf_gz}


        # --- 处理 INDEL VCFs ---
        # 过滤掉不存在的 INDEL VCF 文件，只将存在的写入列表文件
        for f in {input.indel_vcfs}; do
            if [ -f "$f" ]; then
                echo "$f" >> $INDEL_LIST_FILE
            fi
        done

        # 确保至少有一个 INDEL VCF 文件存在以提取头部
        if [ ! -s $INDEL_LIST_FILE ]; then
            echo "Error: No INDEL VCF files found for {wildcards.sample}. This should not happen." >&2
            exit 1
        fi
 
        # 1. 使用 bcftools concat 合并所有存在的 INDEL VCF 文件
        bcftools concat -a -f $INDEL_LIST_FILE > {output.indel_vcf} 2>> {log}
        bgzip -k {output.indel_vcf}
        tabix -p vcf {output.indel_vcf_gz}

        # 清理所有临时文件
        rm -f $SNP_LIST_FILE $INDEL_LIST_FILE
        """

rule varscan2_reformat:
    input:
        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf",
        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf"
    output:
        vcf = "dna/variants/varscan2/{sample}.varscan2.vcf.gz"
    params:
        prefix = "dna/variants/varscan2/{sample}",
        script_dir = workflow.basedir
    log:
        combine = "logs/varscan2/{sample}.combine.log",
        reformat = "logs/varscan2/{sample}.reformat.log"
    conda:
        "../../envs/python.yaml"
    threads: 2
    shell:
        """
        # 1. 合并SNP和INDEL
        (
          grep "^#" {input.snp_vcf} || true
          grep -v "^#" {input.snp_vcf} {input.indel_vcf} | \
            sed 's/^[^:]*://' | \
            sort -k1,1V -k2,2n || true
        ) > {params.prefix}.combined.vcf 2> {log.combine}
        
        # 2. 使用Python脚本重新格式化VCF并压缩
        python {params.script_dir}/scripts/reformat_varscan_vcf.py \
            {params.prefix}.combined.vcf \
            {params.prefix}.reformatted.vcf \
            {output.vcf} \
            2> {log.reformat}
        
        # 3. 清理临时文件
        #rm -f {params.prefix}.combined.vcf {params.prefix}.reformatted.vcf {input.snp_vcf} {input.indel_vcf}
        """

# ========================
# VarScan2 VCF Reformatting (Handles compressed input)
# ========================
#rule varscan2_reformat:
#    input:
#        # **MODIFICATION**: Input files are now the .vcf.gz from the gather step
#        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf.gz",
#        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf.gz"
#    output:
#        vcf = "dna/variants/varscan2/{sample}.varscan2.vcf.gz"
#    params:
#        prefix = "dna/variants/varscan2/{sample}",
#        script_dir = workflow.basedir
#    log:
#        combine = "logs/varscan2/{sample}.combine.log",
#        reformat = "logs/varscan2/{sample}.reformat.log"
#    conda:
#        # **MODIFICATION**: Using a dedicated environment with all tools
#        "../../envs/bcftools.yaml" 
#    threads: 2
#    shell:
#        """
#        # 1. Combine SNP and INDEL VCFs using bcftools for safe decompression
#        (
#          # Get the header from the SNP VCF file
#          bcftools view -h {input.snp_vcf} | \
#          
#          # Get the variant records (no header) from both files, decompressing on the fly
#          bcftools concat -a -H | \
#                -f < (echo -e "{input.snp_vcf}\n{input.indel_vcf}") | \
#            sort -k1,1V -k2,2n ) > {params.prefix}.combined.vcf 2> {log.combine}
#        
#        # 2. Use Python script to reformat the combined VCF and compress the final output.
#        #    This part does not need to change as it operates on the intermediate plain-text file.
#        python {params.script_dir}/scripts/reformat_varscan_vcf.py \
#            {params.prefix}.combined.vcf \
#            {params.prefix}.reformatted.vcf \
#            {output.vcf} \
#            2> {log.reformat}
#        
#        # 3. Clean up temporary files
#        rm -f {params.prefix}.combined.vcf {params.prefix}.reformatted.vcf
#        """

# ========================
#    预规范化VCF (原子化)
# ========================
#ruleorder: normalize_vcf > index_vcf
#rule normalize_vcf:
#    input:
#        vcf = "dna/variants/{caller}/{sample}.{caller}.vcf.gz",
#        vcf_idx = "dna/variants/{caller}/{sample}.{caller}.vcf.gz.tbi",
#        ref = config["reference"]["genome"]
#    output:
#        vcf = "dna/variants/{caller}/{sample}.{caller}.norm.vcf.gz",
#        vcf_idx = "dna/variants/{caller}/{sample}.{caller}.norm.vcf.gz.tbi"
#    log:
#        "logs/bcftools/{caller}/{sample}.{caller}.norm.log"
#    params:
#        # --decompose: 分解 MNP 和复杂 INDEL
#        # -m -any: 分解多等位基因位点
#        # -f ref: 左对齐和标准化
#        norm_opts = "--atomize -m -any -f"
#    conda:
#        "../../envs/bcftools.yaml"  # 确保这个环境有较新版本的bcftools
#    shell:
#        """
#        # 使用纯 bcftools 方案进行分解和标准化
#        bcftools norm {params.norm_opts} {input.ref} {input.vcf} -Oz -o {output.vcf} &> {log}
#
#        # 为新生成的VCF创建索引
#        tabix -p vcf {output.vcf}
#        """

# ========================
# 专门为ASEReadCounter准备的VCF标准化
# ========================
ruleorder: prepare_vcf_for_ase > index_vcf
rule prepare_vcf_for_ase:
    input:
        vcf = "dna/variants/mutect2/{sample}.mutect2.vcf.gz",
        vcf_idx = "dna/variants/mutect2/{sample}.mutect2.vcf.gz.tbi",
        ref = config["reference"]["genome"]
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.ase_ready.vcf.gz",
        vcf_idx = "dna/variants/mutect2/{sample}.mutect2.ase_ready.vcf.gz.tbi"
    log:
        "logs/gatk/{sample}.prepare_ase_vcf.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = 8000,
        time = 60
    shell:
        """
        # 临时文件
        TEMP_VCF1=$(mktemp --suffix='.vcf')
        TEMP_VCF2=$(mktemp --suffix='.vcf')
        TEMP_VCF3=$(mktemp --suffix='.vcf')
        
        # 1. 首先解压VCF文件
        zcat {input.vcf} > $TEMP_VCF1 2>> {log}
        
        # 2. 拆分多等位基因并左对齐
        gatk LeftAlignAndTrimVariants \
            -R {input.ref} \
            -V $TEMP_VCF1 \
            --split-multi-allelics \
            --dont-trim-alleles \
            -O $TEMP_VCF2 \
            2>> {log}
        
        # 3. 选择通过过滤的双等位基因SNV，并去除重复
        gatk SelectVariants \
            -R {input.ref} \
            -V $TEMP_VCF2 \
            --select-type-to-include SNP \
            --restrict-alleles-to BIALLELIC \
            --exclude-filtered true \
            --remove-unused-alternates true \
            --exclude-non-variants true \
            -O $TEMP_VCF3 \
            2>> {log}
        
        # 4. 去除重复位点（保留第一个遇到的）
        awk '
        BEGIN {{OFS="\t"}}
        /^#/ {{print; next}}
        {{
            pos = $1 ":" $2
            if (pos in seen) {{
                next
            }}
            seen[pos] = 1
            print
        }}
        ' $TEMP_VCF3 > $TEMP_VCF3.dedup 2>> {log}
        
        # 5. 压缩和索引
        bgzip -c $TEMP_VCF3.dedup > {output.vcf} 2>> {log}
        
        # 6. 创建索引
        gatk IndexFeatureFile -I {output.vcf} 2>> {log}
        
        # 7. 清理临时文件
        rm -f $TEMP_VCF1 $TEMP_VCF2 $TEMP_VCF3 $TEMP_VCF3.dedup
        """


#在文件顶部或规则前定义辅助函数
#def dynamic_file(caller):
#    """辅助函数：根据配置动态生成文件路径（如果caller启用）"""
#    def wrapper(wildcards):
#        if config["variant_callers"].get(caller, False):
#            return f"dna/variants/{caller}/{wildcards.sample}.{caller}.norm.vcf.gz"
#        return []  # 返回空列表表示忽略该输入
#    return wrapper


#rule merge_variant_callers:
#    input:
#        mutect2 = "dna/variants/mutect2/{sample}.mutect2.norm.vcf.gz",  # 必须存在
#        strelka2 = dynamic_file("strelka2"),  # 动态生成（如果启用）
#        varscan2 = dynamic_file("varscan2")   # 动态生成（如果启用）
#    output:
#        vcf = "dna/variants/merged/{sample}.somatic_merged.vcf.gz",
#        vcf_idx = "dna/variants/merged/{sample}.somatic_merged.vcf.gz.tbi"
#    params:
#        sample_id = "{sample}"
#    log:
#        "logs/merge_variants/{sample}.log"
#    conda:
#        "../../envs/bcftools.yaml"
#    script:
#        "../scripts/merge_variants.py"
