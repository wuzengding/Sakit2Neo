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
        vcf_idx = "dna/variants/mutect2/scatter/{sample}.{chromosome}.withpon.vcf.gz.tbi",
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
        vcf_raw = temp("dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz"),
        vcf_idx_raw = temp("dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz.tbi"),
        stats_raw = temp("dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz.stats"),
        f1r2 = temp("dna/variants/mutect2/scatter/{sample}.{chromosome}.f1r2.tar.gz")
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

#rule get_pileup_summaries_normal:
#    input:
#        tumor_bam = get_normal_bam_scattered,
#        tumor_bai = get_normal_bam_scattered.replace(".bam", ".bam.bai"),
#        intervals = config["reference"]["interval_list"],
#        sites = config["reference"]["gnomad"]
#    output:
#        pileup = "dna/variants/mutect2/{sample}_normal.pileups.table"
#    log:
#        "logs/mutect2/{sample}.normal_pileup.log"
#    conda:
#        "../../envs/gatk.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * 8000,
#        time = lambda wildcards, attempt: attempt * 60,
#        threads = 2
#    shell:
#        """
#        gatk GetPileupSummaries \
#            -I {input.bam} \
#            -V {input.sites} \
#            -L {input.intervals} \
#            -O {output.pileup} \
#            2> {log}
#        """

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
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = lambda wildcards, attempt: attempt * 30,
        threads = 1
    shell:
        """
        TUMOR_ARGS=$(for pileup in {input.tumor_pileups}; do echo -n "-I $pileup "; done)
        NORMAL_ARGS=$(for pileup in {input.normal_pileups}; do echo -n "-matched $pileup "; done)
        gatk CalculateContamination \
             $TUMOR_ARGS \
             $NORMAL_ARGS \
            -O {output.contamination} \
            --tumor-segmentation \
            {output.segments} \
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
        stats = "dna/variants/mutect2/scatter/{sample}.{chromosome}.raw.vcf.gz.stats",
        contamination = "dna/variants/mutect2/{sample}.contamination.table",
        segments = "dna/variants/mutect2/{sample}.segments.table",
        orientation_model = "dna/variants/mutect2/{sample}.read-orientation-model.tar.gz",
        ref = config["reference"]["genome"]
    output:
        vcf = temp("dna/variants/mutect2/scatter/{sample}.{chromosome}.nopon.filtered.vcf.gz"),
        filtering_stats = temp("dna/variants/mutect2/scatter/{sample}.{chromosome}.nopon.filtering.stats")
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
               chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.withpon.vcf.gz",
        vcf_idx = "dna/variants/mutect2/{sample}.mutect2.withpon.vcf.gz.tbi"
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
               chromosome=[c.replace("chr", "") for c in get_chromosomes()])
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.nopon.vcf.gz",
        vcf_idx = "dna/variants/mutect2/{sample}.mutect2.nopon.vcf.gz.tbi"
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
rule mutect2_caller:
    input:
        unpack(get_mutect2_input)
    output:
        vcf = "dna/variants/mutect2/{sample}.mutect2.vcf.gz",
        vcf_idx = "dna/variants/mutect2/{sample}.mutect2.vcf.gz.tbi",
    log:
        "logs/mutect2/{sample}.finalize.log"
    shell:
        """
        # 直接链接或复制最终文件
        ln -sf $(basename {input.vcf}) {output.vcf} 2> {log}
        ln -sf $(basename {input.vcf_idx}) {output.vcf_idx} 2>> {log}
        """

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

rule varscan2_call:
    input:
        tumor_bam = get_tumor_bam,
        normal_bam = get_normal_bam,
        tumor_bai = get_tumor_bai,
        normal_bai = get_normal_bai,
        ref = config["reference"]["genome"],
        intervals = config["reference"]["capture_kit"]
    output:
        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf",
        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf"
    params:
        prefix = "dna/variants/varscan2/{sample}",
        min_coverage = config["varscan2"]["min_coverage"],
        min_var_freq = config["varscan2"]["min_var_freq"],
        p_value = config["varscan2"]["p_value"],
        somatic_p = config["varscan2"]["somatic_p_value"],
        tumor_name = "{sample}_tumor",
        normal_name = "{sample}_normal"
    log:
        mpileup = "logs/varscan2/{sample}.mpileup.log",
        varscan = "logs/varscan2/{sample}.varscan.log"
    conda:
        "../../envs/varscan.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 240,
        threads = 4
    shell:
        """
        # 1. 生成mpileup文件
        samtools mpileup \
            -f {input.ref} \
            -l {input.intervals} \
            -q 1 -Q 13 \
            {input.normal_bam} {input.tumor_bam} \
            > {params.prefix}.mpileup \
            2> {log.mpileup}
        
        # 2. VarScan2体细胞调用
        varscan somatic \
            {params.prefix}.mpileup \
            {params.prefix} \
            --mpileup 1 \
            --min-coverage {params.min_coverage} \
            --min-var-freq {params.min_var_freq} \
            --p-value {params.p_value} \
            --somatic-p-value {params.somatic_p} \
            --strand-filter 1 \
            --output-vcf 1 \
            2> {log.varscan}
        
        # 3. 清理中间文件
        #rm -f {params.prefix}.mpileup
        """
rule varscan2_reformat:
    input:
        snp_vcf = "dna/variants/varscan2/{sample}.snp.vcf",
        indel_vcf = "dna/variants/varscan2/{sample}.indel.vcf",
    output:
        vcf = "dna/variants/varscan2/{sample}.varscan2.vcf.gz",
        vcf_idx = "dna/variants/varscan2/{sample}.varscan2.vcf.gz.tbi"
    params:
        prefix = "dna/variants/varscan2/{sample}",
        script_dir = workflow.basedir
    log:
        combine = "logs/varscan2/{sample}.combine.log",
        reformat = "logs/varscan2/{sample}.reformat.log"
    conda:
        "../../envs/python.yaml"  # 或者创建一个专门用于格式化的轻量级环境
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
#    预规范化VCF (原子化)
# ========================
#rule normalize_vcf:
#    """
#    A general and robust rule to normalize a VCF file from any source.
#    It performs a two-step normalization (split then merge) to guarantee
#    that the output VCF has exactly one canonical record per genomic position,
#    making it compatible with GATK tools.
#    """
#    input:
#        vcf = "dna/variants/{caller}/{sample}.{caller}.vcf.gz",
#        ref = config["reference"]["genome"]
#    output:
#        # **重要**: 这个输出现在是最终的、可用于所有下游分析的规范化文件
#        vcf = "dna/variants/{caller}/{sample}.{caller}.norm.vcf.gz",
#        vcf_idx = "dna/variants/{caller}/{sample}.{caller}.norm.vcf.gz.tbi"
#    log:
#        "logs/bcftools/{caller}/{sample}.{caller}.norm.log"
#    conda:
#        "../../envs/bcftools.yaml"
#    shell:
#        """
#       # 步骤 1: 原子化 - 将所有变异拆分为最简单的双等位基因形式
#        # 使用 bcftools norm -m - 来拆分多等位基因位点，并进行左对齐
#        # 输出到一个临时的管道
#        bcftools norm -m - -f {input.ref} {input.vcf} -O z | \\
#        
#        # 步骤 2: 重组 - 将同一位置的所有简单变异合并回单行
#        # 使用 bcftools norm -m +any 从管道读取，并输出最终文件
#        bcftools norm -m +any -f {input.ref} - -O z -o {output.vcf} > {log} 2>&1
#        
#        # 为最终的规范化文件创建索引
#        tabix -p vcf {output.vcf}
#        """

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
