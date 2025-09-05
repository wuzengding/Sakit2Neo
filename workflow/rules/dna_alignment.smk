# BWA MEM比对DNA序列
rule bwa_mem:
    input:
        r1 = "dna/trimmed/{sample}_{sampletype}_R1.trimmed.fastq.gz",
        r2 = "dna/trimmed/{sample}_{sampletype}_R2.trimmed.fastq.gz",
        ref = config["reference"]["genome"]
    output:
        #temp("dna/aligned/{sample}_{sampletype}.raw.bam")
        "dna/aligned/{sample}_{sampletype}.raw.bam"
    log:
        "logs/bwa/{sample}_{sampletype}.log"
    params:
        read_group = lambda wildcards: f"@RG\\tID:{wildcards.sample}_{wildcards.sampletype}\\tSM:{wildcards.sample}_{wildcards.sampletype}\\tPL:ILLUMINA\\tLB:lib1\\tPU:unit1"
    threads: 40
    conda:
        "../../envs/bwa.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        time = lambda wildcards, attempt: attempt * 240
    shell:
        """
        mkdir -p dna/aligned
        bwa mem -M -t {threads} \
            -R '{params.read_group}' \
            {input.ref} {input.r1} {input.r2} 2> {log} | \
            samtools view -Sb - > {output}
        """

# 比对后进行排序和去重
rule sort_and_fix_tags:
    input:
        "dna/aligned/{sample}_{sampletype}.raw.bam"
    output:
        #temp("dna/aligned/{sample}_{sampletype}.sorted.bam")
        "dna/aligned/{sample}_{sampletype}.sorted.bam"
    log:
        "logs/samtools/{sample}_{sampletype}.sort.log"
    conda:
        "../../envs/bwa.yaml"
    threads: 20
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 180
    shell:
        """
        samtools sort -@ {threads} -O bam -o {output} {input} 2> {log}
        """

# 标记重复序列
rule mark_duplicates:
    input:
        "dna/aligned/{sample}_{sampletype}.sorted.bam"
    output:
        bam = temp("dna/aligned/{sample}_{sampletype}.markdup.bam"),
        #"dna/aligned/{sample}_{sampletype}.markdup.bam",
        metrics = "qc/picard/{sample}_{sampletype}.markdup_metrics.txt"
    log:
        "logs/dedup/{sample}_{sampletype}.markdup.log"
    params:
        tmpdir = temp(directory("dna/aligned/temp/{sample}_{sampletype}_markdup"))
        #extra = "--REMOVE_DUPLICATES false --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT"
    conda:
        "../../envs/sambamba.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        time = lambda wildcards, attempt: attempt * 180
    threads: 4
    shell:
        """
        mkdir -p qc/dedup
        # 使用 sambamba markdup
        # -t: 线程数
        # --tmpdir: 指定临时目录
        # --overflow-list-size: 增加这个值可以处理高深度区域，减少I/O
        # -l 0: 使用无压缩的BAM输出，因为是临时文件，速度更快
        # 输出的 metrics 文件格式与 Picard 兼容
        sambamba markdup -t {threads} --tmpdir={params.tmpdir} \\
            --overflow-list-size 1000000 \\
            -l 0 \\
            {input} {output.bam} 2> {output.metrics}

        # sambamba 的日志(stderr)就是metrics文件，所以我们将它重定向。
        # 如果需要单独的日志，可以这样写: ... {output.bam} 2> {log} 1> {output.metrics}
        # 但通常它的metrics就是最好的日志。
        # 为了与Snakemake的日志管理保持一致，我们复制一份
        cp {output.metrics} {log}
        """

## 碱基质量分数重校准 - 创建重校准表
#rule base_recalibrator:
#    input:
#        bam = "dna/aligned/{sample}_{sampletype}.markdup.bam",
#        ref = config["reference"]["genome"],
#        known_sites1 = config["reference"]["dbsnp"],
#        known_sites2 = config["reference"]["known_indels"]
#    output:
#        recal_table = "dna/aligned/{sample}_{sampletype}.recal_data.table"
#    log:
#        "logs/gatk/{sample}_{sampletype}.bqsr.log"
#    params:
#        extra = "--use-original-qualities"
#    conda:
#        "../../envs/gatk.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * 16000,
#        time = lambda wildcards, attempt: attempt * 120
#    shell:
#        """
#        gatk BaseRecalibrator \
#            -R {input.ref} \
#            -I {input.bam} \
#            --known-sites {input.known_sites1} \
#            --known-sites {input.known_sites2} \
#            -O {output.recal_table} \
#            {params.extra} \
#            2> {log}
#        """

# -----------------------------------------------
#  碱基质量分数重校准 (BQSR) - 最终目录分离版
# -----------------------------------------------

# 步骤 1: 分散 (Scatter) - 按染色体并行计算重校准表
rule bqsr_recalibrator_scatter:
    input:
        bam = "dna/aligned/{sample}_{sampletype}.markdup.bam",
        bai = "dna/aligned/{sample}_{sampletype}.markdup.bam.bai",
        ref = config["reference"]["genome"],
        known_sites1 = config["reference"]["dbsnp"],
        known_sites2 = config["reference"]["known_indels"]
    output:
        # **核心修复**: 每个通配符都用目录分隔
        recal_table = temp("dna/aligned/temp/{sample}_{sampletype}_{chromosome}.recal_data.table")
    log:
        "logs/gatk/{sample}_{sampletype}_{chromosome}.bqsr.log"
    params:
        extra = "--use-original-qualities",
        # **重要**: 我们在shell中使用的染色体名需要带'chr'
        chromosome_name = "chr{chromosome}"
    conda:
        "../../envs/gatk.yaml"
    resources: mem_mb = 8000, time = 60
    shell:
        "mkdir -p $(dirname {output.recal_table}); "
        "gatk BaseRecalibrator -R {input.ref} -I {input.bam} "
        "--known-sites {input.known_sites1} --known-sites {input.known_sites2} "
        "-L {params.chromosome_name} -O {output.recal_table} {params.extra} 2> {log}"

# 步骤 2: 聚集 (Gather) - 合并所有染色体的重校准表
rule bqsr_gather_reports:
    input:
        recal_tables = expand(
            # **核心修复**: 匹配新的、完全无歧义的目录结构
            "dna/aligned/temp/{{sample}}_{{sampletype}}_{chromosome}.recal_data.table",
            chromosome=[c.replace("chr", "") for c in get_chromosomes()]
        )
    output:
        recal_table = "dna/aligned/{sample}_{sampletype}.recal_data.table"
    log:
        "logs/gatk/{sample}_{sampletype}.gather_bqsr.log"
    conda:
        "../../envs/gatk.yaml"
    resources: mem_mb = 4000, time = 30
    shell:
        """
        INPUT_ARGS=$(for table in {input.recal_tables}; do echo -n "-I $table "; done)
        gatk GatherBQSRReports $INPUT_ARGS -O {output.recal_table} 2> {log}
        """



## 碱基质量分数重校准 - 应用重校准表
#rule apply_bqsr:
#    input:
#        bam = "dna/aligned/{sample}_{sampletype}.markdup.bam",
#        ref = config["reference"]["genome"],
#        recal_table = "dna/aligned/{sample}_{sampletype}.recal_data.table"
#    output:
#        bam = "dna/aligned/{sample}_{sampletype}.recal.bam",
#        #bai = "dna/aligned/{sample}_{sampletype}.recal.bai"
#    log:
#        "logs/gatk/{sample}_{sampletype}.apply_bqsr.log"
#    params:
#        extra = "--use-original-qualities"
#    conda:
#        "../../envs/gatk.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: attempt * 16000,
#        time = lambda wildcards, attempt: attempt * 120
#    shell:
#        """
#        gatk ApplyBQSR \
#            -R {input.ref} \
#            -I {input.bam} \
#            --bqsr-recal-file {input.recal_table} \
#            -O {output.bam} \
#            {params.extra} \
#            2> {log}
#      """

# -----------------------------------------------
#  应用碱基质量分数重校准 (ApplyBQSR) - 最终目录分离版
# -----------------------------------------------

# 步骤 1: 分散 (Scatter) - 按染色体并行应用重校准表
rule apply_bqsr_scatter:
    input:
        bam = "dna/aligned/{sample}_{sampletype}.markdup.bam",
        bai = "dna/aligned/{sample}_{sampletype}.markdup.bam.bai",
        ref = config["reference"]["genome"],
        recal_table = "dna/aligned/{sample}_{sampletype}.recal_data.table"
    output:
        # **核心修复**: 每个通配符都用目录分隔
        bam = temp("dna/aligned/temp/{sample}_{sampletype}_{chromosome}.recal.bam")
    log:
        "logs/gatk/{sample}_{sampletype}_{chromosome}.apply_bqsr.log"
    params:
        extra = "--use-original-qualities",
        chromosome_name = "chr{chromosome}"
    conda:
        "../../envs/gatk.yaml"
    resources: mem_mb = 8000, time = 60
    shell:
        "mkdir -p $(dirname {output.bam}); "
        "gatk ApplyBQSR -R {input.ref} -I {input.bam} "
        "--bqsr-recal-file {input.recal_table} "
        "-L {params.chromosome_name} -O {output.bam} {params.extra} 2> {log}"

# 步骤 2: 聚集 (Gather) - 合并所有染色体的BAM文件
rule gather_recal_bams:
    input:
        bams = expand(
            # **核心修复**: 匹配新的、完全无歧义的目录结构
            "dna/aligned/temp/{{sample}}_{{sampletype}}_{chromosome}.recal.bam",
            chromosome=[c.replace("chr", "") for c in get_chromosomes()]
        )
    output:
        bam = "dna/aligned/{sample}_{sampletype}.recal.bam"
        #bai = "dna/aligned/{sample}_{sampletype}.recal.bam.bai"
    log:
        "logs/picard/{sample}_{sampletype}.gather_bams.log"
    conda:
        "../../envs/picard.yaml"
    resources: mem_mb = 8000, time = 60
    shell:
        """
        INPUT_ARGS=$(for bam in {input.bams}; do echo -n "-I $bam "; done)
        picard GatherBamFiles $INPUT_ARGS -O {output.bam} --CREATE_INDEX true 2> {log}
        """


# 创建BAM索引
rule index_bam:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/samtools/index_{prefix}.log"
    conda:
        "../../envs/picard.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = lambda wildcards, attempt: attempt * 60
    threads: 20
    shell:
        "samtools index  -@ {threads} {input} 2> {log}"

# 计算覆盖度统计
rule coverage_stats:
    input:
        bam = "dna/aligned/{sample}_{sampletype}.recal.bam",
        bai = "dna/aligned/{sample}_{sampletype}.recal.bam.bai",
        target_bed = config["reference"]["capture_kit"]
    output:
        stats = "qc/coverage/{sample}_{sampletype}.mosdepth.global.dist.txt",
        hist = "qc/coverage/{sample}_{sampletype}.regions.bed.gz"
    log:
        "logs/bedtools/{sample}_{sampletype}.coverage.log"
    params:
        prefix = "qc/coverage/{sample}_{sampletype}"
    threads: 4
    conda:
        "../../envs/mosdepth.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 30
    shell:
        """
        mkdir -p qc/coverage
        mosdepth -t {threads} -b {input.target_bed} -n {params.prefix} {input.bam} 2> {log}
        """

# 计算插入片段大小分布
rule insert_size_metrics:
    input:
        bam = "dna/aligned/{sample}_{sampletype}.recal.bam"
    output:
        metrics = "qc/picard/{sample}_{sampletype}.insert_size_metrics.txt",
        pdf = "qc/picard/{sample}_{sampletype}.insert_size_histogram.pdf"
    log:
        "logs/picard/{sample}_{sampletype}.insert_size.log"
    conda:
        "../../envs/picard.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    shell:
        """
        mkdir -p qc/picard
        picard CollectInsertSizeMetrics \
            -I {input.bam} \
            -O {output.metrics} \
            -H {output.pdf} \
            -M 0.5 \
            2> {log}
        """

# 生成比对质量指标
rule alignment_metrics:
    input:
        bam = "dna/aligned/{sample}_{sampletype}.recal.bam",
        ref = config["reference"]["genome"]
    output:
        "qc/picard/{sample}_{sampletype}.alignment_metrics.txt"
    log:
        "logs/picard/{sample}_{sampletype}.alignment_metrics.log"
    conda:
        "../../envs/picard.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    shell:
        """
        mkdir -p qc/picard
        picard CollectAlignmentSummaryMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            2> {log}
        """

# 计算靶向区域质量指标
rule hs_metrics:
    input:
        bam = "dna/aligned/{sample}_{sampletype}.sorted.bam",
        ref = config["reference"]["genome"],
        bait = config["reference"]["interval_list"],
        target = config["reference"]["interval_list"]
    output:
        "qc/picard/{sample}_{sampletype}.hs_metrics.txt"
    log:
        "logs/picard/{sample}_{sampletype}.hs_metrics.log"
    conda:
        "../../envs/picard.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    shell:
        """
        mkdir -p qc/picard
        picard CollectHsMetrics \
            -R {input.ref} \
            -I {input.bam} \
            -O {output} \
            -BAIT_INTERVALS {input.bait} \
            -TARGET_INTERVALS {input.target} \
            2> {log}
        """
