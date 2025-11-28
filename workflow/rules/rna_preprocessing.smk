# RNA-seq预处理规则
rule fastp_qc_trim_rna:
    input:
        r1 = get_rna_fq1, # 可以直接传递函数名
        r2 = get_rna_fq2
    output:
        r1 = "rna/trimmed/{sample}_{sampletype}_R1_trimmed.fastq.gz", # 标记为temp，如果下游只用BAM，它会自动清理
        r2 = "rna/trimmed/{sample}_{sampletype}_R2_trimmed.fastq.gz",
        json = "qc/fastp/{sample}_{sampletype}_rna_fastp.json",
        html = "qc/fastp/{sample}_{sampletype}_rna_fastp.html"
    log:
        "logs/fastp/{sample}_{sampletype}.log" # 修正日志名
    params:
        extra = "--detect_adapter_for_pe --trim_poly_g --correction"
    conda:
        "../../envs/fastp.yaml" # 假设在 workflow/rules 目录下
    threads: 4
    resources:
        mem_mb = 8000,
        time = 60
    shell:
        """
        mkdir -p rna/trimmed qc/fastp # 确保目录存在
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --report_title "{wildcards.sample}_{wildcards.sampletype}" \
            --thread {threads} \
            {params.extra} \
            > {log} 2>&1
        """

# STAR比对RNA-seq数据到参考基因组
# ---------------------------------------------------
# 步骤 1: STAR 比对 (增强版)
# - 直接输出按坐标排序的BAM
# - 添加Read Group信息
# - 生成用于融合检测的 Chimeric.out.junction 文件
# ---------------------------------------------------
#rule star_align:
#    input:
#        r1 = "rna/trimmed/{sample}_{sampletype}_R1_trimmed.fastq.gz",
#        r2 = "rna/trimmed/{sample}_{sampletype}_R2_trimmed.fastq.gz",
#        index = config["reference"]["star_index"]
#    output:
#        # 主输出，为后续GATK处理准备
#        bam = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam",
#        # 融合分析所需文件
#        chimeric_junctions = "rna/align/{sample}_{sampletype}/Chimeric.out.junction",
#        # QC日志
#        log_final = "rna/align/{sample}_{sampletype}/Log.final.out"
#    log:
#        "logs/align/{sample}.{sampletype}.log"
#    params:
#        # 输出文件的前缀，STAR会在此目录下生成所有文件
#        prefix = "rna/align/{sample}_{sampletype}/",
#        # 为BAM文件添加Read Group信息
#        read_group = lambda wildcards: f"ID:{wildcards.sample}_{wildcards.sampletype} PL:ILLUMINA SM:{wildcards.sample}_{wildcards.sampletype} LB:{wildcards.sample}"
#    threads: 16
#    resources:
#        mem_mb = 50000, # STAR内存需求较大
#        time = 240
#    conda:
#        "../../envs/star.yaml"
#    shell:
#        """
#        STAR \
#            --runThreadN {threads} \
#            --genomeDir {input.index} \
#            --readFilesIn {input.r1} {input.r2} \
#            --readFilesCommand zcat \
#            --outFileNamePrefix {params.prefix} \
#            --outSAMtype BAM SortedByCoordinate \
#            --outSAMattributes NH HI AS nM NM MD \
#            --outSAMunmapped Within \
#            --outFilterType BySJout \
#            --outFilterMultimapNmax 20 \
#            --alignSJoverhangMin 8 \
#            --alignSJDBoverhangMin 1 \
#            --outFilterMismatchNmax 999 \
#            --outFilterMismatchNoverLmax 0.04 \
#            --alignIntronMin 20 \
#            --alignIntronMax 1000000 \
#            --alignMatesGapMax 1000000 \
#            --chimSegmentMin 15 \
#            --chimJunctionOverhangMin 15 \
#            --outSAMattrRGline {params.read_group} \
#            > {log} 2>&1
#        """

#rule star_align:
#    input:
#        r1 = "rna/trimmed/{sample}_{sampletype}_R1_trimmed.fastq.gz",
#        r2 = "rna/trimmed/{sample}_{sampletype}_R2_trimmed.fastq.gz",
#        index = config["reference"]["star_index"]
#    output:
#        bam = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam",
#        counts = "rna/align/{sample}_{sampletype}/ReadsPerGene.out.tab",
#        chimeric_junctions = "rna/align/{sample}_{sampletype}/Chimeric.out.junction", # Arriba's key input
#        log_final = "rna/align/{sample}_{sampletype}/Log.final.out"
#    log:
#        "logs/align/{sample}.{sampletype}.log"
#    params:
#        prefix = "rna/align/{sample}_{sampletype}/",
#        read_group = lambda wildcards: f"ID:{wildcards.sample}_{wildcards.sampletype}\\tPL:ILLUMINA\\tSM:{wildcards.sample}_{wildcards.sampletype}\\tLB:{wildcards.sample}"
#    threads: 16
#    resources:
#        mem_mb = 50000,
#        time = 480
#    conda:
#        "../../envs/star.yaml"
#    shell:
#        """
#        # STAR参数说明:
#        # - 通用输出与GATK要求: BAM排序输出, 保留未比对reads, 完整SAM属性, 链信息, ReadGroup
#        # - 定量分析: GeneCounts模式
#        # - Arriba优化: 嵌合体检测相关参数
#        mkdir -p {params.prefix}
#        STAR \
#            --#runThreadN {threads} \
#            --genomeDir {input.index} \
#           --readFilesIn {input.r1} {input.r2} \
#            --readFilesCommand zcat \
#            --outFileNamePrefix {params.prefix} \
#            --outSAMtype BAM SortedByCoordinate \
#            --outSAMunmapped Within \
#            --outSAMattributes All \
#            --outSAMstrandField intronMotif \
#            --outSAMattrRGline '{params.read_group}' \
#            --twopassMode Basic \
#            --quantMode GeneCounts \
#            --outSJfilterReads All \
#            --outFilterMultimapNmax 50 \
#            --peOverlapNbasesMin 10 \
#            --alignSplicedMateMapLminOverLmate 0.5 \
#            --alignSJstitchMismatchNmax 5 -1 5 5 \
#            --chimSegmentMin 10 \
#            --chimJunctionOverhangMin 10 \
#            --chimScoreDropMax 30 \
#            --chimScoreJunctionNonGTAG 0 \
#            --chimScoreSeparation 1 \
#            --chimSegmentReadGapMax 3 \
#            --chimMultimapNmax 50 \
#            --chimOutType Junctions WithinBAM HardClip \
#            --outFilterType BySJout \
#            --alignSJoverhangMin 8 \
#            --alignSJDBoverhangMin 1 \
#            --outFilterMismatchNmax 999 \
#            --outFilterMismatchNoverLmax 0.04 \
#            --alignIntronMin 20 \
#            --alignIntronMax 1000000 \
#            --alignMatesGapMax 1000000 \
#            > {log} 2>&1
#        """
# ---------------------------------------------------
# RNA-seq比对 (STAR)
# ---------------------------------------------------
rule star_align:
    input:
        r1 = "rna/trimmed/{sample}_{sampletype}_R1_trimmed.fastq.gz",
        r2 = "rna/trimmed/{sample}_{sampletype}_R2_trimmed.fastq.gz",
        index = config["reference"]["star_index"]
    output:
        bam = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam",
        counts = "rna/align/{sample}_{sampletype}/ReadsPerGene.out.tab",
        chimeric_junctions = "rna/align/{sample}_{sampletype}/Chimeric.out.junction",
        log_final = "rna/align/{sample}_{sampletype}/Log.final.out"
    log:
        "logs/align/{sample}.{sampletype}.log"
    params:
        prefix = "rna/align/{sample}_{sampletype}/",
        # --- THE DEFINITIVE FIX for @RG ---
        # Provide RG info as separate key-value pairs. STAR will format the line correctly.
        rg_id = "{sample}_{sampletype}",
        rg_sm = "{sample}_{sampletype}",
        rg_pl = "ILLUMINA",
        rg_lb = "{sample}"
    threads: 16
    conda:
        "../../envs/star.yaml"
    shell:
        """
        mkdir -p {params.prefix}
        
        STAR \\
            --runThreadN {threads} \\
            --genomeDir {input.index} \\
            --readFilesIn {input.r1} {input.r2} \\
            --readFilesCommand zcat \\
            --outFileNamePrefix {params.prefix} \\
            --outSAMtype BAM SortedByCoordinate \\
            --outSAMunmapped Within \\
            --outSAMattributes All \\
            --outSAMstrandField intronMotif \\
            --outSAMattrRGline ID:{params.rg_id} PL:{params.rg_pl} SM:{params.rg_sm} LB:{params.rg_lb} \\
            --twopassMode Basic \\
            --quantMode GeneCounts \\
            --chimSegmentMin 10 \\
            --chimJunctionOverhangMin 10 \\
            --chimOutType Junctions WithinBAM \\
            --outFilterType BySJout \\
            --alignSJoverhangMin 8 \\
            --alignSJDBoverhangMin 1 \\
            --outFilterMismatchNmax 999 \\
            --outFilterMismatchNoverLmax 0.04 \\
            --alignIntronMin 20 \\
            --alignIntronMax 1000000 \\
            --alignMatesGapMax 1000000 \\
            > {log} 2>&1
        """


# ---------------------------------------------------
# 步骤 2: 标记重复 (GATK最佳实践)
# ---------------------------------------------------
rule mark_duplicates_rna:
    input:
        "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam"
    output:
        bam = "rna/align/{sample}_{sampletype}_dedup.bam",
        metrics = "qc/picard/{sample}_{sampletype}_markdup_metrics.txt"
    log:
        "logs/picard/{sample}_{sampletype}_markdup.log"
    conda:
        "../../envs/picard.yaml"
    shell:
        """
        mkdir -p qc/picard
        picard MarkDuplicates \\
            -I {input} \\
            -O {output.bam} \\
            -M {output.metrics} \\
            --CREATE_INDEX true \\
            --VALIDATION_STRINGENCY SILENT \\
            2> {log}
        """

# RNA-seq质量控制指标
rule rna_qc:
    input:
        bam = "rna/align/{sample}_{sampletype}_dedup.bam",
        refFlat = config["reference"]["refflat"]
        #refgenes = config["reference"]["refgenes"]
    output:
        "qc/rna/{sample}_{sampletype}_RNA_metrics.txt"
    log:
        "logs/picard/{sample}_{sampletype}_metrics.log"
    conda:
        "../../envs/picard.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    shell:
        """
        mkdir -p qc/rna
        picard CollectRnaSeqMetrics \
            -I {input.bam} \
            -O {output} \
            -REF_FLAT {input.refFlat} \
            -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
            2> {log}
        """