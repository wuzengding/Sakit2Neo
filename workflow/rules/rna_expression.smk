# RNA-seq表达定量规则

# 使用featureCounts计数基因表达水平
rule feature_counts:
    input:
        bam = "rna/mapped/{sample}.rna.markdup.bam",
        gtf = config["reference"]["gtf"]
    output:
        counts = "rna/expression/{sample}.rna.counts.tsv",
        summary = "rna/expression/{sample}.rna.counts.tsv.summary"
    log:
        "logs/featurecounts/{sample}.rna.log"
    threads: 8
    conda:
        "../envs/subread.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 120
    shell:
        """
        mkdir -p rna/expression
        featureCounts \
            -T {threads} \
            -a {input.gtf} \
            -o {output.counts} \
            -p \
            -t exon \
            -g gene_id \
            {input.bam} \
            > {log} 2>&1
        """

# 使用Salmon进行转录本定量
rule salmon_quant:
    input:
        r1 = "rna/trimmed/{sample}.rna.R1.trimmed.fastq.gz",
        r2 = "rna/trimmed/{sample}.rna.R2.trimmed.fastq.gz",
        index = config["reference"]["salmon_index"]
    output:
        quant = directory("rna/salmon/{sample}_rna_quant")
    log:
        "logs/salmon/{sample}_rna.log"
    params:
        libtype = "A",
        extra = "--validateMappings --gcBias --seqBias"
    threads: 8
    conda:
        "../envs/salmon.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 120
    shell:
        """
        mkdir -p rna/salmon
        salmon quant \
            -i {input.index} \
            -l {params.libtype} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} \
            -o {output.quant} \
            {params.extra} \
            > {log} 2>&1
        """

# 合并所有样本的counts矩阵
rule merge_counts:
    input:
        counts = lambda wildcards: expand("rna/expression/{sample}_rna.counts.tsv", 
                                        sample=samples["sample_id"].unique(), 
                                        )
    output:
        matrix = "rna/expression/all_samples_rna.counts.matrix.tsv"
    log:
        "logs/merge_counts/all_samples_rna.log"
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = lambda wildcards, attempt: attempt * 60
    script:
        "../scripts/merge_counts.R"

# 合并所有样本的Salmon定量结果
rule merge_salmon:
    input:
        quants = lambda wildcards: expand("rna/salmon/{sample}_rna_quant", 
                                         sample=samples["sample_id"].unique()
                                         )
    output:
        matrix = "rna/salmon/all_samples_rna.salmon.matrix.tsv"
    log:
        "logs/merge_salmon/all_samples_rna.log"
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = lambda wildcards, attempt: attempt * 60
    script:
        "../scripts/merge_salmon.R"

# 基因表达差异分析
rule deseq2_analysis:
    input:
        counts = "rna/expression/all_samples_tumor.counts.matrix.tsv",
        samples = config["samples"],
        annotation = config["reference"]["gene_annotation"]
    output:
        results = "rna/differential_expression/deseq2_results.tsv",
        plots = "reports/differential_expression/deseq2_plots.pdf"
    log:
        "logs/deseq2/analysis.log"
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 180
    script:
        "../scripts/deseq2_analysis.R"

# 转录本表达差异分析
rule dexseq_analysis:
    input:
        quants = "rna/salmon/all_samples_tumor.salmon.matrix.tsv",
        samples = config["samples"],
        tx2gene = config["reference"]["tx2gene"]
    output:
        results = "rna/differential_expression/dexseq_results.tsv",
        plots = "reports/differential_expression/dexseq_plots.pdf"
    log:
        "logs/dexseq/analysis.log"
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 180
    script:
        "../scripts/dexseq_analysis.R"

# 创建表达量热图
rule expression_heatmap:
    input:
        counts = "rna/expression/all_samples_tumor.counts.matrix.tsv",
        deseq2 = "rna/differential_expression/deseq2_results.tsv"
    output:
        heatmap = "reports/plots/expression_heatmap.pdf",
        pca = "reports/plots/expression_pca.pdf"
    log:
        "logs/plots/expression_plots.log"
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    script:
        "../scripts/expression_plots.R"

