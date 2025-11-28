# ---------------------------------------------------
#  处理剪接位点 (GATK最佳实践 for RNA)
# ---------------------------------------------------
rule split_n_cigar_reads:
    input:
        bam = "rna/align/{sample}_{sampletype}_dedup.bam",
        ref = config["reference"]["genome"]
    output:
        # 这是最终用于变异验证和表达定量的BAM文件
        bam = "rna/snv_validation/{sample}_{sampletype}_analysis_ready.bam",
    log:
        "logs/gatk/{sample}_{sampletype}_splitn.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = 16000,
        time = 180
    shell:
        """
        gatk SplitNCigarReads \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.bam} \
            2> {log}
        """

# ================================================
#   使用 GATK ASEReadCounter 验证DNA层面的SNV
# ================================================
rule rna_variant_validation:
    """
    Counts reads supporting REF and ALT alleles at given SNP sites
    in an RNA-seq BAM file. This is used to validate the expression
    of somatic variants found in DNA.
    """
    input:
        # 输入是经过 SplitNCigarReads 处理后的RNA BAM文件
        bam = "rna/snv_validation/{sample}_tumor_analysis_ready.bam",
        #bai = "rna/snv_validation/{sample}_tumor_analysis_ready.bai",
    
        # 输入是DNA层面检测到的、需要验证的SNV位点 (VCF格式)
        sites_vcf = "dna/variants/mutect2/{sample}.mutect2.ase_ready.vcf.gz",
        sites_vcf_idx = "dna/variants/mutect2/{sample}.mutect2.ase_ready.vcf.gz.tbi",
        
        # ASEReadCounter 需要参考基因组
        ref = config["reference"]["genome"]
    output:
        # 输出是一个表格文件，每行是一个位点，包含REF和ALT的read计数
        counts_table = "rna/snv_validation/{sample}_tumor_mutect2_ase_counts.tsv"
    log:
        "logs/gatk/{sample}_tumor_mutect2_ASEReadCounter.log"
    conda:
        # ASEReadCounter 是 GATK 的一部分
        "../../envs/gatk.yaml"
    resources:
        mem_mb = 8000
    shell:
        """
        gatk ASEReadCounter \\
            -R {input.ref} \\
            -I {input.bam} \\
            -V {input.sites_vcf} \\
            -O {output.counts_table} \\
            --verbosity INFO \\
            --output-format TABLE \\
            --min-mapping-quality 20 \\
            --min-base-quality 10 \\
            2> {log}
        """