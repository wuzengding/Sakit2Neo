# 生殖细胞变异检测

# GATK HaplotypeCaller检测变异
rule haplotype_caller:
    input:
        bam = "dna/aligned/{sample}_normal.recal.bam",
        bai = "dna/aligned/{sample}_normal.recal.bam.bai",
        ref = config["reference"]["genome"],
        interval = config["reference"]["capture_kit"]
    output:
        vcf = "dna/variants/germline/{sample}.germline.vcf.gz",
        tbi = "dna/variants/germline/{sample}.germline.vcf.gz.tbi"
    log:
        "logs/gatk/{sample}.haplotypecaller.log"
    params:
        extra = "--standard-min-confidence-threshold-for-calling 20.0"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 180
    shell:
        """
        mkdir -p dna/variants/germline
        gatk HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -L {input.interval} \
            -O {output.vcf} \
            {params.extra} \
            2> {log}
        """

# 过滤生殖细胞变异
rule filter_germline_variants:
    input:
        vcf = "dna/variants/germline/{sample}.germline.vcf.gz",
        ref = config["reference"]["genome"]
    output:
        vcf = "dna/variants/germline/{sample}.germline.filtered.vcf.gz",
        tbi = "dna/variants/germline/{sample}.germline.filtered.vcf.gz.tbi"
    log:
        "logs/gatk/{sample}.filter_germline.log"
    params:
        snp_filter = "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0",
        indel_filter = "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    shell:
        """
        # 对SNPs应用过滤器
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O temp.snps.vcf.gz \
            2> {log}
            
        gatk VariantFiltration \
            -R {input.ref} \
            -V temp.snps.vcf.gz \
            --filter-expression "{params.snp_filter}" \
            --filter-name "SNP_filter" \
            -O temp.filtered.snps.vcf.gz \
            2>> {log}
            
        # 对INDELs应用过滤器
        gatk SelectVariants \
            -R {input.ref} \
            -V {input.vcf} \
            --select-type-to-include INDEL \
            -O temp.indels.vcf.gz \
            2>> {log}
            
        gatk VariantFiltration \
            -R {input.ref} \
            -V temp.indels.vcf.gz \
            --filter-expression "{params.indel_filter}" \
            --filter-name "INDEL_filter" \
            -O temp.filtered.indels.vcf.gz \
            2>> {log}
            
        # 合并过滤后的SNPs和INDELs
        gatk MergeVcfs \
            -I temp.filtered.snps.vcf.gz \
            -I temp.filtered.indels.vcf.gz \
            -O {output.vcf} \
            2>> {log}
            
        rm temp.*.vcf.gz temp.*.vcf.gz.tbi
        """

# VQSR校准生殖细胞变异（如果有足够样本数量）
rule vqsr:
    input:
        vcf = "dna/variants/germline/{sample}.germline.vcf.gz",
        ref = config["reference"]["genome"],
        hapmap = "resources/vqsr/hapmap_3.3.hg38.vcf.gz",
        omni = "resources/vqsr/1000G_omni2.5.hg38.vcf.gz",
        g1k = "resources/vqsr/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp = config["reference"]["dbsnp"],
        mills = config["reference"]["known_indels"]
    output:
        recal = "dna/variants/germline/{sample}.germline.recal",
        tranches = "dna/variants/germline/{sample}.germline.tranches",
        vcf = "dna/variants/germline/{sample}.germline.vqsr.vcf.gz",
        tbi = "dna/variants/germline/{sample}.germline.vqsr.vcf.gz.tbi"
    log:
        "logs/gatk/{sample}.vqsr.log"
    conda:
        "../../envs/gatk.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 180
    shell:
        """
        # SNP变异质量重新校准
        gatk VariantRecalibrator \
            -R {input.ref} \
            -V {input.vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 {input.omni} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.g1k} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
            -an QD -an FS -an MQ -an MQRankSum -an ReadPosRankSum -an SOR \
            --mode SNP \
            -O {output.recal} \
            --tranches-file {output.tranches} \
            2> {log}

        # 应用SNP校准
        gatk ApplyVQSR \
            -R {input.ref} \
            -V {input.vcf} \
            --recal-file {output.recal} \
            --tranches-file {output.tranches} \
            --truth-sensitivity-filter-level 99.0 \
            --mode SNP \
            -O {output.vcf} \
            2>> {log}
        """
