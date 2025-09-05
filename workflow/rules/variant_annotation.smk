# 变异注释规则

# 使用SnpEff注释生殖细胞变异
rule annotate_germline_variants:
    input:
        vcf = "dna/variants/germline/{sample}.germline.filtered.vcf.gz"
    output:
        vcf = "dna/variants/germline/{sample}.germline.annotated.vcf.gz",
        tbi = "dna/variants/germline/{sample}.germline.annotated.vcf.gz.tbi",
        stats = "dna/variants/germline/{sample}.germline.snpeff.stats.html"
    log:
        "logs/snpeff/{sample}.germline.log"
    params:
        genome = config["annotation"]["snpeff_genome"],
        extra = "-Xmx8g"
    conda:
        "../envs/snpeff.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = lambda wildcards, attempt: attempt * 180
    shell:
        """
        snpEff {params.extra} -v {params.genome} -stats {output.stats} \
            {input.vcf} | bgzip > {output.vcf}
        
        tabix -p vcf {output.vcf}
        """

# 使用VEP注释体细胞变异
rule annotate_somatic_variants:
    input:
        vcf = "dna/variants/somatic/{sample}.somatic.consensus.vcf.gz"
    output:
        vcf = "dna/variants/somatic/{sample}.somatic.annotated.vcf.gz",
        tbi = "dna/variants/somatic/{sample}.somatic.annotated.vcf.gz.tbi",
        stats = "dna/variants/somatic/{sample}.somatic.vep.stats.html"
    log:
        "logs/vep/{sample}.somatic.log"
    params:
        cache = config["annotation"]["vep_cache"],
        plugins = "Downstream,REVEL,CADD,dbNSFP",
        assembly = config["annotation"]["assembly"],
        extra = "--everything --hgvs --af --af_1kg --af_gnomad"
    conda:
        "../envs/vep.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 16000,
        time = lambda wildcards, attempt: attempt * 240
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file STDOUT \
            --cache \
            --dir_cache {params.cache} \
            --format vcf \
            --vcf \
            --compress_output bgzip \
            --stats_file {output.stats} \
            --fork 4 \
            --assembly {params.assembly} \
            --species homo_sapiens \
            --buffer_size 5000 \
            --regulatory \
            --force_overwrite \
            --tab \
            --verbose \
            --check_existing \
            --sift b \
            --polyphen b \
            --nearest symbol \
            --plugin {params.plugins} \
            {params.extra} \
            2> {log} | bgzip > {output.vcf}
        
        tabix -p vcf {output.vcf}
        """

# 创建变异摘要报告
rule variant_summary:
    input:
        germline = "dna/variants/germline/{sample}.germline.annotated.vcf.gz",
        somatic = "dna/variants/somatic/{sample}.somatic.annotated.vcf.gz"
    output:
        report = "reports/variants/{sample}.variant_summary.html"
    log:
        "logs/reports/{sample}.variant_summary.log"
    conda:
        "../envs/r.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = lambda wildcards, attempt: attempt * 60
    script:
        "../scripts/variant_summary.R"

