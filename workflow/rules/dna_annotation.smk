# ========================
#  gnomAD人群频率注释 
# ========================

#   步骤 1a: 准备 gnomAD 注释库 (新规则, 全局一次性)
# ========================
rule prepare_gnomad_for_annotation:
    input:
        gnomad = config["reference"]["gnomad"]
    output:
        # 输出一个全局共享的、带chr前缀的、精简的TSV注释文件
        annot_source = "dna/variants/annotated/gnomad.annot_source.chr.tsv.gz",
        annot_source_idx = "dna/variants/annotated/gnomad.annot_source.chr.tsv.gz.tbi",
        contig_map = "dna/variants/annotated/gnomad.contig_map.txt"
    log:
        "logs/annotation/prepare_gnomad.log"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        # ================== 步骤 1: 创建临时的注释源文件 (不处理染色体命名) ==================
        # 提取所需字段，并重命名列头。注意：此时输出的染色体名称是gnomAD原始的 (e.g., '1', '2'...)
        ( \
            echo -e '##fileformat=BCFTOOLS_QUERY_FORMAT' && \
            echo -e '#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]gnomAD_AF\t[6]gnomAD_AF_eas' && \
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_eas\n' {input.gnomad} \
        ) | bgzip -c > {output.annot_source}
        
        tabix -s1 -b2 -e2 {output.annot_source}

        # ================== 步骤 2: 创建染色体名称映射文件 ==================
        # 这个文件将告诉 bcftools annotate 如何将注释源中的 '1' 映射到 VCF 中的 'chr1'
        bcftools view -h {input.gnomad} | grep "^##contig" | sed -e 's/.*ID=\([^,>]*\).*/\1\tchr\1/' > {output.contig_map}
        """

rule annotate_gnomad:
    input:
        vcf = "dna/variants/{caller}/{sample}.{caller}.vcf.gz",
        # 输入变为全局预处理好的注释文件
        annot_source = "dna/variants/annotated/gnomad.annot_source.chr.tsv.gz",
        contig_map = "dna/variants/annotated/gnomad.contig_map.txt"
    output:
        vcf = "dna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz",
        vcf_idx = "dna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz.tbi"
    log:
        "logs/annotation/{caller}/{sample}.gnomad.log"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        # ================== 步骤 3: 使用临时文件和 --rename-chrs 进行注释 ==================
        # --rename-chrs 参数在这里使用是正确的，它会根据映射文件来调整 -a 文件中的染色体名称以匹配输入VCF
        bcftools annotate \
            -a {input.annot_source} \
            -c CHROM,POS,REF,ALT,INFO/gnomAD_AF,INFO/gnomAD_AF_eas \
            --rename-chrs {input.contig_map} \
            -h <(echo -e '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Global Allele Frequency from gnomAD (all populations)">\n##INFO=<ID=gnomAD_AF_eas,Number=A,Type=Float,Description="East Asian Allele Frequency from gnomAD">') \
            {input.vcf} \
            -O z -o {output.vcf} \
            --threads {threads} \
            > {log} 2>&1

        # ================== 步骤 4: 为最终输出创建索引 ==================
        tabix -p vcf {output.vcf} 2>> {log}
        """

#rule annotate_gnomad:
#    input:
#        vcf = "dna/variants/{caller}/{sample}.{caller}.vcf.gz",
#        gnomad = config["reference"]["gnomad"]
#    output:
#        vcf = "dna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz",
#        vcf_idx = "dna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz.tbi"
#    log:
#        "logs/annotation/{caller}/{sample}.gnomad.log"
#    conda:
#        "../../envs/bcftools.yaml"
#    params:
#        annot_source = temp("dna/variants/annotated/{sample}.gnomad_source.tsv.gz"),
#        contig_map = temp("logs/annotation/{caller}/{sample}.contig_map.txt")
#    resources:
#        mem_mb = 8000,
#        time = 180,
#        threads = 6
#    shell:
#        """
#        mkdir -p dna/variants/annotated
#
#        # ================== 步骤 1: 创建临时的注释源文件 (不处理染色体命名) ==================
#        # 提取所需字段，并重命名列头。注意：此时输出的染色体名称是gnomAD原始的 (e.g., '1', '2'...)
#        ( \
#            echo -e '##fileformat=BCFTOOLS_QUERY_FORMAT' && \
#            echo -e '#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]gnomAD_AF\t[6]gnomAD_AF_eas' && \
#            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_eas\n' {input.gnomad} \
#        ) | bgzip -c > {params.annot_source}
#        
#        tabix -s1 -b2 -e2 {params.annot_source}
#
#        # ================== 步骤 2: 创建染色体名称映射文件 ==================
#        # 这个文件将告诉 bcftools annotate 如何将注释源中的 '1' 映射到 VCF 中的 'chr1'
#        bcftools view -h {input.gnomad} | grep "^##contig" | sed -e 's/.*ID=\([^,>]*\).*/\1\tchr\1/' > {params.contig_map}
#
#        # ================== 步骤 3: 使用临时文件和 --rename-chrs 进行注释 ==================
#        # --rename-chrs 参数在这里使用是正确的，它会根据映射文件来调整 -a 文件中的染色体名称以匹配输入VCF
#        bcftools annotate \
#            -a {params.annot_source} \
#            -c CHROM,POS,REF,ALT,INFO/gnomAD_AF,INFO/gnomAD_AF_eas \
#            --rename-chrs {params.contig_map} \
#            -h <(echo -e '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Global Allele Frequency from gnomAD (all populations)">\n##INFO=<ID=gnomAD_AF_eas,Number=A,Type=Float,Description="East Asian Allele Frequency from gnomAD">') \
#            {input.vcf} \
#            -O z -o {output.vcf} \
#            --threads {threads} \
#            > {log} 2>&1
#
#        # ================== 步骤 4: 为最终输出创建索引 ==================
#        tabix -p vcf {output.vcf} 2>> {log}
#        """

# ========================
# 准备COSMIC注释库 (最终强化版：同时处理染色体和字段名)
# ========================
rule prepare_cosmic_for_annotation:
    input:
        cosmic = config["reference"]["cosmic"]
    output:
        # 这个输出是临时的，会被下游规则使用后自动删除
        renamed_cosmic = "dna/variants/annotated/cosmic.renamed.vcf.gz",
        renamed_cosmic_idx = "dna/variants/annotated/cosmic.renamed.vcf.gz.tbi"
    log:
        "logs/annotation/prepare_cosmic.log"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        # 步骤 1: 导出并修改 header
        # 使用管道和 sed 一次性完成 header 的染色体和字段名修改
        bcftools view -h {input.cosmic} 2>{log} | \
            sed -e 's/##contig=<ID=/##contig=<ID=chr/' \
                -e 's/ID=CNT,/ID=COSMIC_CNT,/' \
                -e 's/ID=AA,/ID=COSMIC_AA,/' \
                -e 's/ID=CDS,/ID=COSMIC_CDS,/' \
            > header.txt

        # 步骤 2: 导出数据行，并用 awk 给第一列（染色体）加上 'chr' 前缀，并修改INFO字段
        # 使用管道，避免生成巨大的中间文件
        bcftools view -H {input.cosmic} 2>>{log} | \
            awk 'BEGIN {{OFS="\\t"}} \
                 {{ \
                    $1="chr"$1; \
                    gsub("CNT=", "COSMIC_CNT=", $8); \
                    gsub("AA=", "COSMIC_AA=", $8); \
                    gsub("CDS=", "COSMIC_CDS=", $8); \
                    print \
                 }}' \
            > body.txt
        
        # 步骤 3: 合并修改后的 header 和 body，并压缩
        cat header.txt body.txt | bgzip -c > {output.renamed_cosmic} 2>>{log}

        # 步骤 4: 为这个全新的、完美的注释文件创建索引
        tabix -p vcf {output.renamed_cosmic} 2>>{log}

        # 步骤 5: 清理临时文件
        rm header.txt body.txt
        """

# ========================
# COSMIC数据库注释 (最终简化版：使用完美预处理的文件)
# ========================
rule annotate_cosmic:
    input:
        vcf = "dna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz",
        cosmic = "dna/variants/annotated/cosmic.renamed.vcf.gz"
    output:
        vcf = "dna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz",
        vcf_idx = "dna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz.tbi"
    log:
        "logs/annotation/{sample}.{caller}.cosmic.log"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate \
            -a {input.cosmic} \
            -c 'INFO/COSMIC_CNT,INFO/COSMIC_AA,INFO/COSMIC_CDS' \
            {input.vcf} \
            -O z -o {output.vcf} \
            2> {log}
        
        tabix -p vcf {output.vcf} 2>>{log}
        """

# ========================
# VEP 功能注释 (最终版 - 兼容 VEP 110.1)
# ========================
rule vep_annotate:
    input:
        vcf = "dna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz",
        vcf_idx = "dna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz.tbi",
        fasta = config["reference"]["genome"],
        fai = f'{config["reference"]["genome"]}.fai'
    output:
        vcf = "dna/variants/annotated/{sample}.{caller}.vep.vcf.gz",
        vcf_idx = "dna/variants/annotated/{sample}.{caller}.vep.vcf.gz.tbi",
        html = "dna/variants/annotated/{sample}.{caller}.vep_summary.html"
    params:
        cache_dir = config["vep"]["cache_dir"],
        species = config["vep"]["species"],
        assembly = config["vep"]["assembly"],
        fork = config["vep"]["fork"],
        plugins_dir = config["vep"].get("plugins_dir", ""), # 使用.get()增加健壮性
        plugins = ",".join(config["vep"]["plugins_to_install"])
        
    log:
        "logs/annotation/{sample}.{caller}.vep.log"
    conda:
        "../../envs/vep.yaml"
    resources:
        mem_mb = 16000,
        time = 360,
        threads = 4
    shell:
        """
        # 使用 --everything 参数来开启大部分常用注释功能
        # 这是最稳定、最不容易出错的参数组合
         PLUGIN_ARGS=""
        if [ -n "{params.plugins_dir}" ] && [ -n "{params.plugins}" ]; then
            PLUGIN_ARGS="--dir_plugins {params.plugins_dir} --plugin {params.plugins}"
        fi
        vep   -i {input.vcf}   \
                -o {output.vcf}  \
                --stats_file {output.html}  \
                --vcf  --compress_output bgzip  \
                --force_overwrite   \
                --cache  --dir_cache {params.cache_dir} \
                --offline  \
                --fasta {input.fasta}   \
                --everything   \
                --fork {params.fork}  \
                --protein \
                $PLUGIN_ARGS   2> {log}
        # 为最终输出文件创建索引
        tabix -p vcf {output.vcf}
        """