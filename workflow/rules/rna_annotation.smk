# rna_annotation.smk
# ========================
# RNA 变异注释流程
# 逻辑复用 DNA 注释，但输入输出指向 RNA 目录
# ========================

# 定义辅助函数获取需要注释的样本 (适配 mutect2 和 varscan2)
def get_rna_vcf_inputs(wildcards):
    """
    根据 caller 获取对应的输入文件。
    我们在 rna_somatic_variants.smk 中定义的输出命名格式为:
    rna/variants/{caller}/{sample}.{caller}.vcf.gz
    """
    return f"rna/variants/{wildcards.caller}/{wildcards.sample}.{wildcards.caller}.vcf.gz"

def get_rna_vcf_idx(wildcards):
    return f"rna/variants/{wildcards.caller}/{wildcards.sample}.{wildcards.caller}.vcf.gz.tbi"

# ========================
# 资源准备部分 (Reference Preparation)
# ========================
# 注意：如果你的 pipeline 同时包含 dna_annotation.smk，
# 下面这两个 prepare 规则可能会因为“规则重名”或“产生相同输出”而冲突。
# 最佳实践是将 prepare_gnomad 和 prepare_cosmic 移到一个公共的 common_reference.smk 中。
# 
# 为了保证本脚本的独立性，我在这里保留它们，但使用了 onstart/check 逻辑或假设这些文件
# 由 DNA 流程生成。此处代码复用了 DNA 的输出路径，以节省磁盘空间。
# ========================

# 如果你确定 dna_annotation.smk 会被加载，可以将下面两个 rule (prepare_gnomad_resource, prepare_cosmic_resource) 注释掉。
# 这里为了安全起见，我给规则改了名，但输出文件指向相同路径。Snakemake 会自动处理：如果文件已存在，不会重新运行。

#rule prepare_gnomad_resource:
#    input:
#        gnomad = config["reference"]["gnomad"]
#    output:
#        # 复用 DNA 流程生成的资源文件
#        annot_source = "dna/variants/annotated/gnomad.annot_source.chr.tsv.gz",
#        contig_map = "dna/variants/annotated/gnomad.contig_map.txt"
#    log:
#        "logs/annotation/prepare_gnomad_shared.log"
#    conda:
#        "../../envs/bcftools.yaml"
#    shell:
#        """
#        # 逻辑同 DNA，仅当文件不存在时运行
#        if [ ! -f "{output.annot_source}" ]; then
#            ( \
#                echo -e '##fileformat=BCFTOOLS_QUERY_FORMAT' && \
#                echo -e '#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]gnomAD_AF\t[6]gnomAD_AF_eas' && \
#                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_eas\n' {input.gnomad} \
#            ) | bgzip -c > {output.annot_source}
#            tabix -s1 -b2 -e2 {output.annot_source}
#            
#            bcftools view -h {input.gnomad} | grep "^##contig" | sed -e 's/.*ID=\([^,>]*\).*/\1\tchr\1/' > {output.contig_map}
#        else
#            echo "Shared gnomAD resource already exists. Skipping." > {log}
#        fi
#        """

#rule prepare_cosmic_resource:
#    input:
#        cosmic = config["reference"]["cosmic"]
#    output:
#        # 复用 DNA 流程生成的资源文件
#        renamed_cosmic = "dna/variants/annotated/cosmic.renamed.vcf.gz",
#        renamed_cosmic_idx = "dna/variants/annotated/cosmic.renamed.vcf.gz.tbi"
#    log:
#        "logs/annotation/prepare_cosmic_shared.log"
#    conda:
#        "../../envs/bcftools.yaml"
#    shell:
#        """
#        if [ ! -f "{output.renamed_cosmic}" ]; then
#            bcftools view -h {input.cosmic} 2>{log} | \
#                sed -e 's/##contig=<ID=/##contig=<ID=chr/' \
#                    -e 's/ID=CNT,/ID=COSMIC_CNT,/' \
#                    -e 's/ID=AA,/ID=COSMIC_AA,/' \
#                    -e 's/ID=CDS,/ID=COSMIC_CDS,/' \
#                > header.txt
#
#            bcftools view -H {input.cosmic} 2>>{log} | \
#                awk 'BEGIN {{OFS="\\t"}} \
#                     {{ \
#                        $1="chr"$1; \
#                        gsub("CNT=", "COSMIC_CNT=", $8); \
#                        gsub("AA=", "COSMIC_AA=", $8); \
#                        gsub("CDS=", "COSMIC_CDS=", $8); \
#                        print \
#                     }}' \
#                > body.txt
#            
#            cat header.txt body.txt | bgzip -c > {output.renamed_cosmic} 2>>{log}
#            tabix -p vcf {output.renamed_cosmic} 2>>{log}
#            rm header.txt body.txt
#        else
#             echo "Shared COSMIC resource already exists. Skipping." > {log}
#             # 确保索引也存在，否则重新索引
#             if [ ! -f "{output.renamed_cosmic}.tbi" ]; then
#                 tabix -p vcf {output.renamed_cosmic}
#             fi
#        fi
#        """

# ========================
# 步骤 1: RNA gnomAD 注释
# ========================
rule annotate_gnomad_rna:
    input:
        vcf = get_rna_vcf_inputs,
        vcf_idx = get_rna_vcf_idx,
        annot_source = "dna/variants/annotated/gnomad.annot_source.chr.tsv.gz", # 使用共享资源
        contig_map = "dna/variants/annotated/gnomad.contig_map.txt"
    output:
        vcf = "rna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz",
    log:
        "logs/rna_annotation/{caller}/{sample}.gnomad.log"
    conda:
        "../../envs/bcftools.yaml"
    shell:
        """
        bcftools annotate \
            -a {input.annot_source} \
            -c CHROM,POS,REF,ALT,INFO/gnomAD_AF,INFO/gnomAD_AF_eas \
            --rename-chrs {input.contig_map} \
            -h <(echo -e '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="Global Allele Frequency from gnomAD (all populations)">\n##INFO=<ID=gnomAD_AF_eas,Number=A,Type=Float,Description="East Asian Allele Frequency from gnomAD">') \
            {input.vcf} \
            -O z -o {output.vcf} \
            --threads {threads} \
            > {log} 2>&1
        """

# ========================
# 步骤 2: RNA COSMIC 注释
# ========================
rule annotate_cosmic_rna:
    input:
        vcf = "rna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz",
        vcf_idx = "rna/variants/annotated/{sample}.{caller}.gnomad.vcf.gz.tbi",
        cosmic = "dna/variants/annotated/cosmic.renamed.vcf.gz", # 使用共享资源
        cosmic_idx = "dna/variants/annotated/cosmic.renamed.vcf.gz.tbi"
    output:
        vcf = "rna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz",
    log:
        "logs/rna_annotation/{sample}.{caller}.cosmic.log"
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
# 步骤 3: RNA VEP 功能注释
# ========================
rule vep_annotate_rna:
    input:
        vcf = "rna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz",
        vcf_idx = "rna/variants/annotated/{sample}.{caller}.cosmic.vcf.gz.tbi",
        fasta = config["reference"]["genome"],
    output:
        vcf = "rna/variants/annotated/{sample}.{caller}.vep.vcf.gz",
        vcf_idx = "rna/variants/annotated/{sample}.{caller}.vep.vcf.gz.tbi",
        html = "rna/variants/annotated/{sample}.{caller}.vep_summary.html"
    params:
        cache_dir = config["vep"]["cache_dir"],
        species = config["vep"]["species"],
        assembly = config["vep"]["assembly"],
        fork = config["vep"]["fork"],
        plugins_dir = config["vep"].get("plugins_dir", ""),
        plugins = ",".join(config["vep"]["plugins_to_install"])
    log:
        "logs/rna_annotation/{sample}.{caller}.vep.log"
    conda:
        "../../envs/vep.yaml"
    resources:
        mem_mb = 16000,
        time = 360,
        threads = 4
    shell:
        """
        # RNA VEP 注释参数与 DNA 保持一致，使用 --everything 包含所有信息
        # 如果需要专门针对 RNA 筛选主要转录本，可以考虑添加 --pick (但这里为了与 DNA 结果可比，暂不加)
        
        PLUGIN_ARGS=""
        if [ -n "{params.plugins_dir}" ] && [ -n "{params.plugins}" ]; then
            PLUGIN_ARGS="--dir_plugins {params.plugins_dir} --plugin {params.plugins}"
        fi
        
        vep \
            -i {input.vcf} \
            -o {output.vcf} \
            --stats_file {output.html} \
            --vcf --compress_output bgzip \
            --force_overwrite \
            --cache --dir_cache {params.cache_dir} \
            --offline \
            --fasta {input.fasta} \
            --everything \
            --fork {params.fork} \
            --protein \
            $PLUGIN_ARGS 2> {log}

        tabix -p vcf {output.vcf}
        """