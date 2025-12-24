# 在文件顶部或 Rule 之前添加约束，防止 Snakemake 乱匹配点号
wildcard_constraints:
    sample = "[^/.]+",      # sample 不包含点和斜杠
    type = "dna_normal|dna_tumor|rna_normal|rna_tumor"

def get_ngscheckmate_inputs(wildcards):
    paths = []
    s = wildcards.sample
    # DNA
    if config["modules"].get("wes_processing", True):
        # 只有在 samples 表里存在的才添加
        dna_exists = not samples[(samples['sample_id'] == s) & (samples['data_type'] == 'dna')].empty
        if dna_exists:
            paths.append(f"qc/ngscheckmate/{s}.dna_normal.vcf")
            paths.append(f"qc/ngscheckmate/{s}.dna_tumor.vcf")
    # RNA
    if config["modules"].get("rna_processing", True):
        rna_exists = not samples[(samples['sample_id'] == s) & (samples['data_type'] == 'rna')].empty
        if rna_exists:
            paths.append(f"qc/ngscheckmate/{s}.rna_normal.vcf")
            paths.append(f"qc/ngscheckmate/{s}.rna_tumor.vcf")
    return paths

rule checkmate_pileup:
    input:
        bam = lambda wildcards: {
            "dna_normal": f"dna/aligned/{wildcards.sample}_normal.recal.bam",
            "dna_tumor":  f"dna/aligned/{wildcards.sample}_tumor.recal.bam",
            "rna_normal": f"rna/align/{wildcards.sample}_normal_dedup.bam",
            "rna_tumor":  f"rna/align/{wildcards.sample}_tumor_dedup.bam"
        }[wildcards.type],
        ref = config["reference"]["genome"],
        bed = config["reference"]["ngscheckmate_bed"]
    output:
        vcf = "qc/ngscheckmate/{sample}.{type}.vcf"
    log:
        "logs/ngscheckmate/{sample}.{type}.log"
    conda:
        "../../envs/checkmate.yaml"
    shell:
        """
        exec 2> {log}
        
        # 1. 获取 BAM 文件中存在的染色体列表
        # samtools idxstats 可以快速列出 BAM 中的染色体名称
        samtools idxstats {input.bam} | cut -f1 > {output.vcf}.chroms
        
        # 2. 动态过滤 BED 文件，只保留 BAM 中存在的染色体位点
        # grep -wFf 会匹配整行/整词，确保 chr1 不会误匹配到 chr10
        grep -wFf {output.vcf}.chroms {input.bed} > {output.vcf}.tmp.bed
        
        # 3. 使用过滤后的临时 BED 运行分析
        bcftools mpileup -f {input.ref} -R {output.vcf}.tmp.bed {input.bam} | \
        bcftools call -mv -Ov > {output.vcf}
        
        # 4. 清理临时文件
        rm -f {output.vcf}.chroms {output.vcf}.tmp.bed
        """

rule ngscheckmate:
    input:
        vcfs = get_ngscheckmate_inputs,
        viz_script = workflow.source_path("../scripts/qc_checkmate.R")
    output:
        matrix = "qc/ngscheckmate/{sample}_all.txt",
        matched = "qc/ngscheckmate/{sample}_matched.txt",
        heatmap = "qc/ngscheckmate/{sample}_heatmap.png",
        scatter_dna = "qc/ngscheckmate/{sample}_scatter_dna_normal_vs_dna_tumor.png"
    log:
        "logs/ngscheckmate/{sample}.log"
    params:
        outdir = "qc/ngscheckmate",
        bed = config["reference"]["ngscheckmate_bed"],
        ref = config["reference"]["genome"],
        vcf_list = "qc/ngscheckmate/{sample}.vcf_list.txt"
        
    conda:
        "../../envs/checkmate.yaml"
    shell:
        """
        echo "=== Starting NGSCheckMate ===" > {log}
        export NCM_REF={params.ref}

        # 1. 准备 VCF 列表
        rm -f {params.vcf_list}
        for vcf in {input.vcfs}; do
            readlink -f "$vcf" >> {params.vcf_list}
        done

        # 2. 运行 NGSCheckMate (Python 2.7 环境)
        ncm.py -V \
            -l {params.vcf_list} \
            -bed {params.bed} \
            -O {params.outdir} \
            -N {wildcards.sample} >> {log} 2>&1

        # 3. 运行 R 绘图脚本
        echo "=== Running R Visualization Script ===" >> {log}
        Rscript {input.viz_script} \
            --matched {output.matched} \
            --vcf_dir {params.outdir} \
            --out {params.outdir}/{wildcards.sample} >> {log} 2>&1

        rm -f {params.vcf_list}
        """