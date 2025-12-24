# /workflow/rules/rna_assembly.smk

def get_adaptive_assemble_input(wildcards):
    """自适应判断：返回 Illumina BAM 和可选的 PacBio BAM"""
    sample_row = samples[(samples['sample_id'] == wildcards.sample) & 
                         (samples['sampletype'] == wildcards.sampletype) & 
                         (samples['data_type'] == 'rna')]
    
    inputs = {
        "illumina_bam": f"rna/align/{wildcards.sample}_{wildcards.sampletype}/Aligned.sortedByCoord.out.bam"
    }
    
    # 只有当样本表中有 pacbio 路径时，才添加对应的比对结果作为输入
    if not sample_row.empty and pd.notna(sample_row['pacbio'].iloc[0]) and str(sample_row['pacbio'].iloc[0]).strip() != "":
        inputs["pacbio_bam"] = f"pacbio/aligned/{wildcards.sample}_{wildcards.sampletype}.aligned.sorted.bam"
    
    return inputs

rule stringtie_assemble:
    """第一步：组装 (自适应 Hybrid 或 Illumina-only)"""
    input:
        unpack(get_adaptive_assemble_input),
        ref_gtf = config["reference"]["gtf"]
    output:
        gtf = "rna/assembly/{sample}_{sampletype}.assembly.gtf"
    params:
        # 使用 hasattr 安全判断是否存在 pacbio 输入
        bam_list = lambda wildcards, input: f"{input.illumina_bam} {input.pacbio_bam}" if hasattr(input, "pacbio_bam") else input.illumina_bam,
        extra = lambda wildcards, input: "-L" if hasattr(input, "pacbio_bam") else ""
    log: "logs/stringtie/assemble_{sample}_{sampletype}.log"
    threads: 12
    conda: "../../envs/stringtie.yaml"
    shell:
        """
        stringtie {params.bam_list} \
            {params.extra} \
            -G {input.ref_gtf} \
            -o {output.gtf} \
            -p {threads} \
            -l {wildcards.sample}_{wildcards.sampletype} > {log} 2>&1
        """

rule stringtie_merge:
    """第二步：合并同一个患者的 Tumor 和 Normal 组装结果 (统一 ID)"""
    input:
        t_gtf = "rna/assembly/{sample}_tumor.assembly.gtf",
        n_gtf = "rna/assembly/{sample}_normal.assembly.gtf",
        ref_gtf = config["reference"]["gtf"]
    output:
        merged_gtf = "rna/assembly/merged/{sample}_merged.gtf"
    log: "logs/stringtie/merge_{sample}.log"
    threads: 4
    conda: "../../envs/stringtie.yaml"
    shell:
        "stringtie --merge -G {input.ref_gtf} -o {output.merged_gtf} -p {threads} {input.t_gtf} {input.n_gtf} > {log} 2>&1"

rule stringtie_quant:
    """第三步：最终定量 (基于统一地图，使用 Illumina 数据计数)"""
    input:
        bam = f"rna/align/{{sample}}_{{sampletype}}/Aligned.sortedByCoord.out.bam",
        gtf = "rna/assembly/{sample}_merged.gtf"
    output:
        gtf = "rna/expression/{sample}_{sampletype}.quant.gtf",
        abund = "rna/expression/{sample}_{sampletype}.abundance.txt"
    log: "logs/stringtie/quant_{sample}_{sampletype}.log"
    threads: 8
    conda: "../../envs/stringtie.yaml"
    shell:
        """
        stringtie {input.bam} \
            -G {input.gtf} \
            -e -B \
            -p {threads} \
            -o {output.gtf} \
            -A {output.abund} > {log} 2>&1
        """

rule stringtie_count_matrix:
    """第四步：生成 Count 矩阵用于下游差异表达分析"""
    input:
        t_gtf = "rna/expression/{sample}_tumor.quant.gtf",
        n_gtf = "rna/expression/{sample}_normal.quant.gtf"
    output:
        gene_counts = "rna/expression/{sample}_gene_count_matrix.csv",
        tx_counts = "rna/expression/{sample}_transcript_count_matrix.csv"
    params:
        slist = "rna/expression/{sample}_sample_list.txt"
    log: "logs/stringtie/prepDE_{sample}.log"
    conda: "../../envs/stringtie.yaml"
    shell:
        """
        echo "{wildcards.sample}_tumor {input.t_gtf}" > {params.slist}
        echo "{wildcards.sample}_normal {input.n_gtf}" >> {params.slist}
        prepDE.py -i {params.slist} -g {output.gene_counts} -t {output.tx_counts} > {log} 2>&1
        """