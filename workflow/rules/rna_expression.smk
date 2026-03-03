#RNAseq数据用stringtie二次精确定量

rule stringtie_quant:
    """第三步：最终定量 (基于统一地图，使用 Illumina 数据计数)"""
    input:
        bam = f"rna/align/{{sample}}_{{sampletype}}/Aligned.sortedByCoord.out.bam",
        gtf = "rna/assembly/{sample}_merged.assembly.gtf"
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