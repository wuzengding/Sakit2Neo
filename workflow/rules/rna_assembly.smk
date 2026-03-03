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

# 新增 rule: 对 stringtie_assemble 的输出 GTF 进行排序、bgzip、tabix
rule sort_bgzip_tabix_assembly_gtf:
    """
    对 stringtie_assemble 产生的GTF文件进行排序、bgzip压缩和tabix索引
    """
    input:
        gtf = "rna/assembly/{sample}_{sampletype}.assembly.gtf"
    output:
        gtf_gz = "rna/assembly/{sample}_{sampletype}.assembly.gtf.gz",
        gtf_gz_tbi = "rna/assembly/{sample}_{sampletype}.assembly.gtf.gz.tbi"
    log: "logs/stringtie/{sample}_{sampletype}.index.log"
    threads: 1 # bgzip和tabix通常单线程效率更高
    conda: "../../envs/stringtie.yaml" # 使用更新后的 stringtie.yaml
    shell:
        """
        # 对GTF文件进行排序 (保留注释行)
        (grep '^#' {input.gtf}; grep -v '^#' {input.gtf} | sort -k1,1 -k4,4n) | \
        bgzip -c > {output.gtf_gz} 2> {log} && \
        tabix -p gff {output.gtf_gz} 2>> {log}
        """

rule get_assembly_fasta:
    """
    使用 gffread 提取组装转录本的序列
    -w 输出转录本序列 FASTA
    -g 参考基因组 FASTA
    """
    input:
        # 输入现在是经过排序、bgzip和tabix处理后的 GTF.gz 文件
        gtf = "rna/assembly/{sample}_{sampletype}.assembly.gtf.gz",
        gtf_tbi = "rna/assembly/{sample}_{sampletype}.assembly.gtf.gz.tbi", # 确保索引文件也存在
        ref = config["reference"]["genome"]
    output:
        fasta = "rna/assembly/{sample}_{sampletype}.assembly.transcripts.fa"
    log: "logs/stringtie/{sample}_{sampletype}_get_fasta.log"
    conda: "../../envs/stringtie.yaml"
    shell:
        """
        gffread {input.gtf} \
            -g {input.ref} \
            -w {output.fasta} > {log} 2>&1
        """

rule stringtie_merge:
    """第二步：合并同一个患者的 Tumor 和 Normal 组装结果 (统一 ID)"""
    input:
        t_gtf = "rna/assembly/{sample}_tumor.assembly.gtf.gz",
        t_gtf_tbi = "rna/assembly/{sample}_tumor.assembly.gtf.gz.tbi",
        n_gtf = "rna/assembly/{sample}_normal.assembly.gtf.gz",
        n_gtf_tbi = "rna/assembly/{sample}_normal.assembly.gtf.gz.tbi",
        ref_gtf = config["reference"]["gtf"]
    output:
        merged_gtf = "rna/assembly/{sample}_merged.assembly.gtf"
    log: 
        "logs/stringtie/merge_{sample}.log"
    threads: 4
    conda: 
        "../../envs/stringtie.yaml"
    shell:
        """
        stringtie --merge \
            -G {input.ref_gtf} \
            -o {output.merged_gtf} \
            -p {threads} \
            <(gunzip -c {input.t_gtf}) \
            <(gunzip -c {input.n_gtf}) \
            > {log} 2>&1
        """

# 新增 rule: 对 stringtie_merge 的输出 GTF 进行排序、bgzip、tabix
rule sort_bgzip_tabix_merged_gtf:
    """
    对 stringtie_merge 产生的合并GTF文件进行排序、bgzip压缩和tabix索引
    """
    input:
        gtf = "rna/assembly/{sample}_merged.assembly.gtf"
    output:
        gtf_gz = "rna/assembly/{sample}_merged.assembly.gtf.gz",
        gtf_gz_tbi = "rna/assembly/{sample}_merged.assembly.gtf.gz.tbi"
    log: "logs/stringtie/{sample}.merge.index.log"
    threads: 1
    conda: "../../envs/stringtie.yaml" # 使用更新后的 stringtie.yaml
    shell:
        """
        (grep '^#' {input.gtf}; grep -v '^#' {input.gtf} | sort -k1,1 -k4,4n) | \
        bgzip -c > {output.gtf_gz} 2> {log} && \
        tabix -p gff {output.gtf_gz} 2>> {log}
        """

rule get_fasta_merged:
    """提取合并（统一ID）后的转录本序列，这通常用于后续的新抗原分析"""
    input:
        # 输入现在是经过排序、bgzip和tabix处理后的 GTF.gz 文件
        gtf = "rna/assembly/{sample}_merged.assembly.gtf.gz",
        gtf_tbi = "rna/assembly/{sample}_merged.assembly.gtf.gz.tbi", # 确保索引文件也存在
        ref = config["reference"]["genome"]
    output:
        fasta = "rna/assembly/{sample}_merged.transcripts.fa"
    log: "logs/gffread/{sample}_merged_get_fasta.log" # 增加log文件，保持一致
    conda: "../../envs/stringtie.yaml"
    shell:
        """
        zcat {input.gtf} | \
        gffread -g {input.ref} -w {output.fasta} /dev/stdin > {log} 2>&1
        """
