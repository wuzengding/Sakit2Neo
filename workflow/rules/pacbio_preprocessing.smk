# /workflow/rules/pacbio_preprocessing.smk

def get_raw_pacbio(wildcards):
    subset = samples[
        (samples['sample_id'] == wildcards.sample) & 
        (samples['sampletype'] == wildcards.sampletype) & 
        (samples['data_type'] == 'rna')
    ]
    if subset.empty: return ""
    return str(subset['pacbio'].iloc[0])

def get_flnc_input(wildcards):
    """
    逻辑判断：
    1. 获取原始路径
    2. 如果路径包含 'flnc'，直接返回该路径（跳过 refine）
    3. 如果不包含 'flnc'，返回 refine 规则的输出路径（触发 refine）
    """
    raw_path = get_raw_pacbio(wildcards)
    if "flnc" in raw_path.lower():
        # 如果原始文件已经是 flnc，直接作为比对的输入
        return raw_path
    else:
        # 否则，指向 pacbio_refine 规则的输出，这将强制 Snakemake 运行 refine
        return f"pacbio/refine/{wildcards.sample}_{wildcards.sampletype}.flnc.bam"

rule pacbio_lima:
    input:
        bam = get_raw_pacbio
    output:
        bam = "pacbio/lima/{sample}_{sampletype}.primer_5p--primer_3p.bam",
        xml = "pacbio/lima/{sample}_{sampletype}.consensusreadset.xml",
    log: "logs/pacbio/lima_{sample}_{sampletype}.log"
    params:
        primer = config["reference"]["pacbio_primers"]
    threads: 16
    conda: "../../envs/isoseq.yaml"
    shell:
        """
        lima {input.bam} {params.primer} {output.xml} \
            --hifi-preset ASYMMETRIC \
            --dump-clips \
            --peek-guess \
            --split-named \
            -j {threads} > {log} 2>&1
        """

rule pacbio_refine:
    """产生 FLNC Reads (全长非嵌合)"""
    input:
        #xml = "pacbio/lima/{sample}_{sampletype}.consensusreadset.xml",
        bam = get_raw_pacbio
    output:
        bam = "pacbio/refine/{sample}_{sampletype}.flnc.bam",
        report = "pacbio/refine/{sample}_{sampletype}.flnc.report.csv"
    log: "logs/pacbio/refine_{sample}_{sampletype}.log"
    params:
        primer = config["reference"]["pacbio_primers"]
    threads: 8
    conda: "../../envs/isoseq.yaml"
    shell:
        "isoseq3 refine {input.bam} {params.primer} {output.bam} --require-polya > {log} 2>&1"

rule pacbio_align:
    """
    自动判断输入：
    若是已有的 flnc 则直接比对；
    若是原始数据则等待 pacbio_refine 完成后再比对。
    """
    input:
        bam = get_flnc_input,  # 使用动态判断函数
        ref = config["reference"]["genome"]
    output:
        bam = "pacbio/aligned/{sample}_{sampletype}.aligned.sorted.bam",
        bai = "pacbio/aligned/{sample}_{sampletype}.aligned.sorted.bam.bai"
    threads: 24
    conda: "../../envs/isoseq.yaml"
    shell:
        """
        samtools fastq -@ 4 {input.bam} | \
        minimap2 -ax splice:hq  -s 40 -G 350k -t {threads} --MD -uf --secondary=no {input.ref} - | \
        samtools view -hb -| samtools sort -@ 8 -o {output.bam} -
        samtools index {output.bam}
        """