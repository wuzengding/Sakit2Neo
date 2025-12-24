# /workflow/rules/pacbio_preprocessing.smk

def get_raw_pacbio(wildcards):
    subset = samples[
        (samples['sample_id'] == wildcards.sample) & 
        (samples['sampletype'] == wildcards.sampletype) & 
        (samples['data_type'] == 'rna')
    ]
    if subset.empty: return ""
    return str(subset['pacbio'].iloc[0])


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
        xml = "pacbio/lima/{sample}_{sampletype}.consensusreadset.xml",
    output:
        bam = "pacbio/refine/{sample}_{sampletype}.flnc.bam",
        report = "pacbio/refine/{sample}_{sampletype}.flnc.report.csv"
    log: "logs/pacbio/refine_{sample}_{sampletype}.log"
    params:
        primer = config["reference"]["pacbio_primers"]
    threads: 8
    conda: "../../envs/isoseq.yaml"
    shell:
        "isoseq3 refine {input.xml} {params.primer} {output.bam} --require-polya > {log} 2>&1"

rule pacbio_align:
    """跳过中间 Fastq 存储，直接流式比对"""
    input:
        bam = "pacbio/refine/{sample}_{sampletype}.flnc.bam",
        ref = config["reference"]["genome"]
    output:
        bam = "pacbio/aligned/{sample}_{sampletype}.aligned.sorted.bam",
        bai = "pacbio/aligned/{sample}_{sampletype}.aligned.sorted.bam.bai"
    threads: 24
    conda: "../../envs/isoseq.yaml"
    shell:
        """
        # 使用 samtools fastq 直接通过管道传给 minimap2
        # minimap2 的输入端写 - 代表接收管道数据
        samtools fastq -@ 4 {input.bam} | \
        minimap2 -ax splice:hq -uf --secondary=no -t {threads} {input.ref} - | \
        samtools sort -@ 8 -o {output.bam} -
        samtools index {output.bam}
        """