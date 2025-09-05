rule fastp_qc_trim_dna:
    input:
        r1 = lambda wildcards: get_fq1(wildcards),
        r2 = lambda wildcards: get_fq2(wildcards)
    output:
        r1 = "dna/trimmed/{sample}_{sampletype}_R1.trimmed.fastq.gz",
        r2 = "dna/trimmed/{sample}_{sampletype}_R2.trimmed.fastq.gz",
        json = "qc/fastp/{sample}_{sampletype}_dna.fastp.json",
        html = "qc/fastp/{sample}_{sampletype}_dna.fastp.html"
    log:
        "logs/fastp/{sample}_{sampletype}_dna.log"
    params:
        # 移除显式 adapter 文件，启用自动检测
        extra = "--detect_adapter_for_pe --trim_poly_g --correction"
    conda:
        "../../envs/fastp.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        threads = 20,
        time = lambda wildcards, attempt: attempt * 60
    shell:
        """
        mkdir -p dna/trimmed qc/fastp
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --json {output.json} \
            --html {output.html} \
            --report_title "{wildcards.sample}_{wildcards.sampletype}_DNA" \
            --thread {resources.threads} \
            {params.extra} \
            > {log} 2>&1
        """
