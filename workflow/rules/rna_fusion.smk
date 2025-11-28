# In your Snakefile - REVISED VERSION

rule arriba_fusion_detection:
    """
    Detect gene fusions from STAR alignment outputs using the two-step method,
    enhanced with additional annotation databases.
    """
    input:
        # Arriba 在两步法中不直接使用BAM，而是依赖 Chimeric.out.junction
        chimeric_junctions = "rna/align/{sample}_{sampletype}/Chimeric.out.junction",
        ref_genome = config["reference"]["genome"],
        ref_gtf = config["reference"]["gtf"],
        blacklist = config["reference"]["arriba_blacklist"],
        known_fusions = config["reference"]["arriba_known_fusions"],
        protein_domains = config["reference"]["arriba_protein_domains"]
    output:
        fusions = "rna/fusion/{sample}_{sampletype}.fusion.tsv",
        fusions_discarded = "rna/fusion/{sample}_{sampletype}.fusion.discarded.tsv"
    log:
        "logs/arriba/{sample}_{sampletype}.log"
    params:
        # Arriba v2.4.0+ 推荐使用STAR Aligned.sortedByCoord.out.bam
        # 而不是 Chimeric.out.junction 作为主输入。
        # 我们这里提供一个兼容性强的方案，让用户选择
        # 对于最新的Arriba，直接用bam输入更佳
        bam_input = "rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam"
    threads: 8
    resources:
        mem_mb = 16000
    conda:
        "../../envs/arriba.yaml"
    shell:
        """
        # 注意: Arriba v2.4.0+ 推荐直接使用 -x Aligned.sortedByCoord.out.bam
        # 不再需要 -c Chimeric.out.junction.此命令兼容新旧版本。
        arriba \
            -x {params.bam_input} \
            -a {input.ref_genome} \
            -g {input.ref_gtf} \
            -b {input.blacklist} \
            -k {input.known_fusions} \
            -t {input.known_fusions} \
            -p {input.protein_domains} \
            -o {output.fusions} \
            -O {output.fusions_discarded} \
            > {log} 2>&1
        """