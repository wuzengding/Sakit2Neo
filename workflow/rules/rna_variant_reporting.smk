# workflow/rules/rna_variant_reporting.smk
#
# 此规则生成RNA体细胞变异的综合Excel报告。
# 与DNA报告不同，它不需要外部RNA计数文件，因为Mutect2 VCF本身来源于RNA。

rule report_rna_somatic_variants:
    """
    Parses VEP-annotated RNA VCFs to create a final Excel report.
    Logic adapted for RNA-seq data characteristics (variable depth, no external rescue needed).
    """
    input:
        # --- VCF Files (from RNA pipeline) ---
        somatic_vcf = "rna/variants/annotated/{sample}.mutect2.vep.vcf.gz",

        # --- DYNAMIC INPUT: Conditionally include VarScan2 VCF ---
        varscan_vcf = lambda wildcards:
            f"rna/variants/annotated/{wildcards.sample}.varscan2.vep.vcf.gz"
            if config["variant_callers"].get("varscan2", False) else [],

        # --- Reference Files ---
        gtf = config["reference"]["gtf"],
        fasta = config["reference"]["genome"],
        cancer_genes = config["reference"]["cancer_gene"]
    output:
        # The final, polished Excel report for RNA
        xlsx_report_filter_rna = "reports/{sample}.rna.Filtered_variants.xlsx",
        xlsx_report_Somatic_rna = "reports/{sample}.rna.Somatic_NeoPeptides.xlsx"
    log:
        "logs/reporting/{sample}.rna_somatic_report.log"
    params:
        sample_id = "{sample}",
        # Load tiering parameters (shared with DNA or specific for RNA if needed)
        somatic_tiering_params = config["somatic_tiering"]
    conda:
        "../../envs/python.yaml"
    threads: 1
    resources:
        mem_mb = 8000,
        time = 60
    script:
        "../scripts/generate_rna_somatic_report.py"