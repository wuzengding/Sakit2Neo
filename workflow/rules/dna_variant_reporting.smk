# workflow/rules/variant_reporting.smk
#
# This rule generates a comprehensive, multi-sheet Excel report for somatic variants.
# It integrates VCF annotations, RNA expression support, and cross-caller validation
# to produce a human-readable summary for variant review.

rule report_somatic_variants:
    """
    Parses VEP-annotated VCFs and integrates multiple data sources to create a
    final Excel report for somatic variant prioritization. It uses a sophisticated
    Python script to handle complex logic like exon-aware distance calculation
    and multi-tiered somatic variant classification.
    """
    input:
        # --- VCF Files ---
        somatic_vcf = "dna/variants/annotated/{sample}.mutect2.vep.vcf.gz", # Mutect2 is mandatory

        # --- DYNAMIC INPUT: Conditionally include VarScan2 VCF ---
        # Only require this input if varscan2 is enabled in the config.
        # If disabled, this lambda function returns an empty list, making the input optional.
        varscan_vcf = lambda wildcards:
            f"dna/variants/annotated/{wildcards.sample}.varscan2.vep.vcf.gz"
            if config["variant_callers"]["varscan2"] else [],

        # --- Validation & Annotation Data ---
        rna_counts_tumor = "rna/snv_validation/{sample}_tumor_mutect2_ase_counts.tsv",

        # --- Reference Files for PyEnsembl ---
        gtf = config["reference"]["gtf"],
        fasta = config["reference"]["genome"],
        cancer_genes = config["reference"]["cancer_gene"]
    output:
        # The final, polished Excel report
        xlsx_report_filter = "reports/{sample}.Filtered_variants.xlsx",
        xlsx_report_Somatic = "reports/{sample}.Somatic_NeoPeptides.xlsx"
    log:
        "logs/reporting/{sample}.somatic_report.log"
    params:
        # Pass the sample name to the script
        sample_id = "{sample}",
        
        # Load the tiering parameters directly from the config file
        somatic_tiering_params = config["somatic_tiering"]
    conda:
        # This environment must contain: pandas, openpyxl, pysam, pyensembl
        "../../envs/python.yaml"
    threads: 1 # This script is single-threaded
    resources:
        mem_mb = 8000, # Memory for loading VCFs and pyensembl index
        time = 60      # Expected runtime in minutes
    script:
        # Path to the script relative to the rule's file location
        "../scripts/generate_somatic_report.py"