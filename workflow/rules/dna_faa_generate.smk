# (这个规则应该在 filter_vcf_for_phasing 之前)

rule pause_for_manual_review:
    input:
        report_xlsx = "reports/{sample}.somatic_variants_report.xlsx"
    output:
        manual_check_file = "reports/{sample}.somatic_variants_report_manual_check.xlsx"
    run:
        # This python code runs inside the rule
        import os, sys, shutil

        # If the user has already created the file, this rule is considered complete
        if os.path.exists(output.manual_check_file):
            print(f"INFO: Found manually checked report '{output.manual_check_file}'. Resuming workflow.")
        else:
            # If the file does NOT exist, stop the workflow and print instructions
            print(f"--- WORKFLOW PAUSED FOR MANUAL REVIEW ---", file=sys.stderr)
            print(f"\nACTION REQUIRED:", file=sys.stderr)
            print(f"1. Open the generated report: '{input.report_xlsx}'", file=sys.stderr)
            print(f"2. Review the 'Manual_Select' column. Change 'no' to 'yes' for any variants you wish to include.", file=sys.stderr)
            print(f"3. Save the modified file with a new name:", file=sys.stderr)
            print(f"   -> '{output.manual_check_file}'", file=sys.stderr)
            print(f"\nOnce the file is saved with the correct name, re-run the same Snakemake command to resume the pipeline.", file=sys.stderr)
            
            # Exit with an error code to halt Snakemake
            sys.exit(1)

rule filter_vcf_for_phasing:
    """
    Filters the main VCF to create a smaller VCF containing only the variants
    selected for neoantigen generation in the Excel report.
    """
    input:
        report_xlsx = "reports/{sample}.somatic_variants_report_manual_check.xlsx",
        somatic_vcf = "dna/variants/annotated/{sample}.mutect2.vep.vcf.gz"
    output:
        filtered_vcf = "dna/variants/phasing/{sample}.selected_for_phasing.vcf.gz",
        filtered_vcf_tbi = "dna/variants/phasing/{sample}.selected_for_phasing.vcf.gz.tbi"
    log:
        "logs/phasing/{sample}.filter_vcf.log"
    params:
        # --- THE FIX: Prepare the script path here ---
        # We use workflow.source_path to get the correct path relative to this Snakefile.
        # This path is then safely passed to the shell block as a simple parameter.
        script = workflow.source_path("../scripts/filter_vcf_from_report.py")
    conda:
        "../../envs/python.yaml"
    shell:
        """
        python {params.script} \\
            --report_xlsx {input.report_xlsx} \\
            --input_vcf {input.somatic_vcf} \\
            --output_vcf {output.filtered_vcf} \\
            > {log} 2>&1
        """

rule whatshap_phase:
    """
    Performs phasing on the selected variants using tumor and normal BAMs
    to determine cis/trans relationships.
    """
    input:
        vcf = "dna/variants/phasing/{sample}.selected_for_phasing.vcf.gz",
        tumor_bam = get_tumor_bam,    # e.g., dna/aligned/PATIENT1_tumor.recal.bam
        normal_bam = get_normal_bam,  # e.g., dna/aligned/PATIENT1_normal.recal.bam
        ref = config["reference"]["genome"]
    output:
        phased_vcf = "dna/variants/phasing/{sample}.phased.vcf.gz",
        phased_vcf_tbi = "dna/variants/phasing/{sample}.phased.vcf.gz.tbi"
    log:
        "logs/phasing/{sample}.whatshap.log"
    params:
        # --- NEW: Define the tumor sample name for whatshap ---
        # This resolves the ambiguity when using --ignore-read-groups
        tumor_sample_name = "{sample}_tumor"
    conda:
        "../../envs/whatshap.yaml"
    shell:
        """
        whatshap phase \\
            --reference {input.ref} \\
            --output {output.phased_vcf} \\
            --tag=PS \\
            --indels \\
            --ignore-read-groups \\
            --sample {params.tumor_sample_name} \\
            {input.vcf} \\
            {input.tumor_bam} {input.normal_bam} \\
            > {log} 2>&1

        tabix -p vcf {output.phased_vcf}
        """

rule generate_phase_aware_neoantigen_fasta:
    """
    Generates a phase-aware FASTA file of neoantigen peptides.

    This rule uses the manually curated Excel report to select variants and
    the PHASED VCF from whatshap to extract phase information. It reconstructs
    local haplotypes by incorporating 'in-cis' germline variants before
    generating the final reference and mutant peptide sequences.
    """
    input:
        report_xlsx = "reports/{sample}.somatic_variants_report_manual_check.xlsx",
        phased_vcf = "dna/variants/phasing/{sample}.phased.vcf.gz",
        uniprot_fasta =  config["reference"]["uniport_fasta"]
    output:
        fasta = "reports/{sample}.snv.peptides.faa"
    params:
        # Window size for the peptide sequence around the mutation
        # This can be moved to config.yaml if you want it to be configurable
        window_size = 18,
        script = workflow.source_path("../scripts/generate_neopeptide_fasta.py")
    log:
        "logs/report/{sample}.generate_phase_aware_fasta.log"
    conda:
        # Reuses the same python environment as the reporting rule
        "../../envs/python.yaml"
    threads: 1 # The script is single-threaded
    resources:
        mem_mb = 4000,
        time = 30
    shell:
        """
        python {params.script} \\
            --report_xlsx {input.report_xlsx} \\
            --vcf_file {input.phased_vcf} \\
            --uniprot_fasta {input.uniprot_fasta} \\
            --output_fasta {output.fasta} \\
            --window_size {params.window_size} \\
            > {log} 2>&1
        """