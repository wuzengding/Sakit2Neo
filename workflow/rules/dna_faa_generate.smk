# (这个规则应该在 filter_vcf_for_phasing 之前)

rule pause_for_manual_review:
    input:
        #xlsx_report_filter = "reports/{sample}.Filtered_variants.xlsx",
        xlsx_report_Somatic = "reports/{sample}.Somatic_NeoPeptides.xlsx"
    output:
        manual_check_file = "reports/{sample}.Somatic_NeoPeptides_manual_check.xlsx"
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
            print(f"1. Open the generated report: '{input.xlsx_report_Somatic}'", file=sys.stderr)
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
        xlsx_report_Somatic = "reports/{sample}.Somatic_NeoPeptides_manual_check.xlsx",
        somatic_vcf = "dna/variants/annotated/{sample}.mutect2.vep.vcf.gz"
    output:
        filtered_vcf = "dna/variants/phasing/{sample}.mutect2.selected_for_phasing.vcf.gz",
        filtered_vcf_tbi = "dna/variants/phasing/{sample}.mutect2.selected_for_phasing.vcf.gz.tbi"
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
            --xlsx_report_Somatic {input.xlsx_report_Somatic} \\
            --input_vcf {input.somatic_vcf} \\
            --output_vcf {output.filtered_vcf} \\
            > {log} 2>&1
        """

rule whatshap_phase:
    input:
        vcf = "dna/variants/phasing/{sample}.mutect2.selected_for_phasing.vcf.gz",
        tumor_bam = get_tumor_bam,
        ref = config["reference"]["genome"]
    output:
        phased_vcf = "dna/variants/phasing/{sample}.mutect2.phased.vcf.gz",
        phased_vcf_tbi = "dna/variants/phasing/{sample}.mutect2.phased.vcf.gz.tbi"
    log:
        "logs/phasing/{sample}.whatshap.log"
    params:
        tumor_sample_name = "{sample}_tumor"
    conda:
        "../../envs/whatshap.yaml"
    shell:
        """
        # 1. 安全措施：确保 BAM 索引是最新的
        # 如果 snakemake 已经处理了索引，这一步可能很快，但为了解决卡死问题，强制刷新是明智的
        #samtools index {input.tumor_bam}

        # 2. 运行 Whatshap (不使用管道，分步执行)
        # 警告：请根据"第三步"的检查结果，决定是否保留 --ignore-read-groups
        # 如果 BAM @RG SM tag 等于 {params.tumor_sample_name}，建议删除 --ignore-read-groups
        
        whatshap phase \
            --reference {input.ref} \
            --output {output.phased_vcf}.tmp.vcf \
            --tag=PS \
            --indels \
            --ignore-read-groups \
            --sample {params.tumor_sample_name} \
            {input.vcf} \
            {input.tumor_bam} \
            > {log} 2>&1

        # 3. 检查是否成功生成了文件
        if [ ! -s {output.phased_vcf}.tmp.vcf ]; then
            echo "CRITICAL ERROR: Whatshap produced an empty file. Check logs." >> {log}
            # 打印最后几行日志到屏幕以便调试
            tail -n 20 {log}
            exit 1
        fi

        # 4. 压缩并索引
        bgzip -c {output.phased_vcf}.tmp.vcf > {output.phased_vcf}
        tabix -p vcf {output.phased_vcf}

        # 5. 清理
        rm {output.phased_vcf}.tmp.vcf
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
        xlsx_report_Somatic = "reports/{sample}.Somatic_NeoPeptides_manual_check.xlsx",
        phased_vcf = "dna/variants/phasing/{sample}.mutect2.phased.vcf.gz",
        #phased_vcf = "dna/variants/phasing/{sample}.mutect2.selected_for_phasing.vcf.gz",
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
            --xlsx_report_Somatic {input.xlsx_report_Somatic} \\
            --vcf_file {input.phased_vcf} \\
            --uniprot_fasta {input.uniprot_fasta} \\
            --output_fasta {output.fasta} \\
            --window_size {params.window_size} \\
            > {log} 2>&1
        """