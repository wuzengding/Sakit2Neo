# ====================
# Manta结构变异检测规则
# ====================
rule manta_sv_calling:
    input:
        tumor_bam = "dna/aligned/{sample}_tumor.recal.bam",
        tumor_bai = "dna/aligned/{sample}_tumor.recal.bam.bai",
        normal_bam = "dna/aligned/{sample}_normal.recal.bam",
        normal_bai = "dna/aligned/{sample}_normal.recal.bam.bai",
        ref = config["reference"]["genome"],
        regions = config["reference"]["capture_kit_bedgz"]
    output:
        somaticSV = "dna/sv/{sample}_somaticSV.vcf.gz",
        somaticSV_tbi = "dna/sv/{sample}_somaticSV.vcf.gz.tbi",
        diploidSV = "dna/sv/{sample}_diploidSV.vcf.gz",
        diploidSV_tbi = "dna/sv/{sample}_diploidSV.vcf.gz.tbi",
        candidateSV = "dna/sv/{sample}_candidateSV.vcf.gz",
        candidateSV_tbi = "dna/sv/{sample}_candidateSV.vcf.gz.tbi",
        candidateSmallIndels = "dna/sv/{sample}_candidateSmallIndels.vcf.gz",
        candidateSmallIndels_tbi = "dna/sv/{sample}_candidateSmallIndels.vcf.gz.tbi"
    params:
        run_dir = "dna/sv/",
        exe_dir = "dna/sv/results",
        extra = "--exome" if "capture_kit_bedgz" in config["reference"] else "",
        rename_dict = {  # 修改为使用索引方式访问
            "somaticSV": "{sample}_somaticSV",
            "diploidSV": "{sample}_diploidSV",
            "candidateSV": "{sample}_candidateSV",
            "candidateSmallIndels": "{sample}_candidateSmallIndels"
        }
    log:
        "logs/manta/{sample}.log"
    conda:
        "../../envs/manta.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 32000,
        time = lambda wildcards, attempt: attempt * 240
    threads: 8
    shell:
        """
        # 创建工作目录
        rm -rf {params.run_dir}
        mkdir -p {params.run_dir} $(dirname {log})

        # 配置Manta运行参数
        configManta.py \
            --tumorBam {input.tumor_bam} \
            --normalBam {input.normal_bam} \
            --referenceFasta {input.ref} \
            --runDir {params.run_dir} \
            --callRegions {input.regions} \
            {params.extra} \
            --generateEvidenceBam 2> {log}

        # 执行分析
        {params.run_dir}/runWorkflow.py \
            -m local \
            -j {threads} \
            --quiet >> {log} 2>&1

        # 移动并重命名结果文件 (纯shell版本)
        mv {params.run_dir}/results/variants/somaticSV.vcf.gz {output.somaticSV}
        mv {params.run_dir}/results/variants/somaticSV.vcf.gz.tbi {output.somaticSV_tbi}
        mv {params.run_dir}/results/variants/diploidSV.vcf.gz {output.diploidSV}
        mv {params.run_dir}/results/variants/diploidSV.vcf.gz.tbi {output.diploidSV_tbi}
        mv {params.run_dir}/results/variants/candidateSV.vcf.gz {output.candidateSV}
        mv {params.run_dir}/results/variants/candidateSV.vcf.gz.tbi {output.candidateSV_tbi}
        mv {params.run_dir}/results/variants/candidateSmallIndels.vcf.gz {output.candidateSmallIndels}
        mv {params.run_dir}/results/variants/candidateSmallIndels.vcf.gz.tbi {output.candidateSmallIndels_tbi}
        mv  {params.run_dir}/results/*  {params.run_dir}/

        # 清理临时目录
        rm -rf {params.exe_dir}/variants
        rm -rf {params.run_dir}/workspace
        rm -rf {params.run_dir}/workflow*
        """

# ====================
# Manta结果简化规则 (可选)
# ====================
rule process_manta_results:
    input:
        vcf = "dna/sv/{sample}_somaticSV.vcf.gz"
    output:
        filtered = "dna/sv/{sample}_somaticSV.filtered.vcf.gz",
        report = "dna/sv/{sample}_somaticSV.report.txt"
    log:
        "logs/manta/{sample}_filter.log"
    conda:
        "../../envs/manta.yaml"
    resources:
        mem_mb = 8000,
        time = 60
    shell:
        """
        # 基本过滤 (示例：保留PASS变异)
        bcftools view -f PASS -Oz -o {output.filtered} {input.vcf} 2> {log}

        # 生成统计报告
        bcftools stats {output.filtered} > {output.report} 2>> {log}
        """
