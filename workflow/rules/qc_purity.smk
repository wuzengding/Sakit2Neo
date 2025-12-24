# workflow/rules/qc_purity.smk

# 辅助函数：根据 sample 获取对应的 Normal BAM
def get_normal_bam_for_facets(wildcards):
    # 假设 normal 的 sampletype 命名就是 'normal'
    # 如果你的 csv 中 normal 有其他命名，这里需要通过 pandas 查询
    return f"dna/aligned/{wildcards.sample}_normal.recal.bam"

# 辅助函数：根据 sample 获取对应的 Tumor BAM
def get_tumor_bam_for_facets(wildcards):
    # 假设 tumor 的 sampletype 命名就是 'tumor'
    return f"dna/aligned/{wildcards.sample}_tumor.recal.bam"

# -----------------------------------------------
#  Step 1: 运行 SNP-pileup
#  功能：统计 Normal 和 Tumor 在已知 SNP 位点上的覆盖度
# -----------------------------------------------
rule facets_snp_pileup:
    input:
        normal_bam = get_normal_bam_for_facets,
        tumor_bam = get_tumor_bam_for_facets,
        # 必须确保这两个 BAM 的 index 存在，虽然 snp-pileup 可能不强求，但这是好习惯
        normal_bai = lambda w: get_normal_bam_for_facets(w) + ".bai", 
        tumor_bai = lambda w: get_tumor_bam_for_facets(w) + ".bai",
        vcf = config["reference"]["dbsnp"]
    output:
        # snp-pileup 自动会在输出文件名后加 .gz，所以这里我们在 shell 中指定前缀
        # 但 snakemake 需要监控完整的文件名
        pileup = "dna/purity/{sample}_snp.gz"
    log:
        "logs/facets/{sample}.snp_pileup.log"
    conda:
        "../../envs/facets.yaml"
    params:
        # snp-pileup 的输出前缀 (带.gz)
        out_prefix = "dna/purity/{sample}_snp",
        min_map_quality = 20,
        min_base_quality = 15,
        pseudo_snps = 100,
        min_depth = "25,0" # Normal需大于25X, Tumor无限制(0)以保留LOH信息
    threads: 1
    resources: 
        mem_mb = 8000, 
        time = 120
    shell:
        """
        snp-pileup -g \
            -q {params.min_base_quality} \
            -Q {params.min_map_quality} \
            -P {params.pseudo_snps} \
            -r {params.min_depth} \
            {input.vcf} \
            {params.out_prefix} \
            {input.normal_bam} \
            {input.tumor_bam} \
            > {log} 2>&1
        """

# -----------------------------------------------
#  Step 2: 运行 FACETS R 脚本
#  功能：计算纯度、倍体并绘图
# -----------------------------------------------
rule facets_run:
    input:
        pileup = "dna/purity/{sample}_snp.gz",
        facets_script = workflow.source_path("../scripts/facets_run.R")
    output:
        plot_png = "dna/purity/{sample}.cnv.png",
        plot_pdf = "dna/purity/{sample}.cnv.pdf",
        txt_res  = "dna/purity/{sample}.purity.txt",
        rds      = "dna/purity/{sample}.fit.rds"
    log:
        "logs/facets/{sample}.run.log"
    conda:
        "../../envs/facets.yaml"
    params:
        out_prefix = "dna/purity/{sample}", # 传给 R 脚本的前缀
        cval = 100, # 分割阈值，WES常用100-150
        ndepth = 20 # 正常样本最小深度
    resources: 
        mem_mb = 4000
    shell:
        """
        Rscript {input.facets_script} \
            --input {input.pileup} \
            --out_prefix {params.out_prefix} \
            --cval {params.cval} \
            --ndepth {params.ndepth} \
            > {log} 2>&1
        """