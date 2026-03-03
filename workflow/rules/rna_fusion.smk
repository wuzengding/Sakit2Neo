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

# 辅助函数，如果你的 illumina BAM 文件路径不是简单的 `{sample}_{sampletype}.bam`，
# 并且 Arriba 实际上是处理的 STAR 的 Aligned.sortedByCoord.out.bam。
# 假设你的 illumina bam 路径是 `rna/align/{sample}_{sampletype}/Aligned.sortedByCoord.out.bam`
def get_illumina_aligned_bam(wildcards):
    return f"rna/align/{wildcards.sample}_{wildcards.sampletype}/Aligned.sortedByCoord.out.bam"

# rule for FusionSeeker integrated analysis
#rule fusionseeker_integrated_fusion:
#    """
#    使用 FusionSeeker 整合二代短读长 (STAR BAM) 和三代长读长 (PacBio minimap2 BAM)
#    数据进行融合基因检测。
#    """
#    input:
#        # 二代数据输入：使用 STAR 的比对 BAM 文件
#        illumina_bam = get_illumina_aligned_bam,
#        # 三代数据输入：使用 PacBio 的 Minimap2 比对 BAM 文件
#        pacbio_bam = "pacbio/aligned/{sample}_{sampletype}.aligned.sorted.bam",
#        
#        # 参考基因组 FASTA
#        ref_fasta = config["reference"]["genome"],
#        
#        # FusionSeeker 可能需要的额外参考文件和数据库
#        fusionseeker_db = config["fusionseeker"]["db_path"],
#        # ref_gtf = config["fusionseeker"]["reference_gtf"],
#        # star_index_path = config["fusionseeker"]["star_index"],
#        
#    output:
#        # FusionSeeker 的主要融合结果文件
#        fusions_tsv = "rna/fusion/{sample}_{sampletype}.fusionseeker.tsv",
#        # FusionSeeker 可能有其他输出，例如过滤后的结果、日志文件等，根据实际情况添加
#        # 例如: filtered_fusions = "rna/fusion/{sample}_{sampletype}.fusionseeker.filtered.tsv"
#    
#    params:
#        # FusionSeeker 的额外命令行参数，如线程、内存限制、特定模式等
#        # 这里需要根据 FusionSeeker 的实际命令行参数进行详细配置
#        extra_options = "" 
#    log:
#        "logs/fusionseeker_integrated/{sample}_{sampletype}.log"
#    threads: 16 # FusionSeeker 可能会消耗较多CPU，根据你的资源调整
#    resources:
#        mem_mb = 32000 # 32GB内存，根据数据量和FusionSeeker要求调整
#    conda:
#        "../../envs/fusionseeker.yaml" # 使用为 FusionSeeker 新建的环境
#    shell:
#        """
#        # 这里的 'fusionseeker_command' 需要替换为 FusionSeeker 实际的执行命令 
#        # 例如，如果 FusionSeeker 是一个 Python 脚本，可能是 'python path/to/FusionSeeker.py'
#        # 如果是可执行文件，可能是直接的命令名
#        
#        # 以下是基于常见 FusionSeeker 命令行结构的示例，具体请查阅 FusionSeeker 文档！
#        fusionseeker_command \
#            --genome {input.ref_fasta} \
#            --illumina-bam {input.illumina_bam} \
#            --pacbio-bam {input.pacbio_bam} \
#            --db {input.fusionseeker_db} \
#            --threads {threads} \
#            --output {output.fusions_tsv} \
#            {params.extra_options} \
#            > {log} 2>&1
#        """