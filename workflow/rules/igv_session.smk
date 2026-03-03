import os

# 核心路径转换函数
def convert_path_to_url(local_path, igv_base_url, work_dir):
    """
    将本地绝对路径转换为 HTTP URL。
    替换逻辑：将 "/home/bioinfo/" 替换为 config 中的 igv_url
    """
    if not local_path or local_path == "MISSING":
        return ""
    
    # 1. 获取文件的绝对路径
    # 如果路径已经是绝对路径(如参考基因组)，直接使用
    # 如果是相对路径(Snakemake output)，则拼接到 work_dir (config['outdir'])
    if os.path.isabs(local_path):
        abs_path = local_path
    else:
        abs_path = os.path.join(work_dir, local_path)
    
    # 2. 检查路径是否包含 /home/bioinfo/ (根据你的服务器结构)
    prefix_to_replace = "/home/bioinfo/"
    
    if prefix_to_replace in abs_path:
        # 执行替换
        # 确保 url 以 / 结尾，prefix 也以 / 结尾，直接替换即可
        # 例如: http://ip:port/ + 09.data_CCS_sdfyy/...
        return abs_path.replace(prefix_to_replace, igv_base_url)
    else:
        # 如果路径不在 /home/bioinfo 下，可能无法通过web访问，打印警告并返回原路径
        print(f"Warning: Path {abs_path} does not start with {prefix_to_replace}")
        return abs_path

# 获取输入逻辑保持不变
def get_igv_dna_inputs(wildcards):
    inputs = {}
    if config["modules"]["wes_processing"]:
        inputs["dna_tumor"] = f"dna/aligned/{wildcards.sample}_tumor.recal.bam"
        inputs["dna_normal"] = f"dna/aligned/{wildcards.sample}_normal.recal.bam"
    if config["modules"]["rna_processing"]:
        inputs["rna_tumor"] = f"rna/align/{wildcards.sample}_tumor_dedup.bam"
    return inputs

def get_igv_rna_inputs(wildcards):
    inputs = {}
    if config["modules"]["rna_processing"]:
        inputs["rna_tumor_bam"] = f"rna/align/{wildcards.sample}_tumor_dedup.bam"
        inputs["rna_normal_bam"] = f"rna/align/{wildcards.sample}_normal_dedup.bam"
    if config["modules"]["rna_assembly"]:
        inputs["rna_tumor_gtf"] = f"rna/assembly/{wildcards.sample}_tumor.assembly.gtf.gz"
        inputs["rna_normal_gtf"] = f"rna/assembly/{wildcards.sample}_normal.assembly.gtf.gz"
    if config["modules"]["pacbio_processing"] and has_pacbio_data():
        inputs["pacbio_tumor"] = f"pacbio/aligned/{wildcards.sample}_tumor.aligned.sorted.bam"
        inputs["pacbio_normal"] = f"pacbio/aligned/{wildcards.sample}_normal.aligned.sorted.bam"
    return inputs

rule generate_igv_dna_xml:
    input:
        unpack(get_igv_dna_inputs)
    output:
        "reports/igv/{sample}_DNA.xml"
    params:
        template = config["reference"]["igv_dna_xml"],
        # 从config获取参考文件绝对路径
        genome_path = config["reference"]["genome"],
        gtf_path = config["reference"]["gtfgz"],
        # 从config获取URL前缀
        igv_url = config["reference"]["igv_url"],
        # 获取工作目录绝对路径
        work_dir = config["outdir"]
    run:
        import os
        
        # 1. 转换参考文件路径为URL
        genome_url = convert_path_to_url(params.genome_path, params.igv_url, params.work_dir)
        gtf_url = convert_path_to_url(params.gtf_path, params.igv_url, params.work_dir)

        # 2. 转换样本文件路径为URL
        dna_t_path = input.get("dna_tumor", "MISSING")
        dna_n_path = input.get("dna_normal", "MISSING")
        rna_t_path = input.get("rna_tumor", "MISSING")

        dna_t_url = convert_path_to_url(dna_t_path, params.igv_url, params.work_dir)
        dna_n_url = convert_path_to_url(dna_n_path, params.igv_url, params.work_dir)
        rna_t_url = convert_path_to_url(rna_t_path, params.igv_url, params.work_dir)

        # 3. 读取并替换模板
        if not os.path.exists(params.template):
            raise FileNotFoundError(f"Template not found: {params.template}")
            
        with open(params.template, 'r') as f:
            content = f.read()

        replacements = {
            "__SAMPLE__": wildcards.sample,
            "__GENOME_URL__": genome_url,
            "__GTF_URL__": gtf_url,
            "__DNA_TUMOR_BAM_URL__": dna_t_url,
            "__DNA_NORMAL_BAM_URL__": dna_n_url,
            "__RNA_TUMOR_BAM_URL__": rna_t_url
        }

        for key, value in replacements.items():
            content = content.replace(key, value)

        with open(output[0], 'w') as out:
            out.write(content)

rule generate_igv_rna_xml:
    input:
        unpack(get_igv_rna_inputs)
    output:
        "reports/igv/{sample}_RNA.xml"
    params:
        template = config["reference"]["igv_rna_xml"],
        genome_path = config["reference"]["genome"],
        gtf_path = config["reference"]["gtfgz"],
        igv_url = config["reference"]["igv_url"],
        work_dir = config["outdir"]
    run:
        import os
        
        # 1. 转换参考文件
        genome_url = convert_path_to_url(params.genome_path, params.igv_url, params.work_dir)
        gtf_url = convert_path_to_url(params.gtf_path, params.igv_url, params.work_dir)

        # 2. 转换样本文件
        rna_t_url = convert_path_to_url(input.get("rna_tumor_bam", "MISSING"), params.igv_url, params.work_dir)
        rna_n_url = convert_path_to_url(input.get("rna_normal_bam", "MISSING"), params.igv_url, params.work_dir)
        rna_t_gtf_url = convert_path_to_url(input.get("rna_tumor_gtf", "MISSING"), params.igv_url, params.work_dir)
        rna_n_gtf_url = convert_path_to_url(input.get("rna_normal_gtf", "MISSING"), params.igv_url, params.work_dir)
        
        pb_t_path = input.get("pacbio_tumor", "")
        pb_n_path = input.get("pacbio_normal", "")
        pb_t_url = convert_path_to_url(pb_t_path, params.igv_url, params.work_dir) if pb_t_path else ""
        pb_n_url = convert_path_to_url(pb_n_path, params.igv_url, params.work_dir) if pb_n_path else ""

        # 3. 替换模板
        if not os.path.exists(params.template):
             raise FileNotFoundError(f"Template not found: {params.template}")

        with open(params.template, 'r') as f:
            content = f.read()

        replacements = {
            "__SAMPLE__": wildcards.sample,
            "__GENOME_URL__": genome_url,
            "__GTF_URL__": gtf_url,
            "__RNA_TUMOR_BAM_URL__": rna_t_url,
            "__RNA_NORMAL_BAM_URL__": rna_n_url,
            "__RNA_TUMOR_GTF_URL__": rna_t_gtf_url,
            "__RNA_NORMAL_GTF_URL__": rna_n_gtf_url,
            "__PACBIO_TUMOR_BAM_URL__": pb_t_url,
            "__PACBIO_NORMAL_BAM_URL__": pb_n_url
        }

        for key, value in replacements.items():
            content = content.replace(key, value)
            
        # 清理可能不存在的 PacBio Resource 标签
        # 如果 URL 为空，移除整行 Resource 标签以避免 IGV 报错
        if not pb_t_url:
            content = content.replace('<Resource path="" type="bam"/>', '')
            # 同时也尝试移除带占位符的行（防止某些情况下未替换成功）
            content = content.replace('<Resource path="__PACBIO_TUMOR_BAM_URL__" type="bam"/>', '')
        if not pb_n_url:
            content = content.replace('<Resource path="" type="bam"/>', '')
            content = content.replace('<Resource path="__PACBIO_NORMAL_BAM_URL__" type="bam"/>', '')
            
        with open(output[0], 'w') as out:
            out.write(content)