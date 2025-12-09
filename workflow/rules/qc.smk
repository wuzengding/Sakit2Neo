# ---------------------------------------------------------
# 辅助函数：判断某类样本是否存在
# ---------------------------------------------------------
def has_sample_data(sample_id, sample_type, data_type):
    """
    检查 samples dataframe 中是否存在指定的记录
    """
    subset = samples[
        (samples['sample_id'] == sample_id) & 
        (samples['sampletype'] == sample_type) & 
        (samples['data_type'] == data_type)
    ]
    return not subset.empty

# ---------------------------------------------------------
# Input 函数：只返回 Snakemake 需要构建的依赖文件
# ---------------------------------------------------------
def get_qc_input_files(wildcards):
    input_files = []
    
    # 1. DNA Tumor (通常必须有，但为了保险起见也可以检查)
    if has_sample_data(wildcards.sample, 'tumor', 'dna'):
        input_files.extend([
            f"qc/fastp/{wildcards.sample}_tumor_dna.fastp.json",
            f"qc/coverage/{wildcards.sample}_tumor.mosdepth.summary.txt",
            f"qc/picard/{wildcards.sample}_tumor.alignment_metrics.txt",
            f"qc/picard/{wildcards.sample}_tumor.hs_metrics.txt",
            f"qc/picard/{wildcards.sample}_tumor.insert_size_metrics.txt"
        ])

    # 2. DNA Normal
    if has_sample_data(wildcards.sample, 'normal', 'dna'):
        input_files.extend([
            f"qc/fastp/{wildcards.sample}_normal_dna.fastp.json",
            f"qc/coverage/{wildcards.sample}_normal.mosdepth.summary.txt",
            f"qc/picard/{wildcards.sample}_normal.alignment_metrics.txt",
            f"qc/picard/{wildcards.sample}_normal.hs_metrics.txt",
            f"qc/picard/{wildcards.sample}_normal.insert_size_metrics.txt"
        ])

    # 3. RNA Tumor
    if has_sample_data(wildcards.sample, 'tumor', 'rna'):
        input_files.extend([
            f"qc/fastp/{wildcards.sample}_tumor_rna_fastp.json",
            f"qc/rna/{wildcards.sample}_tumor_RNA_metrics.txt"
        ])

    # 4. RNA Normal
    if has_sample_data(wildcards.sample, 'normal', 'rna'):
        input_files.extend([
            f"qc/fastp/{wildcards.sample}_normal_rna_fastp.json",
            f"qc/rna/{wildcards.sample}_normal_RNA_metrics.txt"
        ])
        
    return input_files

# ---------------------------------------------------------
# Params 函数：返回路径或 "None" 字符串传给 Python
# ---------------------------------------------------------
def get_qc_path_or_none(wildcards, sample_type, data_type, file_pattern):
    """
    如果样本存在，返回格式化后的路径；否则返回字符串 "None"
    """
    if has_sample_data(wildcards.sample, sample_type, data_type):
        return file_pattern.format(sample=wildcards.sample)
    else:
        return "None"

rule collect_qc_report:
    input:
        get_qc_input_files, # 动态依赖：只等待存在的文件
        script = workflow.source_path("../scripts/generate_qc_report.py")
    output:
        xlsx = "reports/{sample}.QC.xlsx"
    log:
        "logs/report/{sample}_qc_aggregation.log"
    conda:
        "../../envs/python.yaml"
    params:
        # --- DNA Tumor Params ---
        dt_fastp = lambda w: get_qc_path_or_none(w, 'tumor', 'dna', "qc/fastp/{sample}_tumor_dna.fastp.json"),
        dt_mos   = lambda w: get_qc_path_or_none(w, 'tumor', 'dna', "qc/coverage/{sample}_tumor.mosdepth.summary.txt"),
        dt_align = lambda w: get_qc_path_or_none(w, 'tumor', 'dna', "qc/picard/{sample}_tumor.alignment_metrics.txt"),
        dt_hs    = lambda w: get_qc_path_or_none(w, 'tumor', 'dna', "qc/picard/{sample}_tumor.hs_metrics.txt"),
        #dt_ins   = lambda w: get_qc_path_or_none(w, 'tumor', 'dna', "qc/picard/{sample}_tumor.insert_size_metrics.txt"),
        
        # --- DNA Normal Params ---
        dn_fastp = lambda w: get_qc_path_or_none(w, 'normal', 'dna', "qc/fastp/{sample}_normal_dna.fastp.json"),
        dn_mos   = lambda w: get_qc_path_or_none(w, 'normal', 'dna', "qc/coverage/{sample}_normal.mosdepth.summary.txt"),
        dn_align = lambda w: get_qc_path_or_none(w, 'normal', 'dna', "qc/picard/{sample}_normal.alignment_metrics.txt"),
        dn_hs    = lambda w: get_qc_path_or_none(w, 'normal', 'dna', "qc/picard/{sample}_normal.hs_metrics.txt"),
        #dn_ins   = lambda w: get_qc_path_or_none(w, 'normal', 'dna', "qc/picard/{sample}_normal.insert_size_metrics.txt"),

        # --- RNA Tumor Params ---
        rt_fastp   = lambda w: get_qc_path_or_none(w, 'tumor', 'rna', "qc/fastp/{sample}_tumor_rna_fastp.json"),
        rt_metrics = lambda w: get_qc_path_or_none(w, 'tumor', 'rna', "qc/rna/{sample}_tumor_RNA_metrics.txt"),
        
        # --- RNA Normal Params ---
        rn_fastp   = lambda w: get_qc_path_or_none(w, 'normal', 'rna', "qc/fastp/{sample}_normal_rna_fastp.json"),
        rn_metrics = lambda w: get_qc_path_or_none(w, 'normal', 'rna', "qc/rna/{sample}_normal_RNA_metrics.txt"),
        
    shell:
        """
        python {input.script} \\
            --dna-tumor-fastp {params.dt_fastp} \\
            --dna-tumor-mosdepth {params.dt_mos} \\
            --dna-tumor-align {params.dt_align} \\
            --dna-tumor-hs {params.dt_hs} \\
            --dna-normal-fastp {params.dn_fastp} \\
            --dna-normal-mosdepth {params.dn_mos} \\
            --dna-normal-align {params.dn_align} \\
            --dna-normal-hs {params.dn_hs} \\
            --rna-tumor-fastp {params.rt_fastp} \\
            --rna-tumor-metrics {params.rt_metrics} \\
            --rna-normal-fastp {params.rn_fastp} \\
            --rna-normal-metrics {params.rn_metrics} \\
            --output {output.xlsx} \\
            > {log} 2>&1
        """