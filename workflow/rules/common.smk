# 常用函数和通用规则
import os
from snakemake.utils import validate
from snakemake.utils import min_version

# 获取染色体列表
def get_chromosomes():
    return [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# 获取样本ID列表
def get_sample_ids():
    return samples['sample_id'].unique().tolist()

# 获取样本的BAM文件路径
def et_bam(wildcards):
    return f"dna/aligned/{wildcards.sample}_{wildcards.sampletype}.recal.bam"

# 获取RNA样本的BAM文件路径
def get_rna_bam(wildcards):
    return f"rna/aligned/{wildcards.sample}_{wildcards.sampletype}.Aligned.sortedByCoord.out.bam"

# 获取PacBio样本的BAM文件路径
def get_pacbio_bam(wildcards):
    return f"pacbio/aligned/{wildcards.sample}_{wildcards.sampletype}.aligned.bam"



# 获取样本的tumor和normal BAM路径
def get_tumor_bam(wildcards):
    return f"dna/aligned/{wildcards.sample}_tumor.recal.bam"

def get_normal_bam(wildcards):
    return f"dna/aligned/{wildcards.sample}_normal.recal.bam"

def get_tumor_bai(wildcards):
    return f"dna/aligned/{wildcards.sample}_tumor.recal.bam.bai"

def get_normal_bai(wildcards):
    return f"dna/aligned/{wildcards.sample}_normal.recal.bam.bai"

def get_tumor_bam_scattered(wildcards):
    return f"dna/aligned/temp/{wildcards.sample}_tumor_{wildcards.chromosome}.recal.bam"

def get_normal_bam_scattered(wildcards):
    return f"dna/aligned/temp/{wildcards.sample}_normal_{wildcards.chromosome}.recal.bam"

def get_tumor_bai_scattered(wildcards):
    return f"dna/aligned/temp/{wildcards.sample}_tumor_{wildcards.chromosome}.recal.bam.bai"

def get_normal_bai_scattered(wildcards):
    return f"dna/aligned/temp/{wildcards.sample}_normal_{wildcards.chromosome}.recal.bam.bai"

# 获取成对的正常/肿瘤样本
def get_paired_samples():
    result = []
    all_samples = set(samples['sample_id'])
    for sample in all_samples:
        types = set(samples[samples['sample_id'] == sample]['sample_type'])
        if 'normal' in types and 'tumor' in types:
            result.append(sample)
    return result


# 创建临时目录
rule create_temp_dir:
    output:
        temp(directory("temp"))
    shell:
        "mkdir -p {output}"

