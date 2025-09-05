#!/bin/bash

# 创建fastq2neoantigen项目目录结构
# 运行方式: bash create_project_structure.sh

# 设置项目根目录
PROJECT_ROOT=$(pwd)
echo "Creating directory structure in: $PROJECT_ROOT"

# 创建主目录
mkdir -p config
mkdir -p workflow/rules
mkdir -p resources/{genome,annotation,databases}
mkdir -p scripts/utils
mkdir -p benchmarks

# 创建日志目录
mkdir -p logs/{qc,alignment,variants}

# 创建结果目录结构
# QC结果
mkdir -p results/00_qc/{fastqc,fastp}

# 比对结果
mkdir -p results/01_alignment/{dna,rna}

# 变异检测结果
mkdir -p results/02_variants/snv_indel/{raw,filtered,annotated}
mkdir -p results/02_variants/{cnv,sv}

# RNA分析结果
mkdir -p results/03_rna_analysis/{expression,fusion,splicing}

# 新抗原预测结果
mkdir -p results/04_neoantigen

# 免疫相关分析
mkdir -p results/05_immune

# 综合报告
mkdir -p results/99_report/{tables,figures}

# 创建空的配置文件和脚本占位文件
touch config/config.yaml
touch config/samples.csv
touch config/resources.yaml
touch workflow/Snakefile
touch scripts/setup.sh
touch scripts/vcf_to_table.py

# 创建规则文件占位符
RULES=("common" "qc" "mapping" "preprocessing" "variant_calling" "variant_filtering" "annotation" "report")
for rule in "${RULES[@]}"; do
  touch "workflow/rules/${rule}.smk"
done

echo "Directory structure created successfully!"
echo "Next steps:"
echo "1. Edit config/config.yaml to customize your pipeline settings"
echo "2. Edit config/samples.csv to add your sample information"
echo "3. Run 'bash scripts/setup.sh' to set up the environment"

# 显示创建的目录结构树
if command -v tree &> /dev/null; then
  echo -e "\nProject structure:"
  tree -L 3 --dirsfirst
else
  echo -e "\nInstall 'tree' command to view the directory structure"
  echo "sudo apt-get install tree   # For Debian/Ubuntu"
  echo "sudo yum install tree       # For CentOS/RHEL"
fi

