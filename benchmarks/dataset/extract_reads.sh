#!/bin/bash

# ====== 输入路径配置 ====== #
# BAM文件路径
tumor_bam="/mnt/user/wzd/05.pipeline_dev/ScanNeo2/results/css006-1/dnaseq/align/dna_tumor_aligned_BWA.bam"
normal_bam="/mnt/user/wzd/05.pipeline_dev/ScanNeo2/results/css006-1/dnaseq/align/dna_normal_aligned_BWA.bam"
rna_tumor_bam="/mnt/user/wzd/05.pipeline_dev/ScanNeo2/results/css006-1/rnaseq/align/rna_tumor_aligned_STAR.bam"
rna_normal_bam="/mnt/user/wzd/09.data_CCS_sdfyy/04.Result.LSTVs/css006_N-1/results/03.Align_NGS/css006_N-1_Aligned.sortedByCoord.out.bam"
rna_pacbio_normal_bam="/mnt/user/wzd/09.data_CCS_sdfyy/04.Result.LSTVs/css006_N-1/results/02.Cluster_Align/css006_N-1_hq_isoforms_mapped_hg38.sorted.bam"
rna_pacbio_tumor_bam="/mnt/user/wzd/09.data_CCS_sdfyy/04.Result.LSTVs/css006_T-1/results/02.Cluster_Align/css006_T-1_hq_isoforms_mapped_hg38.sorted.bam"

# BED文件路径
regionbed="/mnt/user/wzd/05.pipeline_dev/SakitNeo/dataset/regions.bed"

# 输出目录
output_dir="/mnt/user/wzd/05.pipeline_dev/SakitNeo/dataset"

# ====== 多线程配置 ====== #
threads=20  # 总线程数（根据服务器核心数调整）
max_jobs=6  # 最大并行任务数（避免内存爆炸）

# ====== 核心函数 ====== #
extract_to_fastq() {
    local bam=$1
    local prefix=$2
    local sample_name=$3

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 开始处理: $sample_name | BAM: $(basename $bam)"
    
    # Step1: 提取区间并排序（使用50%线程）
    bedtools intersect -a "$bam" -b "$regionbed" | \
    samtools sort -n -@ $((threads/2)) -o "$output_dir/${prefix}_sorted.bam"
    
    # Step2: 转FASTQ（使用100%线程）
    samtools fastq -@ $threads \
      -1 "$output_dir/${prefix}_R1.fastq.gz" \
      -2 "$output_dir/${prefix}_R2.fastq.gz" \
      -0 /dev/null -s /dev/null -n "$output_dir/${prefix}_sorted.bam"
    
    # 清理临时文件
    rm -f "$output_dir/${prefix}_sorted.bam"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] 完成处理: $sample_name"
}

# ====== 主执行逻辑 ====== #
# 创建输出目录
mkdir -p "$output_dir"

# 任务列表（保证最多并行max_jobs个任务）
job_counter=0
for sample_args in \
  "$tumor_bam dna_tumor DNA_Tumor" \
  "$normal_bam dna_normal DNA_Normal" \
  "$rna_tumor_bam rna_tumor RNA_Tumor" \
  "$rna_normal_bam rna_normal RNA_Normal" \
  "$rna_pacbio_normal_bam pacbio_normal PacBio_Normal" \
  "$rna_pacbio_tumor_bam pacbio_tumor PacBio_Tumor"
do
    ((job_counter++))
    echo "启动任务 $job_counter: $sample_args"
    extract_to_fastq $sample_args &
    
    # 控制并行数量
    if (( job_counter % max_jobs == 0 )); then
        wait  # 等待当前批次任务完成
        echo "已完成一批任务，继续下一批..."
    fi
done

# 等待剩余任务
wait
echo "[$(date '+%Y-%m-%d %H:%M:%S')] 所有样本处理完成！输出到: $output_dir"
