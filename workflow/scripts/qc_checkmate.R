#!/usr/bin/env Rscript

# 检查并加载必要的库
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(pheatmap)
  library(reshape2)
})

# 1. 命令行参数设置
option_list <- list(
  make_option(c("-m", "--matched"), type="character", default=NULL, help="Path to *_matched.txt", metavar="character"),
  make_option(c("-d", "--vcf_dir"), type="character", default=NULL, help="Directory containing VCF files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="checkmate_res", help="Output prefix", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$matched) || is.null(opt$vcf_dir)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

# 2. 绘制相关性热图
plot_correlation_heatmap <- function(matched_file, out_prefix) {
  cat("Generating heatmap...\n")
  data <- read.table(matched_file, sep="\t", stringsAsFactors=FALSE)
  
  # 提取样本名并去除后缀
  samples <- sort(unique(c(data$V1, data$V3)))
  display_names <- gsub(".vcf", "", samples)
  
  # 构建对称矩阵
  mat <- matrix(1.0, nrow=length(samples), ncol=length(samples))
  rownames(mat) <- display_names
  colnames(mat) <- display_names
  
  for(i in 1:nrow(data)) {
    s1 <- gsub(".vcf", "", data$V1[i])
    s2 <- gsub(".vcf", "", data$V3[i])
    r  <- data$V4[i]
    mat[s1, s2] <- r
    mat[s2, s1] <- r
  }
  
  # 绘图
  pheatmap(mat, 
           display_numbers = TRUE, 
           number_format = "%.3f",
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           color = colorRampPalette(c("white", "steelblue"))(100),
           main = "Sample Correlation (NGSCheckMate)",
           filename = paste0(out_prefix, "_heatmap.png"),
           width = 8, height = 7)
}

# 3. 解析 VCF 获取 VAF
parse_vcf_vaf <- function(vcf_path) {
  # 读取 VCF (过滤掉 ## 表头)
  vcf <- read.table(vcf_path, sep="\t", comment.char="#", stringsAsFactors=FALSE)
  colnames(vcf) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  
  # 从 INFO 列提取 DP4
  # 逻辑: 匹配 DP4=n,n,n,n
  extract_vaf <- function(info_str) {
    dp4_str <- sub(".*DP4=([0-9,]+).*", "\\1", info_str)
    if (dp4_str == info_str) return(NA) # 没匹配到
    counts <- as.numeric(unlist(strsplit(dp4_str, ",")))
    ref_cnt <- counts[1] + counts[2]
    alt_cnt <- counts[3] + counts[4]
    return(alt_cnt / (ref_cnt + alt_cnt))
  }
  
  vaf <- sapply(vcf$INFO, extract_vaf)
  res <- data.frame(ID = paste(vcf$CHROM, vcf$POS, sep="_"), VAF = vaf)
  return(res[!is.na(res$VAF), ])
}

# 4. 绘制 VAF 散点图
plot_scatter <- function(vcf_dir, out_prefix) {
  pairs <- list(
    c("dna_normal", "dna_tumor"),
    c("rna_normal", "rna_tumor"),
    c("iso_normal", "iso_tumor")
  )
  
  vcf_files <- list.files(vcf_dir, pattern="\\.vcf$")
  
  for (p in pairs) {
    f1_name <- vcf_files[grep(p[1], vcf_files)]
    f2_name <- vcf_files[grep(p[2], vcf_files)]
    
    if (length(f1_name) == 1 && length(f2_name) == 1) {
      cat(paste("Generating scatter plot for", p[1], "vs", p[2], "...\n"))
      vaf1 <- parse_vcf_vaf(file.path(vcf_dir, f1_name))
      vaf2 <- parse_vcf_vaf(file.path(vcf_dir, f2_name))
      
      merged <- merge(vaf1, vaf2, by="ID")
      if (nrow(merged) > 10) {
        p_obj <- ggplot(merged, aes(x=VAF.x, y=VAF.y)) +
          geom_point(alpha=0.4, color="darkblue", size=1) +
          geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
          xlim(0, 1) + ylim(0, 1) +
          labs(title=paste("VAF Correlation:", p[1], "vs", p[2]),
               subtitle=paste("Number of SNPs:", nrow(merged)),
               x=paste("VAF", p[1]), y=paste("VAF", p[2])) +
          theme_minimal()
        
        ggsave(paste0(out_prefix, "_scatter_", p[1], "_vs_", p[2], ".png"), p_obj, width=6, height=6)
      }
    }
  }
}

# 执行
plot_correlation_heatmap(opt$matched, opt$out)
plot_scatter(opt$vcf_dir, opt$out)
cat("All plots generated successfully.\n")