#!/usr/bin/env Rscript

# 加载库
suppressPackageStartupMessages({
    library(argparse)
    library(facets)
    library(ggplot2)
})

# 参数解析
parser <- ArgumentParser(description = "Run FACETS analysis")
parser$add_argument("--input", required=TRUE, help="Path to snp-pileup CSV.gz file")
parser$add_argument("--out_prefix", required=TRUE, help="Output prefix for files")
parser$add_argument("--cval", type="integer", default=100, help="Critical value for segmentation")
parser$add_argument("--min_nhet", type="integer", default=500, help="Minimum number of heterozygous SNPs required")
parser$add_argument("--ndepth", type="integer", default=20, help="Min depth in normal sample")

args <- parser$parse_args()

# 1. 读取数据
cat(sprintf("Reading data from %s...\n", args$input))
# 检查文件是否存在
if (!file.exists(args$input)) {
    stop(paste("Input file not found:", args$input))
}
data <- readSnpMatrix(args$input)

# 2. 预处理
cat("Preprocessing sample...\n")
set.seed(1234)
# ndepth 设为 20 比 35 更温和，适合覆盖度不均的 WES
xx <- preProcSample(data, ndepth=args$ndepth, het.thresh=0.25)

# === 诊断步骤 ===
# 统计 het=1 且 keep=1 的有效杂合位点
num_het_kept <- sum(xx$pmat[, "het"] == 1 & xx$pmat[, "keep"] == 1)
cat(sprintf("Number of Heterozygous SNPs used (het=1 & keep=1): %d\n", num_het_kept))

if (num_het_kept < args$min_nhet) {
    # 如果位点太少，生成一个空的占位文件防止 Snakemake 报错中断，并退出
    warning("Not enough heterozygous SNPs found. Skipping model fitting.")
    
    # 输出一个空的纯度结果
    write.table(data.frame(Purity=NA, Ploidy=NA, Warning="Low_SNP_Count"), 
                file=paste0(args$out_prefix, ".purity.txt"), 
                sep="\t", quote=FALSE, row.names=FALSE)
    
    # 创建空的图作为一个占位符
    png(paste0(args$out_prefix, ".cnv.png"))
    plot(1, type="n", main="Not enough SNPs for FACETS")
    dev.off()
    
    saveRDS(NULL, paste0(args$out_prefix, ".fit.rds"))
    quit(save="no", status=0) 
}

# 3. 分割 (Segmentation)
cat("Running segmentation...\n")
oo <- procSample(xx, cval=args$cval)

cat(sprintf("Number of segments found: %d\n", nrow(oo$out)))

# 4. 拟合模型
cat("Fitting model...\n")
fit <- emcncf(oo)

# 5. 输出文本结果
purity <- fit$purity
ploidy <- fit$ploidy
cat(sprintf("Result: Purity = %s, Ploidy = %s\n", purity, ploidy))

res_df <- data.frame(
    Sample = basename(args$out_prefix),
    Purity = purity,
    Ploidy = ploidy
)
write.table(res_df, file=paste0(args$out_prefix, ".purity.txt"), sep="\t", quote=FALSE, row.names=FALSE)

# 6. 保存 RDS 对象
saveRDS(fit, paste0(args$out_prefix, ".fit.rds"))

# 7. 绘图
# PNG
png(paste0(args$out_prefix, ".cnv.png"), width=1000, height=1000, res=100)
plotSample(x=oo, emfit=fit)
dev.off()

# PDF
pdf(paste0(args$out_prefix, ".cnv.pdf"))
plotSample(x=oo, emfit=fit)
dev.off()

cat("FACETS analysis completed successfully.\n")