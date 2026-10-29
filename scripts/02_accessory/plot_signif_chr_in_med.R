#!/usr/bin/env Rscript
# plot_signif_chr_in_med.R
# 比较每个 Contig 在 med 和 env 组之间的共有 kmer 比例分布差异，绘制带显著性标注的 boxplot 组合图
# 每个染色体一个子图，同一染色体的所有 Contig 合并为一个 boxplot
# 支持按 Species 过滤样本

library(optparse)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(scales)

# 解析命令行参数
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, 
              help="Path to input data file (same as heatmap script)"),
  make_option(c("--sample_anno"), type="character", default=NULL,
              help="Sample annotation file (TSV, first column = sample ID, with MedEnv and Species columns)"),
  make_option(c("--species"), type="character", default=NULL,
              help="Filter samples by Species (e.g., 'solani' or 'falciforme'). If not specified, use all samples."),
  make_option(c("--contig_list"), type="character", default=NULL,
              help="Optional file with list of contigs to analyze (one per line). If not specified, use all contigs."),
  make_option(c("--output_prefix"), type="character", default="chr_med_env_boxplot",
              help="Output file prefix (will output .pdf, .png, and .sig.tsv)"),
  make_option(c("--top_n"), type="numeric", default=20,
              help="Number of most significant chromosomes to plot (default: 20)"),
  make_option(c("--panel_size"), type="numeric", default=1.8,
              help="Size of each subplot panel in inches (default: 1.8)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$sample_anno)) {
  print_help(opt_parser)
  stop("Missing required arguments (--input_file and --sample_anno)", call.=FALSE)
}

# ============================================================
# 1. 读取输入数据
# ============================================================
cat("=== Reading input data ===\n")
d <- read.table(opt$input_file, sep = ",", row.names = 1, header = FALSE)

header_file <- paste0(opt$input_file, ".header")
cat("Reading header file:", header_file, "\n")
header <- read.table(header_file, sep = ",")
colnames(d) <- as.character(header[1, ])

# 移除最后一列
d <- d[, -ncol(d), drop = FALSE]
cat(sprintf("Data matrix: %d Contigs x %d samples\n", nrow(d), ncol(d)))
cat("First 5 column names (samples):", paste(head(colnames(d), 5), collapse = ", "), "\n")

# ============================================================
# 2. 读取染色体分组文件
# ============================================================
cat("\n=== Reading chromosome groups ===\n")
chr_group_file <- paste0(opt$input_file, ".chr_group")
chr_group <- read.table(chr_group_file, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
colnames(chr_group) <- "Chromosome"

chromosomes <- sort(unique(chr_group$Chromosome))
cat(sprintf("Total chromosomes: %d\n", length(chromosomes)))
cat("First 5 chromosomes:", paste(head(chromosomes, 5), collapse = ", "), "\n")

# ============================================================
# 2.5. 可选：读取 contig list 过滤
# ============================================================
if (!is.null(opt$contig_list)) {
  cat("\n=== Reading contig list ===\n")
  cat("Contig list file:", opt$contig_list, "\n")
  contig_list <- readLines(opt$contig_list)
  contig_list <- trimws(contig_list)
  contig_list <- contig_list[nchar(contig_list) > 0]  # 移除空行
  cat(sprintf("Contigs in list: %d\n", length(contig_list)))
  
  # 过滤 chr_group，只保留 list 中的 contigs
  chr_group <- chr_group[rownames(chr_group) %in% contig_list, , drop = FALSE]
  chromosomes <- sort(unique(chr_group$Chromosome))
  cat(sprintf("Chromosomes after contig filter: %d\n", length(chromosomes)))
}

# ============================================================
# 3. 读取样本注释
# ============================================================
cat("\n=== Reading sample annotation ===\n")
anno_df <- read.delim(opt$sample_anno, stringsAsFactors = FALSE)
sample_id_col <- colnames(anno_df)[1]
cat(sprintf("Sample ID column: %s\n", sample_id_col))
cat(sprintf("Annotation columns: %s\n", paste(colnames(anno_df)[-1], collapse = ", ")))

# 检查必需列
if (!"MedEnv" %in% colnames(anno_df)) {
  stop("Error: 'MedEnv' column not found in annotation file", call.=FALSE)
}
if (!"Species" %in% colnames(anno_df)) {
  stop("Error: 'Species' column not found in annotation file", call.=FALSE)
}

# 按 Species 过滤样本
if (!is.null(opt$species)) {
  cat(sprintf("Filtering by Species: %s\n", opt$species))
  anno_df <- anno_df[anno_df$Species == opt$species, ]
  cat(sprintf("Samples after Species filter: %d\n", nrow(anno_df)))
}

# 获取 med 和 env 样本
sample_ids <- anno_df[[sample_id_col]]
med_samples <- sample_ids[anno_df$MedEnv == "med"]
env_samples <- sample_ids[anno_df$MedEnv == "env"]
cat(sprintf("Med samples: %d, Env samples: %d\n", length(med_samples), length(env_samples)))

# 检查样本匹配情况
common_med <- intersect(med_samples, colnames(d))
common_env <- intersect(env_samples, colnames(d))
cat(sprintf("\nCommon med samples in data: %d / %d\n", length(common_med), length(med_samples)))
cat(sprintf("Common env samples in data: %d / %d\n", length(common_env), length(env_samples)))

if (length(common_med) == 0 || length(common_env) == 0) {
  stop("Error: No common samples between annotation and data matrix", call.=FALSE)
}

# ============================================================
# 4. 构建绘图数据框
# ============================================================
cat("\n=== Building plot data ===\n")

global_max <- max(d, na.rm = TRUE)
cat(sprintf("Global max kmer proportion: %.4f\n", global_max))

all_data <- data.frame(
  chromosome = character(),
  contig = character(),
  value = numeric(),
  group = character(),
  stringsAsFactors = FALSE
)

sig_annotations <- data.frame(
  chromosome = character(),
  p_value = numeric(),
  y_position = numeric(),
  stringsAsFactors = FALSE
)

for (chr in chromosomes) {
  # 获取该染色体的所有 Contigs
  contigs <- rownames(chr_group)[chr_group$Chromosome == chr]
  contigs_in_data <- intersect(contigs, rownames(d))
  
  if (length(contigs_in_data) == 0) {
    cat(sprintf("  Skip %s: no contigs in data\n", chr))
    next
  }
  
  # 提取该染色体所有 Contigs 的数据
  chr_data <- d[contigs_in_data, , drop = FALSE]
  
  # 获取 med 和 env 样本的列
  med_cols <- intersect(common_med, colnames(chr_data))
  env_cols <- intersect(common_env, colnames(chr_data))
  
  if (length(med_cols) == 0 || length(env_cols) == 0) {
    cat(sprintf("  Skip %s: no med/env columns\n", chr))
    next
  }
  
  # 提取 med 和 env 数据，转为长格式
  med_vals <- as.vector(as.matrix(chr_data[, med_cols, drop = FALSE]))
  env_vals <- as.vector(as.matrix(chr_data[, env_cols, drop = FALSE]))
  
  chr_plot_data <- data.frame(
    chromosome = chr,
    contig = rep(contigs_in_data, each = length(med_cols) + length(env_cols)),
    value = c(med_vals, env_vals),
    group = c(rep("med", length(med_vals)), rep("env", length(env_vals))),
    stringsAsFactors = FALSE
  )
  
  # 移除 NA
  chr_plot_data <- chr_plot_data[!is.na(chr_plot_data$value), ]
  
  if (nrow(chr_plot_data) == 0) {
    cat(sprintf("  Skip %s: no valid data after NA removal\n", chr))
    next
  }
  
  all_data <- rbind(all_data, chr_plot_data)
  
  # 统计检验
  med_v <- chr_plot_data$value[chr_plot_data$group == "med"]
  env_v <- chr_plot_data$value[chr_plot_data$group == "env"]
  
  test_result <- tryCatch(
    wilcox.test(med_v, env_v),    
    error = function(e) NULL
  )
  
  if (!is.null(test_result)) {
    p_value <- test_result$p.value
    sig_annotations <- rbind(sig_annotations, data.frame(
      chromosome = chr,
      p_value = p_value,
      y_position = global_max * 1.15,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("  %s: %d contigs, %d med vals, %d env vals, p=%.4g\n", 
                chr, length(contigs_in_data), length(med_v), length(env_v), p_value))
  }
}

if (nrow(all_data) == 0) {
  stop("Error: No valid data for plotting", call.=FALSE)
}

cat(sprintf("\nTotal data points: %d, Chromosomes tested: %d\n", 
            nrow(all_data), nrow(sig_annotations)))

# ============================================================
# 5. 只保留最显著的前 N 个染色体（cutoff p < 0.01）
# ============================================================
sig_chrs <- sig_annotations$chromosome[sig_annotations$p_value < 0.01]
cat(sprintf("Significant chromosomes (p < 0.01): %d\n", length(sig_chrs)))

if (length(sig_chrs) == 0) {
  stop("Error: No chromosomes with significant differences (p < 0.01)", call.=FALSE)
}

# 先按 p 值排序，再取最显著的前 N 个
sig_chrs_sorted <- sig_chrs[order(sig_annotations$p_value[match(sig_chrs, sig_annotations$chromosome)])]
if (length(sig_chrs_sorted) > opt$top_n) {
  sig_chrs_sorted <- sig_chrs_sorted[1:opt$top_n]
  cat(sprintf("Limiting to top %d most significant chromosomes\n", opt$top_n))
}

# 过滤数据
all_data <- all_data[all_data$chromosome %in% sig_chrs_sorted, ]
sig_annotations_plot <- sig_annotations[sig_annotations$chromosome %in% sig_chrs_sorted, ]

# 设置因子顺序（按显著性 p 值排序）
sig_order <- sig_annotations_plot$chromosome[order(sig_annotations_plot$p_value)]
# 显示层换名：图里统一显示 clinical / non-clinical；源数据与内部统计仍用 med / env
all_data$group_disp <- ifelse(all_data$group == "med", "clinical", "non-clinical")
all_data$group_disp <- factor(all_data$group_disp, levels = c("clinical", "non-clinical"))
all_data$chromosome <- factor(all_data$chromosome, levels = sig_order)
sig_annotations_plot$chromosome <- factor(sig_annotations_plot$chromosome, levels = sig_order)

n_plot <- length(sig_order)
ncols <- ceiling(sqrt(n_plot))
nrows <- ceiling(n_plot / ncols)
cat(sprintf("Plotting %d chromosomes in %d x %d grid\n", n_plot, ncols, nrows))

cat("Chromosomes in plot (ordered by p-value):\n")
for (i in seq_along(sig_order)) {
  chr <- sig_order[i]
  p_val <- sig_annotations$p_value[sig_annotations$chromosome == chr]
  cat(sprintf("  %d. %s (p=%.4g)\n", i, chr, p_val))
}

# ============================================================
# 6. 绘图
# ============================================================
cat("\n=== Generating plot ===\n")

p <- ggplot(all_data, aes(x = group_disp, y = value, fill = group_disp)) +
  geom_boxplot(outlier.size = 0.3, width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 0.2, alpha = 0.3) +
  geom_signif(
    data = sig_annotations_plot,
    aes(xmin = "clinical", xmax = "non-clinical", annotations = ifelse(p_value < 0.001, "***",
                                                   ifelse(p_value < 0.01, "**",
                                                   ifelse(p_value < 0.05, "*", "ns"))),
        y_position = y_position),
    manual = TRUE,
    size = 0.4,
    textsize = 3.5,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("clinical" = "#E41A1C", "non-clinical" = "#377EB8")) +
  scale_y_continuous(
    limits = c(0, global_max * 1.3),
    labels = percent_format(accuracy = 1),
    breaks = pretty(c(0, global_max), n = 4)
  ) +
  facet_wrap(~ chromosome, ncol = ncols) +
  labs(x = NULL, y = "Kmer Proportion", fill = "Clinical status") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_rect(fill = "white", color = "grey80"),
    axis.text.x = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 6),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    plot.margin = margin(5, 5, 5, 5, unit = "pt")
  )

# ============================================================
# 7. 输出显著性结果 TSV
# ============================================================
sig_tsv <- paste0(opt$output_prefix, ".sig.tsv")
sig_output <- sig_annotations[, c("chromosome", "p_value")]
sig_output$significance <- ifelse(sig_output$p_value < 0.001, "***",
                           ifelse(sig_output$p_value < 0.01, "**",
                           ifelse(sig_output$p_value < 0.05, "*", "ns")))
sig_output <- sig_output[order(sig_output$p_value), ]
write.table(sig_output, file = sig_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Significance results (all %d chromosomes) saved to: %s\n", nrow(sig_output), sig_tsv))

# ============================================================
# 8. 计算输出尺寸并保存绘图
# ============================================================
# 每个子图固定为正方形，总图尺寸根据子图数量和网格计算
panel_size <- opt$panel_size
# 估算边距（左侧y轴标签 + 右侧图例 + 顶部标题）
margin_x <- 1.5  # 英寸
margin_y <- 1.2  # 英寸

out_width <- ncols * panel_size + margin_x
out_height <- nrows * panel_size + margin_y

output_pdf <- paste0(opt$output_prefix, ".pdf")
output_png <- paste0(opt$output_prefix, ".png")

cat(sprintf("Output dimensions: %.1f x %.1f inches (%d panels in %d x %d grid)\n", 
            out_width, out_height, n_plot, ncols, nrows))

cat(sprintf("Saving PDF to %s...\n", output_pdf))
pdf(output_pdf, width = out_width, height = out_height)
print(p)
dev.off()

cat(sprintf("Saving PNG to %s...\n", output_png))
png(output_png, width = out_width, height = out_height, units = "in", res = 300)
print(p)
dev.off()

cat("Done!\n")
