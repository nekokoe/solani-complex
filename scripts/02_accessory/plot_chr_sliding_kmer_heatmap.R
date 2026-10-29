#!/usr/bin/env Rscript

options(error = function() traceback(sys.calls(), max.lines = 200))
options(warn = 1)

library(optparse)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL,
              help="Path to input data file (CSV, rows=contigs, cols=samples)", metavar="character"),
  make_option(c("--chromosome_list"), type="character", default=NULL,
              help="File with list of chromosomes to plot (one per line). If not provided, all chromosomes will be plotted.", metavar="character"),
  make_option(c("--exclude_sample"), type="character", default=NULL,
              help="Sample name to exclude", metavar="character"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Output file prefix (will output .pdf and .png for each chromosome)", metavar="character"),
  make_option(c("--width"), type="numeric", default=14,
              help="Plot width in inches (default: 14)", metavar="numeric"),
  make_option(c("--height"), type="numeric", default=4,
              help="Plot height in inches (default: 4)", metavar="numeric"),
  make_option(c("--sample_order"), type="character", default=NULL,
              help="File with sample order (tsv, single column). If not provided, samples will be clustered.", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Missing required argument: --input_file", call.=FALSE)
}

# ============================================================
# 1. 读取数据
# ============================================================
cat("=== Reading input data ===\n")
cat("Input file:", opt$input_file, "\n")

d <- read.table(opt$input_file, sep = ",", row.names = 1, header = FALSE)

header_file <- paste0(opt$input_file, ".header")
cat("Header file:", header_file, "\n")
header <- read.table(header_file, sep = ",")
colnames(d) <- as.character(header[1, ])

# 移除最后一列
d <- d[, -ncol(d)]
cat(sprintf("Data dimensions: %d contigs x %d samples\n", nrow(d), ncol(d)))

# ============================================================
# 2. 读取染色体分组
# ============================================================
chr_group_file <- paste0(opt$input_file, ".chr_group")
cat("Chromosome group file:", chr_group_file, "\n")
chr_group <- read.table(chr_group_file, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
colnames(chr_group) <- "Chromosome"

# 确保 chr_group 与数据矩阵行名一致
common_contigs <- intersect(rownames(d), rownames(chr_group))
d <- d[common_contigs, , drop = FALSE]
chr_group <- chr_group[common_contigs, , drop = FALSE]

# ============================================================
# 3. 排除样本
# ============================================================
if (!is.null(opt$exclude_sample)) {
  cat(sprintf("Excluding sample: %s\n", opt$exclude_sample))
  d <- d[, !colnames(d) %in% opt$exclude_sample, drop = FALSE]
  cat(sprintf("Data after exclusion: %d x %d\n", nrow(d), ncol(d)))
}

# ============================================================
# 4. 样本排序
# ============================================================
if (!is.null(opt$sample_order)) {
  cat("\n=== Using custom sample order ===\n")
  sample_order_df <- read.table(opt$sample_order, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  sample_order <- sample_order_df[[1]]
  sample_order <- intersect(sample_order, colnames(d))
  cat("Sample order:", paste(sample_order, collapse = ", "), "\n")
} else {
  cat("\n=== Clustering samples (ward.D2) ===\n")
  sample_dist <- dist(t(as.matrix(d)), method = "euclidean")
  sample_hclust <- hclust(sample_dist, method = "ward.D2")
  sample_order <- colnames(d)[sample_hclust$order]
  cat("Sample order:", paste(sample_order, collapse = ", "), "\n")
}

# ============================================================
# 5. 确定要绘制的染色体列表
# ============================================================
all_chromosomes <- unique(chr_group$Chromosome)

if (!is.null(opt$chromosome_list)) {
  cat("\n=== Reading chromosome list ===\n")
  cat("Chromosome list file:", opt$chromosome_list, "\n")
  chr_list <- readLines(opt$chromosome_list)
  chr_list <- trimws(chr_list)
  chr_list <- chr_list[nchar(chr_list) > 0]
  
  # 只保留数据中存在的染色体
  chr_list <- intersect(chr_list, all_chromosomes)
  
  if (length(chr_list) == 0) {
    stop("Error: No valid chromosomes found in the list", call.=FALSE)
  }
  
  cat(sprintf("Chromosomes to plot: %s\n", paste(chr_list, collapse = ", ")))
} else {
  chr_list <- all_chromosomes
  cat(sprintf("\nNo chromosome list provided, plotting all %d chromosomes\n", length(chr_list)))
}

# ============================================================
# 6. 为每条染色体单独绘图
# ============================================================
output_prefix <- ifelse(is.null(opt$output_prefix), opt$input_file, opt$output_prefix)

breaks <- seq(0, 1.0, by = 0.001)
color_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))

for (chr_id in chr_list) {
  cat(sprintf("\n=== Plotting chromosome: %s ===\n", chr_id))
  
  # 提取该染色体的 contigs
  chr_contigs <- rownames(chr_group[chr_group$Chromosome == chr_id, , drop = FALSE])
  chr_contigs <- intersect(chr_contigs, rownames(d))
  
  if (length(chr_contigs) == 0) {
    cat(sprintf("  Skip %s: no contigs found\n", chr_id))
    next
  }
  
  cat(sprintf("  Contigs: %d\n", length(chr_contigs)))
  
  # 检查 contig 前缀是否与样本名重复，如有重复则从该热图中去除对应样本
  contig_prefixes <- unique(gsub("^([^_]+)_.*$", "\\1", chr_contigs))
  samples_to_remove <- intersect(contig_prefixes, colnames(d))
  if (length(samples_to_remove) > 0) {
    cat(sprintf("  Removing %d overlapping sample(s): %s\n", 
                length(samples_to_remove), paste(samples_to_remove, collapse = ", ")))
    chr_data <- d[chr_contigs, setdiff(colnames(d), samples_to_remove), drop = FALSE]
  } else {
    chr_data <- d[chr_contigs, , drop = FALSE]
  }
  
  # 提取数据
  chr_data$contig <- rownames(chr_data)
  chr_data$contig <- factor(chr_data$contig, levels = chr_contigs)
  
  # 更新样本顺序（去除被移除的样本）
  current_sample_order <- intersect(sample_order, colnames(chr_data))
  
  # 转换为长格式
  chr_long <- melt(chr_data,
                   id.vars = "contig",
                   variable.name = "sample",
                   value.name = "kmer_prop")
  
  chr_long$sample <- factor(chr_long$sample, levels = current_sample_order)
  
  # 动态决定是否显示 Y 轴样本 ID
  n_samples <- length(current_sample_order)
  show_y_labels <- n_samples <= 20
  
  # 绘图
  p <- ggplot(chr_long, aes(x = contig, y = sample, fill = kmer_prop)) +
    geom_tile() +
    scale_fill_gradientn(
      colors = color_palette,
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      labels = scales::percent_format(accuracy = 1),
      name = "Kmer\nProportion"
    ) +
    labs(x = NULL, y = "Sample") +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = if (show_y_labels) element_text(size = 7) else element_blank(),
      axis.ticks.y = if (show_y_labels) element_line() else element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7)
    )
  
  # 计算输出尺寸：固定 cell 宽度为 2mm，高度为 2mm
  cell_width_mm <- 3
  cell_height_mm <- 3
  plot_area_width_mm <- length(chr_contigs) * cell_width_mm
  plot_area_width_in <- plot_area_width_mm / 25.4
  plot_area_height_mm <- n_samples * cell_height_mm
  plot_area_height_in <- plot_area_height_mm / 25.4
  margin_x_in <- 2.0  # 左侧y轴标签 + 右侧图例
  margin_y_in <- 0.5  # 顶部和底部边距
  chr_width <- plot_area_width_in + margin_x_in
  chr_height <- plot_area_height_in + margin_y_in
  
  # 输出
  chr_output_prefix <- paste0(output_prefix, "_", chr_id)
  chr_pdf <- paste0(chr_output_prefix, ".pdf")
  chr_png <- paste0(chr_output_prefix, ".png")
  
  cat(sprintf("  Contigs: %d, Samples: %d, Cell size: %dx%dmm, Plot size: %.1f x %.1f inches\n", 
              length(chr_contigs), n_samples, cell_width_mm, cell_height_mm, chr_width, chr_height))
  cat(sprintf("  Saving PDF to %s\n", chr_pdf))
  pdf(chr_pdf, width = chr_width, height = chr_height)
  print(p)
  dev.off()
  
  cat(sprintf("  Saving PNG to %s\n", chr_png))
  png(chr_png, width = chr_width, height = chr_height, units = "in", res = 300)
  print(p)
  dev.off()
}

cat("\nDone!\n")
