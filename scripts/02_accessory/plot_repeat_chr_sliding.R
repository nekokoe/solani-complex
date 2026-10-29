#!/usr/bin/env Rscript
# plot_repeat_chr_sliding.R
# 绘制染色体重复序列滑窗比例堆积图（类似 MCscan/ALLMAPS 风格）

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales)

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript plot_repeat_chr_sliding.R -i <input.tsv> -f <genome.fasta> [-o output.pdf] [-m minlength] [-c color1,color2,...]\n")
  cat("  -i  stats_repeat_chr_sliding.py 输出的 TSV 文件\n")
  cat("  -f  基因组 fasta 文件，用于获取染色体长度\n")
  cat("  -o  输出 PDF 文件 (默认: repeat_chr_sliding.pdf)\n")
  cat("  -m  最小染色体长度，低于此值不绘制 (默认: 100000)\n")
  cat("  -c  逗号分隔的颜色列表，用于 repeat_family\n")
  quit(status = 1)
}

input_file <- NULL
fasta_file <- NULL
output_file <- "repeat_chr_sliding.pdf"
custom_cols <- NULL
min_length <- 100000  # 默认 100kb

i <- 1
while (i <= length(args)) {
  switch(args[i],
    "-i" = { input_file <- args[i + 1]; i <- i + 2 },
    "-f" = { fasta_file <- args[i + 1]; i <- i + 2 },
    "-o" = { output_file <- args[i + 1]; i <- i + 2 },
    "-c" = { custom_cols <- args[i + 1]; i <- i + 2 },
    "-m" = { min_length <- as.integer(args[i + 1]); i <- i + 2 },
    { cat(sprintf("Unknown argument: %s\n", args[i])); i <- i + 1 }
  )
}

if (is.null(input_file) || is.null(fasta_file)) {
  cat("Error: -i and -f are required.\n")
  quit(status = 1)
}

# 从 fasta 文件获取染色体长度
get_chrom_lengths <- function(fasta_path) {
  chrom_lengths <- list()
  current_chrom <- NULL
  cur_len <- 0
  con <- file(fasta_path, "r")
  while (length(line <- readLines(con, n = 1)) > 0) {
    line <- trimws(line)
    if (startsWith(line, ">")) {
      if (!is.null(current_chrom)) {
        chrom_lengths[[current_chrom]] <- cur_len
      }
      current_chrom <- trimws(strsplit(sub("^>", "", line), "\\s+")[[1]][1])
      cur_len <- 0
    } else {
      cur_len <- cur_len + nchar(line)
    }
  }
  close(con)
  if (!is.null(current_chrom)) {
    chrom_lengths[[current_chrom]] <- cur_len
  }
  return(chrom_lengths)
}

# 读取 fasta 获取染色体长度
cat("读取 fasta 文件获取染色体长度...\n")
chrom_lengths <- get_chrom_lengths(fasta_file)
total_chr <- length(chrom_lengths)

# 过滤掉长度小于 min_length 的染色体
chrom_lengths <- chrom_lengths[sapply(chrom_lengths, function(x) x >= min_length)]
cat(sprintf("共 %d 条染色体，过滤后保留 %d 条 (>= %d bp)\n", total_chr, length(chrom_lengths), min_length))

# 按长度从长到短排序（上长下短）
chrom_order <- names(sort(unlist(chrom_lengths), decreasing = TRUE))

if (length(chrom_order) == 0) {
  cat("Error: 没有符合条件的染色体，请检查 -m 参数或输入文件。\n")
  quit(status = 1)
}

# 读取重复序列统计结果
cat("读取重复序列滑窗统计结果...\n")
repeat_data <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# 只保留在 chrom_order 中的染色体数据
repeat_data <- repeat_data[repeat_data$chr_id %in% chrom_order, , drop = FALSE]
cat(sprintf("共 %d 行数据\n", nrow(repeat_data)))

# 自动获取窗口大小和步长
window_size <- unique(repeat_data$end - repeat_data$start)
if (length(window_size) > 1) {
  window_size <- window_size[1]
  cat(sprintf("Warning: 检测到多种窗口大小，使用 %d\n", window_size))
}
# 从同一染色体内相邻窗口的 start 差值推断步长
if (nrow(repeat_data) > 1) {
  first_chr <- unique(repeat_data$chr_id)[1]
  chr_starts <- sort(unique(repeat_data$start[repeat_data$chr_id == first_chr]))
  step_size <- if (length(chr_starts) > 1) diff(chr_starts)[1] else window_size
} else {
  step_size <- window_size
}
cat(sprintf("窗口大小: %d, 步长: %d\n", window_size, step_size))

# 获取所有 repeat_family 类型
families <- sort(unique(repeat_data$repeat_family))
cat(sprintf("共 %d 种 repeat family: %s\n", length(families), paste(families, collapse = ", ")))

# 设置颜色（参考 plot_chromosome_repeat_content.R 的 Spectral 渐变配色）
if (!is.null(custom_cols)) {
  my_colors <- strsplit(custom_cols, ",")[[1]]
} else {
  n_fam <- length(families)
  spectral_colors <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(n_fam)
  my_colors <- spectral_colors
}
names(my_colors) <- families

# 为每条染色体生成完整的窗口列表，使用 list + bind_rows 避免 rbind 列数不匹配
cat("构建完整的窗口列表...\n")
plot_list <- list()

for (idx in seq_along(chrom_order)) {
  chr <- chrom_order[idx]
  chr_length <- chrom_lengths[[chr]]
  # 生成该染色体的所有窗口
  starts <- seq(0, chr_length - window_size, by = step_size)
  if (length(starts) == 0) {
    starts <- 0
  }
  chr_windows <- data.frame(chr_id = chr, start = starts, end = starts + window_size, stringsAsFactors = FALSE)

  # 获取该染色体的 repeat 数据
  chr_repeat <- repeat_data[repeat_data$chr_id == chr, ]

  if (nrow(chr_repeat) == 0) {
    # 该染色体没有 repeat，所有 window 的 proportion = 0
    chr_windows$repeat_family <- "none"
    chr_windows$proportion <- 0
  } else {
    chr_plot <- merge(chr_windows,
                      chr_repeat[, c("start", "end", "repeat_family", "proportion")],
                      by = c("start", "end"), all.x = TRUE)
    chr_plot$proportion[is.na(chr_plot$proportion)] <- 0
    chr_plot$repeat_family[is.na(chr_plot$repeat_family)] <- "none"
    chr_windows <- chr_plot
  }

  plot_list[[idx]] <- chr_windows
}

# 使用 bind_rows 安全合并（兼容不同列数/类型）
plot_data <- dplyr::bind_rows(plot_list)

# 计算非重复序列部分
plot_data <- plot_data %>%
  group_by(chr_id, start) %>%
  mutate(repeat_sum = sum(proportion)) %>%
  ungroup() %>%
  mutate(non_repeat = pmax(1 - repeat_sum, 0))

# 为每个窗口添加一行 non_repeat
window_unique <- plot_data %>%
  distinct(chr_id, start, non_repeat) %>%
  mutate(repeat_family = "non_repeat", proportion = non_repeat) %>%
  select(chr_id, start, repeat_family, proportion)

# 保留原始的 repeat 行
repeat_rows <- plot_data %>%
  filter(repeat_family != "none") %>%
  select(chr_id, start, repeat_family, proportion)

# 合并：有 repeat 的窗口用 repeat 行 + non_repeat 行；无 repeat 的窗口只用 non_repeat 行
plot_long <- dplyr::bind_rows(repeat_rows, window_unique) %>%
  select(chr_id, start, repeat_family, proportion)

# 设置因子顺序（按长度从长到短，与 chrom_order 一致）
plot_long$chr_id <- factor(plot_long$chr_id, levels = chrom_order)
all_families <- c("non_repeat", families)
my_colors["non_repeat"] <- "#E0E0E0"
plot_long$repeat_family <- factor(plot_long$repeat_family, levels = all_families)

# 绘制
cat("绘制染色体重复序列图...\n")

# 准备染色体名称标签（放在每个子图右侧）
chr_labels <- data.frame(
  chr_id = chrom_order,
  label = chrom_order,
  x_max = sapply(chrom_order, function(ch) {
    chr_data <- plot_long[plot_long$chr_id == ch, ]
    if (nrow(chr_data) == 0) return(0)
    max(chr_data$start) / 1e6 + window_size / 1e6
  }),
  stringsAsFactors = FALSE
)

p <- ggplot(plot_long, aes(x = start / 1e6, y = proportion, fill = repeat_family)) +
  geom_col(position = "stack", width = window_size / 1e6 * 0.95, show.legend = TRUE) +
  geom_text(data = chr_labels, aes(x = x_max, y = 0.5, label = label),
            hjust = 0, vjust = 0.5, angle = 0, size = 3, fontface = "bold",
            inherit.aes = FALSE) +
  facet_grid(rows = vars(chr_id), scales = "free_x") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.15))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = my_colors, drop = FALSE) +
  labs(
    x = "Position (Mb)",
    y = "Proportion",
    fill = "Repeat family"
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10, margin = margin(t = 5)),
    strip.text = element_blank(),
    strip.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(2, "mm"),
    panel.grid = element_blank(),
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    plot.margin = margin(r = 2, unit = "cm")
  )

# 保存 PDF
num_chr <- length(chrom_order)
chr_height <- max(8, num_chr * 0.6)
cat(sprintf("保存输出到 %s (尺寸: 12 x %.1f 英寸, %d 条染色体)...\n", output_file, chr_height, num_chr))
pdf(output_file, width = 12, height = chr_height)
print(p)
dev.off()
cat("Done!\n")
