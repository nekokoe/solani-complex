#!/usr/bin/env Rscript

# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(optparse)

# 命令行参数解析
option_list <- list(
  make_option(c("-g", "--gene_chrom_file"), type="character", default=NULL,
              help="基因与染色体对应关系文件路径 (tsv格式)", metavar="character"),
  make_option(c("-o", "--gene_hog_file"), type="character", default=NULL,
              help="基因与HOG对应关系文件路径 (tsv格式)", metavar="character"),
  make_option(c("-p", "--hog_presence_files"), type="character", default=NULL,
              help="HOG与泛基因组区间对应关系文件路径 (逗号分隔多个文件)", metavar="character"),
  make_option(c("-f", "--output_format"), type="character", default="png",
              help="输出图片格式 (png, pdf, svg等)", metavar="character"),
  make_option(c("-w", "--width"), type="numeric", default=12,
              help="图片宽度", metavar="numeric"),
  make_option(c("--height_per_plot"), type="numeric", default=6,
              help="每个子图的高度", metavar="numeric"),
  make_option(c("-t", "--title"), type="character", default="Chromosome Gene Interval Distribution",
              help="图片标题", metavar="character"),
  make_option(c("--output_file"), type="character", default="chromosome_gene_interval_plot",
              help="输出文件名前缀", metavar="character")
)

# 创建OptionParser对象，禁用自动添加的帮助选项
option_parser <- OptionParser(option_list=option_list, add_help_option=TRUE)

oopt <- parse_args(option_parser)

# 参数验证
if (is.null(oopt$gene_chrom_file) || is.null(oopt$gene_hog_file) || is.null(oopt$hog_presence_files)) {
  cat("错误: 必须提供所有必要的输入文件路径\n\n")
  print_help(option_parser)
  quit(status=1)
}

# 读取文件函数
read_tsv_file <- function(file_path, col_names, ...) {
  if (!file.exists(file_path)) {
    stop(paste("文件不存在:", file_path))
  }
  tryCatch({
    read.table(file_path, sep = "\t", header = FALSE, col.names = col_names, stringsAsFactors = FALSE, ...)
  }, error = function(e) {
    stop(paste("读取文件失败", file_path, ":", e$message))
  })
}

# 读取基因与染色体对应关系
message("读取基因与染色体对应关系...")
gene_chrom <- read_tsv_file(oopt$gene_chrom_file, col_names = c("TranscriptID", "GeneID", "ContigID"))

# 使用TranscriptID作为基因ID
gene_chrom <- gene_chrom %>% 
  select(GeneID = TranscriptID, ContigID)

# 读取基因与HOG对应关系
message("读取基因与HOG对应关系...")
gene_hog <- read_tsv_file(oopt$gene_hog_file, col_names = c("GeneID", "SampleID", "HOG_ID", "OG_ID"))

# 处理HOG存在率文件列表
hog_presence_files <- str_split(oopt$hog_presence_files, ",")[[1]]

# 定义处理单个HOG存在率文件的函数
process_hog_presence_file <- function(hog_presence_file) {
  message(paste("处理HOG存在率文件:", hog_presence_file))
  hog_presence <- read_tsv_file(hog_presence_file, col_names = c("HOG", "Presence_Rate", "Presence_Rate_Interval"))
  
  # 合并数据
  merged_data <- gene_chrom %>%
    left_join(gene_hog %>% select(GeneID, HOG_ID), by = "GeneID") %>%
    left_join(hog_presence, by = c("HOG_ID" = "HOG"))
  
  # 处理未找到的HOG，归为0-20%区间
  merged_data <- merged_data %>%
    mutate(Presence_Rate_Interval = ifelse(is.na(Presence_Rate_Interval), "0-20%", Presence_Rate_Interval))
  
  # 统计每个染色体上各区间的基因数量
  chrom_interval_counts <- merged_data %>%
    group_by(ContigID, Presence_Rate_Interval) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(ContigID, Presence_Rate_Interval)
  
  # 计算每个染色体的基因总数
  chrom_total_counts <- merged_data %>%
    group_by(ContigID) %>%
    summarise(total = n(), .groups = "drop")
  
  # 计算比例
  chrom_interval_ratios <- chrom_interval_counts %>%
    left_join(chrom_total_counts, by = "ContigID") %>%
    mutate(ratio = count / total)
  
  # 为文件名添加标识符
  file_label <- basename(hog_presence_file)
  file_label <- str_replace_all(file_label, "\\..*", "")
  
  list(data = chrom_interval_ratios, label = file_label)
}

# 处理所有HOG存在率文件
processed_data_list <- map(hog_presence_files, process_hog_presence_file)

# 确定染色体的排序（基于第一个文件的数据进行聚类）
message("计算染色体聚类顺序...")

# 为第一个数据集创建比例矩阵用于聚类
first_data <- processed_data_list[[1]]$data

# 创建宽格式矩阵
ratio_matrix <- first_data %>%
  pivot_wider(names_from = Presence_Rate_Interval, values_from = ratio, values_fill = 0)

# 提取数值部分进行聚类
numeric_data <- ratio_matrix %>% select(-ContigID)

# 执行层次聚类
if (nrow(numeric_data) > 1) {
  # 计算欧氏距离
  dist_matrix <- dist(numeric_data)
  # 层次聚类
  hclust_result <- hclust(dist_matrix, method = "ward.D2")
  # 获取聚类后的染色体顺序并确保唯一性
  chromosome_order <- ratio_matrix$ContigID[hclust_result$order]
  chromosome_order <- unique(chromosome_order)
} else {
  # 如果只有一个染色体，则直接使用并确保唯一性
  chromosome_order <- unique(ratio_matrix$ContigID)
}

# 为所有数据集设置染色体顺序
processed_data_list <- map(processed_data_list, function(item) {
  # 确保只使用数据中存在的染色体ID
  common_contigs <- intersect(item$data$ContigID, chromosome_order)
  # 如果没有共同的染色体，使用数据中所有唯一的染色体
  if (length(common_contigs) == 0) {
    common_contigs <- unique(item$data$ContigID)
  }
  item$data$ContigID <- factor(item$data$ContigID, levels = common_contigs)
  item
})

# 定义颜色方案
interval_colors <- c(
  "Core" = "#E64B35",  # 红色
  "Softcore" = "#4DBBD5",  # 蓝色
  "Shell" = "#00A087",  # 绿色
  "Cloud" = "#3C5488"  # 深蓝色
)

# 自定义颜色函数（处理未预定义的区间）
get_interval_colors <- function(intervals) {
  default_colors <- scales::hue_pal()(length(intervals))
  names(default_colors) <- intervals
  # 替换预定义的颜色
  for (interval in names(interval_colors)) {
    if (interval %in% intervals) {
      default_colors[interval] <- interval_colors[interval]
    }
  }
  default_colors
}

# 生成所有区间的颜色映射
all_intervals <- unique(unlist(map(processed_data_list, function(item) item$data$Presence_Rate_Interval)))
all_colors <- get_interval_colors(all_intervals)

# 为每个数据集创建ggplot对象
message("生成图形...")
plot_list <- map(processed_data_list, function(item) {
  data <- item$data
  label <- item$label
  
  # 计算标签位置（堆叠条形图的中间位置）
  # 确保Presence_Rate_Interval有固定顺序，以保证标签位置正确
  interval_order <- c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%", "Core", "Softcore", "Shell", "Cloud")
  
  label_positions <- data %>%
    group_by(ContigID) %>%
    # 确保按区间顺序排列，使堆叠顺序一致
    arrange(factor(Presence_Rate_Interval, levels = interval_order)) %>%
    mutate(cum_ratio = cumsum(ratio),
           label_pos = cum_ratio - ratio/2,
           # 只在比例大于0.05时显示标签
           display_label = ifelse(ratio > 0.05, paste0(round(ratio*100, 1), "%"), "")) %>%
    ungroup()
  
  # 创建ggplot
  p <- ggplot(data, aes(x = ContigID, y = ratio, fill = Presence_Rate_Interval)) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    geom_text(data = label_positions, 
              aes(x = ContigID, y = label_pos, label = display_label),
              size = 3, color = "white", fontface = "bold") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = all_colors, name = "Genome Interval") +
    labs(title = label, 
         x = "Chromosome", 
         y = "Gene Proportion") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5)
    )
  
  p
})

# 组合多个子图
library(ggpubr)
message("组合图形...")

# 计算总高度
total_height <- oopt$height_per_plot * length(plot_list)

# 组合图形
combined_plot <- ggarrange(plotlist = plot_list, 
                          nrow = length(plot_list),
                          ncol = 1,
                          common.legend = TRUE,
                          legend = "right",
                          align = "v")

# 添加主标题
if (!is.null(oopt$title)) {
  combined_plot <- annotate_figure(combined_plot, 
                                  top = text_grob(oopt$title, 
                                                 face = "bold", 
                                                 size = 16))
}

# 输出文件路径
output_file <- paste0(oopt$output_file, ".", oopt$output_format)

# 保存图形
message(paste("保存图形到:", output_file))
ggsave(output_file, 
       combined_plot, 
       width = oopt$width, 
       height = total_height, 
       dpi = 300, 
       device = oopt$output_format)

# 输出统计信息
message("\n统计信息:")
message(paste("总染色体数量:", length(chromosome_order)))
message(paste("总泛基因组区间类别:", length(all_intervals)))
message(paste("处理的HOG存在率文件数量:", length(hog_presence_files)))
message("\n程序执行完成!")