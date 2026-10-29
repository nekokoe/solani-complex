#!/usr/bin/env Rscript

# 加载必要的库
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(optparse)

# 检查并安装ggnewscale包（用于创建多个fill scale）
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

# 解析命令行参数
parse_args <- function() {
  option_list <- list(
    make_option(c("-q", "--qfile"), type = "character", default = NULL, 
                help = "Input admixture .Q file path", metavar = "character"),
    make_option(c("-i", "--idfile"), type = "character", default = NULL, 
                help = "Input sample ID file path", metavar = "character"),
    make_option(c("-g", "--groups"), type = "character", default = NULL, 
                help = "Sample group information TSV file", metavar = "character"),
    make_option(c("-s", "--species"), type = "character", default = NULL, 
                help = "Sample species information TSV file", metavar = "character"),
    make_option(c("-o", "--output"), type = "character", default = NULL, 
                help = "Output image file path", metavar = "character"),
    make_option(c("--width"), type = "double", default = 10, 
                help = "Image width (inches)", metavar = "double"),
    make_option(c("--height"), type = "double", default = 4, 
                help = "Image height (inches)", metavar = "double"),
    make_option(c("--dpi"), type = "integer", default = 300, 
                help = "Image resolution (DPI)", metavar = "integer"),
    make_option(c("--cluster_method"), type = "character", default = "ward.D2", 
                help = "Clustering method (default: ward.D2)", metavar = "character"),
    make_option(c("--no_cluster"), action = "store_true", default = FALSE, 
                help = "Do not perform clustering"),
    make_option(c("--sample_order"), type = "character", default = NULL, 
                help = "Sample order file: one sample ID per line, in desired order", metavar = "character"),
    make_option(c("--legend_title"), type = "character", default = "Populations", 
                help = "Legend title", metavar = "character")
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  opt <- optparse::parse_args(opt_parser)
  
  # 检查必需参数
  if (is.null(opt$qfile) || is.null(opt$idfile) || is.null(opt$groups) || 
      is.null(opt$species) || is.null(opt$output)) {
    print_help(opt_parser)
    stop("All input/output files are required.", call. = FALSE)
  }
  
  return(opt)
}

# 读取并预处理数据
read_data <- function(qfile_path, idfile_path, groups_path, species_path) {
  # 读取.Q文件
  q_df <- read.table(qfile_path, header = FALSE)
  num_pops <- ncol(q_df)
  colnames(q_df) <- paste0("Pop", 1:num_pops)
  
  # 读取样本ID文件（只读取第一列）
  id_df <- read.table(idfile_path, header = FALSE, colClasses = "character")
  sample_ids <- id_df[[1]]
  
  # 检查样本数量是否匹配
  if (length(sample_ids) != nrow(q_df)) {
    stop(sprintf("Sample count mismatch: %d in ID file, %d in Q file", 
                 length(sample_ids), nrow(q_df)))
  }
  
  # 设置样本ID为行名
  rownames(q_df) <- sample_ids
  
  # 读取分组信息（TSV文件，使用\t分隔）
  groups_df <- read.table(groups_path, sep = "\t", header = FALSE, 
                         col.names = c("SampleID", "Group"), 
                         colClasses = c("character", "character"))
  
  # 读取物种信息（TSV文件，使用\t分隔）
  species_df <- read.table(species_path, sep = "\t", header = FALSE, 
                          col.names = c("SampleID", "Species"), 
                          colClasses = c("character", "character"))
  
  # 合并分组和物种信息
  sample_info <- merge(groups_df, species_df, by = "SampleID", all = TRUE)
  
  # 填充缺失的物种信息
  sample_info$Species[is.na(sample_info$Species)] <- "Unknown"
  
  # 填充缺失的分组信息
  sample_info$Group[is.na(sample_info$Group)] <- "Unknown"
  
  # 确保所有样本都包含在sample_info中
  all_samples <- data.frame(SampleID = sample_ids, stringsAsFactors = FALSE)
  sample_info <- merge(all_samples, sample_info, by = "SampleID", all.x = TRUE)
  
  # 填充合并后可能产生的新缺失值
  sample_info$Group[is.na(sample_info$Group)] <- "Unknown"
  sample_info$Species[is.na(sample_info$Species)] <- "Unknown"
  
  # 检查是否有样本被标记为Unknown
  unknown_groups <- sample_info$SampleID[sample_info$Group == "Unknown"]
  if (length(unknown_groups) > 0) {
    cat("Warning: The following samples have unknown group: ", 
        paste(unknown_groups, collapse = ", "), "\n")
  }
  
  unknown_species <- sample_info$SampleID[sample_info$Species == "Unknown"]
  if (length(unknown_species) > 0) {
    cat("Warning: The following samples have unknown species: ", 
        paste(unknown_species, collapse = ", "), "\n")
  }
  
  # 将分组和物种信息添加到q_df中
  q_df$SampleID <- rownames(q_df)
  q_df <- merge(q_df, sample_info, by = "SampleID")
  rownames(q_df) <- q_df$SampleID
  
  return(q_df)
}

# 对样本进行聚类并返回排序后的索引
cluster_samples <- function(q_df, method = "ward.D2") {
  # 只选择数值型列进行聚类
  numeric_cols <- which(sapply(q_df, is.numeric))
  numeric_data <- q_df[, numeric_cols]
  
  # 执行层次聚类
  dist_matrix <- dist(numeric_data, method = "euclidean")
  hclust_result <- hclust(dist_matrix, method = method)
  
  # 获取聚类后的样本顺序
  ordered_samples <- hclust_result$labels[hclust_result$order]
  
  # 根据物种信息调整样本顺序，确保相同物种的样本相邻
  if ("Species" %in% colnames(q_df)) {
    # 创建一个包含样本ID和物种信息的数据框
    sample_species <- q_df[ordered_samples, "Species", drop = FALSE]
    sample_species$SampleID <- rownames(sample_species)
    
    # 按物种排序，然后在每个物种内部保持聚类顺序
    # 首先按物种分组，然后在每个组内保持原有顺序
    ordered_by_species <- character(0)
    species_levels <- unique(sample_species$Species)
    
    for (species in species_levels) {
      # 获取该物种的样本
      species_samples <- sample_species$SampleID[sample_species$Species == species]
      # 按照它们在聚类结果中的顺序添加
      for (sample in ordered_samples) {
        if (sample %in% species_samples) {
          ordered_by_species <- c(ordered_by_species, sample)
        }
      }
    }
    
    ordered_samples <- ordered_by_species
  }
  
  return(ordered_samples)
}

# 定义物种颜色映射函数（深色系）
get_species_colors <- function(species) {
  # 对物种进行排序，确保每次运行时顺序一致
  species <- sort(unique(species))
  
  # 深色系颜色映射（参考ggpubr的nature系列配色）
  color_map <- c()
  
  # 重点物种的固定颜色
  key_species_colors <- c(
    "solani" = "#E15759",  # 红色系（nature配色）
    "falciformis" = "#F28E2B",  # 橙色系（nature配色）
    "keratoplastica" = "#4E79A7",  # 蓝色系（nature配色）
    "unknown" = "#B0B0B0"  # 灰色系
  )
  
  # 其他物种的颜色库（使用spectral配色表的深色部分，增加取色密度）
  library(RColorBrewer)
  spectral_colors <- brewer.pal(11, "Spectral")
  # 过滤掉与固定颜色相似的颜色，只保留深色部分
  base_colors <- spectral_colors[c(1,2,4,5,7,8,10,11)]
  # 使用colorRampPalette生成更多颜色，避免重复
  color_generator <- colorRampPalette(base_colors)
  
  # 分离重点物种和其他物种
  key_species_list <- intersect(tolower(species), names(key_species_colors))
  other_species_list <- setdiff(species, species[tolower(species) %in% key_species_list])
  
  # 为重点物种分配颜色
  for (sp in species) {
    sp_lower <- tolower(sp)
    if (sp_lower %in% key_species_list) {
      color_map[sp] <- key_species_colors[sp_lower]
    }
  }
  
  # 为其他物种分配颜色
  if (length(other_species_list) > 0) {
    available_colors <- color_generator(length(other_species_list))
    for (i in seq_along(other_species_list)) {
      color_map[other_species_list[i]] <- available_colors[i]
    }
  }
  
  return(color_map)
}

# 定义来源颜色映射函数（浅色系）
get_group_colors <- function(groups) {
  # 对分组进行排序，确保每次运行时顺序一致
  groups <- sort(unique(groups))
  
  # 浅色系颜色映射（参考ggpubr的nature系列配色）
  color_map <- c()
  
  # 固定来源的颜色（进一步加深的浅色系）
  key_group_colors <- c(
    "human" = "#F5B7B1",  # 浅红色系（进一步加深）
    "med" = "#F5B7B1",  # 浅红色系（进一步加深）
    "plant" = "#A9DFBF",  # 浅绿色系（进一步加深）
    "soil" = "#F9E79F",  # 浅土黄色系（进一步加深）
    "unknown" = "#D0D0D0"  # 浅灰色系（进一步加深）
  )
  
  # 其他来源的颜色库（使用spectral配色表的中色系部分，增加取色密度）
  library(RColorBrewer)
  spectral_colors <- brewer.pal(11, "Spectral")
  # 过滤掉与固定颜色相似的颜色，只保留中色系部分
  base_colors <- spectral_colors[c(1,3,5,8,10,11)]
  # 使用colorRampPalette生成更多颜色，避免重复
  color_generator <- colorRampPalette(base_colors)
  
  # 分离固定来源和其他来源
  key_group_list <- intersect(tolower(groups), names(key_group_colors))
  other_group_list <- setdiff(groups, groups[tolower(groups) %in% key_group_list])
  
  # 为固定来源分配颜色
  for (group in groups) {
    group_lower <- tolower(group)
    if (group_lower %in% key_group_list) {
      color_map[group] <- key_group_colors[group_lower]
    }
  }
  
  # 为其他来源分配颜色
  if (length(other_group_list) > 0) {
    available_colors <- color_generator(length(other_group_list))
    for (i in seq_along(other_group_list)) {
      color_map[other_group_list[i]] <- available_colors[i]
    }
  }
  
  return(color_map)
}

# 绘制admixture堆积柱状图
plot_admixture <- function(q_df, output_path, width = 12, height = 8, dpi = 300, legend_title = "Populations") {
  # 准备数据用于ggplot
  # 获取数值型列（群体数据）
  numeric_cols <- which(sapply(q_df, is.numeric))
  pop_data <- q_df[, numeric_cols]
  pop_names <- colnames(pop_data)

  # 转换为长格式
  long_data <- pop_data %>% 
    as.data.frame() %>% 
    mutate(SampleID = rownames(pop_data)) %>% 
    pivot_longer(cols = all_of(pop_names), 
                 names_to = "Population", 
                 values_to = "Proportion")
  
  # 添加分组和物种信息
  long_data <- long_data %>% 
    left_join(q_df[, c("SampleID", "Group", "Species")], by = "SampleID")
  
  # 创建样本顺序的因子变量，保持聚类顺序
  long_data$SampleID <- factor(long_data$SampleID, levels = q_df$SampleID)
  
  # 使用类似npg文献的配色方案（用于群体）
  npg_colors <- c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
  )
  
  # 如果群体数量超过配色数量，循环使用颜色
  colors <- rep(npg_colors, times = ceiling(length(pop_names) / length(npg_colors)))
  
  # 准备样本标签数据
  sample_labels_data <- q_df[, c("SampleID", "Group", "Species")]
  sample_labels_data$SampleID <- factor(sample_labels_data$SampleID, levels = q_df$SampleID)
  
  # 定义颜色映射
  species_colors <- get_species_colors(sample_labels_data$Species)
  group_colors <- get_group_colors(sample_labels_data$Group)
  
  # 重新设计绘图逻辑，使用ggnewscale来创建多个fill scale
  # 创建主图表，只包含堆积图
  final_plot <- ggplot() +
    # 绘制堆积图
    geom_bar(
      data = long_data,
      aes(x = SampleID, y = Proportion, fill = Population),
      stat = "identity", 
      position = "stack", 
      width = 1.0
    ) +
    # 为Population设置颜色映射
  scale_fill_manual(
    name = legend_title,
    values = colors,
    breaks = pop_names,
    labels = pop_names
  ) +
  # 新的scale用于Group
  new_scale_fill() +
  # 来源（分组）的颜色注释条（上方，浅色系）
  geom_rect(
    data = sample_labels_data,
    aes(
      xmin = as.numeric(SampleID) - 0.5,
      xmax = as.numeric(SampleID) + 0.5,
      ymin = -0.08,
      ymax = -0.02,
      fill = Group
    ),
    inherit.aes = FALSE
  ) +
  # 为Group设置颜色映射
  scale_fill_manual(
    name = "Group",
    values = group_colors,
    breaks = names(group_colors),
    labels = names(group_colors)
  ) +
  # 新的scale用于Species
  new_scale_fill() +
  # 物种的颜色注释条（下方，深色系）
  geom_rect(
    data = sample_labels_data,
    aes(
      xmin = as.numeric(SampleID) - 0.5,
      xmax = as.numeric(SampleID) + 0.5,
      ymin = -0.15,
      ymax = -0.09,
      fill = Species
    ),
    inherit.aes = FALSE
  ) +
  # 为Species设置颜色映射
  scale_fill_manual(
    name = "Species",
    values = species_colors,
    breaks = names(species_colors),
    labels = names(species_colors)
  ) +
    # 设置标签和主题
    labs(y = "Ancestry Proportion", x = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 8),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.2, "cm")
    )
  
  # 统一设置所有图例的大小和排列方向
  final_plot <- final_plot +
    theme(
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.direction = "horizontal",
      legend.position = "bottom",
      legend.box = "vertical"
    )
  
  # 调整坐标轴范围，为颜色注释留出空间
  final_plot <- final_plot +
    coord_cartesian(ylim = c(-0.18, 1), clip = "off") +
    # 调整图边距，为颜色注释留出空间
    theme(
      plot.margin = margin(t = 5, r = 10, b = 25, l = 10)
    )
  
  # 提取基础路径（不含扩展名）
  base_path <- tools::file_path_sans_ext(output_path)
  png_path <- paste0(base_path, ".png")
  pdf_path <- paste0(base_path, ".pdf")
  
  # 保存为PNG格式
  png(png_path, width = width, height = height, units = "in", res = dpi)
  print(final_plot)
  dev.off()
  
  # 保存为PDF格式
  pdf(pdf_path, width = width, height = height)
  print(final_plot)
  dev.off()
  
  cat(sprintf("Plots saved to %s and %s\n", png_path, pdf_path))
}

# 主函数
main <- function() {
  # 解析命令行参数
  opt <- parse_args()
  
  tryCatch({
    # 读取数据
    cat("Reading data...\n")
    q_df <- read_data(opt$qfile, opt$idfile, opt$groups, opt$species)
    cat(sprintf("Successfully read data: %d samples, %d populations\n", 
                nrow(q_df), length(which(sapply(q_df, is.numeric)))))
    
    # 处理样本顺序
    if (!is.null(opt$sample_order)) {
      cat("Using manual sample order from file...\n")
      manual_order <- readLines(opt$sample_order)
      # 检查所有手动指定的样本是否都在数据中
      missing_samples <- setdiff(manual_order, q_df$SampleID)
      if (length(missing_samples) > 0) {
        cat(sprintf("Warning: %d samples in order file not found in data: %s\n", 
                    length(missing_samples), paste(missing_samples, collapse = ", ")))
      }
      # 只保留在数据中存在的样本
      manual_order <- intersect(manual_order, q_df$SampleID)
      # 处理没有在手动顺序中指定的样本
      additional_samples <- setdiff(q_df$SampleID, manual_order)
      if (length(additional_samples) > 0) {
        cat(sprintf("Warning: %d samples not in order file, appending to end: %s\n", 
                    length(additional_samples), paste(additional_samples, collapse = ", ")))
        manual_order <- c(manual_order, additional_samples)
      }
      q_df <- q_df[manual_order, ]
    } else if (!opt$no_cluster) {
      # 进行聚类分析
      cat(sprintf("Performing %s clustering...\n", opt$cluster_method))
      ordered_samples <- cluster_samples(q_df, method = opt$cluster_method)
      q_df <- q_df[ordered_samples, ]
      cat("Clustering completed\n")
    }
    
    # 绘制图形
    cat("Generating plot...\n")
    plot_admixture(q_df, opt$output, width = opt$width, height = opt$height, dpi = opt$dpi, legend_title = opt$legend_title)
    
  }, error = function(e) {
    cat(sprintf("Error: %s\n", e$message))
    quit(status = 1)
  })
}

# 执行主函数
if (!interactive()) {
  main()
}
