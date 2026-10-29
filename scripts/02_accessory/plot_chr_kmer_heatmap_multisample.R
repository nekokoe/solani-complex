#!/usr/bin/env Rscript

# 设置错误处理选项，显示完整的调用栈信息
options(error = function() traceback(sys.calls(), max.lines = 200))
options(warn = 1) # 实时显示警告信息

library(optparse)
library(pheatmap)
library(RColorBrewer)

# 定义命令行参数
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, 
              help="Path to input data file", metavar="character"),
  make_option(c("--exclude_sample"), type="character", default=NULL, 
              help="Sample name to exclude", metavar="character"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Path to output heatmap file prefix (will output both png and pdf)", metavar="character"),
  make_option(c("--width"), type="numeric", default=16,
              help="Heatmap width in inches (default: 16)", metavar="numeric"),
  make_option(c("--height"), type="numeric", default=9,
              help="Heatmap height in inches (default: 9)", metavar="numeric"),
  make_option(c("--annotation_file"), type="character", default=NULL, 
              help="Path to sample annotation file (TSV with header)", metavar="character"),
  make_option(c("--sample_order"), type="character", default=NULL,
              help="File with sample order (tsv, single column). If provided, columns will not be clustered and will be ordered by this list.", metavar="character")
)

# 解析命令行参数
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必需参数
if (is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

# 确保输出目录存在
if (!is.null(opt$output_file)) {
  output_dir <- dirname(opt$output_file)
  if (!is.null(output_dir) && output_dir != "" && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive=TRUE)
    cat(paste("Creating output directory:", output_dir, "\n"))
  }
}

# 参考plot_admixture.R的配色方案
# 创建配色方案 - 与plot_admixture.R保持一致的颜色方案
create_group_colors <- function(group_names) {
  # 定义各色系颜色列表 - 基于plot_admixture.R的配色风格
  # 高饱和度红色系（用于med）- 使用更鲜艳的红色
  med_red_colors <- c('#FF0000', '#FF3333', '#FF5555', '#FF6666', '#FF4444', '#FF2222')
  # 较低饱和度红色系（用于Human）
  human_red_colors <- c('#E41A1C', '#CC0000', '#B22222', '#A52A2A', '#8B0000', '#CD5C5C')
  # 蓝色系（用于env）
  blue_colors <- c('#1E90FF', '#6495ED', '#87CEEB', '#4682B4', '#5F9EA0', '#708090')
  # 橙色系（用于Mammal）
  orange_colors <- c('#FF8C00', '#FFA726', '#FFB74D', '#FF9800', '#F57C00', '#EF6C00')
  # 土黄色系（用于Soil）
  soil_colors <- c('#DAA520', '#E6BE8A', '#CCA01D', '#CD853F', '#DEB887', '#BC8F8F')
  # 草绿色系（用于Plant）
  plant_colors <- c('#228B22', '#3CB371', '#2E8B57', '#4CAF50', '#66BB6A', '#81C784')
  # 深灰色（用于Unknown）
  unknown_color <- '#555555'
  
  # 用于其他标签的备用颜色列表
  other_colors <- c('#377EB8', '#984EA3', '#A65628', '#F781BF',
                    '#999999', '#7570B3', '#E7298A', '#E6AB02', '#1E88E5', 
                    '#42A5F5', '#90CAF9', '#1976D2', '#1565C0', '#0D47A1', 
                    '#795548', '#A1887F', '#D7CCC8', '#6D4C41', '#5D4037')
  
  # 创建颜色映射
  color_map <- character(length(group_names))
  names(color_map) <- group_names
  
  # 为已知类型的组分配特定色系
  color_index <- 1
  
  for (group_name in group_names) {
    # 转为小写进行匹配
    group_lower <- tolower(group_name)
    
    if (grepl("med", group_lower)) {
      # med使用高饱和度红色系
      color_map[group_name] <- med_red_colors[color_index %% length(med_red_colors) + 1]
    } else if (grepl("human", group_lower)) {
      # human使用较低饱和度红色系
      color_map[group_name] <- human_red_colors[color_index %% length(human_red_colors) + 1]
    } else if (grepl("env", group_lower)) {
      # env使用蓝色系
      color_map[group_name] <- blue_colors[color_index %% length(blue_colors) + 1]
    } else if (grepl("mammal", group_lower)) {
      # 使用橙色系
      color_map[group_name] <- orange_colors[color_index %% length(orange_colors) + 1]
    } else if (grepl("soil", group_lower)) {
      # 使用土黄色系
      color_map[group_name] <- soil_colors[color_index %% length(soil_colors) + 1]
    } else if (grepl("plant", group_lower)) {
      # 使用草绿色系
      color_map[group_name] <- plant_colors[color_index %% length(plant_colors) + 1]
    } else if (grepl("unknown", group_lower)) {
      # 使用深灰色
      color_map[group_name] <- unknown_color
    } else {
      # 使用备用颜色
      color_map[group_name] <- other_colors[color_index %% length(other_colors) + 1]
    }
    
    color_index <- color_index + 1
  }
  
  return(color_map)
}

# 加载样本注释
load_sample_annotations <- function(annotation_file) {
  cat(paste("Reading sample annotation file:", annotation_file, "\n"))
  # 读取注释文件
  annotation_df <- read.delim(annotation_file, stringsAsFactors=FALSE)
  
  # 检查文件是否有足够的列
  if (ncol(annotation_df) < 2) {
    stop("Annotation file must have at least 2 columns: SampleID and at least one annotation column")
  }
  
  # 获取样本ID列（第一列）
  sample_id_col <- colnames(annotation_df)[1]
  # 获取注释列名（除第一列外的所有列）
  annotation_columns <- colnames(annotation_df)[-1]
  
  cat(paste("Loaded annotations for", nrow(annotation_df), "samples with", 
            length(annotation_columns), "annotation types:", 
            paste(annotation_columns, collapse=", "), "\n"))
  
  # 为pheatmap准备annotation_col数据框
  # 创建行名向量
  sample_ids <- annotation_df[[sample_id_col]]
  
  # 提取注释列
  annotation_col <- as.data.frame(annotation_df[, annotation_columns, drop=FALSE])
  
  # 设置行名
  rownames(annotation_col) <- sample_ids
  
  # 将所有列转换为因子，以确保正确的颜色映射
  for (col in colnames(annotation_col)) {
    annotation_col[[col]] <- as.factor(annotation_col[[col]])
  }
  
  return(list(annotation_col=annotation_col, annotation_columns=annotation_columns))
}

# 生成样本注释的函数（只生成注释，不修改矩阵）
generate_sample_annotations <- function(matrix_df, sample_annotations) {
  # 初始化返回值
  result <- list(
    annotation_col = NULL,
    annotation_colors = list()
  )
  
  # 处理样本注释（列注释）
  if (!is.null(sample_annotations)) {
    # 确保annotation_col与矩阵列名匹配
    annotation_col <- sample_annotations$annotation_col
    common_samples <- intersect(colnames(matrix_df), rownames(annotation_col))
    
    if (length(common_samples) > 0) {
      cat(paste("Number of common samples:", length(common_samples), "\n"))
      # 只保留共有样本
      annotation_col <- annotation_col[common_samples, , drop=FALSE]
      
      # 检查是否只有一个物种，如果是则记录物种名称
      for (col in colnames(annotation_col)) {
        if (tolower(col) == "species" && length(unique(annotation_col[[col]])) == 1) {
          # 记录唯一的物种名称
          single_species_name <- unique(annotation_col[[col]])[1]
          break
        }
      }
    } else {
      cat("Warning: No common samples between data matrix and annotations. No filtering applied.\n")
    }
  }
  
  # 如果还有注释列，则添加到结果中
  if (!is.null(annotation_col) && ncol(annotation_col) > 0) {
    result$annotation_col <- annotation_col
    
    # 为每个注释列创建配色方案
    for (col in colnames(annotation_col)) {
      group_names <- unique(annotation_col[[col]])
      result$annotation_colors[[col]] <- create_group_colors(group_names)
    }
  }
  
  return(result)
}

# 主函数
main <- function() {
  tryCatch({
    # 读取数据文件
    cat(paste("Reading input file:", opt$input_file, "\n"))
    d <- read.table(opt$input_file, sep = ",", row.names = 1, header = F)

    # 读取header文件
    header_file <- paste(opt$input_file, ".header", sep="")
    cat(paste("Reading header file:", header_file, "\n"))
    header <- read.table(header_file, sep=",")
    header <- as.matrix(header)
    colnames(d) <- header[1,]

    # 移除最后一列
    d <- d[,-dim(d)[2]]
    cat(paste("Original data dimensions:", dim(d)[1], "x", dim(d)[2], "\n"))

    # 读取chr_group文件，从contig ID提取样本信息用于列注释
    anno_col_file <- paste(opt$input_file, ".chr_group", sep="")
    cat(paste("Reading chromosome annotation file:", anno_col_file, "\n"))
    anno_col_df <- read.table(anno_col_file, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
    colnames(anno_col_df) <- "Chromosomes"

    # 从contig ID提取样本ID（第一个下划线之前）
    contig_ids <- rownames(anno_col_df)
    sample_ids_from_contig <- gsub("^([^_]+)_.*$", "\\1", contig_ids)
    anno_col_df$SampleID <- sample_ids_from_contig
    cat(paste("Extracted sample IDs from contig names\n"))

    # 创建样本颜色映射（用于列注释）
    unique_samples <- unique(sample_ids_from_contig)
    if (length(unique_samples) > 1) {
      sample_colors <- brewer.pal(min(length(unique_samples), 9), "Set1")
      if (length(unique_samples) > 9) {
        sample_colors <- colorRampPalette(sample_colors)(length(unique_samples))
      }
      names(sample_colors) <- sort(unique_samples)
      cat(paste("Detected", length(unique_samples), "samples, using sample-based coloring for columns\n"))
    } else {
      sample_colors <- NULL
      cat("Single sample detected\n")
    }

    # 处理样本注释文件（用于行注释）
    anno_row <- NULL
    annotation_colors <- list()

    if (!is.null(opt$annotation_file)) {
      cat("Loading row annotations from:", opt$annotation_file, "\n")
      anno_row_df <- read.delim(opt$annotation_file, stringsAsFactors = FALSE)
      sample_id_col <- colnames(anno_row_df)[1]
      rownames(anno_row_df) <- anno_row_df[[sample_id_col]]
      anno_row <- anno_row_df[, -1, drop = FALSE]

      for (col in colnames(anno_row)) {
        anno_row[[col]] <- as.factor(anno_row[[col]])
      }

      for (col in colnames(anno_row)) {
        group_names <- unique(anno_row[[col]])
        annotation_colors[[col]] <- create_group_colors(group_names)
        names(annotation_colors[[col]]) <- group_names
      }
      cat(paste("Row annotations loaded:", paste(colnames(anno_row), collapse = ", "), "\n"))
    }

    # 处理排除样本
    if (!is.null(opt$exclude_sample)) {
      cat(paste("Excluding sample:", opt$exclude_sample, "\n"))
      d <- d[, !colnames(d) %in% c(opt$exclude_sample), drop = FALSE]
      cat(paste("Data dimensions after exclusion:", dim(d)[1], "x", dim(d)[2], "\n"))
    }

    # 行列转置：样本到行，contig到列
    cat("Transposing matrix: samples to rows, contigs to columns\n")
    d <- t(d)
    cat(paste("Transposed data dimensions:", dim(d)[1], "x", dim(d)[2], "\n"))

    # 获取转置后的contig名（列名），用于计算gap
    current_contigs <- colnames(d)
    sample_ids_for_gaps <- gsub("^([^_]+)_.*$", "\\1", current_contigs)

    # 获取转置后的样本顺序（行名）
    current_samples <- rownames(d)

    # 只保留有注释的样本
    if (!is.null(anno_row)) {
      common_samples <- intersect(rownames(anno_row), current_samples)
      if (length(common_samples) > 0) {
        anno_row <- anno_row[common_samples, , drop = FALSE]
        d <- d[common_samples, , drop = FALSE]
        cat(paste("Filtered to", length(common_samples), "samples with annotations\n"))
      } else {
        cat("Warning: No common samples between data and annotations\n")
        anno_row <- NULL
      }
    }

    # 根据sample_order重新排序
    cluster_rows <- TRUE
    if (!is.null(opt$sample_order)) {
      cat(paste("Reading sample order from:", opt$sample_order, "\n"))
      sample_order_df <- read.table(opt$sample_order, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      sample_order_list <- sample_order_df[[1]]
      common_samples <- intersect(rownames(d), sample_order_list)
      if (length(common_samples) > 0) {
        sample_order_list <- sample_order_list[sample_order_list %in% rownames(d)]
        d <- d[sample_order_list, , drop = FALSE]
        if (!is.null(anno_row)) {
          anno_row <- anno_row[sample_order_list[sample_order_list %in% rownames(anno_row)], , drop = FALSE]
        }
        cluster_rows <- FALSE
        cat(paste("Rows reordered according to sample list,", length(sample_order_list), "samples kept\n"))
      } else {
        cat("Warning: No common samples between data and sample order list\n")
      }
    }

    # 设置颜色断点
    breaks <- seq(0,1.0,by=0.001)

    # 设置输出文件路径
    output_prefix <- ifelse(is.null(opt$output_prefix), opt$input_file, opt$output_prefix)
    output_pdf <- paste0(output_prefix, ".pdf")
    output_png <- paste0(output_prefix, ".png")
    cat(paste("Output files:", output_pdf, "and", output_png, "\n"))

    # 计算列方向的gap（不同样本之间）
    sample_changes <- which(sample_ids_for_gaps[-1] != sample_ids_for_gaps[-length(sample_ids_for_gaps)])
    gaps_col <- sample_changes
    cat(paste("Adding column gaps at positions:", paste(gaps_col, collapse = ", "), "\n"))

    # 为列注释准备数据框
    anno_col_final <- data.frame(
      SampleID = sample_ids_for_gaps,
      stringsAsFactors = TRUE
    )
    rownames(anno_col_final) <- current_contigs

    anno_col_colors <- list(SampleID = sample_colors)

    cat("Generating heatmap...\n")
    pheatmap(d,
             scale = "none",
             cluster_rows = cluster_rows,
             cluster_cols = FALSE,
             clustering_method = "ward.D2",
             annotation_row = anno_row,
             annotation_col = anno_col_final,
             annotation_colors = c(annotation_colors, anno_col_colors),
             color = colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(length(breaks)),
             legend_breaks = seq(0,1.0,by=0.2),
             breaks = breaks,
             show_rownames = F,
             show_colnames = T,
             gaps_col = gaps_col,
             filename = output_pdf,
             width = opt$width,
             height = opt$height)

    pheatmap(d,
             scale = "none",
             cluster_rows = cluster_rows,
             cluster_cols = FALSE,
             clustering_method = "ward.D2",
             annotation_row = anno_row,
             annotation_col = anno_col_final,
             annotation_colors = c(annotation_colors, anno_col_colors),
             color = colorRampPalette(c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))(length(breaks)),
             legend_breaks = seq(0,1.0,by=0.2),
             breaks = breaks,
             show_rownames = F,
             show_colnames = T,
             gaps_col = gaps_col,
             filename = output_png,
             width = opt$width,
             height = opt$height)

    cat(paste("Heatmaps saved to:", output_pdf, "and", output_png, "\n"))

  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    traceback()
    stop(e)
  })
}

# 运行主函数
main()

