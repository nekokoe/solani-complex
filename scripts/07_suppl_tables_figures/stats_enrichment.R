#!/usr/bin/env Rscript
# 功能富集分析计算脚本
# 分析子集HOG与背景HOG之间的功能注释富集或耗竭情况
# 直接按照给定的子集HOG列表和背景集HOG列表进行统计分析
#
# 参数说明：
# --annotation_file: HOG注释文件路径 (hog_annotation_integrator_v2.py的输出)
# --annotation_column: 要分析的注释列名，默认为'KEGG_Pathway'
# --subset_hog_list: 子集HOG列表文件路径，每行一个HOG_ID
# --background_hog_list: 背景HOG列表文件路径，每行一个HOG_ID
# --output_prefix: 输出文件前缀，默认为'functional_enrichment'
#
# 示例用法：
# Rscript stats_enrichment.R \
#   --annotation_file hog_annotation.tsv \
#   --annotation_column KEGG_Pathway \
#   --subset_hog_list subset_hogs.txt \
#   --background_hog_list background_hogs.txt \
#   --output_prefix functional_enrichment

library(optparse)
library(tidyverse)

parse_arguments <- function() {
  option_list <- list(
    make_option(c("--annotation_file"), type="character", default=NULL,
                help="HOG注释文件路径 (hog_annotation_integrator_v2.py的输出)", metavar="character"),
    make_option(c("--annotation_column"), type="character", default="KEGG_Pathway",
                help="要分析的注释列名，默认为'KEGG_Pathway'", metavar="character"),
    make_option(c("--subset_hog_list"), type="character", default=NULL,
                help="子集HOG列表文件路径，每行一个HOG_ID", metavar="character"),
    make_option(c("--background_hog_list"), type="character", default=NULL,
                help="背景HOG列表文件路径，每行一个HOG_ID", metavar="character"),
    make_option(c("--output_prefix"), type="character", default="functional_enrichment",
                help="输出文件前缀，默认为'functional_enrichment'", metavar="character")
  )
  
  parser <- OptionParser(option_list=option_list, 
                        description="功能富集分析工具，直接按照给定的子集HOG列表和背景集HOG列表分析功能富集情况")
  
  args <- parse_args(parser)
  
  if (is.null(args$annotation_file) || is.null(args$subset_hog_list) || is.null(args$background_hog_list)) {
    print_help(parser)
    stop("缺少必需参数", call.=FALSE)
  }
  
  return(args)
}

read_and_process_annotation <- function(annotation_file, annotation_column) {
  cat("读取并处理注释文件...\n")
  
  annotation_df <- read_tsv(annotation_file, col_types = cols())
  
  if (!(annotation_column %in% colnames(annotation_df))) {
    stop(paste("错误：注释文件中不存在列", annotation_column), call.=FALSE)
  }
  
  cat(paste("注释文件总行数:", nrow(annotation_df), "\n"))
  cat(paste("注释文件中唯一HOG_ID数量:", n_distinct(annotation_df$HOG_ID), "\n"))
  
  annotation_data <- annotation_df %>%
    select(HOG_ID, !!sym(annotation_column))
  
  cat(paste("过滤前的HOG-注释对数量:", nrow(annotation_data), "\n"))
  
  annotation_data <- annotation_data %>%
    filter(!!sym(annotation_column) != "-", !is.na(!!sym(annotation_column)), !!sym(annotation_column) != "")
  
  cat(paste("过滤无效注释后的HOG-注释对数量:", nrow(annotation_data), "\n"))
  cat(paste("过滤后唯一HOG_ID数量:", n_distinct(annotation_data$HOG_ID), "\n"))
  
  result <- annotation_data %>%
    mutate(annotations = str_split(!!sym(annotation_column), "[,;]", simplify = FALSE)) %>%
    unnest(annotations) %>%
    mutate(annotations = str_trim(annotations)) %>%
    filter(annotations != "") %>%
    rename(Annotation = annotations) %>%
    distinct(HOG_ID, Annotation)
  
  cat(paste("处理完成，共得到", nrow(result), "条HOG-注释映射\n"))
  cat(paste("处理后唯一HOG_ID数量:", n_distinct(result$HOG_ID), "\n"))
  return(result)
}

read_hog_list <- function(hog_list_file) {
  cat("读取HOG列表文件...\n")
  
  hogs <- read_lines(hog_list_file)
  hogs <- str_trim(hogs)
  hogs <- hogs[hogs != ""]
  
  cat(paste("读取完成，共", length(hogs), "个HOG\n"))
  return(hogs)
}

calculate_functional_stats <- function(annotation_data, subset_hog_list, background_hog_list) {
  subset_total_hogs <- length(subset_hog_list)
  background_total_hogs <- length(background_hog_list)
  
  cat(paste("子集总HOG数:", subset_total_hogs, "\n"))
  cat(paste("背景总HOG数:", background_total_hogs, "\n"))
  
  background_functional_hogs <- annotation_data %>%
    filter(HOG_ID %in% background_hog_list) %>%
    group_by(Annotation) %>%
    summarise(Background_Functional_Hogs = n_distinct(HOG_ID), .groups = "drop")
  
  subset_functional_hogs <- annotation_data %>%
    filter(HOG_ID %in% subset_hog_list) %>%
    group_by(Annotation) %>%
    summarise(Subset_Functional_Hogs = n_distinct(HOG_ID), .groups = "drop")
  
  functional_stats <- subset_functional_hogs %>%
    full_join(background_functional_hogs, by = "Annotation") %>%
    mutate(
      Background_Total_Hogs = background_total_hogs,
      Subset_Total_Hogs = subset_total_hogs,
      Subset_Functional_Hogs = replace_na(Subset_Functional_Hogs, 0),
      Background_Functional_Hogs = replace_na(Background_Functional_Hogs, 0)
    ) %>%
    filter(!(Subset_Functional_Hogs == 0 & Background_Functional_Hogs == 0))
  
  cat("\n功能统计数据摘要：\n")
  cat("- 功能注释总数:", nrow(functional_stats), "\n")
  cat("- 子集功能HOG数范围:", min(functional_stats$Subset_Functional_Hogs), "-", max(functional_stats$Subset_Functional_Hogs), "\n")
  cat("- 背景功能HOG数范围:", min(functional_stats$Background_Functional_Hogs), "-", max(functional_stats$Background_Functional_Hogs), "\n")
  
  cat(paste("计算完成，共", nrow(functional_stats), "个功能\n"))
  return(functional_stats)
}

perform_hypergeometric_test <- function(functional_stats) {
  cat("执行超几何分布检验...\n")
  
  result <- functional_stats %>%
    mutate(
      k = pmax(0, round(Subset_Functional_Hogs)),
      m = pmax(0, round(Background_Functional_Hogs)),
      n = pmax(0, round(Background_Total_Hogs - Background_Functional_Hogs)),
      K = pmax(1, round(Subset_Total_Hogs)),
      
      is_valid_parameters = m >= 0 & n >= 0 & K > 0 & (m + n) > 0,
      
      Enrichment_Pvalue = ifelse(is_valid_parameters,
                               phyper(q = pmax(0, k - 1), m = m, n = n, k = pmin(K, m + n), lower.tail = FALSE),
                               NA),
      
      Enrichment_Fold = ifelse(is_valid_parameters,
                             (k / K) / (m / (m + n)),
                             NA),
      
      Pvalue = ifelse(is_valid_parameters, Enrichment_Pvalue, NA),
      
      Type = ifelse(is_valid_parameters, "Enrichment", NA)
    ) %>%
    mutate(
      Enrichment_FDR = p.adjust(Enrichment_Pvalue, method = "fdr"),
      FDR_Adjusted_Pvalue = Enrichment_FDR,
      log10_FDR = -log10(FDR_Adjusted_Pvalue)
    ) %>%
    arrange(Pvalue)
  
  valid_params_count = sum(result$is_valid_parameters, na.rm = TRUE)
  invalid_params_count = nrow(result) - valid_params_count
  
  if (invalid_params_count > 0) {
    cat(paste("警告：有", invalid_params_count, "个功能注释的超几何分布参数无效，已标记为NA\n"))
  }
  cat(paste("超几何分布检验完成，", valid_params_count, "个功能注释通过参数有效性检查\n"))
  
  result = result %>%
    select(-starts_with("is_valid_"))
  
  return(result)
}

save_results <- function(test_results, output_prefix) {
  cat("保存分析结果...\n")
  
  selected_columns <- test_results %>% 
    select(
      Annotation,
      Subset_Functional_Hogs,
      Background_Functional_Hogs,
      Subset_Total_Hogs,
      Background_Total_Hogs,
      k, m, n, K,
      Enrichment_Fold,
      Pvalue,
      Enrichment_Pvalue,
      Type,
      FDR_Adjusted_Pvalue,
      log10_FDR
    )
  
  output_file <- paste0(output_prefix, ".full_results.tsv")
  write_tsv(selected_columns, output_file)
  cat(paste("完整结果已保存到", output_file, "\n"))
  
  significant_results <- selected_columns %>% filter(Pvalue < 0.05 & Type == "Enrichment")
  cat("仅保存显著富集(Enrichment)结果\n")
  
  if (nrow(significant_results) > 0) {
    sig_output_file <- paste0(output_prefix, ".significant_results.tsv")
    write_tsv(significant_results, sig_output_file)
    cat(paste("显著结果已保存到", sig_output_file, "(共", nrow(significant_results), "条)\n"))
  } else {
    cat("未发现符合条件的显著结果 (Pvalue < 0.05)\n")
  }
}

main <- function() {
  start_time <- Sys.time()
  
  cat("===== 功能富集分析开始 =====\n")
  
  args <- parse_arguments()
  
  annotation_data <- read_and_process_annotation(args$annotation_file, args$annotation_column)
  subset_hog_list <- read_hog_list(args$subset_hog_list)
  background_hog_list <- read_hog_list(args$background_hog_list)
  
  functional_stats <- calculate_functional_stats(annotation_data, subset_hog_list, background_hog_list)
  
  test_results <- perform_hypergeometric_test(functional_stats)
  
  save_results(test_results, args$output_prefix)
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  cat(paste("\n===== 功能富集分析完成 =====\n"))
  cat(paste("总运行时间:", round(runtime, 2), "分钟\n"))
}

main()