#!/usr/bin/env Rscript
# 功能富集结果气泡图绘制脚本
# 读取功能富集分析结果，绘制气泡图
#
# 参数说明：
# --enrichment_results: 富集分析结果文件路径 (stats_enrichment.R的输出)
# --output_prefix: 输出文件前缀，默认为'functional_enrichment'
# --top_n: 绘制top N显著功能，默认为20
# --annotation_label_file: 可选的注释标签文件路径，TSV格式，两列：功能ID和注释文字
# --sort_by: 排序方式，可选值：Gene_Ratio, P-value, Enrichment_Fold，默认为Gene_Ratio
# --title: 图表标题，默认为空
#
# 示例用法：
# Rscript plot_enrichment.R \
#   --enrichment_results functional_enrichment.full_results.tsv \
#   --output_prefix functional_enrichment \
#   --top_n 30 \
#   --sort_by Gene_Ratio

library(optparse)
library(tidyverse)
library(ggplot2)
library(ggrepel)

parse_arguments <- function() {
  option_list <- list(
    make_option(c("--enrichment_results"), type="character", default=NULL,
                help="富集分析结果文件路径 (stats_enrichment.R的输出)", metavar="character"),
    make_option(c("--output_prefix"), type="character", default="functional_enrichment",
                help="输出文件前缀，默认为'functional_enrichment'", metavar="character"),
    make_option(c("--top_n"), type="integer", default=20,
                help="绘制top N显著功能，默认为20", metavar="integer"),
    make_option(c("--annotation_label_file"), type="character", default=NULL,
                help="可选的注释标签文件路径，TSV格式，两列：功能ID和注释文字", metavar="character"),
    make_option(c("--sort_by"), type="character", default="Gene_Ratio",
                help="排序方式，可选值：Gene_Ratio, P-value, Enrichment_Fold，默认为Gene_Ratio", metavar="character"),
    make_option(c("--title"), type="character", default="",
                help="图表标题，默认为空", metavar="character")
  )
  
  parser <- OptionParser(option_list=option_list, 
                        description="功能富集结果气泡图绘制工具")
  
  args <- parse_args(parser)
  
  if (is.null(args$enrichment_results)) {
    print_help(parser)
    stop("缺少必需参数 --enrichment_results", call.=FALSE)
  }
  
  return(args)
}

wrap_label <- function(label, max_chars = 30) {
  if (nchar(label) <= max_chars) {
    return(label)
  }
  # 从第max_chars个位置向前找第一个空格
  substr_30 <- substr(label, 1, max_chars)
  space_pos <- regexpr(" ", substr(label, max_chars, nchar(label)), fixed = TRUE)
  if (space_pos > 0) {
    wrap_pos <- max_chars + space_pos - 1
    return(paste0(substr(label, 1, wrap_pos - 1), "\n", substr(label, wrap_pos + 1, nchar(label))))
  }
  return(label)
}

plot_bubble_chart <- function(test_results, output_prefix, top_n, annotation_label_file = NULL, sort_by = "Gene_Ratio", title = "") {
  cat("绘制气泡图...\n")
  
  filtered_results <- test_results %>% filter(Type == "Enrichment")
  cat("仅显示富集(Enrichment)结果\n")
  
  filtered_results <- filtered_results %>%
    mutate(
      log10_Pvalue = -log10(Pvalue),
      Gene_Ratio = (Subset_Functional_Hogs / Subset_Total_Hogs) * 100,
      Pvalue_Category = case_when(
        log10_Pvalue >= 3 ~ "3+",
        log10_Pvalue >= 2 ~ "2-3",
        log10_Pvalue >= 1 ~ "1-2",
        TRUE ~ "0-1"
      ),
      Pvalue_Category = factor(Pvalue_Category, levels = c("3+", "2-3", "1-2", "0-1"))
    )
  
  if (nrow(filtered_results) < top_n) {
    top_n <- nrow(filtered_results)
    cat(paste("注意：功能数量不足", top_n, "个，将绘制所有功能\n"))
  }
  
  if (!sort_by %in% c("Gene_Ratio", "P-value", "Enrichment_Fold")) {
    stop("排序参数必须是 Gene_Ratio, P-value 或 Enrichment_Fold")
  }
  
  if (sort_by == "Gene_Ratio") {
    top_results <- filtered_results %>%
      arrange(desc(Gene_Ratio)) %>%
      slice_head(n = top_n)
    cat(paste("按Gene_Ratio降序选择前", top_n, "个功能\n"))
  } else if (sort_by == "P-value") {
    top_results <- filtered_results %>%
      arrange(Pvalue) %>%
      slice_head(n = top_n)
    cat(paste("按P-value升序选择前", top_n, "个功能\n"))
  } else if (sort_by == "Enrichment_Fold") {
    top_results <- filtered_results %>%
      arrange(desc(Enrichment_Fold)) %>%
      slice_head(n = top_n)
    cat(paste("按Enrichment_Fold降序选择前", top_n, "个功能\n"))
  }
  
  if (!is.null(annotation_label_file)) {
    cat(paste("读取注释标签文件:", annotation_label_file, "\n"))
    annotation_labels <- read_tsv(annotation_label_file, col_names = c("Annotation", "Annotation_Text"), col_types = cols())
    
    top_results <- top_results %>%
      left_join(annotation_labels, by = "Annotation") %>%
      mutate(
        full_label = ifelse(
          !is.na(Annotation_Text),
          paste0(Annotation, " ", Annotation_Text),
          Annotation
        ),
        Annotation_Label = sapply(full_label, wrap_label)
      )
    
    cat(paste("成功合并", nrow(annotation_labels), "个注释标签\n"))
  } else {
    top_results <- top_results %>%
      mutate(Annotation_Label = sapply(Annotation, wrap_label))
  }
  
  if (nrow(top_results) == 0) {
    cat("没有符合条件的功能可绘制\n")
    return(NULL)
  }
  
  min_enrichment <- min(top_results$Enrichment_Fold)
  max_enrichment <- max(top_results$Enrichment_Fold)
  cat(paste("Enrichment Score范围：", round(min_enrichment, 2), "-" , round(max_enrichment, 2), "\n"))
  
  if (sort_by == "Gene_Ratio") {
    y_order <- reorder(top_results$Annotation_Label, -top_results$Gene_Ratio)
  } else if (sort_by == "P-value") {
    y_order <- reorder(top_results$Annotation_Label, top_results$Pvalue)
  } else if (sort_by == "Enrichment_Fold") {
    y_order <- reorder(top_results$Annotation_Label, -top_results$Enrichment_Fold)
  }
  
  if (title == "") {
    title <- "Functional Enrichment Analysis"
  }
  
  p <- ggplot(top_results, 
              aes(x = Gene_Ratio, 
                  y = y_order,
                  size = Enrichment_Fold,
                  color = Pvalue_Category)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("3+" = "#E64B35", "2-3" = "#F39B7F", "1-2" = "#4DBBD5", "0-1" = "#B09C85")) +
    scale_size(range = c(3, 10)) +
    scale_y_discrete(limits = rev) +
    labs(
      title = title,
      x = "Gene Ratio (%)",
      y = "",
      size = "Enrichment Fold",
      color = "-log10(P-value)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.y = element_text(hjust = 1, size = 7),
      axis.text.x = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 1, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
    ) +
    guides(size = guide_legend(order = 2), color = guide_legend(order = 1))
  
  output_file <- paste0(output_prefix, ".bubble_plot.png")
  ggsave(output_file, p, width = 10, height = max(10, top_n * 0.35 + 2), dpi = 300)
  cat(paste("气泡图已保存到", output_file, "\n"))
  
  pdf_output_file <- paste0(output_prefix, ".bubble_plot.pdf")
  ggsave(pdf_output_file, p, width = 10, height = max(10, top_n * 0.35 + 2))
  cat(paste("PDF已保存到", pdf_output_file, "\n"))
  
  return(p)
}

main <- function() {
  start_time <- Sys.time()
  
  cat("===== 功能富集气泡图绘制开始 =====\n")
  
  args <- parse_arguments()
  
  cat(paste("读取富集分析结果:", args$enrichment_results, "\n"))
  test_results <- read_tsv(args$enrichment_results, col_types = cols())
  cat(paste("读取完成，共", nrow(test_results), "条记录\n"))
  
  plot_bubble_chart(test_results, args$output_prefix, args$top_n, args$annotation_label_file, args$sort_by, args$title)
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  cat(paste("\n===== 气泡图绘制完成 =====\n"))
  cat(paste("总运行时间:", round(runtime, 2), "分钟\n"))
}

main()