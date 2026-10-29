setwd("/share/pasteur/zhuqh/fusarium/acc_chr_enrich")

#加载包
library(MicrobiotaProcess)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(patchwork)

# 获取当前目录下所有 *.ko.list 文件
ko_list_files <- list.files(pattern = "\\.ko\\.list$")

cat(sprintf("Found %d KO list files\n", length(ko_list_files)))

for (ko_file in ko_list_files) {
  cat(sprintf("\nProcessing: %s\n", ko_file))

  # 读取 ko 列表
  kolist <- read.table(ko_file)

  # 构建输出文件名（去掉 .ko.list 后缀）
  base_name <- sub("\\.ko\\.list$", "", ko_file)

  # KEGG 富集分析
  enrich_result <- enrichKEGG(
    gene = kolist$V1,
    organism = "ko",
    pvalueCutoff = 0.05,
    pAdjustMethod = "none"
  )

  # 保存富集分析结果
  enrich_tsv <- paste0(base_name, ".ko.enrich_result.tsv")
  write.table(enrich_result, file = enrich_tsv, row.names = F, quote = F, sep = "\t")

  # 检查是否有富集结果
  if (is.null(enrich_result) || nrow(enrich_result) == 0) {
    cat(sprintf("No significant enrichment results for %s\n", ko_file))
    next
  }

  # 绘制 dotplot 并保存为 PDF
  pdf_file <- paste0(base_name, ".ko.enrichment.pdf")
  pdf(pdf_file, width = 6, height = 6)

  p <- dotplot(
    enrich_result,
    x = "RichFactor",
    orderBy = "x",
    size = "GeneRatio",
    color = "pvalue",
    showCategory = 20,
    label_format = 40
  ) +
    ggtitle(base_name) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )

  print(p)
  dev.off()

  cat(sprintf("Enrichment result saved to: %s\n", enrich_tsv))
  cat(sprintf("Plot saved to: %s\n", pdf_file))
}

cat(sprintf("\nAll done. Processed %d files.\n", length(ko_list_files)))
