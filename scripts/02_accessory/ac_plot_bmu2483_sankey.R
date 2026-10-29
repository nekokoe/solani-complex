#!/usr/bin/env Rscript
# S21b (候选): BMU02483 从临床 AC 获得基因的 sankey 图 (nature-figure, backend = R)
# 核心: 展示 BMU02483 基因组中共享基因的 AC 来源流向 (Top ACs + Other),
#       宽度 ∝ 共享基因数, 按 block/mosaic 类型着色。
suppressMessages({library(ggalluvial); library(ggplot2)})

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid = element_blank()))

RED <- "#C44E52"; BLUE <- "#3B7DD8"; GREY <- "#BBBBBB"; DARK <- "#444444"

sd <- read.delim("06.accessory-chrs/BMU02483_gene_content/sankey_ac_shared.tsv")
# AC 类型
v3 <- read.delim("06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv")
tmap <- setNames(v3$kmer_type, v3$AC)
sd$type <- ifelse(sd$AC == "Other ACs", "Other",
                  ifelse(tmap[sd$AC] == "type1_stable", "Block", "Mosaic"))
sd$type <- factor(sd$type, levels = c("Mosaic", "Block", "Other"))
# AC 顺序 (按 shared 降序)
sd$AC <- factor(sd$AC, levels = sd$AC[order(sd$shared)])

# alluvial 数据: 两列 (source, target) + weight
allu <- data.frame(
  source = as.character(sd$AC),
  target = "BMU02483",
  weight = sd$shared,
  type = as.character(sd$type),
  stringsAsFactors = FALSE)

# 左侧标签: AC 名 + 共享基因数; 右侧: BMU02483
sd$src_lab <- paste0(sd$AC, " (", sd$shared, ")")
src_labs <- setNames(sd$src_lab, sd$AC)

p <- ggplot(allu, aes(axis1 = source, axis2 = target, y = weight)) +
  geom_alluvium(aes(fill = type), width = 1/6, alpha = 0.85, curve_type = "linear") +
  geom_stratum(width = 1/6, fill = "grey92", colour = "grey40", linewidth = 0.3) +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 2.3, nudge_x = -0.05, hjust = 1) +
  scale_fill_manual(values = c("Mosaic" = RED, "Block" = BLUE, "Other" = GREY),
                    name = "AC type") +
  scale_x_discrete(expand = c(0.15, 0.15)) +
  labs(title = "Origin of genes shared between BMU02483 and clinical accessory chromosomes") +
  theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"),
        axis.text.x = element_text(size = 7, face = "bold"))

out <- "06.accessory-chrs/BMU02483_gene_content/figure_BMU2483_sankey"
ggsave(paste0(out, ".pdf"), p, width = 7, height = 5.5, device = cairo_pdf)
ggsave(paste0(out, ".png"), p, width = 7, height = 5.5, dpi = 300)
cat("saved ->", out, ".pdf/.png\n")
