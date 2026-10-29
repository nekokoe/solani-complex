#!/usr/bin/env Rscript
# S21: BMU02483 特有基因是否来自 clinical AC? (nature-figure, backend = R)
# 核心结论:
#   (1) 170 个 BMU02483 特有 HOG (vs tree-order 相邻 17 个 F. falciformis 环境株)
#       中 12.4% (21) 落在临床 mosaic AC (BMU03267/BMU01933) 的基因集,
#       是相邻环境株对照 (4.2%) 的 ~3 倍 (OR 4.3, P = 1.6e-7) —
#       "邻居没有 + clinical 有"的基因在 2483 中显著富集, 支持 AC 获得。
#   (2) 170 特有中 40% (68) 落在 13 TGS 的 accessory chromosomes 基因集
#       (背景 13.8%, OR 4.2), 高于相邻 env 对照 (16%)。
#   (3) 21 个 clinical-AC 特有基因中 12 个为低频 (presence <=30%),
#       含 MFS 转运蛋白、HET、CFEM、RTA1 等候选临床适应性基因。
# 面板:
#   a: 170 特有 vs 相邻 env 对照: clinical-AC 基因占比 (核心)
#   b: 21 个 clinical-AC 特有基因中低频候选基因的功能类
#   c: 170 特有中 AC vs 非 AC 占比 (vs env 对照)
#   d: AC -> BMU02483 基因流向 sankey
suppressMessages({library(ggplot2); library(patchwork); library(ggalluvial)})

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid = element_blank()))

RED <- "#C44E52"; BLUE <- "#3B7DD8"; GOLD <- "#E0B74E"; GREY <- "#BBBBBB"; DARK <- "#444444"

# ---- a: clinical-AC fraction in 2483-specific vs env control ----
a_df <- data.frame(
  group = factor(c("BMU02483-specific\n(170 HOGs)", "Nearby env strains\n(20,133 HOGs)",
                   "Genome background\n(29,080 HOGs)"),
                 levels = c("BMU02483-specific\n(170 HOGs)", "Nearby env strains\n(20,133 HOGs)",
                            "Genome background\n(29,080 HOGs)")),
  pct = c(12.4, 4.2, 3.2),
  n_in = c(21, 838, 939),
  n_all = c(170, 20133, 29080))
p_a <- ggplot(a_df, aes(x = group, y = pct, fill = group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(RED, BLUE, GREY)) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct, n_in, n_all)),
            vjust = -0.3, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL,
       y = "HOGs in clinical-AC gene set (%)",
       title = "BMU02483-specific HOGs enriched for clinical-AC content") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6.5, lineheight = 0.8))

# ---- b: low-frequency clinical-AC candidate genes ----
b_df <- data.frame(
  Class = c("HET domain", "MFS transporter", "CFEM", "RTA1",
            "Ankyrin", "NACHT", "Other"),
  n = c(2, 1, 1, 1, 1, 1, 5))
b_df$Class <- factor(b_df$Class, levels = b_df$Class[order(b_df$n)])
p_b <- ggplot(b_df, aes(x = Class, y = n, fill = Class)) +
  geom_col(show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "n low-frequency genes\nin clinical-AC set",
       title = "Clinical-adaptive candidates among the 21 specific HOGs") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5))

# ---- c: AC vs non-AC in 2483-specific vs env control ----
c_df <- data.frame(
  group = factor(c("BMU02483-specific\n(170 HOGs)", "Nearby env strains\n(20,133 HOGs)",
                   "Genome background\n(29,080 HOGs)"),
                 levels = c("BMU02483-specific\n(170 HOGs)", "Nearby env strains\n(20,133 HOGs)",
                            "Genome background\n(29,080 HOGs)")),
  pct_ac = c(40.0, 16.0, 13.8),
  n_ac = c(68, 3229, 4002),
  n_all = c(170, 20133, 29080))
p_c <- ggplot(c_df, aes(x = group, y = pct_ac, fill = group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c(RED, BLUE, GREY)) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct_ac, n_ac, n_all)),
            vjust = -0.3, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL,
       y = "HOGs on accessory chromosomes (%)",
       title = "Specific HOGs are enriched on accessory chromosomes") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6.5, lineheight = 0.8))

# ---- d: sankey (AC -> BMU02483 gene flow) ----
sd <- read.delim("06.accessory-chrs/BMU02483_gene_content/sankey_ac_shared.tsv")
v3 <- read.delim("06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv")
tmap <- setNames(v3$kmer_type, v3$AC)
sd$type <- ifelse(sd$AC == "Other ACs", "Other",
                  ifelse(tmap[sd$AC] == "type1_stable", "Block", "Mosaic"))
sd$type <- factor(sd$type, levels = c("Mosaic", "Block", "Other"))
sd$AC <- factor(sd$AC, levels = sd$AC[order(sd$shared)])
allu <- data.frame(source = as.character(sd$AC), target = "BMU02483",
                   weight = sd$shared, type = as.character(sd$type))
p_d <- ggplot(allu, aes(axis1 = source, axis2 = target, y = weight)) +
  geom_alluvium(aes(fill = type), width = 1/8, alpha = 0.85, curve_type = "linear") +
  geom_stratum(width = 1/8, fill = "grey92", colour = "grey40", linewidth = 0.3) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.2, nudge_x = -0.12, hjust = 1) +
  scale_fill_manual(values = c("Mosaic" = RED, "Block" = BLUE, "Other" = GREY),
                    name = "AC type") +
  scale_x_discrete(expand = c(0.22, 0.1)) +
  labs(title = "Origin of genes shared with BMU02483") +
  theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"),
        axis.text.x = element_text(size = 7, face = "bold"),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank())

# ---- combine: a|b / c / d ----
fig <- (p_a | p_b) / p_c / p_d +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 10))

out <- "06.accessory-chrs/BMU02483_gene_content/figure_BMU2483_specific_origin"
ggsave(paste0(out, ".pdf"), fig, width = 9, height = 9, device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 9, height = 9, dpi = 300)
cat("saved ->", out, ".pdf/.png (4 panels)\n")
