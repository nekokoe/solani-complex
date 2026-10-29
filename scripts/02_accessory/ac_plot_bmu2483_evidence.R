#!/usr/bin/env Rscript
# AC-T1: BMU02483 基因成分证据图 (nature-figure, backend = R)
# 科学问题: BMU02483 与相邻 falciformis 环境株相比特异携带的 HOG,
# 是否显著富集在 clinical AC 的稀有基因集上?
# 面板:
#   a: BMU02483 vs falciformis env 株的 HOG 携带 (特异 145 / 共享 / 缺失)
#   b: 特异 HOG 在 clinical AC 稀有基因集 vs 背景的富集 (Fisher OR + 阈值扫描)
#   c: 落在 clinical AC 基因集的特异 HOG 功能注释 (Pfam 类)
#   d: 与 kmer 证据对照 (5.9x RF 富集) 的整合视图
suppressMessages({library(ggplot2); library(patchwork)})

BASE <- "06.accessory-chrs/BMU02483_gene_content"

# ---- 数据 ----
stats <- read.delim(file.path(BASE, "bm2483_gene_content_stats.tsv"), header = FALSE,
                    col.names = c("key", "value"))
gv <- function(k) stats$value[stats$key == k]
n_spec <- as.numeric(gv("bm2483_specific_hogs"))
n_carried <- as.numeric(gv("bm2483_carried_hogs"))
n_env <- as.numeric(gv("falc_env_count"))
n_clin_target <- 156

enr <- read.delim(file.path(BASE, "bm2483_enrichment_results.tsv"))
# 主检验 (freq<=30% 最优)
main <- enr[enr$analysis == "specific_vs_rare_clinicalAC_30pct", ]
or_main <- main$odds_ratio
p_main <- main$fisher_p

spec <- read.delim(file.path(BASE, "bm2483_specific_hogs_annotated.tsv"))
spec$in_clinical_AC <- as.logical(spec$in_clinical_AC)

# ---- a: HOG carriage composition ----
frac_shared <- n_carried - n_spec
a_df <- data.frame(
  category = c("BMU02483-specific\n(absent in all 18 env)", 
               "Shared with env\nisolates"),
  n = c(n_spec, frac_shared),
  stringsAsFactors = FALSE)

# ---- b: enrichment across frequency thresholds ----
b_df <- enr[grep("specific_vs_rare_clinicalAC_", enr$analysis), ]
b_df$threshold <- c(20, 30, 40, 50)
b_df$sig <- b_df$fisher_p < 0.05

# ---- c: functional classes of the 8 specific HOGs in clinical AC ----
spec_in_ac <- spec[spec$in_clinical_AC, ]
func_class <- function(pfam) {
  pf <- tolower(pfam)
  if (grepl("pkinase|kinase", pf)) return("Protein kinase")
  if (grepl("ank", pf)) return("Ankyrin")
  if (grepl("het", pf)) return("HET")
  if (grepl("transp|permease|mfs", pf)) return("Transporter")
  if (grepl("oxidoreduct|dehydrogen|adh_|fad_|nad_|p450|monooxygen|reductase", pf)) return("Oxidoreductase")
  if (grepl("glyco_hydro|cazy|hydrolase", pf)) return("CAZyme/hydrolase")
  "Other"
}
c_df <- data.frame(HOG = spec_in_ac$HOG_ID, 
                   Class = vapply(spec_in_ac$Pfam_Domain, func_class, character(1)))

# ---- d: evidence integration ----
d_df <- data.frame(
  evidence = c("RF clinical HOG enrichment (kmer)", 
               "Specific HOGs in clinical-AC rare genes (this study)"),
  enrichment = c(5.9, round(or_main, 1)),
  p = c("P < 0.001", sprintf("P = %.1e", p_main)))

# ---- plot ----
red <- "#C44E52"; blue <- "#3B7DD8"; grey <- "#BBBBBB"; dark <- "#444444"

p_a <- ggplot(a_df, aes(x = category, y = n, fill = category)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c(red, grey)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  geom_text(aes(label = format(n, big.mark = ",")), vjust = -0.4, size = 3.5) +
  labs(x = NULL, y = "HOGs carried by BMU02483",
       title = "Gene-content comparison vs 18 F. falciformis environmental isolates") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(size = 7.5), plot.title = element_text(size = 8.5))

p_b <- ggplot(b_df, aes(x = threshold, y = odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = grey) +
  geom_line(colour = blue, linewidth = 0.7) +
  geom_point(aes(colour = sig), size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = red, "FALSE" = grey)) +
  scale_x_continuous(breaks = c(20, 30, 40, 50), limits = c(18, 52)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  geom_text(aes(label = sprintf("P=%.0e", fisher_p)), hjust = -0.05, vjust = -0.8, size = 2.4) +
  labs(x = "Max presence rate of clinical-AC gene set (%)",
       y = "Odds ratio (BMU02483-specific HOGs)",
       title = "Enrichment of BMU02483-specific HOGs in clinical-AC rare genes") +
  theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 8.5))

p_c <- ggplot(c_df, aes(x = Class, fill = Class)) +
  geom_bar(show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "n specific HOGs in clinical AC",
       title = "Functional classes of the 8 BMU02483-specific HOGs on clinical ACs") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 7),
        plot.title = element_text(size = 8.5))

p_d <- ggplot(d_df, aes(x = evidence, y = enrichment, fill = evidence)) +
  geom_col(width = 0.55, show.legend = FALSE) +
  scale_fill_manual(values = c(blue, red)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  geom_text(aes(label = sprintf("%.1f\u00d7  %s", enrichment, p)),
            vjust = -0.5, size = 2.8) +
  labs(x = NULL, y = "Enrichment (fold)",
       title = "Independent evidence for clinical-AC acquisition") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(size = 6.5), plot.title = element_text(size = 8.5))

combo <- (p_a / p_b) | (p_c / p_d)
fig <- combo + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 11))

out <- file.path(BASE, "figure_BMU2483_gene_evidence")
ggsave(paste0(out, ".pdf"), fig, width = 10, height = 7.5, device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 10, height = 7.5, dpi = 300)
cat("saved ->", out, ".pdf/.png\n")
cat("main: OR =", or_main, "P =", p_main, "\n")
