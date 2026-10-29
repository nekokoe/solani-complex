#!/usr/bin/env Rscript
# ac_plot_transfer_mechanism_combined_v2.R — S13 (a-e), rebuilt around the new
# panel a.
#
# Only panel a changes: it now breaks the transfer counts out by species and
# labels each C-A bar with its share of that group's accessory-chromosome genes.
# Panels b-e are copied verbatim from ac_plot_transfer_mechanism_combined.R, so
# their appearance is untouched. The layout is re-tuned because panel a is now
# square-ish rather than wide.

suppressMessages({library(ggplot2); library(patchwork); library(reshape2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s13a_panel.R"))

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid = element_blank()))

RED <- S13_RED; BLUE <- S13_BLUE

# ---- a: AC-related gene transfer, by species (new) ----
p_a <- s13a_panel(s13a_table("06.accessory-chrs/gene_transfer_211/suppl_table.gene_transfer_211.tsv"),
                  fontsize = 8)

# ---- b: TE flanking enrichment ----
b_df <- read.delim("06.accessory-chrs/gene_transfer_211/suppl_table.repeat_enrichment_211.tsv")
b_df$flank <- ifelse(b_df$Metric == "flank_1kb", "1 kb", "5 kb")
b_df$Family <- factor(b_df$Family, levels = rev(c("DNA", "LTR", "LINE", "RC", "Unknown")))
p_b <- ggplot(b_df, aes(x = Odds_ratio, y = Family, fill = flank)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("1 kb" = RED, "5 kb" = BLUE)) +
  geom_text(aes(label = ifelse(Odds_ratio >= 4, sprintf("%.1f", Odds_ratio), "")),
            position = position_dodge(width = 0.7), hjust = -0.1, size = 2.0) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Odds ratio (transfer vs whole-genome flanking)",
       y = NULL, fill = "Flank", title = "TE-flanking enrichment of transferred genes") +
  guides(fill = guide_legend(title = "Flank"))

# ---- c: fine-family ranking (1kb + 5kb dual dimension) ----
c_df <- read.delim("06.accessory-chrs/gene_transfer_211/suppl_table.repeat_enrichment_211.fine_both_flanks.tsv")
c_top <- c_df[order(c_df$OR_1kb, decreasing = TRUE), ][1:10, ]
c_top <- c_top[order(c_top$OR_1kb), ]
c_top$Fine_family <- factor(c_top$Fine_family, levels = c_top$Fine_family)
c_long <- data.frame(
  Fine_family = rep(c_top$Fine_family, 2),
  flank = factor(rep(c("1 kb", "5 kb"), each = nrow(c_top)),
                 levels = c("1 kb", "5 kb")),
  OR = c(c_top$OR_1kb, c_top$OR_5kb))
p_c <- ggplot(c_long, aes(x = OR, y = Fine_family, colour = flank)) +
  geom_point(size = 2.0, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("1 kb" = RED, "5 kb" = BLUE)) +
  labs(x = "Odds ratio (flank, top families)", y = NULL, colour = "Flank",
       title = "Repeat families most associated with transfer") +
  theme(legend.position = "bottom", legend.key.size = unit(0.25, "cm"))

# ---- d: cassette test ----
d_df <- read.delim("06.accessory-chrs/gene_transfer/suppl_table.cassette_test.tsv")
d_df <- d_df[d_df$Cluster == "cluster100", ]
d_df$lab <- ifelse(d_df$Specific_TE_only, "Specific TE pairs", "All TE pairs")
d_df$lab <- factor(d_df$lab, levels = c("All TE pairs", "Specific TE pairs"))
p_d <- ggplot(d_df, aes(x = lab)) +
  geom_col(aes(y = observed_frac * 100, fill = "Observed"), width = 0.5) +
  geom_point(aes(y = null_mean_frac * 100, colour = "Permutation null"), size = 2.4) +
  geom_errorbar(aes(ymin = (null_mean_frac - null_sd_frac) * 100,
                    ymax = (null_mean_frac + null_sd_frac) * 100,
                    colour = "Permutation null"),
                width = 0.15) +
  scale_fill_manual(values = c("Observed" = RED)) +
  scale_colour_manual(values = c("Permutation null" = BLUE)) +
  geom_text(aes(y = observed_frac * 100, label = sprintf("%.1f%%", observed_frac * 100)),
            vjust = -0.6, size = 2.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Clusters sharing the same TE bracket (%)",
       title = "TE-bracketed mobile cassettes (100% identity clusters)") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 5.5),
        legend.key.size = unit(0.25, "cm"))

# ---- e: TE-pair heatmap (upstream x downstream flanking families) ----
pair <- as.matrix(read.delim("06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.tsv",
                             row.names = 1, check.names = FALSE))
# Denominator = the number of TE-BRACKETED genes (both flanks a repeat family
# other than Unknown/Simple_repeat/Low_complexity), read from the companion file.
# It is NOT sum(pair): the 12 x 12 grid keeps only pairs whose two families are both
# inside the top 12, so it covers 501 of the 581 bracketed genes (86%).
pair_total <- as.numeric(read.delim(
  "06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.total.tsv",
  check.names = FALSE)$te_bracketed_genes[1])
pair_pct <- 100 * pair / pair_total
pair_floor <- min(pair_pct[pair_pct > 0]) / 2
pair_l <- reshape2::melt(pmax(pair_pct, pair_floor),
                         varnames = c("up", "down"), value.name = "pct")
short <- function(s) gsub("DNA/|RC/|LINE/|LTR/", "", s)
p_e <- ggplot(pair_l, aes(x = down, y = up, fill = pct)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = RED, trans = "log10",
                      name = "% of genes\n(log)") +
  scale_x_discrete(labels = short) + scale_y_discrete(labels = short) +
  labs(x = "Downstream flanking TE family", y = "Upstream flanking TE family",
       title = "Bracketing TE pairs, as % of all bracketed genes (5 kb)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4.5),
        axis.text.y = element_text(size = 4.5),
        legend.position = "right", legend.key.height = unit(0.4, "cm"),
        legend.title = element_text(size = 5.5))

# ---- combine: 2x2 + full-width heatmap; row 1 gets the extra height panel a needs
fig <- (p_a | p_b) / (p_c | p_d) / p_e +
  plot_layout(heights = c(1.25, 0.95, 0.95)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))

out <- "06.accessory-chrs/gene_transfer_211/figure_combined_transfer_mechanism_v2"
ggsave(paste0(out, ".pdf"), fig, width = 9.5, height = 11, device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 9.5, height = 11, dpi = 300)
cat("saved ->", out, ".pdf/.png\n")
