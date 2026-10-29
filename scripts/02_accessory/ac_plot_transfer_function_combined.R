#!/usr/bin/env Rscript
# 组2: 转移功能组合图 (原 S16 + S17 + S19 -> 新 S16)
# nature-figure, backend = R
# 核心结论: AC 相关转移移动不同功能库 ——
#   a A-A (AC-AC) vs C-A (AC-core) 交换基因功能 (激酶 A-A 富集 / 转运·氧化还原·CAZyme C-A 富集)
#   b bulk vs scatter 转移 AC 货物功能 (bulk 富集转运/TF/effector/甲基转移酶)
#   c 种内 (solani/falciforme/keratoplastica) vs 跨种转移货物功能
suppressMessages({library(ggplot2); library(patchwork)})

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid = element_blank()))

RED <- "#C44E52"; BLUE <- "#3B7DD8"; GREY <- "#BBBBBB"; DARK <- "#444444"

# ---- a: A-A vs C-A function ----
a_df <- read.delim("06.accessory-chrs/gene_transfer_211/suppl_table.transfer_function.tsv")
a_df$Category <- gsub("Pfam_", "", a_df$Category)
a_df <- a_df[a_df$FDR < 0.05, ]
a_df$sig <- ifelse(a_df$Odds_ratio > 1, "A-A enriched", "C-A enriched")
a_df$Category <- factor(a_df$Category, levels = a_df$Category[order(a_df$Odds_ratio)])
p_a <- ggplot(a_df, aes(x = Odds_ratio, y = Category, fill = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = GREY) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("A-A enriched" = RED, "C-A enriched" = BLUE)) +
  geom_text(aes(label = sprintf("%.2f", Odds_ratio)), hjust = -0.1, size = 2.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Odds ratio (A-A vs C-A exchanged genes)", y = NULL,
       title = "Functional divergence of exchange modes (FDR < 0.05)") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.size = unit(0.25, "cm"))

# ---- b: bulk vs scatter ----
b_df <- read.delim("06.accessory-chrs/gene_transfer_211/suppl_table.ac_class_function.tsv")
b_df$Category <- gsub("Pfam_", "", b_df$Category)
b_df <- b_df[b_df$FDR < 0.05, ]
b_df$Category <- factor(b_df$Category, levels = rev(b_df$Category[order(b_df$Odds_ratio)]))
p_b <- ggplot(b_df, aes(x = Odds_ratio, y = Category)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = GREY) +
  geom_segment(aes(xend = 1, yend = Category), linewidth = 0.4, colour = "#888888") +
  geom_point(size = 2.8, colour = "#222222") +
  geom_text(aes(label = sprintf("%.1f", Odds_ratio)), hjust = -0.4, size = 2.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Odds ratio (bulk vs scatter-transfer ACs)",
       y = NULL, title = "Cargo of block-transfer (bulk) ACs (FDR < 0.05)")

# ---- c: intra-species vs cross-species function ----
c_full <- read.delim("06.accessory-chrs/gene_transfer_211/suppl_table.intra_species_function.tsv")
c_full$Category <- gsub("Pfam_", "", c_full$Category)
c_full$Mode_lab <- c("intra_solani" = "Intra F. solani",
                     "intra_falciformis" = "Intra F. falciforme",
                     "intra_keratoplastica" = "Intra F. keratoplastica",
                     "cross_species" = "Cross-species")[c_full$Mode]
c_full$Mode_lab <- factor(c_full$Mode_lab, levels = c("Intra F. solani", "Intra F. falciforme",
                                                      "Intra F. keratoplastica", "Cross-species"))
cats <- unique(c_full$Category)
modes <- levels(c_full$Mode_lab)
grid <- expand.grid(Category = cats, Mode_lab = modes, stringsAsFactors = FALSE)
c_df <- merge(grid, c_full[, c("Category", "Mode_lab", "Odds_ratio", "FDR")],
              by = c("Category", "Mode_lab"), all.x = TRUE)
c_df$sig <- !is.na(c_df$FDR) & c_df$FDR < 0.05
c_df$na <- is.na(c_df$Odds_ratio)
p_c <- ggplot(c_df, aes(x = Category, y = Mode_lab, fill = Odds_ratio)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_tile(data = c_df[c_df$na, ], fill = "#F5F5F5", colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 1,
                       na.value = "#F5F5F5", name = "Odds ratio") +
  geom_text(data = c_df[!c_df$na, ],
            aes(label = ifelse(sig, sprintf("%.2f*", Odds_ratio),
                               sprintf("%.2f", Odds_ratio))), size = 2.0) +
  labs(x = NULL, y = NULL,
       title = "Cargo of intra-species and cross-species transfers") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
        legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.8, "cm"),
        legend.title = element_text(size = 5.5))

# ---- combine ----
fig <- (p_a | p_b) / p_c +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))

out <- "06.accessory-chrs/gene_transfer_211/figure_combined_transfer_function"
ggsave(paste0(out, ".pdf"), fig, width = 8.5, height = 7, device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 8.5, height = 7, dpi = 300)
cat("saved ->", out, ".pdf/.png\n")
