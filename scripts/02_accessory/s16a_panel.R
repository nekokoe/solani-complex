#!/usr/bin/env Rscript
# s16a_panel.R — S16 panel a as a standalone square panel: functional divergence
# between the two exchange modes (A-A vs C-A exchanged genes).
#
# Same conventions as s13a_panel.R / s13e_panel.R / s7a_panel.R: square canvas,
# font sizes driven by one number, NimbusSans (Arial resolves to it on this host),
# and the colours the S16 assembly already uses (RED = A-A enriched, BLUE = C-A
# enriched).
#
# The message: A-A exchange moves kinases, C-A exchange moves transporters,
# oxidoreductases and CAZymes (10 of 24 categories survive FDR < 0.05).
#
# Source: 06.accessory-chrs/gene_transfer_211/suppl_table.transfer_function.tsv

suppressMessages({library(ggplot2)})

S16_RED <- "#C44E52"; S16_BLUE <- "#3B7DD8"; S16_GREY <- "#BBBBBB"

read_transfer_function <- function(path, fdr = 0.05) {
  df <- read.delim(path, stringsAsFactors = FALSE)
  need <- c("Category", "Odds_ratio", "FDR")
  if (!all(need %in% names(df))) stop("missing columns: ",
                                      paste(setdiff(need, names(df)), collapse = ", "))
  df <- df[df$FDR < fdr, ]
  df$Category <- gsub("Pfam_", "", df$Category)
  df$sig <- ifelse(df$Odds_ratio > 1, "A-A enriched", "C-A enriched")
  # keep the published ordering (ascending odds ratio) so the panel matches S16
  df$Category <- factor(df$Category, levels = df$Category[order(df$Odds_ratio)])
  df
}

s16a_panel <- function(df, family = "Arial", fontsize = 9,
                       # the one-line form needs ~118 mm at 12 pt, wider than the
                       # panel's text column, so it is wrapped
                       title = "Functional divergence of exchange modes\n(FDR < 0.05)") {
  ggplot(df, aes(x = Odds_ratio, y = Category, fill = sig)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = S16_GREY) +
    geom_col(width = 0.62) +
    scale_fill_manual(values = c("A-A enriched" = S16_RED,
                                 "C-A enriched" = S16_BLUE)) +
    geom_text(aes(label = sprintf("%.2f", Odds_ratio)), hjust = -0.12,
              size = (fontsize - 0.6) / 2.845) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = "Odds ratio (A-A vs C-A exchanged genes)", y = NULL, title = title) +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(size = fontsize, colour = "black"),
      axis.text.y = element_text(size = fontsize, colour = "black"),
      axis.title.x = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = fontsize),
      legend.key.size = unit(4, "mm"),
      legend.margin = margin(t = -1),
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 3, 4)
    )
}
