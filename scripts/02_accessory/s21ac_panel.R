#!/usr/bin/env Rscript
# s21ac_panel.R — S21 panels a and c as standalone square panels.
#
#   a = share of HOGs that sit in the clinical accessory-chromosome gene set,
#       for BMU02483-specific HOGs versus two controls
#   c = share of HOGs that sit on accessory chromosomes, same three groups
#
# Both are three-bar panels sharing the same group labels and colours as the S21
# assembly, so the two read as a pair.
#
# The numbers are the ones hard-coded in scripts/ac_plot_bmu2483_origin.R (the S21
# script loads no table for these two panels); they are repeated here rather than
# recomputed, and the denominators travel with each bar.

suppressMessages({library(ggplot2)})

S21_RED <- "#C44E52"; S21_BLUE <- "#3B7DD8"; S21_GREY <- "#BBBBBB"

GROUPS <- c("BMU02483-specific\n(170 HOGs)", "Nearby env strains\n(20,133 HOGs)",
            "Genome background\n(29,080 HOGs)")

S21_PANELS <- list(
  a = list(pct = c(12.4, 4.2, 3.2), n_in = c(21, 838, 939),
           ylab = "HOGs in clinical-AC gene set (%)",
           title = "BMU02483-specific HOGs enriched\nfor clinical-AC content"),
  c = list(pct = c(40.0, 16.0, 13.8), n_in = c(68, 3229, 4002),
           ylab = "HOGs on accessory chromosomes (%)",
           title = "Specific HOGs are enriched\non accessory chromosomes")
)

s21_data <- function(which) {
  p <- S21_PANELS[[which]]
  if (is.null(p)) stop("which must be 'a' or 'c'")
  data.frame(group = factor(GROUPS, levels = GROUPS), pct = p$pct,
             n_in = p$n_in, n_all = c(170, 20133, 29080),
             lab = sprintf("%.1f%%\n(%s/%s)", p$pct,
                           format(p$n_in, big.mark = ",", trim = TRUE),
                           format(c(170, 20133, 29080), big.mark = ",", trim = TRUE)))
}

s21_panel <- function(df, which, family = "Arial", fontsize = 10) {
  p <- S21_PANELS[[which]]
  ggplot(df, aes(x = group, y = pct, fill = group)) +
    geom_col(width = 0.62) +
    scale_fill_manual(values = c(S21_RED, S21_BLUE, S21_GREY)) +
    geom_text(aes(label = lab), vjust = -0.28, lineheight = 0.95,
              size = (fontsize - 0.8) / 2.845) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(x = NULL, y = p$ylab, title = p$title) +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(size = fontsize, colour = "black", lineheight = 0.85),
      axis.text.y = element_text(size = fontsize, colour = "black"),
      axis.title.y = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 4, 4)
    )
}
