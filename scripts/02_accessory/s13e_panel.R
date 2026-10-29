#!/usr/bin/env Rscript
# s13e_panel.R — S13 panel e as a small square panel: the TE-pair heatmap
# (upstream flanking family x downstream flanking family, 5 kb), for pasting into
# the figure.
#
# Same conventions as s13a_panel.R: square canvas, legible font sizes, NimbusSans
# (Arial resolves to it on this host), and the S13 palette.
#
# The fill is the SHARE OF ALL TE-BRACKETED GENES that use that upstream x downstream
# pair, not the raw count. The denominator is therefore the number of TE-bracketed
# genes (both flanks a repeat family other than Unknown / Simple_repeat /
# Low_complexity), read from the companion total file written by
# plot_s5_s6_boxplots_s15_heatmap.py. Note sum(matrix) is NOT that denominator: the
# 12 x 12 grid only keeps pairs whose two families are both inside the top 12, so it
# covers 501 of the 581 bracketed genes (86%). Both numbers go in the caption.
#
# The log colour scale needs a floor because 44 of the 144 cells are exactly zero;
# they are floored just below the smallest non-empty share.
#
# The x tick labels are rotated 45 degrees. An earlier note here claimed 45 degrees
# was impossible because a label needs ~44 pt of HORIZONTAL footprint; that is the
# wrong criterion. Parallel rotated labels are separated along the axis PERPENDICULAR
# to their baseline, so the test is
#     slot * sin(theta)  >=  ink thickness          (ink ~ fontsize * 0.95)
# At 120 mm the tile pitch is 15.32 pt, so 15.32 * sin 45 = 10.83 pt against 8.55 pt
# of ink: clearance +2.28 pt, and 45 degrees also uses ~3 mm less vertical space than
# 90 degrees. The axis-aligned bounding boxes of parallel 45-degree labels always
# intersect, so a bbox-overlap count must NOT be used to reject the angle.

suppressMessages({library(ggplot2); library(reshape2)})

S13_RED <- "#C44E52"

read_te_pair_matrix <- function(path) {
  m <- as.matrix(read.delim(path, row.names = 1, check.names = FALSE))
  if (nrow(m) != ncol(m)) stop("expected a square matrix")
  m
}

# the number of TE-bracketed genes, i.e. the denominator of every cell
read_te_pair_total <- function(path =
    "06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.total.tsv") {
  if (!file.exists(path)) stop("missing denominator file: ", path,
                               " (re-run scripts/plot_s5_s6_boxplots_s15_heatmap.py)")
  t <- read.delim(path, check.names = FALSE)
  as.numeric(t$te_bracketed_genes[1])
}

s13e_panel <- function(m, denom, family = "Arial",
                       title = "Bracketing TE pairs\n(% of TE-bracketed genes, 5 kb)",
                       fontsize = 10, xangle = 45) {
  font_axis <- fontsize
  font_title <- fontsize + 1
  font_legend <- fontsize - 1
  short <- function(s) gsub("DNA/|RC/|LINE/|LTR/", "", s)
  lv <- rownames(m)
  pct <- 100 * m / denom
  floor_pct <- min(pct[pct > 0]) / 2          # 0 cells need a positive value on a log scale
  df <- reshape2::melt(pmax(pct, floor_pct), varnames = c("up", "down"),
                       value.name = "pct")
  df$up <- factor(df$up, levels = lv)
  df$down <- factor(df$down, levels = lv)
  # keep the true count alongside, so cells can be labelled with either
  df$n <- as.vector(m)
  df$pct_true <- as.vector(pct)
  bmax <- max(pct)
  brks <- c(0.1, 0.5, 1, 2, 5)
  brks <- brks[brks <= bmax]

  ggplot(df, aes(x = down, y = up, fill = pct)) +
    geom_tile(colour = "grey60", linewidth = 0.25) +
    scale_fill_gradient(low = "white", high = S13_RED, trans = "log10",
                        # breaks outside the data range are silently dropped by
                        # ggplot, which leaves the top of the bar unlabelled
                        breaks = brks,
                        labels = formatC(brks, format = "fg"),
                        name = "% of TE-\nbracketed genes") +
    scale_x_discrete(labels = function(x) short(x), expand = c(0, 0)) +
    scale_y_discrete(labels = function(x) short(x), expand = c(0, 0)) +
    coord_fixed() +
    labs(x = "Downstream flanking TE family", y = "Upstream flanking TE family",
         title = title) +
    theme_minimal(base_size = 9, base_family = family) +
    theme(
      panel.grid = element_blank(),
      # vjust matters at 45 degrees: vjust = 0 hangs the labels downward past the
      # space ggplot reserved, pushing them off the bottom of the page (and the axis
      # title ends up above them). hjust = 1, vjust = 1 is the correct 45-degree form.
      axis.text.x = element_text(angle = xangle, hjust = 1,
                                 vjust = if (xangle == 90) 0.5 else 1,
                                 size = font_axis, colour = "black"),
      axis.text.y = element_text(size = font_axis, colour = "black"),
      axis.title.x = element_text(size = font_axis + 0.5),
      axis.title.y = element_text(size = font_axis + 0.5),
      plot.title = element_text(size = font_title, face = "bold", hjust = 0),
      legend.position = "right",
      legend.title = element_text(size = font_legend),
      legend.text = element_text(size = font_legend - 0.5),
      # do NOT pin legend.key.height: a fixed value made the legend taller than
      # the panel, and the bar + its title were clipped off the top of the page
      legend.key.width = unit(4, "mm"),
      plot.margin = margin(5, 4, 4, 5)
    )
}
