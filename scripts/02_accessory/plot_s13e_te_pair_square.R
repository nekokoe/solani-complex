#!/usr/bin/env Rscript
# plot_s13e_te_pair_square.R — S13 panel e as a standalone small square panel.
#
# Usage:  Rscript scripts/plot_s13e_te_pair_square.R [--size_mm 90]

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s13e_panel.R"))

option_list <- list(
  make_option(c("-m", "--matrix"), type = "character",
              default = "06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.tsv"),
  make_option(c("-t", "--total"), type = "character",
              default = "06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.total.tsv"),
  make_option(c("-o", "--out"), type = "character",
              default = "06.accessory-chrs/gene_transfer_211/fig_s13e_te_pair_heatmap"),
  make_option("--size_mm", type = "numeric", default = 120),
  make_option("--dpi", type = "numeric", default = 600),
  make_option("--fontsize", type = "numeric", default = 10),
  make_option("--xangle", type = "numeric", default = 45,
              help = "90 fits at 90 mm; 45 needs a ~187 mm heatmap"),
  make_option("--family", type = "character", default = "Arial")
)
args <- parse_args(OptionParser(option_list = option_list))

m <- read_te_pair_matrix(args$matrix)
denom <- read_te_pair_total(args$total)
cat(sprintf("matrix %d x %d, counts %d..%d, %d zeros; denominator %d TE-bracketed genes (grid covers %.1f%%)\n",
            nrow(m), ncol(m), min(m), max(m), sum(m == 0), denom, 100 * sum(m) / denom))
p <- s13e_panel(m, denom, family = args$family, fontsize = args$fontsize,
                xangle = args$xangle)

W <- args$size_mm / 25.4
cat(sprintf("panel %.1f x %.1f mm\n", args$size_mm, args$size_mm))
cairo_pdf(paste0(args$out, ".pdf"), width = W, height = W, family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(args$out, ".png"), width = W, height = W, units = "in",
              res = args$dpi, background = "white")
print(p); dev.off()
cat("wrote", paste0(args$out, ".pdf"), "and", paste0(args$out, ".png"), "\n")
