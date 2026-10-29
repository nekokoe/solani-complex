#!/usr/bin/env Rscript
# plot_s16a_transfer_function.R — S16 panel a as a standalone square panel.
#
# Usage:  Rscript scripts/plot_s16a_transfer_function.R [--size_mm 120] [--fontsize 9]

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s16a_panel.R"))

option_list <- list(
  make_option(c("-t", "--table"), type = "character",
              default = "06.accessory-chrs/gene_transfer_211/suppl_table.transfer_function.tsv"),
  make_option(c("-o", "--out"), type = "character",
              default = "06.accessory-chrs/gene_transfer_211/fig_s16a_transfer_function"),
  make_option("--size_mm", type = "numeric", default = 120),
  make_option("--fontsize", type = "numeric", default = 11),
  make_option("--dpi", type = "numeric", default = 600),
  make_option("--fdr", type = "numeric", default = 0.05),
  make_option("--family", type = "character", default = "Arial")
)
args <- parse_args(OptionParser(option_list = option_list))

df <- read_transfer_function(args$table, args$fdr)
cat(sprintf("categories at FDR < %g: %d of %d\n", args$fdr, nrow(df),
            nrow(read.delim(args$table, stringsAsFactors = FALSE))))
print(df[order(df$Odds_ratio), c("Category", "Odds_ratio", "sig")], row.names = FALSE)

p <- s16a_panel(df, family = args$family, fontsize = args$fontsize)
W <- args$size_mm / 25.4
cat(sprintf("panel %.1f x %.1f mm\n", args$size_mm, args$size_mm))
cairo_pdf(paste0(args$out, ".pdf"), width = W, height = W, family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(args$out, ".png"), width = W, height = W, units = "in",
              res = args$dpi, background = "white")
print(p); dev.off()
cat("wrote", paste0(args$out, ".pdf"), "and", paste0(args$out, ".png"), "\n")
