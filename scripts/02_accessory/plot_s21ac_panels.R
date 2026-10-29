#!/usr/bin/env Rscript
# plot_s21ac_panels.R — S21 panels a and c as standalone square panels.
#
# Usage:
#   Rscript scripts/plot_s21ac_panels.R --panel a
#   Rscript scripts/plot_s21ac_panels.R --panel c

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s21ac_panel.R"))

option_list <- list(
  make_option("--panel", type = "character", default = "a", help = "a | c"),
  make_option(c("-o", "--out"), type = "character", default = NULL),
  make_option("--size_mm", type = "numeric", default = 120),
  make_option("--fontsize", type = "numeric", default = 10),
  make_option("--dpi", type = "numeric", default = 600),
  make_option("--family", type = "character", default = "Arial")
)
args <- parse_args(OptionParser(option_list = option_list))

df <- s21_data(args$panel)
out <- if (is.null(args$out))
  file.path("06.accessory-chrs/BMU02483_gene_content",
            paste0("fig_s21", args$panel, "_bmu2483_origin")) else args$out
p <- s21_panel(df, args$panel, family = args$family, fontsize = args$fontsize)

cat(sprintf("panel %s: %s groups, values %s\n", args$panel, nrow(df),
            paste(df$pct, collapse = " / ")))
W <- args$size_mm / 25.4
cat(sprintf("panel %.1f x %.1f mm\n", args$size_mm, args$size_mm))
cairo_pdf(paste0(out, ".pdf"), width = W, height = W, family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(out, ".png"), width = W, height = W, units = "in",
              res = args$dpi, background = "white")
print(p); dev.off()
cat("wrote", paste0(out, ".pdf"), "and", paste0(out, ".png"), "\n")
