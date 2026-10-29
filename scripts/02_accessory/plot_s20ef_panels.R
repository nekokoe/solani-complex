#!/usr/bin/env Rscript
# plot_s20ef_panels.R — S20 panels e and f as standalone files.
#
#   --panel e : 240 x 120 mm (2:1, as requested)
#   --panel f : 120 x 120 mm (square, same row height so it lines up with e)
#
# Usage:
#   Rscript scripts/plot_s20ef_panels.R --panel e
#   Rscript scripts/plot_s20ef_panels.R --panel f [--no_ylabels]

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s20ef_panel.R"))

option_list <- list(
  make_option("--panel", type = "character", default = "e", help = "e | f"),
  make_option("--sliding", type = "character",
              default = "06.accessory-chrs/BMU02483_gene_content/sliding/sliding_2483_prop.tsv"),
  make_option("--matrix", type = "character",
              default = "06.accessory-chrs/kmer_analysis/kmer_matrix.ACx211.tsv"),
  make_option("--class", type = "character",
              default = "06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv"),
  make_option(c("-o", "--out"), type = "character", default = NULL),
  make_option("--width_mm", type = "numeric", default = NA),
  make_option("--height_mm", type = "numeric", default = NA),
  make_option("--fontsize", type = "numeric", default = 9),
  make_option("--dpi", type = "numeric", default = 600),
  make_option("--no_ylabels", type = "logical", default = FALSE,
              help = "panel f only: drop the chromosome names (for placing it next to e)"),
  make_option("--family", type = "character", default = "Arial")
)
args <- parse_args(OptionParser(option_list = option_list))

sl <- s20e_data(args$sliding)

if (args$panel == "e") {
  W <- if (is.na(args$width_mm)) 240 else args$width_mm
  H <- if (is.na(args$height_mm)) W / 2 else args$height_mm
  out <- if (is.null(args$out)) "06.accessory-chrs/BMU02483_gene_content/fig_s20e_bmu2483_sliding" else args$out
  p <- s20e_panel(sl, family = args$family, fontsize = args$fontsize)
  cat(sprintf("panel e: %d chromosomes x %d windows\n", length(levels(sl$AC)), nrow(sl)))
} else if (args$panel == "f") {
  W <- if (is.na(args$width_mm)) 120 else args$width_mm
  H <- if (is.na(args$height_mm)) 120 else args$height_mm
  out <- if (is.null(args$out)) "06.accessory-chrs/BMU02483_gene_content/fig_s20f_whole_ac_sharing" else args$out
  b <- s20f_data(args$matrix, args$class, sl)
  p <- s20f_panel(b, family = args$family, fontsize = args$fontsize,
                  ylabels = !args$no_ylabels)
  print(table(b$type))
} else stop("--panel must be e or f")

cat(sprintf("panel %s: %g x %g mm\n", args$panel, W, H))
cairo_pdf(paste0(out, ".pdf"), width = W / 25.4, height = H / 25.4,
          family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(out, ".png"), width = W / 25.4, height = H / 25.4,
              units = "in", res = args$dpi, background = "white")
print(p); dev.off()
cat("wrote", paste0(out, ".pdf"), "and", paste0(out, ".png"), "\n")
