#!/usr/bin/env Rscript
# plot_s13a_transfer_by_species.R — S13 panel a, standalone square panel.
#
# The published panel a showed only the cohort totals (A-A 8,988 / C-A 2,907)
# under a title reading "(211 genomes, 100% identity)", which reads as if the 211
# genomes were the unit of counting. They are not: the counted genes are the 13
# chromosome-level (TGS) assemblies' genes and the 211 genomes only supply the
# mmseqs clusters (see s13a_panel.R and stats_gene_transfer_211.py). This redraw
# states that split, breaks the counts out by species, and labels each C-A bar
# with its share of that group's accessory-chromosome genes.
#
# Usage:  Rscript scripts/plot_s13a_transfer_by_species.R [--size_mm 90]

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s13a_panel.R"))

option_list <- list(
  make_option(c("-t", "--table"), type = "character",
              default = "06.accessory-chrs/gene_transfer_211_99/suppl_table.gene_transfer_211.tsv"),
  make_option(c("-o", "--out"), type = "character",
              default = "06.accessory-chrs/gene_transfer_211_99/fig_s13a_transfer_by_species"),
  make_option("--size_mm", type = "numeric", default = 120),
  make_option("--dpi", type = "numeric", default = 600),
  make_option("--with_all", type = "logical", default = TRUE),
  make_option("--fontsize", type = "numeric", default = 11),
  make_option("--subtitle", type = "logical", default = FALSE),
  make_option("--pct", type = "logical", default = FALSE),
  make_option("--family", type = "character", default = "Arial",
              help = "Arial resolves to Nimbus Sans on this host, as in the sibling panels")
)
args <- parse_args(OptionParser(option_list = option_list))

rows <- s13a_table(args$table, args$with_all)
p <- s13a_panel(rows, family = args$family, fontsize = args$fontsize,
                show_subtitle = args$subtitle, show_pct = args$pct)

cat(sprintf("panel %.1f x %.1f mm\n", args$size_mm, args$size_mm))
print(as.data.frame.matrix(tapply(rows$n, list(rows$grp, rows$type), identity)))

W <- args$size_mm / 25.4
cairo_pdf(paste0(args$out, ".pdf"), width = W, height = W, family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(args$out, ".png"), width = W, height = W, units = "in",
              res = args$dpi, background = "white")
print(p); dev.off()
cat("wrote", paste0(args$out, ".pdf"), "and", paste0(args$out, ".png"), "\n")
