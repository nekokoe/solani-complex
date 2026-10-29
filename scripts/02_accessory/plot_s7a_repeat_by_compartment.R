#!/usr/bin/env Rscript
# plot_s7a_repeat_by_compartment.R — S7 panel a as a standalone square panel.
#
# Usage:  Rscript scripts/plot_s7a_repeat_by_compartment.R [--size_mm 120] [--fontsize 9]

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s7a_panel.R"))

option_list <- list(
  make_option(c("-t", "--table"), type = "character",
              default = "06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_content.tsv"),
  make_option(c("-o", "--out"), type = "character",
              default = "06.accessory-chrs/repeat_terminal/fig_s7a_repeat_by_compartment"),
  make_option("--size_mm", type = "numeric", default = 120),
  make_option("--fontsize", type = "numeric", default = 10),
  make_option("--dpi", type = "numeric", default = 600),
  make_option("--tests", type = "character",
              default = "06.accessory-chrs/repeat_terminal/suppl_table.repeat_terminal_tests.tsv",
              help = "Mann-Whitney table the significance brackets are read from"),
  make_option("--signif", type = "logical", default = TRUE),
  make_option("--family", type = "character", default = "Arial")
)
args <- parse_args(OptionParser(option_list = option_list))

cc <- read_repeat_content(args$table)
df <- s7a_data(cc)
cat(sprintf("chromosomes: %d core, %d accessory\n",
            sum(cc$Type == "core"), sum(cc$Type == "accessory")))
print(round(tapply(df$repeat_pct, df$compartment, median), 2))

sg <- if (args$signif) s7a_signif(args$tests, levels(df$compartment)) else NULL
if (!is.null(sg)) {
  tt <- read.delim(args$tests, stringsAsFactors = FALSE, check.names = FALSE)
  cat("significance brackets (pooled Mann-Whitney):\n")
  for (i in seq_along(sg$annotations)) {
    cat(sprintf("   %-34s %s\n", paste(sg$comparisons[[i]], collapse = " vs "),
                sg$annotations[i]))
  }
}
p <- s7a_panel(df, family = args$family, fontsize = args$fontsize, signif = sg)
W <- args$size_mm / 25.4
cat(sprintf("panel %.1f x %.1f mm\n", args$size_mm, args$size_mm))
cairo_pdf(paste0(args$out, ".pdf"), width = W, height = W, family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(args$out, ".png"), width = W, height = W, units = "in",
              res = args$dpi, background = "white")
print(p); dev.off()
cat("wrote", paste0(args$out, ".pdf"), "and", paste0(args$out, ".png"), "\n")
