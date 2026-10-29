#!/usr/bin/env Rscript
# plot_s20a_classification.R — standalone k-mer panels built from s20a_panel.R:
#   --mode landscape          S20a  the k-mer classification landscape (78 ACs)
#   --mode clinenv --species falciformis   S20g  sharing with clinical vs non-clinical
#   --mode clinenv --species solani        S20h  isolates of the same species
# Usage:  Rscript scripts/plot_s20a_classification.R [--mode ...] [--size_mm 120]

suppressMessages({library(optparse); library(ggplot2)})
HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(HERE) || HERE == "") HERE <- "scripts"
source(file.path(HERE, "s20a_panel.R"))

option_list <- list(
  make_option("--dual", type="character", default="06.accessory-chrs/kmer_analysis/kmer_type.dual_domain.tsv"),
  make_option("--class", type="character", default="06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv"),
  make_option("--matrix", type="character", default="06.accessory-chrs/kmer_analysis/kmer_matrix.ACx211.tsv"),
  make_option("--meta", type="character", default="05.pan-genome/bmu_ref_groups.211.txt",
              help="MedEnv per isolate -> the AC's clinical/environmental source"),
  make_option(c("-o","--out"), type="character", default="",
              help="default: fig_s20a_classification, or fig_s20<letter>_clinenv_<species> in clinenv mode"),
  make_option("--size_mm", type="numeric", default=120),
  make_option("--fontsize", type="numeric", default=11),
  make_option("--dpi", type="numeric", default=600),
  make_option("--legend", type="character", default="bottom", help="bottom | right | none"),
  make_option("--mode", type="character", default="landscape",
              help="landscape (default) | clinenv: clinical vs non-clinical sharing WITHIN one species"),
  make_option("--species", type="character", default="all",
              help="for --mode clinenv: all (pool solani+falciformis+keratoplastica, panel S20g) or a comma-separated list"),
  make_option("--wide", type="logical", default=FALSE,
              help="2:1 panel for the main figure: plot on the left half, legend on the right half"),
  make_option("--family", type="character", default="Arial")
)
args <- parse_args(OptionParser(option_list=option_list))

# The S20 composite already uses panel letters a-f, and the standalone panels are
# named after the composite panel they were extracted from (S20a = panel a, etc.).
# The clinical/non-clinical panels are NOT panels of that composite, so they take the
# next free letters instead of reusing b/c/d, which would read as "extracted from
# panel b".
CLINENV_LETTER <- c(all = "g")   # only the pooled panel has an assigned letter
if (args$out == "") {
  args$out <- if (args$mode == "clinenv") {
    lt <- CLINENV_LETTER[[args$species]]
    if (is.null(lt))
      stop("no panel letter for '", args$species,
           "'; the assigned letter is for '-\u2013species all'. Pass -o explicitly.")
    "06.accessory-chrs/kmer_analysis/fig_s20g_clinenv_pooled_wide"
  } else "06.accessory-chrs/kmer_analysis/fig_s20a_classification"
}

if (args$mode == "clinenv") {
  spvec <- if (args$species == "all")
    c("solani", "falciformis", "keratoplastica") else strsplit(args$species, ",")[[1]]
  a <- s20a_clinenv_data(args$matrix, args$meta, spvec, args$class)
  nc <- attr(a, "n_clin"); ne <- attr(a, "n_env")
  cat(sprintf("species pooled: %s\n", paste(spvec, collapse = " + ")))
  cat(sprintf("  %d ACs (%s) | POOLED clinical isolates %d, non-clinical %d\n",
              nrow(a), paste(sprintf("%s %d", names(table(a$species)), table(a$species)),
                             collapse = ", "), nc, ne))
  cat("  species-blind: every AC is averaged over the SAME pooled isolate sets; the AC's own donor is excluded\n")
  cat(sprintf("  mean sharing: clinical %.1f%%  non-clinical %.1f%%  |  clinical > non-clinical for %d/%d ACs\n",
              mean(a$share_clin), mean(a$share_env), sum(a$share_clin > a$share_env), nrow(a)))
  cat(sprintf("  Wilcoxon signed-rank p = %.4g\n",
              wilcox.test(a$share_clin, a$share_env, paired = TRUE)$p.value))
  # S20g is drawn in the 2:1 form (plot left half, legend right half), like S20a-wide
  p <- s20a_panel_wide(a, family=args$family, fontsize=args$fontsize,
                       title="AC sharing with clinical vs non-clinical isolates\n(78 ACs of three species, species boundaries ignored)",
                       xcol="share_clin", ycol="share_env", xlab=CLIN_LAB, ylab=ENV_LAB,
                       ellipses=FALSE, diagonal=TRUE,
                       diagonal_labels=c("shared more with\nclinical isolates",
                                         "shared more with\nnon-clinical isolates"))
  W <- 2 * args$size_mm/25.4; H <- args$size_mm/25.4
  cat(sprintf("wide panel %.0f x %.0f mm (plot left half, legend right half)\n", args$size_mm*2, args$size_mm))
} else {
a <- s20a_data(args$dual, args$class, args$matrix, args$meta)
cat(sprintf("ACs: %d  (Block %d / Mosaic %d; Clinical %d / Environmental %d)\n",
            nrow(a), sum(a$cls == "Block"), sum(a$cls == "Mosaic"),
            sum(a$source == "Clinical"), sum(a$source == "Environmental")))
cat(sprintf("carriers %d..%d | same-species sharing %.1f-%.1f%% | other-species %.1f-%.1f%%\n",
            min(a$n_carriers), max(a$n_carriers), min(a$share_in), max(a$share_in),
            min(a$share_out), max(a$share_out)))

if (isTRUE(args$wide)) {
  p <- s20a_panel_wide(a, family=args$family, fontsize=args$fontsize)
  W <- 2 * args$size_mm / 25.4; H <- args$size_mm / 25.4
  cat(sprintf("wide panel %.0f x %.0f mm (plot left half, legend right half)\n",
              args$size_mm * 2, args$size_mm))
} else {
  lp <- if (args$legend == "none") "none" else args$legend
  p <- s20a_panel(a, family=args$family, fontsize=args$fontsize, legend_position=lp)
  W <- args$size_mm / 25.4; H <- W
}
}
cairo_pdf(paste0(args$out, ".pdf"), width=W, height=H, family=args$family); print(p); dev.off()
ragg::agg_png(paste0(args$out, ".png"), width=W, height=H, units="in", res=args$dpi, background="white"); print(p); dev.off()
cat("wrote", paste0(args$out, ".pdf"), "and", paste0(args$out, ".png"), "\n")
