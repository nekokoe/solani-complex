#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# =============================================================================
# Candidate spok figures (C-2) -- statistical + structural schematic candidates
# for the supplementary material.  Five figures, all from the C-1/C-2 products:
#
#   fig_spok_A_proportion   observed share on AC vs the length-implied null
#   fig_spok_B_by_isolate   per-isolate AC/core split for both gene families
#   fig_spok_C_forest       density ratio with 95% CI, across identity thresholds
#   fig_spok_D_genemap      gene-map schematic: where spok and Captain actually sit
#   fig_spok_E_model        conceptual model, FOSC (paper) vs FSSC (this study)
#
# Conventions inherited from this project's panel builders:
#   * font family = "Arial"            (resolves to NimbusSans; do NOT use "sans")
#   * grid::unit() must be qualified   (an unqualified unit() yields a blank PDF)
#   * ggforce has no geom_rrect -- geom_shape(radius = ...) instead
#   * one colour scale per ggplot      (ggnewscale for a second)
#   * every render is measured afterwards (>= 1.2 mm of ink-free margin per side)
#
# Usage: Rscript scripts/spok_figures.R [outdir]
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ggforce)
  library(scales)
})

# ---- resolve the repo root from this script's own location (never from cwd) --
.argv <- commandArgs(trailingOnly = FALSE)
.f <- sub("^--file=", "", .argv[grep("^--file=", .argv)])
BASE <- if (length(.f)) normalizePath(dirname(dirname(.f[[1]]))) else getwd()
DATA <- file.path(BASE, "06.accessory-chrs", "starship_ac")
OUT <- if (length(commandArgs(trailingOnly = TRUE)))
  commandArgs(trailingOnly = TRUE)[[1]] else file.path(DATA, "figs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
cat("BASE:", BASE, "\nOUT :", OUT, "\n")

FONT <- "Arial"
theme_spok <- function(base = 12) {
  theme_bw(base_size = base, base_family = FONT) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", linewidth = 0.5),
          plot.title = element_text(size = base + 1, face = "bold"),
          plot.subtitle = element_text(size = base - 2, colour = "grey25"),
          axis.title = element_text(size = base),
          axis.text = element_text(size = base - 1, colour = "black"),
          legend.title = element_text(size = base - 1),
          legend.text = element_text(size = base - 2),
          # reserve space so text can never sit flush against the canvas edge
          plot.margin = margin(5.5, 9, 5.5, 5.5))
}
# colour-blind-safe
C_AC   <- "#D1495B"   # accessory
C_CORE <- "#4A6FA5"   # core
C_CAP  <- "#8A7B4F"   # Captain / DUF3435
C_SPOK <- "#2E8B7A"   # spok

# ------------------------------------------------------------------ read data
cap <- read.delim(file.path(DATA, "suppl_table.starship_gene_positions.tsv"),
                  stringsAsFactors = FALSE)
spok <- read.delim(file.path(DATA, "spok_unique_genes.tsv"), stringsAsFactors = FALSE)
lens <- read.delim(file.path(DATA, "suppl_table.starship_lengths.tsv"),
                   stringsAsFactors = FALSE)
sens_spok <- read.delim(file.path(DATA, "spok_threshold_sensitivity.tsv"),
                        stringsAsFactors = FALSE)

AC_MB   <- lens$length_Mb[lens$contig_class == "accessory"]
CORE_MB <- lens$length_Mb[lens$contig_class == "core"]
P_AC    <- AC_MB / (AC_MB + CORE_MB)
cat(sprintf("AC share of length = %.4f  (AC %.1f Mb / core %.1f Mb)\n", P_AC, AC_MB, CORE_MB))

summ <- rbind(
  data.frame(family = "Starship Captain\n(DUF3435 = PF11917)", n_ac = sum(cap$contig_class == "accessory"),
             n_core = sum(cap$contig_class == "core")),
  data.frame(family = "spok homologs\n(Spok3 homology)", n_ac = sum(spok$contig_class == "accessory"),
             n_core = sum(spok$contig_class == "core"))
)
summ$n_tot <- summ$n_ac + summ$n_core
summ$p_obs <- summ$n_ac / summ$n_tot
summ$ratio <- summ$p_obs / P_AC
summ$pval  <- mapply(function(k, n) binom.test(k, n, P_AC, alternative = "greater")$p.value,
                     summ$n_ac, summ$n_tot)
cat("\n", paste(capture.output(print(summ, row.names = FALSE)), collapse = "\n"), "\n")

# =============================================================================
# A -- observed share on AC vs the length-implied null
# =============================================================================
pa <- ggplot(summ, aes(x = family, y = p_obs * 100)) +
  geom_col(aes(fill = family), width = 0.55, colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = P_AC * 100, linetype = "dashed", linewidth = 0.5, colour = "grey20") +
  geom_text(aes(label = sprintf("%d/%d\n%.1f%%", n_ac, n_tot, p_obs * 100)),
            vjust = -0.25, size = 3.4, family = FONT, lineheight = 1.05) +
  # significance bracket
  annotate("segment", x = 1, xend = 2, y = max(summ$p_obs * 100) * 1.32,
           yend = max(summ$p_obs * 100) * 1.32, linewidth = 0.4) +
  annotate("text", x = 1.5, y = max(summ$p_obs * 100) * 1.35, size = 3.0, family = FONT,
           lineheight = 1.1,
           label = sprintf("binomial p\nCaptain %.3f  |  spok %.4f",
                           summ$pval[1], summ$pval[2])) +
  scale_fill_manual(values = c(C_CAP, C_SPOK), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.30)),
                     labels = function(x) paste0(x, "%")) +
  labs(title = "Which gene family sits on\naccessory chromosomes?",
       subtitle = sprintf("bars = share on AC\nnull = %.1f%% from AC genome length", P_AC * 100),
       x = NULL, y = "genes on accessory chromosomes") +
  theme_spok()

# =============================================================================
# B -- per-isolate split
# =============================================================================
iso_order <- c("BMU01933","BMU02028","BMU02474","BMU03267","BMU06673","BMU12683",
               "BMU12684","BMU12685","BMU12686","Fu3.1","Fu6.1","SB1","SX13")
tab <- rbind(
  data.frame(isolate = cap$isolate,  family = "Captain (DUF3435)", cls = cap$contig_class),
  data.frame(isolate = spok$isolate, family = "spok homolog",      cls = spok$contig_class)
)
tab <- tab[tab$cls %in% c("accessory", "core"), ]
tab$family <- factor(tab$family, levels = c("Captain (DUF3435)", "spok homolog"))
tab$cls <- factor(tab$cls, levels = c("core", "accessory"))
tab$isolate <- factor(tab$isolate, levels = rev(iso_order))

mk_pb <- function(base = 11) ggplot(tab, aes(x = isolate, fill = cls)) +
  geom_bar(position = position_dodge(width = 0.78), width = 0.7,
           colour = "black", linewidth = 0.3) +
  facet_wrap(~family, nrow = 1) +
  coord_flip() +
  scale_fill_manual(values = c(core = C_CORE, accessory = C_AC), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = pretty_breaks()) +
  labs(title = "Per-isolate gene counts by chromosome class",
       subtitle = "13 chromosome-level assemblies; accessory = AC contigs",
       x = NULL, y = "number of genes") +
  theme_spok(base = base) +
  theme(legend.position = "top",
        strip.text = element_text(size = base * 0.9, face = "bold"))
pb <- mk_pb()

# =============================================================================
# C -- density ratio with 95% CI across identity thresholds
# =============================================================================
wilson <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- qnorm(1 - (1 - conf) / 2)
  p <- k / n
  d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(ctr - half, ctr + half)
}
fr <- do.call(rbind, lapply(seq_len(nrow(sens_spok)), function(i) {
  r <- sens_spok[i, ]
  n <- r$on_accessory + r$on_core
  ci <- wilson(r$on_accessory, n)
  data.frame(family = "spok homolog", threshold = r$pident_threshold,
             ratio = if (is.na(r$density_ratio_AC_over_core)) NA else r$density_ratio_AC_over_core,
             lo = ci[1] / P_AC, hi = ci[2] / P_AC, k = r$on_accessory, n = n)
}))
fr <- rbind(fr, data.frame(family = "Captain (DUF3435)", threshold = NA,
                           ratio = summ$ratio[1], lo = wilson(summ$n_ac[1], summ$n_tot[1])[1] / P_AC,
                           hi = wilson(summ$n_ac[1], summ$n_tot[1])[2] / P_AC,
                           k = summ$n_ac[1], n = summ$n_tot[1]))
fr$lab <- ifelse(is.na(fr$threshold), "Captain (all)", sprintf("spok, id >= %g%%", fr$threshold))
# a threshold with no AC hit cannot be placed on a log axis at all; drop it and
# say so, rather than letting ggplot silently remove the row.
drop_lab <- fr$lab[!is.na(fr$k) & fr$k == 0]
fr <- fr[is.na(fr$k) | fr$k > 0, ]
fr$lab <- factor(fr$lab, levels = fr$lab[order(fr$family, fr$threshold)])
fr$family <- factor(fr$family, levels = c("Captain (DUF3435)", "spok homolog"))

mk_pc <- function(base = 11) ggplot(fr, aes(x = ratio, y = lab, colour = family)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5, colour = "grey30") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.16, linewidth = 0.5) +
  geom_point(aes(size = n), shape = 18) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3, 5), labels = c("0.5x", "1x", "2x", "3x", "5x")) +
  scale_colour_manual(values = c(C_CAP, C_SPOK), guide = "none") +
  scale_size_continuous(name = "genes", range = c(2.2, 4.6), breaks = c(20, 44, 68)) +
  labs(title = "Enrichment on accessory chromosomes, with 95% CI",
       subtitle = paste0("observed / expected share from AC length",
                         if (length(drop_lab))
                           sprintf("\n0 AC hits, not shown: id >= %s%%",
                                   paste(unique(sens_spok$pident_threshold[
                                     sens_spok$on_accessory == 0]), collapse = "/")) else ""),
       x = "observed / expected density on AC", y = NULL) +
  theme_spok(base = base) +
  theme(legend.position = "right")
pc <- mk_pc()

# =============================================================================
# D -- gene-map schematic: where spok and Captain actually sit
# =============================================================================
read_coords <- function(iso) {
  gff <- file.path(BASE, "02.annotation", "per_isolates", iso,
                   "funannotate", "predict_results", "fusarium.gff3")
  if (!file.exists(gff)) return(NULL)
  ln <- grep("\tmRNA\t", readLines(gff, warn = FALSE), value = TRUE, fixed = TRUE)
  if (!length(ln)) return(NULL)
  f <- strsplit(ln, "\t", fixed = TRUE)
  tid <- vapply(f, function(x) {
    a <- strsplit(x[[9]], ";", fixed = TRUE)[[1]]
    a <- a[grepl("^ID=", a)][1]
    if (is.na(a)) NA_character_ else sub("^ID=", "", a)
  }, character(1))
  data.frame(isolate = iso, contig = vapply(f, `[[`, character(1), 1),
             start = as.numeric(vapply(f, `[[`, character(1), 4)),
             end   = as.numeric(vapply(f, `[[`, character(1), 5)),
             gene  = tid, stringsAsFactors = FALSE)
}
co <- do.call(rbind, lapply(unique(c(cap$isolate, spok$isolate)), read_coords))
co$len <- NA_real_
cl <- do.call(rbind, lapply(unique(c(cap$isolate, spok$isolate)), function(iso) {
  p <- file.path(BASE, "02.annotation", "per_isolates", iso, "repeatmasker",
                 sprintf("%s.assembly.fasta.masked.length", iso))
  if (!file.exists(p)) return(NULL)
  d <- read.table(p, stringsAsFactors = FALSE)
  data.frame(isolate = iso, contig = d[[1]], clen = d[[2]], stringsAsFactors = FALSE)
}))

mark <- rbind(
  data.frame(co[co$gene %in% spok$protein, c("isolate","contig","start","end","gene")],
             kind = "spok"),
  data.frame(co[co$gene %in% cap$tid, c("isolate","contig","start","end","gene")],
             kind = "Captain"))
mark <- unique(mark)
# the AC contigs that carry at least one spok gene
ac_spok <- unique(spok[spok$contig_class == "accessory", c("isolate", "contig")])
ac_spok$clen <- cl$clen[match(paste(ac_spok$isolate, ac_spok$contig),
                              paste(cl$isolate, cl$contig))]
ac_spok <- ac_spok[!is.na(ac_spok$clen), ]
ac_spok <- ac_spok[order(ac_spok$isolate, ac_spok$contig), ]
ac_spok$row <- seq_len(nrow(ac_spok))

mk <- merge(mark, ac_spok[, c("isolate","contig","row")], by = c("isolate","contig"))
mk$pos <- (mk$start + mk$end) / 2

# The question is how many *spok* genes share a contig with a Captain gene --
# NOT how many Captain genes sit on those contigs (an easy and wrong substitution).
co_loc <- unique(mk$contig[mk$kind == "Captain"])
n_spok_ac <- length(unique(paste(spok$isolate[spok$contig_class == "accessory"],
                                 spok$contig[spok$contig_class == "accessory"])))
pd <- ggplot(ac_spok, aes(y = row)) +
  geom_segment(aes(x = 0, xend = clen / 1000, y = row, yend = row),
               linewidth = 1.6, colour = "grey82", lineend = "round") +
  geom_segment(data = mk, aes(x = pos / 1000, xend = pos / 1000, y = row - 0.30,
                              yend = row + 0.30, colour = kind),
               linewidth = 1.1, inherit.aes = FALSE) +
  scale_colour_manual(values = c(spok = C_SPOK, Captain = C_CAP), name = NULL) +
  scale_y_continuous(breaks = ac_spok$row,
                     labels = sprintf("%s  %s  (%.1f Mb)", ac_spok$isolate,
                                      ac_spok$contig, ac_spok$clen / 1e6),
                     expand = expansion(add = 0.8)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  labs(title = "Structural context of the accessory-chromosome spok genes",
       subtitle = paste0("every AC contig carrying a spok homolog; ticks mark gene positions",
                         "\n", length(co_loc), " of ", n_spok_ac,
                         " AC contigs with spok also carry a Captain gene"),
       x = "position on the contig (kb)", y = NULL) +
  theme_spok(base = 10) +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 7.6, family = FONT),
        panel.grid.major.y = element_blank())

# =============================================================================
# E -- conceptual model: FOSC (paper) vs FSSC (this study)
# =============================================================================
box <- function(xmin, ymin, xmax, ymax, fill, col = "black", lw = 0.5, r = 0.02)
  geom_shape(data = data.frame(x = c(xmin, xmax, xmax, xmin), y = c(ymin, ymin, ymax, ymax)),
             aes(x = x, y = y), fill = fill, colour = col, linewidth = lw,
             radius = grid::unit(r, "npc"), inherit.aes = FALSE)
txt <- function(x, y, label, size = 3.1, ...)
  annotate("text", x = x, y = y, label = label, size = size, family = FONT, ...)
seg <- function(x, xend, y, yend, ...) annotate("segment", x = x, xend = xend, y = y, yend = yend, ...)

pe <- ggplot() +
  xlim(0, 10) + ylim(0, 10) +
  # ---- left: the FOSC model from the paper
  box(0.25, 5.4, 4.35, 9.4, "#F3F4F6") +
  txt(2.3, 9.05, "F. oxysporum (Puebla Planas et al.)", size = 3.4, fontface = "bold") +
  seg(0.6, 4.0, 6.15, 6.15, linewidth = 1.8, colour = "grey80", lineend = "round") +
  box(1.55, 5.75, 3.05, 6.55, "#E7D9B0") +
  txt(2.3, 6.15, "Starship", size = 2.9) +
  seg(1.55, 1.55, 5.75, 6.55, linewidth = 1.0, colour = C_CAP) +
  txt(1.30, 6.15, "Captain", size = 2.3, angle = 90, colour = C_CAP) +
  seg(2.6, 2.6, 5.72, 6.58, linewidth = 1.5, colour = C_SPOK) +
  txt(2.80, 7.05, "spok", size = 2.9, colour = C_SPOK, fontface = "bold") +
  txt(2.3, 5.05, "dispensable chromosome (DC)", size = 2.9, colour = "grey25") +
  txt(2.3, 4.35, "spok kills cells that lose\nthe DC  ->  DC maintained", size = 2.9, lineheight = 1.05) +
  # ---- right: this study
  box(5.65, 5.4, 9.75, 9.4, "#F3F4F6") +
  txt(7.7, 9.05, "F. solani complex (this study)", size = 3.4, fontface = "bold") +
  seg(6.0, 9.4, 6.15, 6.15, linewidth = 1.8, colour = "grey80", lineend = "round") +
  seg(7.2, 7.2, 6.0, 6.3, linewidth = 2.0, colour = C_SPOK) +
  txt(7.35, 6.15, "spok", size = 3.0, colour = C_SPOK, fontface = "bold", hjust = 0) +
  txt(6.0, 5.05, "accessory chromosome (AC)", size = 2.9, colour = "grey25", hjust = 0) +
  seg(6.0, 9.4, 3.75, 3.75, linewidth = 1.8, colour = "grey80", lineend = "round") +
  seg(6.6, 6.6, 3.6, 3.9, linewidth = 1.6, colour = C_CAP) +
  seg(8.0, 8.0, 3.6, 3.9, linewidth = 1.6, colour = C_CAP) +
  txt(6.0, 3.22, "core chromosome", size = 2.9, colour = "grey25", hjust = 0) +
  txt(6.6, 4.28, "Captain", size = 2.4, colour = C_CAP) +
  txt(8.0, 4.28, "Captain", size = 2.4, colour = C_CAP) +
  txt(7.7, 2.10,
      sprintf("spok IS enriched on AC (%.2fx, p = %.4f)\nbut Captain is NOT (%.2fx, p = %.3f)\n=> the two are decoupled",
              summ$ratio[2], summ$pval[2], summ$ratio[1], summ$pval[1]),
      size = 3.0, lineheight = 1.15) +
  # ---- shared legend
  seg(3.6, 4.05, 0.92, 0.92, linewidth = 1.6, colour = C_SPOK) +
  txt(4.15, 0.92, "spok (spore-killer driver)", size = 2.9, hjust = 0) +
  seg(3.6, 4.05, 0.42, 0.42, linewidth = 1.6, colour = C_CAP) +
  txt(4.15, 0.42, "Captain (Starship marker)", size = 2.9, hjust = 0) +
  labs(title = "Two different questions: is the driver there,\nand is its marker there?",
       subtitle = "left = the mechanism demonstrated in F. oxysporum;  right = what the FSSC data show") +
  theme_spok(base = 11) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank())

# =============================================================================
# render
# =============================================================================
specs <- list(
  list(name = "fig_spok_A_proportion", plot = pa, w = 132, h = 108),
  list(name = "fig_spok_B_by_isolate", plot = pb, w = 130, h = 120),
  list(name = "fig_spok_C_forest",     plot = pc, w = 145, h = 118),
  list(name = "fig_spok_D_genemap",    plot = pd, w = 185, h = 125),
  list(name = "fig_spok_E_model",      plot = pe, w = 175, h = 105)
)
for (s in specs) {
  pdf_f <- file.path(OUT, paste0(s$name, ".pdf"))
  png_f <- file.path(OUT, paste0(s$name, ".png"))
  ggsave(pdf_f, s$plot, width = s$w, height = s$h, units = "mm", device = cairo_pdf)
  ggsave(png_f, s$plot, width = s$w, height = s$h, units = "mm", dpi = 300)
  cat(sprintf("wrote %-26s %4.0f x %4.0f mm  (%.1f kB pdf)\n",
              s$name, s$w, s$h, file.info(pdf_f)$size / 1024))
}
cat("\ndone\n")

# =============================================================================
# composites for the supplementary material
#   S29 = a proportion | b per-isolate | c forest          (statistics)
#   S30 = a gene map   | b conceptual model                (structure)
#
# Composing is NOT free: a panel designed at 130-145 mm and then squeezed into
# an ~89 mm cell keeps its fixed point sizes, so its text collides (that is
# exactly what happened on the first attempt -- 25 overlaps in S29). The two
# shrunk panels are therefore re-rendered with a font scale reduced in
# proportion (7/11), and the canvas height is derived from the natural aspect
# ratios so patchwork never distorts or compresses a panel.
# =============================================================================
ar <- function(p, w, h)
  p + theme(aspect.ratio = h / w,
            plot.tag = element_text(family = FONT, face = "bold", size = 11))

W29 <- 180
S29_R2 <- (W29 - 2) / 2                       # two panels per row, ~1 mm gutter
h29_r1 <- W29 / (132 / 108)                   # A at full width
h29_r2 <- S29_R2 / (130 / 120)                # B is the taller of the pair

s29 <- ar(pa, 132, 108) / (ar(mk_pb(7), 130, 120) | ar(mk_pc(7), 145, 118)) +
  plot_layout(heights = c(h29_r1, h29_r2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = 11))

W30 <- 180
s30 <- ar(pd, 185, 125) / ar(pe, 175, 105) +
  plot_layout(heights = c(W30 / (185 / 125), W30 / (175 / 105))) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(family = FONT, face = "bold", size = 11))

for (s in list(list(n = "Supplementary-Fig-S29-spok-statistics", p = s29,
                    w = W29, h = h29_r1 + h29_r2 + 6),
               list(n = "Supplementary-Fig-S30-spok-structure",  p = s30,
                    w = W30, h = W30 / (185 / 125) + W30 / (175 / 105) + 6))) {
  for (dev in c("pdf", "png")) {
    ggsave(file.path(OUT, paste0(s$n, ".", dev)), s$p, width = s$w, height = s$h,
           units = "mm", dpi = 300,
           device = if (dev == "pdf") cairo_pdf else "png")
  }
  cat(sprintf("wrote %-46s %3.0f x %3.0f mm  (%.1f kB)\n", s$n, s$w, s$h,
              file.info(file.path(OUT, paste0(s$n, ".pdf")))$size / 1024))
}
cat("\ncomposites done\n")
