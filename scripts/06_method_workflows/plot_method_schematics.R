#!/usr/bin/env Rscript
# plot_method_schematics.R — the four method schematics (S2, S8, S9, S11), redrawn.
#
# WHY THIS SCRIPT EXISTS
#   The previous versions were produced with the BASE `pdf()` device and
#   `fontfamily = "sans"` (-> Helvetica), which has no glyphs for ³ ⁴ ⁵ β α ρ, so every
#   such character was rendered as "." (e.g. "10. / 10. conidia per larva",
#   "IL–1. / TNF–."). They also carried pre-99%-identity and pre-21-strain numbers, and
#   their layout left large empty gaps with dashed frames that the arrows crossed.
#
#   Here: cairo_pdf + Arial (-> NimbusSans, full Greek/superscript coverage), one shared
#   visual language across the four figures, and the current numbers:
#     99% identity  -> 12,902 transfer-involved genes, A-A 9,344 / C-A 1,236 on ACs
#     21 strains    -> 17 clinical + 4 environmental, 4,484 variable HOGs,
#                      40,356 tests modelled / 31,388 in the FDR set, P >= 0.136
#
#   Every figure is written as PDF *and* 300-dpi PNG from the same object (project rule).
#
# Usage: Rscript scripts/plot_method_schematics.R [--only S26 S27 S28 S11]

suppressMessages({library(optparse); library(ggplot2); library(grid)})

option_list <- list(
  make_option("--out_dir", type = "character", default = "00.manuscript/supplementary"),
  make_option("--dpi", type = "numeric", default = 300),
  make_option("--family", type = "character", default = "Arial"),
  make_option("--only", type = "character", default = "S26,S27,S28,S11")
)
args <- parse_args(OptionParser(option_list = option_list))

# ---- shared visual language -----------------------------------------------------
FILL_BOX  <- "#E8F2F9"
EDGE_BOX  <- "#3D7EA6"
FILL_RAIL <- "#F2F5F7"
COL_TEXT  <- "#1A1A1A"
COL_ARROW <- "#4D4D4D"
FS_BASE   <- 7.6            # box text
FS_HEAD   <- 9.2            # section header
FS_TITLE  <- 10.5           # figure title

BOX_H <- 0.78               # box height in layout units
GAP   <- 0.30               # gap between boxes (room for the arrow)
SECT_GAP <- 0.62            # gap between sections
WRAP  <- 76                 # characters per line inside a full-width box

wrap_txt <- function(s, width = WRAP) {
  paste(strwrap(s, width = width), collapse = "\n")
}

# ---- layout engine --------------------------------------------------------------
# sections: list(list(title = , steps = c(...)))
build_schematic <- function(sections, title) {
  rows <- list(); heads <- list(); y <- 0; bands <- list()
  for (si in seq_along(sections)) {
    sec <- sections[[si]]
    y_head <- y
    y <- y - 0.55                                  # header consumes 0.55
    # the band must end ABOVE the first box top (y_head - 0.55); otherwise the lower
    # arrow piece starts below its own target and is drawn pointing upwards.
    bands[[length(bands) + 1]] <- data.frame(y0 = y_head + 0.12, y1 = y_head - 0.28)
    heads <- c(heads, list(data.frame(y = y_head - 0.09, label = sec$title)))
    for (k in seq_along(sec$steps)) {
      rows[[length(rows) + 1]] <- data.frame(
        ymin = y - BOX_H, ymax = y, label = wrap_txt(sec$steps[k]),
        first = k == 1, last = k == length(sec$steps), y = y - BOX_H / 2)
      y <- y - BOX_H - GAP
    }
    y <- y - SECT_GAP + GAP
  }
  boxes <- do.call(rbind, rows)
  heads <- do.call(rbind, heads)
  bands <- do.call(rbind, bands)

  # arrows: inside a section, box k -> k+1; between sections, last -> first
  # An arrow that crosses a section band would run straight through the section
  # header text (the band sits between the two boxes). Split such an arrow at the
  # band so the header stays uncontested.
  pieces <- list()
  for (i in seq_len(nrow(boxes) - 1)) {
    y0 <- boxes$ymin[i]; y1 <- boxes$ymax[i + 1]
    cut <- bands$y0 > y1 & bands$y1 < y0          # band lying between the two boxes
    if (any(cut)) {
      bd <- bands[cut, ][1, ]
      pieces[[length(pieces) + 1]] <- data.frame(x = 0.5, y0 = y0, y1 = bd$y0)
      pieces[[length(pieces) + 1]] <- data.frame(x = 0.5, y0 = bd$y1, y1 = y1 + 0.06)
    } else {
      pieces[[length(pieces) + 1]] <- data.frame(x = 0.5, y0 = y0, y1 = y1 + 0.06)
    }
  }
  seg <- do.call(rbind, pieces)

  ggplot() +
    geom_rect(data = bands, aes(xmin = 0.015, xmax = 0.985, ymin = y0, ymax = y1),
              fill = FILL_RAIL, colour = NA) +
    geom_rect(data = bands, aes(xmin = 0.015, xmax = 0.030, ymin = y0, ymax = y1),
              fill = EDGE_BOX, colour = NA) +
    geom_segment(data = seg, aes(x = x, xend = x, y = y0, yend = y1 + 0.06),
                 linewidth = 0.32, colour = COL_ARROW,
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    geom_rect(aes(xmin = 0.055, xmax = 0.945, ymin = boxes$ymin, ymax = boxes$ymax),
              fill = FILL_BOX, colour = EDGE_BOX, linewidth = 0.28) +
    geom_text(data = boxes, aes(x = 0.5, y = y, label = label),
              size = FS_BASE / 2.845, lineheight = 1.0, colour = COL_TEXT, family = args$family) +
    geom_text(data = heads, aes(x = 0.045, y = y, label = label),
              hjust = 0, size = FS_HEAD / 2.845, fontface = "bold",
              colour = EDGE_BOX, family = args$family) +
    annotate("text", x = 0.5, y = max(bands$y0) + 0.62, label = title,
             size = FS_TITLE / 2.845, fontface = "bold", colour = COL_TEXT, family = args$family) +
    coord_cartesian(xlim = c(0, 1),
                    ylim = c(min(boxes$ymin) - 0.10, max(bands$y0) + 0.85), clip = "off") +
    theme_void(base_size = FS_BASE, base_family = args$family) +
    theme(plot.margin = margin(3, 3, 3, 3))
}

render <- function(p, name, width_mm) {
  # height follows the canvas we laid out: 1 y-unit ~ 12 mm
  yr <- p$coordinates$limits$y
  h_mm <- (yr[2] - yr[1]) * 12.2
  w <- width_mm / 25.4; h <- h_mm / 25.4
  cairo_pdf(file.path(args$out_dir, paste0(name, ".pdf")),
            width = w, height = h, family = args$family)
  print(p); dev.off()
  ragg::agg_png(file.path(args$out_dir, paste0(name, ".png")),
                width = w, height = h, units = "in", res = args$dpi, background = "white")
  print(p); dev.off()
  cat(sprintf("wrote %s.pdf/.png  (%.0f x %.0f mm)\n", name, width_mm, h_mm))
}

# ---- the four figures -----------------------------------------------------------
FIG <- list()

FIG$S26 <- list(
  title = "Supplementary Figure S26: Pan-genome analysis workflow",
  sections = list(
    list(title = "1  Sample collection & genome assembly",
         steps = c("211 FSSC isolates: 140 BMU (132 clinical + 8 environmental) + 71 public",
                   "Illumina 2\u00d7150 bp | PacBio HiFi (9 strains) | ONT",
                   "SPAdes / hifiasm / Flye / nextPolish \u2192 BUSCO \u2265 95% complete")),
    list(title = "2  Repeat masking, gene prediction & annotation",
         steps = c("RepeatModeler2 \u2192 CD-HIT \u2192 RepeatMasker (species-specific library)",
                   "Funannotate gene prediction (UniProt + BUSCO proteins)",
                   "eggNOG-mapper 5.0 / dbCAN / Pfam / DFVF / PHI-base / SignalP 6")),
    list(title = "3  Pan-genome construction & clinical association",
         steps = c("OrthoFinder (DIAMOND + MCL) \u2192 29,080 HOGs, 5 frequency bins",
                   "PCoA of core CNV versus accessory PAV + PERMANOVA (Bray\u2013Curtis)",
                   "Random forest on PAV (70/30 split, balanced weights, 5-fold CV)",
                   "512 clinical-associated HOGs (Bonferroni P < 0.05)",
                   "KEGG / Pfam enrichment (clusterProfiler)")),
    list(title = "4  Species identification",
         steps = c("5-barcode phylogeny (TEF1 / RPB1 / RPB2 / CAM / TUB2; MAFFT + IQ-TREE)",
                   "FastANI (95% cutoff) + SNP phylogeny (ADMIXTURE, FST, AMOVA)"))
  ))

FIG$S27 <- list(
  title = "Supplementary Figure S27: Accessory-chromosome analysis workflow",
  sections = list(
    list(title = "1  AC identification (13 chromosome-level assemblies)",
         steps = c("Core-gene proportion \u00d7 repeat content \u2192 accessory / core boundary",
                   "Terminal repeat quantification (100-kb windows, 80\u2013120 kb scan)")),
    list(title = "2  Gene transfer & TE mediation (99% identity)",
         steps = c("mmseqs2 99%-identity clusters over 211 isolates (12,902 transfer genes)",
                   "A\u2013A / C\u2013A classification (9,344 / 1,236 among the 13,825 genes on ACs)",
                   "TE flanking enrichment (1 kb / 5 kb) + cassette permutation test",
                   "Hub chromosomes & intra- / cross-species cargo functions")),
    list(title = "3  k-mer classification (KMC + kmer-db)",
         steps = c("Shared k-mer matrix: 78 ACs \u00d7 211 isolates (19-mers)",
                   "block (species-stable) versus mosaic (cross-species)",
                   "Mosaic subtypes 1\u20133 (hierarchical clustering of sharing breadth)")),
    list(title = "4  Clinical association & case study",
         steps = c("Species-mix-corrected clinical enrichment \u2192 22 clinical ACs",
                   "RF clinical HOG location: mosaic AC 44.9% versus 33.3% (OR 1.63)",
                   "BMU02483: 5.9\u00d7 RF enrichment; 170 specific HOGs"))
  ))

FIG$S28 <- list(
  title = "Supplementary Figure S28: Pathogenicity assessment workflow",
  sections = list(
    list(title = "1  Galleria mellonella virulence model",
         steps = c("21 strains (17 clinical + 4 environmental) \u00d7 6 conditions",
                   "28 \u00b0C / 37 \u00b0C \u00d7 10\u00b3 / 10\u2074 / 10\u2075 conidia per larva",
                   "Kaplan\u2013Meier + log-rank + Cox HR; day-10 mortality and survival AUC",
                   "Temperature preference = mean mortality difference (28 \u00b0C \u2212 37 \u00b0C)")),
    list(title = "2  Virulence\u2013genotype association",
         steps = c("16 literature pathogenicity genes (mmseqs easy-search, fident \u2265 0.3)",
                   "Copy number and best-hit identity versus 6 phenotypes (Spearman, MW, FDR)",
                   "Genome-wide: 4,484 HOGs \u00d7 9 traits (40,356 tests; 31,388 in FDR set)",
                   "Permutation (\u2265 200\u00d7) + FDR + leave-one-out; candidates |\u03c1| \u2265 0.55",
                   "28 \u00b0C / 37 \u00b0C modules \u2192 Pfam hypergeometric enrichment")),
    list(title = "3  THP-1 macrophage co-culture (immunostimulation)",
         steps = c("PMA-differentiated THP-1 + spores (10% / 1% cell ratio, 24 h)",
                   "qPCR: IL-1\u03b2 / TNF-\u03b1 / CXCL1 / CSF2;   ELISA: TNF-\u03b1 / CSF2")),
    list(title = "4  Clinical versus environmental comparison",
         steps = c("Mann\u2013Whitney U across 13 phenotypes (P \u2265 0.136)"))
  ))

FIG$S11 <- list(
  title = "Supplementary Figure S11: Method \u2014 are transferred genes mediated by specific repeat families?",
  sections = list(
    list(title = "1  Transfer detection (211-genome context)",
         steps = c("All 211 FSSC genomes \u2192 mmseqs protein clusters at 99% identity",
                   "Cluster spans \u2265 2 chromosomes AND involves \u2265 1 accessory chromosome",
                   "AC-related transfer candidates: A\u2013A (accessory\u2013accessory), C\u2013A (core\u2013accessory)")),
    list(title = "2  Filter candidates on the 13 chromosome-level assemblies",
         steps = c("12,902 transfer-involved genes (A\u2013A 9,344 + C\u2013A 1,236 on ACs)",
                   "GFF3 genes + RepeatMasker (non-nested hits) + accessory-chromosome calls",
                   "Versus whole-genome background: 198,987 other genes")),
    list(title = "3  Flanking test (1 kb and 5 kb)",
         steps = c("For each gene: repeat families overlapping or flanking the gene",
                   "Per-family Fisher exact test: transfer candidates versus genome background",
                   "Conclusion: LINE \u00d717.2, RC \u00d715.6, DNA \u00d713.6\u201316.5, LTR \u00d712.1\u201313.5",
                   "Simple repeats and low-complexity sequence depleted (\u00d70.29\u20130.60)")),
    list(title = "4  Robustness: core-resident matched control",
         steps = c("Accessory-derived genes on core chromosomes versus other core genes",
                   "DNA \u00d710.6, LINE \u00d78.2, LTR \u00d77.9 (all P < 0.001)"))
  ))

# ---- run ------------------------------------------------------------------------
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
want <- strsplit(args$only, ",")[[1]]
for (s in names(FIG)) {
  if (!s %in% want) next
  p <- build_schematic(FIG[[s]]$sections, FIG[[s]]$title)
  render(p, paste0("Supplementary-Fig-", s, "-",
                   switch(s, S26 = "pangenome-method", S27 = "AC-method",
                          S28 = "pathogenicity-method", S11 = "method-abstract")),
         180)
}
