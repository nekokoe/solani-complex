#!/usr/bin/env Rscript
# Method workflow figures (nature-figure, backend = R, grid-based)
# S2 / S8 / S9 ProcessOn-style flowcharts, Electron-UI palette
suppressMessages({library(grid)})

ELEC_BLUE  <- "#47848F"; ELEC_LIGHT <- "#9FEAF9"; ELEC_MID <- "#74B1BE"
ELEC_DARK  <- "#2C3E50"; ELEC_GREY <- "#6B7280"; ARROW_COL <- "#4A5568"

draw_box <- function(x, y, w, h, text, fill = ELEC_LIGHT, border = ELEC_BLUE,
                     text_col = ELEC_DARK, size = 8, bold = FALSE) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(2.5, "mm"), gp = gpar(fill = fill, col = border, lwd = 1.2))
  grid.text(text, x = unit(x, "npc"), y = unit(y, "npc"),
            gp = gpar(col = text_col, fontsize = size,
                      fontface = if (bold) "bold" else "plain", fontfamily = "sans"))
}
# 箭头: 从 (x, y1) 到 (x, y2), y1 > y2
draw_arrow <- function(x, y1, y2) {
  grid.lines(x = unit(c(x, x), "npc"), y = unit(c(y1, y2), "npc"),
             arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
             gp = gpar(col = ARROW_COL, lwd = 1.2, fill = ARROW_COL))
}
# 分组虚线框: 中心 (x,y) 高度 h, 标题在框内顶部
draw_group <- function(x, y, w, h, title) {
  grid.roundrect(x = unit(x, "npc"), y = unit(y, "npc"),
                 width = unit(w, "npc"), height = unit(h, "npc"),
                 r = unit(3, "mm"), gp = gpar(fill = "white", col = ELEC_GREY,
                                              lwd = 0.8, lty = "dashed"))
  grid.text(title, x = unit(x, "npc"), y = unit(y + h/2 - 0.028, "npc"),
            gp = gpar(col = ELEC_DARK, fontsize = 8.5, fontface = "bold", fontfamily = "sans"))
}
# 组内链式流程: boxes = 文本向量, y0 = 首框中心, step = 框间距
calc_chain <- function(h, boxes, y0, step) {
  n <- length(boxes)
  c(top = unname(y0) + h/2, bottom = unname(y0) - (n - 1) * step - h/2)
}
render_chain <- function(x, w, h, boxes, y0, step) {
  n <- length(boxes)
  for (k in seq_len(n)) {
    yk <- y0 - (k - 1) * step
    draw_box(x, yk, w, h, boxes[k])
    if (k < n) draw_arrow(x, yk - h/2 - 0.003, yk - step + h/2 + 0.003)
  }
}

W <- 0.60; CX <- 0.5; BH <- 0.040; STEP <- 0.058

# =============== S2 ===============
pdf("00.manuscript/supplementary/Supplementary-Fig-S2-pangenome-method.pdf", width = 7.5, height = 9)
grid.newpage(); pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.rect(gp = gpar(fill = "white", col = NA))
g1 <- calc_chain(BH, c("211 FSSC isolates: 140 BMU (132 clinical + 8 environmental) + 71 public",
                         "Illumina 2\u00d7150 bp  |  PacBio HiFi (9 strains)  |  ONT",
                         "SPAdes / hifiasm / Flye / nextPolish  \u2192  BUSCO \u2265 95% complete"), 0.93, STEP)
draw_group(CX, mean(c(g1["top"], g1["bottom"])), 0.86, g1["top"]-g1["bottom"]+0.085, "Sample collection & genome assembly")
g2 <- calc_chain(BH, c("RepeatModeler2 \u2192 CD-HIT \u2192 RepeatMasker (species-specific library)",
                         "Funannotate gene prediction (UniProt + BUSCO proteins)",
                         "eggNOG-mapper 5.0 / dbCAN / Pfam / DFVF / PHI-base / SignalP 6"), g1["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g2["top"], g2["bottom"])), 0.86, g2["top"]-g2["bottom"]+0.085, "Repeat masking, gene prediction & annotation")
g3 <- calc_chain(BH, c("OrthoFinder (DIAMOND + MCL) \u2192 29,080 HOGs, 5 frequency bins",
                         "PCoA: core CNV vs accessory PAV  +  PERMANOVA (Bray\u2013Curtis)",
                         "Random forest (PAV, 70/30 split, balanced weights, 5-fold CV)",
                         "512 clinical-associated HOGs (Bonferroni P < 0.05)",
                         "KEGG / Pfam enrichment (clusterProfiler)"), g2["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g3["top"], g3["bottom"])), 0.86, g3["top"]-g3["bottom"]+0.085, "Pan-genome construction & clinical association")
g4 <- calc_chain(BH, c("5-barcode phylogeny (TEF1/RPB1/RPB2/CAM/TUB2; MAFFT + IQ-TREE)",
                         "FastANI (95% cutoff) + SNP phylogeny (ADMIXTURE, FST, AMOVA)"), g3["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g4["top"], g4["bottom"])), 0.86, g4["top"]-g4["bottom"]+0.085, "Species identification")
render_chain(CX, W, BH, c("211 FSSC isolates: 140 BMU (132 clinical + 8 environmental) + 71 public",
                          "Illumina 2\u00d7150 bp  |  PacBio HiFi (9 strains)  |  ONT",
                          "SPAdes / hifiasm / Flye / nextPolish  \u2192  BUSCO \u2265 95% complete"), 0.93, STEP)
render_chain(CX, W, BH, c("RepeatModeler2 \u2192 CD-HIT \u2192 RepeatMasker (species-specific library)",
                          "Funannotate gene prediction (UniProt + BUSCO proteins)",
                          "eggNOG-mapper 5.0 / dbCAN / Pfam / DFVF / PHI-base / SignalP 6"), g1["bottom"]-0.10, STEP)
render_chain(CX, W, BH, c("OrthoFinder (DIAMOND + MCL) \u2192 29,080 HOGs, 5 frequency bins",
                          "PCoA: core CNV vs accessory PAV  +  PERMANOVA (Bray\u2013Curtis)",
                          "Random forest (PAV, 70/30 split, balanced weights, 5-fold CV)",
                          "512 clinical-associated HOGs (Bonferroni P < 0.05)",
                          "KEGG / Pfam enrichment (clusterProfiler)"), g2["bottom"]-0.10, STEP)
render_chain(CX, W, BH, c("5-barcode phylogeny (TEF1/RPB1/RPB2/CAM/TUB2; MAFFT + IQ-TREE)",
                          "FastANI (95% cutoff) + SNP phylogeny (ADMIXTURE, FST, AMOVA)"), g3["bottom"]-0.10, STEP)
# 组间箭头
draw_arrow(CX, g1["bottom"]-0.015, g2["top"]+0.015)
draw_arrow(CX, g2["bottom"]-0.015, g3["top"]+0.015)
draw_arrow(CX, g3["bottom"]-0.015, g4["top"]+0.015)
grid.text("Supplementary Figure S2: Pan-genome analysis workflow", x = 0.04, y = 0.988, just = "left",
          gp = gpar(col = ELEC_DARK, fontsize = 11, fontface = "bold", fontfamily = "sans"))
popViewport(); dev.off()

# =============== S8 ===============
pdf("00.manuscript/supplementary/Supplementary-Fig-S8-AC-method.pdf", width = 7.5, height = 9.5)
grid.newpage(); pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.rect(gp = gpar(fill = "white", col = NA))
g1 <- calc_chain(BH, c("Core-gene proportion \u00d7 repeat content \u2192 accessory / core boundary",
                         "Terminal repeat quantification (100-kb windows, 80\u2013120 kb scan)"), 0.93, STEP)
draw_group(CX, mean(c(g1["top"], g1["bottom"])), 0.86, g1["top"]-g1["bottom"]+0.085, "AC identification (13 chromosome-level assemblies)")
g2 <- calc_chain(BH, c("mmseqs2 100% identity clusters across 211 isolates (11,895 genes)",
                         "A\u2013A (8,988) / C\u2013A (2,907) transfer classification",
                         "TE flanking enrichment (1 kb / 5 kb) + cassette permutation test",
                         "Hub chromosomes & intra-/cross-species cargo functions"), g1["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g2["top"], g2["bottom"])), 0.86, g2["top"]-g2["bottom"]+0.085, "Gene transfer & TE mediation")
g3 <- calc_chain(BH, c("Shared k-mer matrix: 78 ACs \u00d7 211 isolates (19-mers)",
                         "block (species-stable) vs mosaic (cross-species)",
                         "Mosaic subtypes 1\u20133 (hierarchical clustering of sharing breadth)"), g2["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g3["top"], g3["bottom"])), 0.86, g3["top"]-g3["bottom"]+0.085, "k-mer classification (KMC + kmer-db)")
g4 <- calc_chain(BH, c("Species-mix\u2013corrected clinical enrichment \u2192 22 clinical ACs",
                         "RF clinical HOG location: mosaic AC 44.9% vs 33.3% (OR 1.63)",
                         "BMU02483: 5.9\u00d7 RF enrichment; 170 specific HOGs"), g3["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g4["top"], g4["bottom"])), 0.86, g4["top"]-g4["bottom"]+0.085, "Clinical association & case study")
render_chain(CX, W, BH, c("Core-gene proportion \u00d7 repeat content \u2192 accessory / core boundary",
                          "Terminal repeat quantification (100-kb windows, 80\u2013120 kb scan)"), 0.93, STEP)
render_chain(CX, W, BH, c("mmseqs2 100% identity clusters across 211 isolates (11,895 genes)",
                          "A\u2013A (8,988) / C\u2013A (2,907) transfer classification",
                          "TE flanking enrichment (1 kb / 5 kb) + cassette permutation test",
                          "Hub chromosomes & intra-/cross-species cargo functions"), g1["bottom"]-0.10, STEP)
render_chain(CX, W, BH, c("Shared k-mer matrix: 78 ACs \u00d7 211 isolates (19-mers)",
                          "block (species-stable) vs mosaic (cross-species)",
                          "Mosaic subtypes 1\u20133 (hierarchical clustering of sharing breadth)"), g2["bottom"]-0.10, STEP)
render_chain(CX, W, BH, c("Species-mix\u2013corrected clinical enrichment \u2192 22 clinical ACs",
                          "RF clinical HOG location: mosaic AC 44.9% vs 33.3% (OR 1.63)",
                          "BMU02483: 5.9\u00d7 RF enrichment; 170 specific HOGs"), g3["bottom"]-0.10, STEP)
draw_arrow(CX, g1["bottom"]-0.015, g2["top"]+0.015)
draw_arrow(CX, g2["bottom"]-0.015, g3["top"]+0.015)
draw_arrow(CX, g3["bottom"]-0.015, g4["top"]+0.015)
grid.text("Supplementary Figure S8: Accessory-chromosome analysis workflow", x = 0.04, y = 0.988, just = "left",
          gp = gpar(col = ELEC_DARK, fontsize = 11, fontface = "bold", fontfamily = "sans"))
popViewport(); dev.off()

# =============== S9 ===============
pdf("00.manuscript/supplementary/Supplementary-Fig-S9-pathogenicity-method.pdf", width = 7.5, height = 9.2)
grid.newpage(); pushViewport(viewport(width = unit(1, "npc"), height = unit(1, "npc")))
grid.rect(gp = gpar(fill = "white", col = NA))
g1 <- calc_chain(BH, c("15 strains (11 clinical + 4 environmental) \u00d7 6 conditions",
                         "28 \u00b0C / 37 \u00b0C \u00d7 10\u00b3 / 10\u2074 / 10\u2075 conidia per larva",
                         "Kaplan\u2013Meier + Log-rank + Cox HR; day-10 mortality & survival AUC",
                         "Temperature preference = mean mortality difference (28 \u00b0C \u2212 37 \u00b0C)"), 0.93, STEP)
draw_group(CX, mean(c(g1["top"], g1["bottom"])), 0.86, g1["top"]-g1["bottom"]+0.085, "Galleria mellonella virulence model")
g2 <- calc_chain(BH, c("16 literature pathogenicity genes (mmseqs easy-search, fident \u2265 0.3)",
                         "Copy number & best-hit identity vs 6 phenotypes (Spearman, MW, FDR)",
                         "Genome-wide: 3,630 variable HOGs \u00d7 9 traits (39,930 tests)",
                         "Permutation (100\u00d7) + FDR + leave-one-out; candidates |\u03c1| \u2265 0.55",
                         "28 \u00b0C / 37 \u00b0C modules \u2192 Pfam hypergeometric enrichment"), g1["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g2["top"], g2["bottom"])), 0.86, g2["top"]-g2["bottom"]+0.085, "Virulence\u2013genotype association")
g3 <- calc_chain(BH, c("PMA-differentiated THP-1 + spores (10% / 1% cell ratio, 24 h)",
                         "qPCR: IL-1\u03b2 / TNF-\u03b1 / CXCL1 / CSF2   ELISA: TNF-\u03b1 / CSF2"), g2["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g3["top"], g3["bottom"])), 0.86, g3["top"]-g3["bottom"]+0.085, "THP-1 macrophage co-culture (immunostimulation)")
g4 <- calc_chain(BH, c("Mann\u2013Whitney U across all phenotypes (P \u2265 0.10)"), g3["bottom"]-0.10, STEP)
draw_group(CX, mean(c(g4["top"], g4["bottom"])), 0.86, g4["top"]-g4["bottom"]+0.085, "Clinical vs environmental comparison")
render_chain(CX, W, BH, c("15 strains (11 clinical + 4 environmental) \u00d7 6 conditions",
                          "28 \u00b0C / 37 \u00b0C \u00d7 10\u00b3 / 10\u2074 / 10\u2075 conidia per larva",
                          "Kaplan\u2013Meier + Log-rank + Cox HR; day-10 mortality & survival AUC",
                          "Temperature preference = mean mortality difference (28 \u00b0C \u2212 37 \u00b0C)"), 0.93, STEP)
render_chain(CX, W, BH, c("16 literature pathogenicity genes (mmseqs easy-search, fident \u2265 0.3)",
                          "Copy number & best-hit identity vs 6 phenotypes (Spearman, MW, FDR)",
                          "Genome-wide: 3,630 variable HOGs \u00d7 9 traits (39,930 tests)",
                          "Permutation (100\u00d7) + FDR + leave-one-out; candidates |\u03c1| \u2265 0.55",
                          "28 \u00b0C / 37 \u00b0C modules \u2192 Pfam hypergeometric enrichment"), g1["bottom"]-0.10, STEP)
render_chain(CX, W, BH, c("PMA-differentiated THP-1 + spores (10% / 1% cell ratio, 24 h)",
                          "qPCR: IL-1\u03b2 / TNF-\u03b1 / CXCL1 / CSF2   ELISA: TNF-\u03b1 / CSF2"), g2["bottom"]-0.10, STEP)
render_chain(CX, W, BH, c("Mann\u2013Whitney U across all phenotypes (P \u2265 0.10)"), g3["bottom"]-0.10, STEP)
draw_arrow(CX, g1["bottom"]-0.015, g2["top"]+0.015)
draw_arrow(CX, g2["bottom"]-0.015, g3["top"]+0.015)
draw_arrow(CX, g3["bottom"]-0.015, g4["top"]+0.015)
grid.text("Supplementary Figure S9: Pathogenicity assessment workflow", x = 0.04, y = 0.988, just = "left",
          gp = gpar(col = ELEC_DARK, fontsize = 11, fontface = "bold", fontfamily = "sans"))
popViewport(); dev.off()
cat("S2/S8/S9 method flowcharts saved\n")
