#!/usr/bin/env Rscript
# S27 v4 (nature-figure, R + ComplexHeatmap) — distribution of
# temperature-bias-enriched HOGs across the 15 strains.
# Col annotation = strain group (37C-biased red / intermediate grey /
# 28C-biased blue). Row annotation = chitin-related (red) vs
# infection-relevant (orange) vs other (grey). Rows split into the two
# toolkits with block titles. Copy numbers shown in cells.

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
})

# Paths are resolved from this script's own location so the script runs from any
# working directory. It previously mixed cwd-relative (tmp/...) and root-relative
# (07.lethal_groups/...) paths, so it only worked from genomewide_assoc/.
.HERE <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                         value = TRUE)[1]))
BASE <- normalizePath(file.path(.HERE, "..", ".."))
if (!dir.exists(file.path(BASE, "07.lethal_groups"))) BASE <- getwd()

# strain sets read from the per-strain CMH classification (P < 0.05, stratified by
# inoculum dose) instead of a hard cutoff on the continuous effect size
.ph <- read.delim(file.path(BASE, "07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv"),
                  row.names = 1, check.names = FALSE)
.bias <- read.delim(file.path(BASE, "07.lethal_groups/phenotype_analysis/temperature_bias_cmh.tsv"),
                    stringsAsFactors = FALSE)
B37 <- sort(.bias$strain[.bias$class == "37C_biased"])
B28 <- sort(.bias$strain[.bias$class == "28C_biased"])
MID <- sort(setdiff(rownames(.ph), c(B37, B28)))
STRAINS <- c(B37, MID, B28)

d <- read.delim(file.path(BASE, "07.lethal_groups/genomewide_assoc/tmp/s27_long.tsv"), stringsAsFactors = FALSE)
m <- read.delim(file.path(BASE, "07.lethal_groups/genomewide_assoc/tmp/s27_meta.tsv"), stringsAsFactors = FALSE)
ann <- read.delim(file.path(BASE, "07.lethal_groups/genomewide_assoc/tmp/s27_ann.tsv"), stringsAsFactors = FALSE,
                  check.names = FALSE)

# single-line row labels: HOG ID followed by (gene name · function)
# — gene name and Pfam domain come only from the real annotation file
CHIT <- c("Glyco_hydro_18", "Glyco_hydro_19", "Glyco_hydro_20",
          "Glyco_hydro_12", "Chitin_synth", "chitin", "Glyco_hydro_28")
INF <- c("Hce2", "CFEM", "LysM", "HET", "Peptidase_S", "Effector")

# representative Pfam domain: chitin domain first for chitin HOGs, then
# infection-relevant domains, otherwise the first domain in the list
pick_repr <- function(pfs, chit, inf) {
  pfs <- pfs[pfs != "" & pfs != "-"]
  if (length(pfs) == 0) return("")
  if (chit) {
    m2 <- grep("Chitin_synth_2", pfs)
    if (length(m2)) return(pfs[m2[1]])
    for (k in CHIT) {
      hit <- grep(k, pfs, fixed = TRUE)
      if (length(hit)) return(pfs[hit[1]])
    }
  }
  if (inf) {
    for (k in INF) {
      hit <- grep(k, pfs, fixed = TRUE)
      if (length(hit)) return(pfs[hit[1]])
    }
  }
  pfs[1]
}

ord <- m[order(m$block), ]
aix <- match(ord$HOG, ann$HOG)
label2 <- vapply(seq_len(nrow(ord)), function(i) {
  gn <- strsplit(ann$Gene_Name[aix[i]], ";")[[1]][1]
  if (is.na(gn) || gn == "-") gn <- ""
  repr <- pick_repr(strsplit(ann$Pfam_Domain[aix[i]], ";")[[1]],
                    ord$chit[i], ord$inf[i])
  if (gn != "") gn else repr
}, character(1))
ord$label <- ifelse(is.na(label2) | label2 == "",
                    ord$HOG, paste0(ord$HOG, " (", label2, ")"))
mat <- matrix(0, nrow = nrow(ord), ncol = length(STRAINS),
              dimnames = list(ord$label, STRAINS))
for (i in seq_len(nrow(ord))) {
  h <- ord$HOG[i]
  mat[i, ] <- d$copy[d$HOG == h][match(STRAINS, d$strain[d$HOG == h])]
}
block <- factor(ord$block, levels = c("37C", "28C"))

# column annotation: strain temperature bias group
col_ann <- HeatmapAnnotation(
  bias = anno_simple(rep("37C", 2) -> g37,
                     col = c("37C" = "#C44E52")),
  show_annotation_name = FALSE)
grp <- c(rep("37C-biased", length(B37)), rep("intermediate", length(MID)),
         rep("28C-biased", length(B28)))
names(grp) <- STRAINS
col_ann <- HeatmapAnnotation(
  group = grp,
  col = list(group = c("37C-biased" = "#C44E52",
                       "intermediate" = "#BBBBBB",
                       "28C-biased" = "#3B7DD8")),
  annotation_legend_param = list(title = "Strain"),
  simple_anno_size = unit(0.35, "cm"), show_annotation_name = FALSE)

# row annotation: chitin / infection category
row_cat <- ifelse(ord$chit, "Chitin", ifelse(ord$inf, "Infection", "Other"))
row_ann <- rowAnnotation(
  cat = row_cat,
  col = list(cat = c("Chitin" = "#C44E52", "Infection" = "#E08E0B",
                     "Other" = "#D9D9D9")),
  annotation_legend_param = list(title = "HOG category"),
  width = unit(0.35, "cm"))

# copy-number colour scale: max observed = 2 (0: 560 cells, 1: 271, 2: 9)
col_fun <- colorRamp2(c(0, 1, 2), c("#FFFFFF", "#FDBE85", "#A63603"))
# Presence is shown by the cell fill alone: 660 of the 693 cells hold a copy number of
# exactly 1, so printing every value only added ink. Numerals are kept for the few
# multi-copy cells (>=2), where the copy number carries information.
cell_fun <- function(j, i, x, y, w, h, fill) {
  grid.rect(x, y, w, h, gp = gpar(fill = fill, col = NA))
  v <- mat[i, j]
  if (v >= 2)
    grid.text(v, x, y, gp = gpar(fontsize = 6, col = ifelse(v >= 3, "white", "grey15")))
}

ht <- Heatmap(
  mat, name = "Copy number",
  col = col_fun,
  rect_gp = gpar(type = "none"), cell_fun = cell_fun,
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = block,
  row_title = c(sprintf("37 \u00b0C-biased toolkit (n = %d)", sum(block == "37C")),
                sprintf("28 \u00b0C-biased toolkit (n = %d)", sum(block == "28C"))),
  row_title_gp = gpar(fontsize = 11, fontface = "bold",
                      col = c("#A63603", "#08306B")),
  row_gap = unit(6, "mm"),
  row_names_gp = gpar(fontsize = 8,
                      col = ifelse(ord$chit, "#C44E52",
                                   ifelse(ord$inf, "#E08E0B", "#444444")),
                      fontface = ifelse(ord$chit, "bold", "plain")),
  row_names_max_width = unit(10, "cm"),
  column_names_gp = gpar(fontsize = 9, fontface = "bold",
                         col = c(rep("#C44E52", length(B37)),
                                 rep("grey40", length(MID)),
                                 rep("#3B7DD8", length(B28)))),
  column_names_rot = 45,
  top_annotation = col_ann,
  left_annotation = row_ann,
  heatmap_legend_param = list(
    title = "Copy number", title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 8)),
  border = TRUE)

pdf(file.path(BASE, "07.lethal_groups/genomewide_assoc/figure_hog_bias_S27v4.pdf"),
    width = 8, height = 12)
draw(ht, merge_legend = TRUE,
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()
png(file.path(BASE, "07.lethal_groups/genomewide_assoc/figure_hog_bias_S27v4.png"),
    width = 1600, height = 2400, res = 200)
draw(ht, merge_legend = TRUE,
     annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()
cat("S27 v4 saved (ComplexHeatmap)\n")
