#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# Plot PAV heatmaps of the RF-selected clinical-associated HOGs for
# F. solani and F. falciformis, in the style of Figure 3E/H of the manuscript
# (pheatmap, ward.D2, viridis endpoints for 0/1 PAV), with the sample-side
# MedEnv annotation on the left and the gene-side annotation on top replaced
# by the EFFECTOR status of each HOG (EffectorP-3; HOG with >=1 predicted
# effector gene = 'Yes').
#
# Inputs:
#   --feature-matrix   long-format HOG count matrix (SampleID, HOG_ID, GeneCount)
#   --top-features     RF-selected HOG list for the species (one per line)
#   --sample-list      cohort sample list for the species (one per line)
#   --group-file       sample annotation (Sample, MedEnv, ...) bmu_ref_groups.211.txt
#   --effector-hog     effector_master.hog_level.tsv (HOG_ID, EffectorHOG)
#   --species          solani | falciformis
#   --output-prefix    output .pdf/.png prefix
#
# Usage:
#   Rscript plot_rf_heatmap.R \
#     --feature-matrix 05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv \
#     --top-features 05.pan-genome/classify_med_env_shorter/top_significant_hogs.solani.list \
#     --sample-list 05.pan-genome/classify_med_env_shorter/all.solani.sample.list \
#     --group-file bmu_ref_groups.211.txt \
#     --effector-hog 08.pan_020_genes/effector_master.hog_level.tsv \
#     --species solani \
#     --output-prefix 09.rf_clinical_genes/heatmap_rf.solani

suppressMessages({
  library(optparse)
  library(pheatmap)
  library(dplyr)
  library(tidyr)
  library(RColorBrewer)
  library(grid)
})

option_list <- list(
  make_option(c("--feature-matrix"), type = "character", default = NULL, dest = "feature_matrix"),
  make_option(c("--top-features"), type = "character", default = NULL, dest = "top_features"),
  make_option(c("--sample-list"), type = "character", default = NULL, dest = "sample_list"),
  make_option(c("--group-file"), type = "character", default = NULL, dest = "group_file"),
  make_option(c("--effector-hog"), type = "character", default = NULL, dest = "effector_hog"),
  make_option(c("--species"), type = "character", default = NULL),
  make_option(c("--output-prefix"), type = "character", default = "./heatmap", dest = "output_prefix"),
  make_option(c("--fig-w"), type = "double", default = 10.0, dest = "fig_w"),
  make_option(c("--fig-h"), type = "double", default = 9.34, dest = "fig_h"),
  make_option(c("--cellwidth"), type = "double", default = NULL, dest = "cellwidth"),
  make_option(c("--cellheight"), type = "double", default = NULL, dest = "cellheight"),
  make_option(c("--mini"), action = "store_true", default = FALSE, dest = "mini",
              help = "mini panel for the main figure: no sample or HOG labels, all annotation legends kept"),
  make_option(c("--fontsize"), type = "double", default = NULL, dest = "fontsize",
              help = "override the pheatmap fontsize (mini default 4)"),
  make_option(c("--drop-main"), action = "store_true", default = FALSE, dest = "drop_main",
              help = "mini only: omit the internal title (the manuscript caption carries it)")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$feature_matrix), !is.null(opt$top_features),
          !is.null(opt$sample_list), !is.null(opt$group_file),
          !is.null(opt$effector_hog), !is.null(opt$species))

# ---------- 1. top features (in importance order) ----------
top_feats <- readLines(opt$top_features)
top_feats <- top_feats[nzchar(trimws(top_feats))]
# the per-species classify output ships a header row and (optionally) an Importance column
if (grepl("\t", top_feats[1]) || grepl("^(Feature|HOG_ID)\b", top_feats[1], ignore.case = TRUE)) {
  cat("Dropping header row from", opt$top_features, "\n")
  top_feats <- top_feats[-1]
}
top_feats <- sub("\t.*$", "", top_feats)          # keep only the HOG id if columns follow
top_feats <- top_feats[nzchar(trimws(top_feats))]
cat("Top features:", length(top_feats), "\n")

# ---------- 2. sample list + MedEnv + Source annotation ----------
samples <- readLines(opt$sample_list)
samples <- samples[nzchar(trimws(samples))]
groups <- read.delim(opt$group_file, stringsAsFactors = FALSE)
colnames(groups)[1] <- "Sample"
med_env <- setNames(groups$MedEnv, groups$Sample)
source_map <- setNames(groups$Host, groups$Sample)
source_levels <- c("Human", "Plant", "Soil", "Reptile", "Insect", "Mammal",
                   "Stingray", "Nematode", "Arthropod", "Unknown")
source_vals <- source_map[samples]
source_vals[is.na(source_vals)] <- "Unknown"
# Displayed labels are clinical / non-clinical. The SOURCE DATA keeps its med / env values
# (bmu_ref_groups.*.txt is untouched); the relabelling happens here, at annotation-build time.
med_env_disp <- c(med = "clinical", env = "non-clinical", unknown = "unknown")
annot_row <- data.frame(
  `Clinical status` = factor(med_env_disp[med_env[samples]],
                             levels = c("clinical", "non-clinical", "unknown")),
  Source = factor(source_vals, levels = source_levels),
  row.names = samples,
  check.names = FALSE
)

# ---------- 3. effector annotation (gene-side, replaces Presence Rate) ----------
# EffectorP-3 class subdivided into Cytoplasmic / Apoplastic / Both /
# Non-effector based on the HOG-level counts of each class.
eff <- read.delim(opt$effector_hog, stringsAsFactors = FALSE)
eff_cyto <- setNames(eff$n_cytoplasmic, eff$HOG_ID)
eff_apo <- setNames(eff$n_apoplastic, eff$HOG_ID)
eff_class_of <- function(h) {
  cy <- eff_cyto[h]
  ap <- eff_apo[h]
  if (is.na(cy) || is.na(ap)) return("Non-effector")   # absent from effector table
  if (cy > 0 && ap > 0) return("Both")
  if (cy > 0) return("Cytoplasmic")
  if (ap > 0) return("Apoplastic")
  return("Non-effector")
}
eff_cat <- vapply(top_feats, eff_class_of, character(1))
eff_annot <- data.frame(
  Effector = factor(eff_cat, levels = c("Cytoplasmic", "Apoplastic", "Both",
                                        "Non-effector")),
  row.names = top_feats)

# ---------- 4. build PAV matrix (samples x features) ----------
mat <- read.delim(opt$feature_matrix, stringsAsFactors = FALSE)
colnames(mat)[1:3] <- c("SampleID", "HOG_ID", "GeneCount")
mat <- mat[mat$SampleID %in% samples & mat$HOG_ID %in% top_feats, ]
mat$Present <- as.numeric(mat$GeneCount > 0)
mat_pav <- tidyr::pivot_wider(mat[, c("SampleID", "HOG_ID", "Present")],
                              id_cols = "SampleID", names_from = "HOG_ID",
                              values_from = "Present") %>%
  as.data.frame()
rownames(mat_pav) <- mat_pav$SampleID
mat_pav$SampleID <- NULL
mat_pav[is.na(mat_pav)] <- 0
mat_pav <- mat_pav[, top_feats]              # importance order
mat_pav <- as.matrix(mat_pav)
# remove features whose PAV is constant across the species cohort (present in
# all or absent in all isolates): they carry no clinical/non-clinical signal
const_feats <- colnames(mat_pav)[apply(mat_pav, 2, function(x) length(unique(x)) == 1)]
if (length(const_feats) > 0) {
  cat("Removing", length(const_feats), "constant-PAV features:", paste(const_feats, collapse=", "), "\n")
  mat_pav <- mat_pav[, !colnames(mat_pav) %in% const_feats, drop = FALSE]
}
cat("Matrix:", nrow(mat_pav), "samples x", ncol(mat_pav), "features\n")
stopifnot(max(mat_pav) == 1)                  # PAV binary

# ---------- 5. colors ----------
annot_colors <- list(
  # The clinical-status bar sits directly above the Source bar, so it must not reuse the
  # saturated Source hues: med used to be literally the same hex as Human (#E41A1C) and env
  # clashed with the Source blues. Pastel1 red / blue / grey keeps the bars visually separate.
  `Clinical status` = c("clinical" = "#FBB4AE", "non-clinical" = "#B3CDE3",
                        "unknown" = "#CCCCCC"),
  Source = c("Human" = "#E41A1C", "Plant" = "#228B22", "Soil" = "#DAA520",
             "Reptile" = "#984EA3", "Insect" = "#FF8C00", "Mammal" = "#FFA726",
             "Stingray" = "#5F9EA0", "Nematode" = "#377EB8", "Arthropod" = "#8B0000",
             "Unknown" = "#555555"),
  Effector = c("Cytoplasmic" = "#E64B35", "Apoplastic" = "#4C72B0",
             "Both" = "#984EA3", "Non-effector" = "#F2F2F2")
)
pav_colors <- c("#440154", "#FDE725")         # viridis endpoints: 0 absent / 1 present

n_feat <- ncol(mat_pav)
n_samp <- nrow(mat_pav)
# cell sizes drive the rendered heatmap-body size (independent of the device).
# Default: square natural body (cellheight = cellwidth * n_cols/n_rows); the
# calibration loop then adjusts cell sizes until the RENDERED body is square.
# mini panels get a small fixed figure and no axis labels; the annotation legends
# (MedEnv, Source, Effector) and the PAV colour key are all retained.
mini <- isTRUE(opt$mini)
cellwidth_mm <- ifelse(is.null(opt$cellwidth), if (mini) 0.12 else 0.9, opt$cellwidth)
cellheight_mm <- ifelse(is.null(opt$cellheight), cellwidth_mm * n_feat / n_samp, opt$cellheight)
# The mini panels are rendered on a deliberately generous device and then trimmed to their
# ink box by 09.rf_clinical_genes/mini/crop_mini.py. Sizing the device to the gtable instead
# squeezed the annotation legends: pheatmap then silently dropped legend rows (the MedEnv
# legend and most of the Source legend never made it into the file).
fig_w <- if (mini && opt$fig_w == 10.0) 8.0 else opt$fig_w
fig_h <- if (mini && opt$fig_h == 9.34) 5.0 else opt$fig_h
title <- if (mini && isTRUE(opt$drop_main)) NA else paste("PAV of RF-selected HOGs (P < 0.05) -", opt$species)
fontsize <- if (!is.null(opt$fontsize)) opt$fontsize else if (mini) 4 else 11
fs_row <- if (mini) 4 else ifelse(n_samp > 80, 5, 6)
tree_h <- if (mini) 7 else 20

# ---- render -------------------------------------------------------------------------
# pheatmap fills only as much of the device as the cell sizes need, so a small figure is
# left with a large empty border. Building the gtable with silent = TRUE lets us measure
# its natural size and draw it tightly for the mini panels.
ph <- pheatmap(mat_pav,
               scale = "none",
               clustering_method = "ward.D2",
               cluster_rows = TRUE, cluster_cols = TRUE,
               treeheight_row = tree_h, treeheight_col = tree_h,
               cellwidth = cellwidth_mm, cellheight = cellheight_mm,
               color = pav_colors,
               annotation_row = annot_row,
               annotation_col = eff_annot,
               annotation_colors = annot_colors,
               show_rownames = !mini,
               show_colnames = FALSE,
               fontsize = fontsize, fontsize_row = fs_row, fontsize_col = fs_row,
               border_color = if (mini) NA else "grey90",
               main = title,
               silent = TRUE)
g <- ph$gtable

pdf(paste0(opt$output_prefix, ".pdf"), width = fig_w, height = fig_h)
grid::grid.newpage(); grid::grid.draw(g); dev.off()
png(paste0(opt$output_prefix, ".png"), width = fig_w, height = fig_h, units = "in", res = 300)
grid::grid.newpage(); grid::grid.draw(g); dev.off()
if (mini) {
  grDevices::svg(paste0(opt$output_prefix, ".svg"), width = fig_w, height = fig_h)
  grid::grid.newpage(); grid::grid.draw(g); dev.off()
  ragg::agg_tiff(paste0(opt$output_prefix, ".tif"), width = fig_w, height = fig_h, units = "in",
                 res = 600, background = "white")
  grid::grid.newpage(); grid::grid.draw(g); dev.off()
}

cat("Saved:", paste0(opt$output_prefix, ".pdf/.png"), "\n")
