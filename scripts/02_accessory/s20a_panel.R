#!/usr/bin/env Rscript
# s20a_panel.R — the standalone k-mer panels of the S20 family.
#
#   S20a  the k-mer classification landscape of the 78 accessory chromosomes
#   S20g  F. falciformis: AC sharing with clinical vs non-clinical isolates
#   S20h  F. solani:      AC sharing with clinical vs non-clinical isolates
#
# S20a is panel a of the delivered S20 composite, so it reuses that letter. S20g/S20h
# are NOT panels of that composite (which already has a-f) and therefore take the next
# free letters rather than b/c/d, which would read as extractions of panels b/c/d.
#
# Panel a uses one data function and two render modes:
#
#   s20a_panel()       standalone square panel (120 mm)
#   s20a_panel_wide()  2:1 panel for the main figure: the plot occupies the LEFT
#                      half and the whole legend the RIGHT half, so the plot stays
#                      square while the three guides get a full half to themselves
#
# ENCODING (revised)
#   x      mean k-mer sharing with isolates of the SAME species (%)
#   y      mean k-mer sharing with isolates of OTHER species (%)
#   colour AC source: clinical vs environmental host isolate
#   shape  AC class: block vs mosaic  (>=95% of the carriers from one species)
#   size   number of carrier isolates (k-mer sharing >= 0.4)
#
# Why the axes changed: the previous version plotted intra_/inter_breadth_0.4 —
# "the fraction of isolates whose k-mer sharing with this AC is >= 0.4". That is a
# THRESHOLDED CARRIER PREVALENCE, and the threshold is invisible in the plot, which
# is what makes the axis hard to read. The table already carries the continuous
# quantity (intra_share / inter_share = the MEAN sharing), so the axes now use that:
# a plain 0-100% similarity with no hidden cut-off. The carrier count keeps the
# >=0.4 definition, and it is named explicitly in the size legend.
#
# Cost of the switch, measured: the old thresholded quantity separated block from
# mosaic perfectly on the cross-species axis (P(block<mosaic) = 1.00); the mean
# separates them at 0.90. Both are strong; the mean buys readability with that.
#
# The transfer-mode encoding (block transfer / scattered transfer) was dropped: it
# was read from suppl_table.ac_transfer_profile.tsv and duplicated the block/mosaic
# distinction without adding information.
#
# Sources
#   06.accessory-chrs/kmer_analysis/kmer_type.dual_domain.tsv        intra_/inter_share
#   06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv  frac_top_species
#   06.accessory-chrs/kmer_analysis/kmer_matrix.ACx211.tsv          carrier counts
#   05.pan-genome/bmu_ref_groups.211.txt                             MedEnv of the host isolate

suppressMessages({library(ggplot2); library(patchwork); library(cowplot)})

S20_CLIN <- "#C44E52"    # clinical
S20_ENV  <- "#3B7DD8"    # environmental

s20a_data <- function(dual_path, class_path, matrix_path, meta_path) {
  dd <- read.delim(dual_path, stringsAsFactors = FALSE)
  v3 <- read.delim(class_path, stringsAsFactors = FALSE)
  mat <- read.delim(matrix_path, row.names = 1, check.names = FALSE)
  meta <- read.delim(meta_path, stringsAsFactors = FALSE)

  n_carriers <- apply(mat >= 0.4, 1, sum)
  a <- merge(dd, v3[, c("AC", "kmer_type", "frac_top_species")], by = "AC")
  a$n_carriers <- n_carriers[a$AC]

  # AC id -> host isolate -> clinical/environmental source
  a$iso <- sub("_[^_]+$", "", a$AC)
  a <- merge(a, meta[, c("BMU_ID", "MedEnv")], by.x = "iso", by.y = "BMU_ID",
             all.x = TRUE)
  if (any(is.na(a$MedEnv)))
    warning("no MedEnv for: ", paste(unique(a$AC[is.na(a$MedEnv)]), collapse = ", "))

  # block = at least 95% of the carriers come from a single species
  a$cls <- factor(ifelse(a$frac_top_species >= 0.95, "Block", "Mosaic"),
                  levels = c("Block", "Mosaic"))
  a$source <- factor(ifelse(a$MedEnv == "med", "Clinical", "Environmental"),
                     levels = c("Clinical", "Environmental"))
  # continuous similarity, as a percentage
  a$share_in <- 100 * a$intra_share
  a$share_out <- 100 * a$inter_share
  a
}

XLAB <- "Mean k-mer sharing with the same species (%)"
YLAB <- "Mean k-mer sharing with other species (%)"
CLIN_LAB <- "Mean k-mer sharing with the clinical isolates (%)"
ENV_LAB  <- "Mean k-mer sharing with the non-clinical isolates (%)"

# Per-species view: does an AC share more with the CLINICAL or the NON-CLINICAL
# isolates of its own species?
#
# The AC's own donor isolate is EXCLUDED from both means, and that is not cosmetic:
# the matrix gives the donor a sharing of exactly 1.0, so including it would inflate
# whichever group the donor belongs to -- a bias that is collinear with the AC's own
# source and would fabricate a clinical preference. For F. solani it also changes the
# result materially (26/42 ACs with the donor, 30/42 without).
s20a_clinenv_data <- function(matrix_path, meta_path, species, class_path) {
  # `species` may be a vector, and the pooling is deliberately SPECIES-BLIND: the ACs
  # of every listed species go into one panel and every AC is averaged over the SAME
  # pooled isolate sets (all the listed species' clinical isolates, and all their
  # non-clinical ones). Species boundaries are not respected anywhere, by request.
  # (Comparing each AC only against its OWN species is a different question, and it
  # is what makes F. keratoplastica look like it goes the other way: that species
  # contributes only 4 non-clinical isolates, so its own-species means are unstable.)
  mat <- read.delim(matrix_path, row.names = 1, check.names = FALSE)
  meta <- read.delim(meta_path, stringsAsFactors = FALSE)
  v3 <- read.delim(class_path, stringsAsFactors = FALSE)
  spp <- setNames(as.character(meta$Species), as.character(meta$BMU_ID))
  med <- setNames(as.character(meta$MedEnv), as.character(meta$BMU_ID))

  acs <- rownames(mat)[spp[sub("_[^_]+$", "", rownames(mat))] %in% species]
  clin_all <- names(spp)[!is.na(spp) & spp %in% species & med == "med"]
  env_all  <- names(spp)[!is.na(spp) & spp %in% species & med == "env"]

  n_carriers <- apply(mat >= 0.4, 1, sum)
  rows <- lapply(acs, function(a) {
    own <- sub("_[^_]+$", "", a)
    spc <- spp[own]                      # recorded for reporting only
    clin <- clin_all                     # the POOLED set: no species restriction
    env  <- env_all
    data.frame(AC = a, species = spc,
               n_carriers = n_carriers[a],
               # unlist: indexing a one-row data.frame returns a data.frame, and
               # mean() of that is NA, not an error -- a silent failure
               share_clin = 100 * mean(unlist(mat[a, setdiff(clin, own)])),
               share_env  = 100 * mean(unlist(mat[a, setdiff(env, own)])))
  })
  a <- do.call(rbind, rows)
  a <- merge(a, v3[, c("AC", "frac_top_species")], by = "AC")
  a$cls <- factor(ifelse(a$frac_top_species >= 0.95, "Block", "Mosaic"),
                  levels = c("Block", "Mosaic"))
  a$source <- factor(ifelse(med[sub("_[^_]+$", "", a$AC)] == "med",
                            "Clinical", "Environmental"),
                     levels = c("Clinical", "Environmental"))
  attr(a, "n_clin") <- length(clin_all); attr(a, "n_env") <- length(env_all)
  a
}

s20a_panel <- function(a, family = "Arial", fontsize = 11,
                       legend_position = "bottom",
                       title = "K-mer classification landscape (78 ACs)",
                       xcol = "share_in", ycol = "share_out",
                       xlab = XLAB, ylab = YLAB, ellipses = TRUE,
                       diagonal = FALSE,
                       diagonal_labels = c("shared more with the x group",
                                           "shared more with the y group")) {
  base <- ggplot(a, aes(x = .data[[xcol]], y = .data[[ycol]]))

  # y = x is the null line for a PAIRED comparison: if an AC is shared equally with
  # the two groups being compared, its two means coincide and the point sits on the
  # diagonal. Points BELOW it are shared more with the x group, points ABOVE it more
  # with the y group, so the side of the line is the direction of any preference.
  # Both axes are forced to a common range so the diagonal is visually 45 degrees;
  # with different limits the line is still y = x but the quadrants read wrong.
  if (isTRUE(diagonal)) {
    rng <- range(c(a[[xcol]], a[[ycol]]), na.rm = TRUE)
    pad <- diff(rng) * 0.06
    lim <- c(rng[1] - pad, rng[2] + pad)
    base <- base +
      geom_abline(intercept = 0, slope = 1, linetype = "22",
                  colour = "grey35", linewidth = 0.4) +
      # coord_fixed, not just equal limits: the panel is ~5% wider than tall, so
      # with equal limits alone y = x renders at ~43.6 degrees and the quadrants
      # look skewed. Fixing the aspect makes the null line exactly 45 degrees.
      coord_fixed(xlim = lim, ylim = lim)
  }

  # 95% data ellipses for the two classes. They go on the BASE so they are drawn
  # UNDER the points. Outlines in one neutral grey, separated by linetype rather
  # than by colour: the colour channel already encodes the AC source, and a second
  # colour scale would make the two legends look interchangeable when they are not.
  # type = "norm" is the bivariate-normal DATA ellipse (~95% of the class's ACs),
  # not a confidence ellipse of the centroid, which would be far smaller and could
  # not show the overlap. NOTE the layer must be appended to `base` in its own
  # statement: writing `ggplot(...) + if (ellipses) stat_ellipse(...) + geom_point()`
  # makes the IF swallow the whole remaining chain, and ellipses = FALSE then returns
  # NULL for every layer.
  if (isTRUE(ellipses)) {
    base <- base +
      stat_ellipse(aes(group = cls, linetype = cls), type = "norm", level = 0.95,
                   colour = "grey25", linewidth = 0.4, show.legend = TRUE) +
      scale_linetype_manual(values = c(Block = "solid", Mosaic = "22"),
                            name = "95% ellipse")
  }

  p <- base +
    geom_point(aes(size = n_carriers, shape = cls, colour = source), alpha = 0.85) +
    scale_colour_manual(values = c(Clinical = S20_CLIN, Environmental = S20_ENV),
                        name = "AC source") +
    scale_shape_manual(values = c(Block = 17, Mosaic = 16), name = "AC class") +
    scale_size_continuous(range = c(1.2, 4.5),
                          name = "Carriers (k-mer \u2265 0.4)") +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = fontsize - 1, colour = "black"),
      axis.title = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold", hjust = 0),
      legend.position = legend_position,
      legend.text = element_text(size = fontsize - 1),
      legend.title = element_text(size = fontsize - 1),
      legend.key.size = unit(4, "mm"),
      legend.box.spacing = unit(1, "mm"),
      legend.margin = margin(t = -1),
      legend.spacing.x = unit(1.5, "mm"),
      legend.box = "vertical",
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 3, 4)
    )

  # guides, in a fixed order; the ellipse guide is only requested when it exists
  if (isTRUE(diagonal)) {
    p <- p + annotate("text", x = lim[2] - diff(lim) * 0.02, y = lim[1] + diff(lim) * 0.04,
                      label = diagonal_labels[1], hjust = 1, vjust = 0,
                      size = fontsize / 2.845, colour = "grey35") +
             annotate("text", x = lim[1] + diff(lim) * 0.02, y = lim[2] - diff(lim) * 0.04,
                      label = diagonal_labels[2], hjust = 0, vjust = 1,
                      size = fontsize / 2.845, colour = "grey35")
  }

  g <- list(colour = guide_legend(order = 1),
            shape  = guide_legend(order = 2),
            size   = guide_legend(order = 4))
  if (isTRUE(ellipses))
    g$linetype <- guide_legend(order = 3, override.aes = list(linewidth = 0.6))
  p + do.call(guides, g)
}

# 2:1 panel: plot on the left half, the entire legend on the right half.
s20a_panel_wide <- function(a, family = "Arial", fontsize = 11,
                            title = "K-mer classification landscape (78 ACs)",
                            xcol = "share_in", ycol = "share_out",
                            xlab = XLAB, ylab = YLAB, ellipses = TRUE,
                            diagonal = FALSE,
                            diagonal_labels = c("shared more with the x group",
                                                "shared more with the y group")) {
  p <- s20a_panel(a, family = family, fontsize = fontsize,
                  legend_position = "right", title = title,
                  xcol = xcol, ycol = ycol, xlab = xlab, ylab = ylab,
                  ellipses = ellipses, diagonal = diagonal,
                  diagonal_labels = diagonal_labels)
  leg <- get_legend(p + theme(legend.box = "vertical",
                              legend.box.spacing = unit(2, "mm"),
                              legend.key.size = unit(4.5, "mm")))
  if (is.null(leg)) stop("could not extract the legend")
  (p + theme(legend.position = "none")) + wrap_elements(leg) +
    plot_layout(widths = c(1, 1))
}
