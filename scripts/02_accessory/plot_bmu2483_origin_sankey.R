#!/usr/bin/env Rscript
# plot_bmu2483_origin_sankey.R — where BMU02483's genes could have come from.
#
# FOUR columns, coloured by the donor's ISOLATION BACKGROUND:
#
#   BMU02483 genes -> isolate source -> species source -> known chromosome source
#
#   species     falciformis / inter-species / unknown
#   isolate     from nearby env isolates / from clinical isolates / both / unknown
#   chromosome  ACs / CCs / unknown
#
#   species "inter-species" = at least one donor belongs to a species other than
#   F. falciformis, including genes that still carry a falciformis donor as well;
#   "unknown" = the gene clusters with no other isolate at all. isolate "unknown" =
#   no donor, or donors that are environmental isolates outside the +/-10 window.
#
# Read from suppl_table.bmu2483_gene_origin_flows_99.tsv, which
# scripts/bmu2483_gene_origin_sankey.py writes off the >=99%-identity clustering.
#
# DESIGN NOTES
#   * Portrait, 3.2:4 (width:height). The height is the 3:4 one and only the width grew,
#     so the inter-column gutters (the flows) are ~15% wider than at 3:4 while the
#     labels keep the vertical room they had. Four columns in a narrow canvas cannot take
#     horizontal stratum labels, so every label is rotated 90 degrees.
#   * No axes: the value axis carries no information the labels do not (each label
#     states its own gene count), so theme_void removes the axis, its ticks and the
#     axis title. The stage names survive as axis.text.x on top.
#   * KEEPING THE TEXT OFF THE RIBBONS. A label is put INSIDE its lode whenever the
#     lode is tall enough, which is the only truly ribbon-free place in a middle
#     column. Whether it fits is decided with REAL text metrics (strwidth on a
#     scratch device at the final size), not a per-character estimate: measured
#     widths run from 0.25 to 0.57 x font size depending on the characters, so any
#     single factor misplaces labels in both directions. A label that still does not
#     fit goes into the inter-column gutter, but the LAST column's gutters are
#     ribbon-free (the flows end at that lode), so its labels are drawn as plain
#     text there with no backing box; only a middle-column overflow gets a white box.
#   * Label placement is DETERMINISTIC and never guessed: the stratum rectangles are
#     read back from ggplot_build(), so a label is anchored to the lode ggplot actually
#     drawn rather than to a re-derived cumsum. A node too small to label in place gets
#     a BOLD, larger, one-line call-out in the gutter, so that it is not overlooked.
#
# Usage: Rscript scripts/plot_bmu2483_origin_sankey.R [--size_w 100] [--size_h 133.33]

suppressMessages({library(optparse); library(ggplot2); library(ggalluvial)})

option_list <- list(
  make_option(c("-i","--in_tsv"), type="character",
              default="06.accessory-chrs/BMU02483_gene_content/suppl_table.bmu2483_gene_origin_flows_99.tsv"),
  make_option(c("-o","--out"), type="character",
              default="06.accessory-chrs/BMU02483_gene_content/fig_s35_bmu2483_origin_sankey"),
  make_option("--size_w", type="numeric", default=106.67),  # 3.2:4 with the 133.33 mm height
  make_option("--size_h", type="numeric", default=133.33),
  # Stratum width in x units. Narrower strata = wider inter-column gutters (the flows);
  # the two trade against each other one-for-one because the canvas width is fixed.
  make_option("--lode_w", type="numeric", default=0.46),
  make_option("--dpi", type="numeric", default=600),
  make_option("--fontsize", type="numeric", default=8),
  # Labels that cannot fit inside their node go into the inter-column gutter, where they
  # belong to a SMALL node that is otherwise easy to miss: draw them bold and larger so
  # they read as a call-out rather than as a stray label.
  make_option("--gutter_fontsize", type="numeric", default=10),
  # Ribbon width floor, OFF by default: flows are drawn at their true proportion.
  # A non-zero value (e.g. 0.5) lifts sub-pixel flows into view, at the cost of
  # over-representing them; it is a DISPLAY-only change, the counts never move.
  make_option("--min_ribbon_pt", type="numeric", default=0),
  make_option("--family", type="character", default="Arial")
)
args <- parse_args(OptionParser(option_list=option_list))

d <- read.delim(args$in_tsv, stringsAsFactors = FALSE)
names(d) <- c("spec", "iso", "chrom", "n")

SPEC_LEV  <- c("falciformis", "inter-species", "unknown")
ISO_LEV   <- c("from nearby env isolates", "from clinical isolates", "both", "unknown")
CHROM_LEV <- c("ACs", "CCs", "unknown")

for (v in c("spec","iso","chrom")) {
  lv <- get(c(spec= "SPEC_LEV", iso = "ISO_LEV", chrom = "CHROM_LEV")[[v]])
  miss <- setdiff(unique(d[[v]]), lv)
  if (length(miss)) stop("unmapped category in ", v, ": ", paste(miss, collapse=", "))
  d[[v]] <- factor(d[[v]], levels = lv)
}

# short names for the in-figure labels; the full wording goes in the caption
SHORT <- c("falciformis" = "falciformis", "inter-species" = "inter-species",
           "unknown" = "unknown", "both" = "both",
           "from nearby env isolates" = "nearby env",
           "from clinical isolates" = "clinical",
           "ACs" = "ACs", "CCs" = "CCs")

d$flow <- paste(d$spec, d$iso, d$chrom, sep = " | ")

# stratum names MUST be unique across the whole plot: ggalluvial identifies a stratum
# by name, so the bare "unknown" labels of three columns would be merged into one node.
sid <- function(axis, v) paste0(axis, "\u2502", v)

# Column order: genes -> isolate source -> species source -> chromosome source.
# Each row is built by NAME, never by AX[i]: an index-based build silently pairs the
# wrong stratum with an axis as soon as this order is changed.
AX <- c("genes", "isolate", "species", "chromosome")
long <- rbind(
  data.frame(axis = "genes",      stratum = "BMU02483",                value = d$n, iso = d$iso, flow = d$flow),
  data.frame(axis = "isolate",    stratum = sid("isolate", d$iso),     value = d$n, iso = d$iso, flow = d$flow),
  data.frame(axis = "species",    stratum = sid("species", d$spec),    value = d$n, iso = d$iso, flow = d$flow),
  data.frame(axis = "chromosome", stratum = sid("chromosome", d$chrom),value = d$n, iso = d$iso, flow = d$flow)
)
long$axis    <- factor(long$axis, levels = AX)
long$stratum <- factor(long$stratum,
                       levels = c("BMU02483",
                                  sid("isolate", ISO_LEV), sid("species", SPEC_LEV),
                                  sid("chromosome", CHROM_LEV)))
# label text per stratum id
LAB <- c("BMU02483" = "16,864")
for (v in SPEC_LEV)  LAB[sid("species", v)]    <- SHORT[[v]]
for (v in ISO_LEV)   LAB[sid("isolate", v)]    <- SHORT[[v]]
for (v in CHROM_LEV) LAB[sid("chromosome", v)] <- SHORT[[v]]

ISO_COL <- c("from nearby env isolates" = "#3B7DD8",
             "from clinical isolates"    = "#C44E52",
             "both"                      = "#8E6BB5",
             "other"                     = "#9A9A9A",
             "unknown"                   = "#9A9A9A")

LODE_W <- args$lode_w   # stratum width, in x units (half-width = LODE_W / 2)
GAP_X  <- 0.40  # offset of an overflow label from the column centre

W <- args$size_w / 25.4; H <- args$size_h / 25.4
mm_per_pt  <- 25.4 / 72
panel_h_pt <- (args$size_h - 20 - 8) / mm_per_pt      # ~20 mm header + ~8 mm margins

tot <- sum(long$value[long$axis == "genes"])          # true gene total
gpp0 <- tot * 1.02 / panel_h_pt                       # display units per point, unfloored
floor_g <- if (args$min_ribbon_pt > 0) ceiling(args$min_ribbon_pt * gpp0) else 0
long$w <- pmax(long$value, floor_g)                   # same floor on every axis -> consistent
n_floor <- sum(long$value < floor_g & long$axis == "genes")

p0 <- ggplot(long, aes(x = axis, stratum = stratum, alluvium = flow, y = w)) +
  geom_alluvium(aes(fill = iso), width = LODE_W, alpha = 0.72, curve_type = "cubic") +
  geom_stratum(width = LODE_W, fill = "grey94", colour = "grey35", linewidth = 0.3) +
  scale_fill_manual(values = ISO_COL, name = NULL) +
  scale_x_discrete(position = "top",
                   labels = c(genes = "BMU02483\ngenes", species = "species\nsource",
                              isolate = "isolate\nsource",
                              chromosome = "known chromosome\nsource")) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.01))) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# ---- read the lode geometry back out of ggplot ---------------------------------
b <- ggplot_build(p0)
lod <- b$data[[2]][, c("x", "ymin", "ymax", "y", "stratum")]
lod$stratum <- as.character(lod$stratum)
tot_w <- max(lod$ymax)                                # display total (>= true total)

# ---- REAL text metrics, not a per-character guess ------------------------------
# strwidth needs an open device with the final font, so measure before drawing.
fmt <- function(x) formatC(as.integer(round(x)), format = "d", big.mark = ",")
st_tot <- tapply(long$value, long$stratum, sum)
# The gene-pool node carries its own name too, in the same two-line form as the
# category nodes, so the first column reads like the rest.
lod$label <- ifelse(lod$stratum == "BMU02483",
                    sprintf("BMU02483\n(%s)", fmt(tot)),
                    sprintf("%s\n%s", LAB[lod$stratum], fmt(st_tot[lod$stratum])))

scratch <- tempfile(fileext = ".pdf")
cairo_pdf(scratch, width = W, height = H, family = args$family, pointsize = args$fontsize)
line_pt <- function(s) {
  vapply(strsplit(s, "\n")[[1]],
         function(x) strwidth(x, units = "inches") * 72, numeric(1))
}
lod$len_pt <- vapply(lod$label, function(s) max(line_pt(s)), numeric(1))
lod$nline  <- vapply(lod$label, function(s) length(strsplit(s, "\n")[[1]]), integer(1))
# the same label folded onto ONE line ("ACs 495"), which is what a gutter call-out
# uses; measured in BOLD at the gutter size, because that is how it will be DRAWN.
# Measuring it at the in-node size and weight under-estimated it by ~30% and let the
# label run past the panel edge, where the device silently clipped it.
lod$lab1   <- gsub("\n", " ", lod$label)
par(font = 2)
lod$glen_pt <- vapply(lod$lab1, function(s)
  strwidth(s, units = "inches", cex = args$gutter_fontsize / args$fontsize) * 72, numeric(1))
par(font = 1)
dev.off(); unlink(scratch)

# pt -> data units, from the panel's own height
genes_per_pt <- tot_w * 1.02 / panel_h_pt

lod$h    <- lod$ymax - lod$ymin
lod$h_pt <- lod$h / genes_per_pt

# lode width in points, from the panel the plot actually got
panel_w_pt <- args$size_w / 25.4 * 72 - 10            # minus the 5 pt plot margins
pt_per_x   <- panel_w_pt / diff(b$layout$panel_params[[1]]$x.range)
lode_w_pt  <- LODE_W * pt_per_x

# Three placements, tried in order. A wide-but-short lode (e.g. isolate "clinical",
# 32 x 20 pt) cannot hold a ROTATED label but holds a two-line HORIZONTAL one, so it
# still lands inside the node and never touches a ribbon.
lod$rot_fit   <- lod$h_pt >= lod$len_pt + 3
lod$horiz_h   <- lod$nline * args$fontsize * 0.92
lod$horiz_fit <- !lod$rot_fit & lod$h_pt >= lod$horiz_h + 2 &
                 lode_w_pt >= lod$len_pt + 6

# Everything that fits neither in-node form goes into the gutter as a BOLD, larger,
# one-line call-out: a node too small to be labelled in place is exactly the one a
# reader would otherwise miss, so it is worth marking out (chromosome "ACs", 8.6 pt).
inn   <- lod[lod$rot_fit, ]                       # rotated, inside the lode
hor   <- lod[lod$horiz_fit, ]                     # two lines, inside the lode
out   <- lod[!lod$rot_fit & !lod$horiz_fit, ]     # call-out in the gutter

if (nrow(out)) {
  # Clamp against the panel using each label's OWN half length, measured at the size
  # and weight it is actually drawn at: clamping the centre alone let a top-most label
  # run past the panel edge, where the device silently clipped it.
  g   <- genes_per_pt
  pad <- 2 * g
  half <- out$glen_pt / 2 * g
  out$y_lab <- pmin(pmax(out$y, half + pad), tot_w - half - pad)
  for (ax in unique(out$x)) {                    # only labels in one gutter collide
    i <- which(out$x == ax)
    o <- i[order(-out$y_lab[i])]                 # top to bottom
    if (length(o) > 1) {
      for (k in seq_along(o)[-1]) {              # push down
        need <- (out$glen_pt[o[k - 1]] + out$glen_pt[o[k]]) / 2 * g + pad
        if (out$y_lab[o[k - 1]] - out$y_lab[o[k]] < need)
          out$y_lab[o[k]] <- out$y_lab[o[k - 1]] - need
      }
      for (k in rev(seq_along(o))[-1]) {         # lift a stack that overflowed
        lo_k <- out$glen_pt[o[k]] / 2 * g + pad
        if (out$y_lab[o[k]] < lo_k) {
          out$y_lab[o[k]] <- lo_k
          need <- (out$glen_pt[o[k - 1]] + out$glen_pt[o[k]]) / 2 * g + pad
          if (out$y_lab[o[k]] - out$y_lab[o[k - 1]] < need)
            out$y_lab[o[k - 1]] <- out$y_lab[o[k]] + need
        }
      }
    }
  }
  # The flows END at the last column, so the strip to the right of its lodes carries
  # no ribbon: those overflow labels are plain text. Every other gutter is crossed by
  # ribbons, so those labels get a tight white box to stay readable.
  out$boxed <- out$x < max(lod$x)
}

p <- p0 +
  { if (nrow(out)) geom_segment(data = out,
      aes(x = x + LODE_W / 2, xend = x + GAP_X - 0.06, y = y, yend = y_lab),
      inherit.aes = FALSE, linewidth = 0.25, colour = "grey25") } +
  { if (nrow(inn)) geom_text(data = inn,
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      angle = 90, size = args$fontsize / 2.845, lineheight = 0.92, colour = "black") } +
  { if (nrow(hor)) geom_text(data = hor,
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      angle = 0, size = args$fontsize / 2.845, lineheight = 0.92, colour = "black") } +
  { if (nrow(out) && any(!out$boxed)) geom_text(data = out[!out$boxed, ],
      aes(x = x + GAP_X, y = y_lab, label = lab1), inherit.aes = FALSE,
      angle = 90, size = args$gutter_fontsize / 2.845, fontface = "bold",
      colour = "black") } +
  { if (nrow(out) && any(out$boxed)) geom_label(data = out[out$boxed, ],
      aes(x = x + GAP_X, y = y_lab, label = lab1), inherit.aes = FALSE,
      angle = 90, size = args$gutter_fontsize / 2.845, fontface = "bold",
      colour = "black", fill = "white", linewidth = 0,
      label.padding = unit(0.5, "mm"), alpha = 0.95) } +
  theme_void(base_size = args$fontsize, base_family = args$family) +
  theme(
    axis.text.x.top   = element_text(size = args$fontsize, colour = "black",
                                     lineheight = 0.92, margin = margin(b = 3)),
    legend.position   = "bottom",
    legend.text       = element_text(size = args$fontsize - 0.5),
    legend.key.size   = unit(3.4, "mm"),
    legend.margin     = margin(t = -3, b = 0),
    plot.margin       = margin(5, 5, 5, 5)
  )

cairo_pdf(paste0(args$out, ".pdf"), width = W, height = H, family = args$family)
print(p); dev.off()
ragg::agg_png(paste0(args$out, ".png"), width = W, height = H, units = "in",
              res = args$dpi, background = "white")
print(p); dev.off()

cat(sprintf("wrote %s.pdf and .png  (%.2f x %.2f mm, ratio W:H = 1:%.3f)\n",
            args$out, args$size_w, args$size_h, args$size_h / args$size_w))
cat(sprintf("labels: %d rotated inside, %d two-line inside, %d call-outs in the gutter (%d plain, %d boxed)\n",
            nrow(inn), nrow(hor), nrow(out),
            if (nrow(out)) sum(!out$boxed) else 0,
            if (nrow(out)) sum(out$boxed) else 0))
if (nrow(out)) cat(sprintf("  call-outs (bold %.1f pt, one line): %s\n",
                           args$gutter_fontsize, paste(out$lab1, collapse = ", ")))
cat(sprintf("lode = %.1f x %.1f pt (width computed from the panel)\n", lode_w_pt, lod$h_pt[1]))
cat(sprintf("min ribbon = %.2f pt -> floor %d genes; %d flows floored; display total %s (+%d = +%.2f%%)\n",
            args$min_ribbon_pt, floor_g, n_floor, fmt(tot_w), tot_w - tot,
            100 * (tot_w - tot) / tot))
for (i in seq_len(nrow(lod)))
  cat(sprintf("  %-30s h=%5.1fpt len=%5.1fpt  %s\n", lod$stratum[i], lod$h_pt[i], lod$len_pt[i],
              ifelse(lod$rot_fit[i], "rotated inside",
                     ifelse(lod$horiz_fit[i], "two-line inside", "call-out (gutter)"))))
