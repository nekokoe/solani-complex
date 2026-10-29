#!/usr/bin/env Rscript
# =====================================================================================
# plot_fig2_HIJ.R — Figure 2 panels H, I, J  (small square panels, for manual assembly)
#
#   H : nucleotide diversity (pi) boxplot, 9 groups, faceted by species; outliers shown;
#       pastel species backdrop matching Figure 2A (solani light blue, falciformis light
#       yellow, keratoplastica light purple)
#   I : dot matrix, pairwise significance of the panel-H pi distributions
#       size ~ -log10(q);  red intensity ~ significance;  non-significant = dark grey
#   J : dot matrix, pairwise Weir & Cockerham Fst between 11 groups
#       size ~ Fst;        red intensity ~ permutation z;  non-significant = dark grey
#
# For each panel two files are produced:
#   <name>.{pdf,svg,tif,png}            panel including its legend
#   <name>_nolegend.{pdf,svg,tif,png}   62 x 62 mm, no legend (for tight assembly)
# Half-matrix versions (right-angled bottom-right, labels bottom + right) are produced for
# the two symmetric matrices: fig2I_half, fig2J_half.
#
# Usage: Rscript plot_fig2_HIJ.R [datdir] [outdir]
# =====================================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(scales); library(grid); library(ggnewscale)
})

args   <- commandArgs(trailingOnly = TRUE)
datdir <- if (length(args) >= 1) args[1] else "data"
outdir <- if (length(args) >= 2) args[2] else "panels"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- palette ------------------------------------------------------------------------
COL_INK   <- "#1A1A1A"
COL_GRID  <- "#EDEDED"
COL_SIG   <- "#C0392B"          # significant dots (binary key, panel J)
COL_NS    <- "#4D4D4D"          # non-significant dots: dark grey
RED_LOW   <- "#F2C6BE"          # barely significant
RED_HIGH  <- "#7B1010"          # most significant

# Figure 2A species backdrop (pastel). EDIT THESE THREE VALUES to match Figure 2A exactly.
SP_BG <- c(solani = "#D6E4F0", falciformis = "#FFF3CC", keratoplastica = "#E6DAF2")
# box fill for the three source strata
SI_BREW <- c(all = "#4C72B0", env = "#55A868", clinical = "#C44E52")

PANEL_MM <- as.numeric(Sys.getenv("FIG_CANVAS", "62"))
W_IJ  <- as.numeric(Sys.getenv("FIG_W_IJ", as.character(PANEL_MM)))  # matrix canvas width
H_NL  <- as.numeric(Sys.getenv("FIG_H_NL", as.character(PANEL_MM)))  # matrix height, no legend
H_LG  <- as.numeric(Sys.getenv("FIG_H_LG", "76"))                    # matrix height, with legend
ONLY  <- Sys.getenv("FIG_ONLY", "")
want  <- function(x) ONLY == "" || x %in% strsplit(ONLY, ",")[[1]]
BASE_PT  <- 5.5

# Rotated axis labels overflow the space ggplot allocates for them (ggplot accounts only
# for the UNrotated text height), so the matrix panels need explicit margins. At 45 degrees
# the long group names ("kerato.-clinical") need ~10 mm of extra top+right room and still
# clipped at a 62 mm canvas; at 90 degrees the overflow is only ~2 mm, so 4 mm margins on
# all four sides keep every label inside. Verified by the layout sweep: all gaps >= 3.9 mm.
MX_L <- as.numeric(Sys.getenv("FIG_MX_L", "4"))
MX_R <- as.numeric(Sys.getenv("FIG_MX_R", "4"))
MX_T <- as.numeric(Sys.getenv("FIG_MX_T", "4"))
MX_B <- as.numeric(Sys.getenv("FIG_MX_B", "4"))
XANG <- as.numeric(Sys.getenv("FIG_XANG", "45"))

theme_panel <- function() {
  theme_minimal(base_size = BASE_PT, base_family = "sans") +
  theme(
    text              = element_text(colour = COL_INK, size = 5.5),
    axis.text         = element_text(colour = COL_INK, size = 5.2),
    axis.title        = element_text(size = 5.5),
    axis.ticks        = element_line(colour = COL_INK, linewidth = 0.25),
    axis.ticks.length = unit(0.6, "mm"),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_rect(colour = COL_INK, fill = NA, linewidth = 0.3),
    plot.margin       = margin(1.5, 1.5, 1.5, 1.5, "mm"),
    legend.position   = "none"
  )
}

# stack the guides vertically: side by side they were wider than the 62 mm canvas and got
# clipped at the image edges.
legend_bottom <- theme(
  legend.position    = "bottom",
  legend.box         = "vertical",
  legend.box.just    = "center",
  legend.key.size    = unit(2.8, "mm"),
  legend.text        = element_text(size = 5.0),
  legend.title       = element_text(size = 5.0),
  legend.margin      = margin(0, 2, 0, 2),
  legend.box.spacing = unit(0.8, "mm"))

## ---- export ------------------------------------------------------------------------
save_panel <- function(p, name, with_legend = TRUE, w_mm = PANEL_MM, h_mm = PANEL_MM) {
  if (!with_legend) name <- paste0(name, "_nolegend")
  q <- if (with_legend) p + legend_bottom else p + theme(legend.position = "none")
  grDevices::cairo_pdf(file.path(outdir, paste0(name, ".pdf")),
                       width = w_mm / 25.4, height = h_mm / 25.4, family = "sans")
  print(q); grDevices::dev.off()

  f_svg <- file.path(outdir, paste0(name, ".svg"))
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(f_svg, width = w_mm / 25.4, height = h_mm / 25.4)
  } else {
    grDevices::svg(f_svg, width = w_mm / 25.4, height = h_mm / 25.4)
  }
  print(q); grDevices::dev.off()

  ragg::agg_tiff(file.path(outdir, paste0(name, ".tif")),
                 width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 600,
                 background = "white")
  print(q); grDevices::dev.off()

  if (!with_legend && isTRUE(all.equal(w_mm, 62)) && isTRUE(all.equal(h_mm, 62))) {
    # canonical-size assembly panel: literal so the physical size is unambiguous
    ggsave(file.path(outdir, paste0(name, ".png")), q,
           width = 62, height = 62, units = "mm", dpi = 600, bg = "white")
  } else {
    ragg::agg_png(file.path(outdir, paste0(name, ".png")),
                  width = w_mm / 25.4, height = h_mm / 25.4, units = "in", res = 600,
                  background = "white")
    print(q); grDevices::dev.off()
  }
  cat("  wrote", name, "\n")
}

## ---- group ordering / labels --------------------------------------------------------
SP_ORDER  <- c("solani", "falciformis", "keratoplastica")
# facet strips: species epithet only (no genus initial), set in italic.
# label_parsed needs the plotmath source in the values, so the stripping variable is a
# dedicated column holding e.g. italic("solani").
SP_LABEL  <- c(solani = "solani", falciformis = "falciformis",
               keratoplastica = "keratoplastica")
sp_plot  <- function(x) sprintf('italic("%s")', x)
STR_LABEL <- c(all = "all", env = "env", clinical = "clinical")

pi_lv <- c("solani-all", "solani-env", "solani-clinical",
           "falciformis-all", "falciformis-env", "falciformis-clinical",
           "keratoplastica-all", "keratoplastica-env", "keratoplastica-clinical")

fj_lv <- c("falciformis-all", "falciformis-pop1", "falciformis-pop2", "falciformis-pop3",
           "solani-all", "solani-pop1", "solani-pop2", "solani-pop3",
           "keratoplastica-all", "keratoplastica-pop1", "keratoplastica-pop2")

# readable labels: abbreviate the species epithet only, never the stratum
# Matrix tick labels: species epithet written out in full, bold + italic; the stratum /
# population suffix stays upright so the species reads as a name.
#
# Implemented with base-R plotmath rather than ggtext markdown: the installed ggtext is 0.2.0
# and does not register an element_grob method for ggplot2 4.0.3, so element_markdown passes
# the markup through literally. Installing a newer ggtext is blocked here because its
# dependencies curl and xml2 fail to build (missing system libraries). plotmath renders a true
# bold-italic face (verified: NimbusSans-BoldItalic embedded in the PDF).
TICK_PT <- as.numeric(Sys.getenv("FIG_TICK_PT", "4.6"))
# Right-hand axis of the half matrices.
#
# ggplot2 4.0 defaults a RIGHT-positioned y axis to hjust = 1, i.e. the labels are
# right-aligned: their right edges line up and the shorter names float away from the axis
# (measured spread of left edges: 8.1 mm vs 0.08 mm with hjust = 0). Setting hjust = 0 is
# what makes the labels left-aligned against the axis.
#
# The residual gap between the panel and the text is 1.0 mm and is fixed by the ggplot2 4.0
# axis guide: axis.text.y margin (left or right, positive or negative), axis.ticks.length and
# axis.ticks.length.y all leave it unchanged (verified experimentally). TICK_LEN therefore
# only controls how long the tick marks are drawn.
TICK_LEN <- as.numeric(Sys.getenv("FIG_TICK_LEN", "0.4"))

SP_RE <- "^(keratoplastica|falciformis|solani)"
pm_lab <- function(x) {
  x <- as.character(x)
  sp <- regmatches(x, regexpr(SP_RE, x))
  rest <- sub(SP_RE, "", x)
  parse(text = sprintf('bolditalic("%s")*"%s"', sp, rest))
}
lab_pi <- pm_lab
lab_fj <- pm_lab

# score_col is the continuous significance measure; it MUST be carried through the
# mirror, otherwise the colour aesthetic resolves to NA and every dot renders grey50.
build_matrix <- function(dt, lv, size_col, q_col, score_col) {
  m  <- CJ(a = lv, b = lv, sorted = FALSE, unique = TRUE)
  mk <- function(A, B) dt[, .(a = A, b = B, size = get(size_col),
                              q = get(q_col), score = get(score_col))]
  d  <- rbind(mk(dt$group_a, dt$group_b), mk(dt$group_b, dt$group_a))
  m  <- merge(m, d, by = c("a", "b"), all.x = TRUE)
  m[, a := factor(a, levels = rev(lv))]
  m[, b := factor(b, levels = lv)]
  m[, sig := ifelse(!is.na(q) & q < 0.05, "sig", "ns")]
  m[sig == "ns", score := NA_real_]
  m
}

# =====================================================================================
# PANEL H — pi boxplot
# =====================================================================================
message("Panel H")
w <- fread(file.path(datdir, "panelH_pi_windows.tsv"))
w[, group_id := factor(group_id, levels = pi_lv)]
w[, species  := factor(species, levels = SP_ORDER)]
w[, stratum  := factor(stratum, levels = c("all", "env", "clinical"))]
w[, sp_plot  := factor(sp_plot(as.character(species)),
                       levels = sp_plot(SP_ORDER))]

# the backdrop layer must key off the SAME variable the panel is faceted by, otherwise the
# same rectangles are drawn in every facet and the per-species tints are lost
sp_keys <- sp_plot(SP_ORDER)
bg <- data.table(sp_plot = factor(sp_keys, levels = sp_keys))

# y upper limit: values above this are clipped on purpose (a few extreme windows), so the
# boxes and the bulk of the outliers stay readable.
y_max <- as.numeric(Sys.getenv("FIG_H_YMAX", "0.025"))

pH <- ggplot(w, aes(x = stratum, y = pi)) +
  geom_rect(data = bg, aes(fill = sp_plot), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.55, inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(SP_BG, sp_keys), guide = "none") +
  ggnewscale::new_scale_fill() +
  geom_boxplot(aes(fill = stratum), outlier.shape = 16, outlier.size = 0.30,
               outlier.alpha = 0.45, outlier.stroke = 0,
               linewidth = 0.22, width = 0.60, colour = COL_INK) +
  scale_fill_manual(values = SI_BREW, name = NULL) +
  scale_x_discrete(labels = STR_LABEL, drop = FALSE) +
  scale_y_continuous(labels = function(v) format(v, scientific = FALSE),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = NULL, y = expression("Nucleotide diversity ("*pi*")")) +
  coord_cartesian(ylim = c(0, y_max)) +
  facet_grid(. ~ sp_plot, scales = "free_x", space = "free_x",
             labeller = labeller(sp_plot = label_parsed)) +
  theme_panel() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.0),
        axis.text.y = element_text(size = 6.0),
        axis.title.y = element_text(size = 6.3),
        panel.spacing = unit(1.1, "mm"),
        strip.text = element_text(size = 5.9, colour = COL_INK),
        strip.background = element_rect(fill = "white", colour = NA),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.size = unit(2.8, "mm"),
        legend.text = element_text(size = 5.8),
        legend.title = element_text(size = 5.8),
        legend.box.spacing = unit(0.8, "mm"))
if (want("H")) {
  # height reduced to 80% of the previous 62 / 68 mm; width unchanged
  save_panel(pH, "fig2H", with_legend = TRUE,  w_mm = PANEL_MM, h_mm = 54)
  save_panel(pH, "fig2H", with_legend = FALSE, w_mm = PANEL_MM, h_mm = 50)
}

# =====================================================================================
# DOT MATRIX HELPERS
# =====================================================================================
# non-significant cells get a fixed-size dark grey dot; only significant cells are mapped
# to the size and colour channels, so both legends describe exactly the same population.
# mode = "gradient": colour is a continuous score (used by panel I, where -log10 q is
#   monotone in the effect and spans a usable range).
# mode = "binary": colour is the significance decision itself (used by panel J, where every
#   pair sits on the permutation floor so a continuous significance score carries no
#   information -- and a standardised score z is uncorrelated with Fst, Spearman -0.16).
dot_layers <- function(m, size_range = c(0.10, 3.00), ns_size = 0.75, mode = "gradient") {
  if (mode == "binary") {
    return(list(
      geom_point(data = m[!is.na(size)], aes(size = size, colour = sig)),
      scale_colour_manual(values = c(sig = COL_SIG, ns = COL_NS), name = NULL,
                          labels = c(sig = "P < 0.05", ns = "n.s.")),
      scale_size_continuous(range = size_range)
    ))
  }
  list(
    geom_point(data = m[sig == "ns" & !is.na(size)], size = ns_size, colour = COL_NS),
    ggnewscale::new_scale_colour(),
    geom_point(data = m[sig == "sig" & !is.na(size)], aes(size = size, colour = score)),
    scale_colour_gradient(low = RED_LOW, high = RED_HIGH),
    scale_size_continuous(range = size_range)
  )
}

# =====================================================================================
# PANEL I — pi significance dot matrix
# =====================================================================================
message("Panel I")
pi_pairs <- fread(file.path(datdir, "panelI_pi_pairwise.tsv"))
pi_pairs[, sigscore := -log10(pmax(q_block, 1e-12))]

# size = effect size (|delta pi| between the two groups), colour = strength of evidence
# (-log10 of the FDR-adjusted permutation p). Same split of channels as panel J, where size
# is Fst and colour is significance.
pi_pairs[, dpi := abs(as.numeric(delta_pi))]
mI <- build_matrix(pi_pairs[, .(group_a, group_b, size = dpi, q_block, score = sigscore)],
                   pi_lv, "size", "q_block", "score")

mk_full_matrix <- function(m, lv, xlab_fn, mode = "gradient") {
  ggplot(m, aes(x = b, y = a)) +
    geom_tile(fill = "white", colour = COL_GRID, linewidth = 0.15) +
    dot_layers(m, mode = mode) +
    scale_x_discrete(labels = xlab_fn, position = "top", limits = lv, drop = FALSE) +
    scale_y_discrete(labels = xlab_fn, limits = rev(lv), drop = FALSE) +
    labs(x = NULL, y = NULL) +
    coord_equal() +
    theme_panel() +
    theme(axis.text.x = ggtext::element_markdown(angle = XANG, hjust = 0, vjust = 0, size = TICK_PT),
          axis.text.y = ggtext::element_markdown(size = TICK_PT),
          plot.margin = margin(MX_T, MX_R, MX_B, MX_L, "mm"))
}

pI <- mk_full_matrix(mI, pi_lv, lab_pi)
pI <- pI +   scale_size_continuous(range = c(0.10, 3.00), breaks = c(0.001, 0.002, 0.003),
                        labels = c("0.001", "0.002", "0.003")) +
  labs(colour = expression(-log[10]~"BH-FDR adjusted "~italic(P)), size = expression(italic(Delta)*pi))
if (want("I")) {
  save_panel(pI, "fig2I", with_legend = TRUE,  w_mm = W_IJ, h_mm = H_LG)
  save_panel(pI, "fig2I", with_legend = FALSE, w_mm = W_IJ, h_mm = H_NL)
}

# ---- half matrix: right angle at the bottom-right, labels bottom + right -------------
half_cells <- function(m, lv) {
  x <- copy(m)
  x[, pa := match(as.character(a), lv)]
  x[, pb := match(as.character(b), lv)]
  x[pb >= pa]
}

mk_half_matrix <- function(m, lv, xlab_fn, mode = "gradient") {
  h <- half_cells(m, lv)
  h[, a := factor(a, levels = lv)]        # ascending -> diagonal runs bottom-left to top-right
  h[, b := factor(b, levels = lv)]
  g <- ggplot(h, aes(x = b, y = a)) +
    geom_tile(fill = "white", colour = COL_GRID, linewidth = 0.15) +
    dot_layers(h, mode = mode) +
    scale_x_discrete(labels = xlab_fn, position = "bottom", limits = lv, drop = FALSE) +
    scale_y_discrete(labels = xlab_fn, position = "right", limits = lv, drop = FALSE) +
    labs(x = NULL, y = NULL) +
    coord_equal() +
    theme_panel() +
    theme(axis.text.x = element_text(angle = XANG, hjust = 1, vjust = 1, size = TICK_PT),
          # right-hand axis: left-aligned so every label starts at the axis, with the default
          # tick length and text margin trimmed so the text hugs the axis instead of floating
          # ~1.7 mm away from it
          axis.text.y = element_text(size = TICK_PT, hjust = 0),
          axis.ticks.length.y = unit(TICK_LEN, "mm"),
          plot.margin = margin(MX_T, MX_R, MX_B, MX_L, "mm"))
  g
}

pI_half <- mk_half_matrix(mI, pi_lv, lab_pi) +
    scale_size_continuous(range = c(0.10, 3.00), breaks = c(0.001, 0.002, 0.003),
                        labels = c("0.001", "0.002", "0.003")) +
  labs(colour = expression(-log[10]~"BH-FDR adjusted "~italic(P)), size = expression(italic(Delta)*pi))
if (want("Ihalf")) {
  save_panel(pI_half, "fig2I_half", with_legend = TRUE,  w_mm = W_IJ, h_mm = H_LG)
  save_panel(pI_half, "fig2I_half", with_legend = FALSE, w_mm = W_IJ, h_mm = H_NL)
}

# =====================================================================================
# PANEL J — Fst dot matrix
# =====================================================================================
message("Panel J")
fstf <- file.path(datdir, "panelJ_fst_matrix_v2.tsv")
if (!file.exists(fstf)) fstf <- file.path(datdir, "panelJ_fst_matrix.tsv")
fp <- fread(fstf)

# ---- colour = p-value gradient ------------------------------------------------------
# The exact permutation p-value is 0.001 for all 55 pairs (1000 permutations put the floor
# at 1/(1000+1)), so an empirical p-value gradient is flat by construction. The p-value is
# therefore taken from the normal approximation to the permutation null of each pair:
#   P = 2 * Phi(-|z|),  z = (Fst_obs - mean(null)) / sd(null)
# Colour then encodes strength of evidence and size encodes effect size (Fst), which are
# deliberately different quantities.
fp[, z := as.numeric(z)]
fp[, plog := -log10(2 * pnorm(-abs(z)))]
fp[!is.finite(plog), plog := NA_real_]
PCAP <- as.numeric(Sys.getenv("FIG_PCAP", "300"))
fp[plog > PCAP, plog := PCAP]
message(sprintf("  panel J colour: -log10 P, range %.1f to %.1f (capped at %.0f)",
                min(fp$plog, na.rm = TRUE), max(fp$plog, na.rm = TRUE), PCAP))

mJ <- build_matrix(fp[, .(group_a, group_b, size = fst, p_perm, score = plog)],
                   fj_lv, "size", "p_perm", "score")

pJ <- mk_full_matrix(mJ, fj_lv, lab_fj, mode = "gradient") +
  scale_colour_gradient(low = RED_LOW, high = RED_HIGH, breaks = c(20, 100, 200, 300),
                        labels = c("20", "100", "200", "\u2265300")) +
  labs(colour = expression(-log[10]~italic(P)), size = expression(italic(F)[ST]))
if (want("J")) {
  save_panel(pJ, "fig2J", with_legend = TRUE,  w_mm = W_IJ, h_mm = H_LG)
  save_panel(pJ, "fig2J", with_legend = FALSE, w_mm = W_IJ, h_mm = H_NL)
}

pJ_half <- mk_half_matrix(mJ, fj_lv, lab_fj, mode = "gradient") +
  scale_colour_gradient(low = RED_LOW, high = RED_HIGH, breaks = c(20, 100, 200, 300),
                        labels = c("20", "100", "200", "\u2265300")) +
  labs(colour = expression(-log[10]~italic(P)), size = expression(italic(F)[ST]))
if (nzchar(Sys.getenv("FIG_DIAG"))) {
  g <- ggplotGrob(pJ_half)
  lay <- g$layout
  wmm <- as.numeric(grid::convertWidth(g$widths, "mm"))
  cat("\n--- half-matrix layout (mm) ---\n")
  for (i in seq_len(nrow(lay))) {
    nm <- as.character(lay$name[i]); l <- lay$l[i]; r <- lay$r[i]
    x0 <- if (l > 1) sum(wmm[1:(l-1)]) else 0
    x1 <- sum(wmm[1:r])
    if (nm %in% c("panel", "axis-r", "axis-l", "spacer")) cat(sprintf("  %-10s cols %d-%d  x %.2f -> %.2f mm (w %.2f)\n", nm, l, r, x0, x1, x1-x0))
  }
  cat(sprintf("  tick length set to %.2f mm\n", TICK_LEN))
}

if (want("Jhalf")) {
  save_panel(pJ_half, "fig2J_half", with_legend = TRUE,  w_mm = W_IJ, h_mm = H_LG)
  save_panel(pJ_half, "fig2J_half", with_legend = FALSE, w_mm = W_IJ, h_mm = H_NL)
}

# =====================================================================================
# standalone legend strips (for flexible assembly)
# =====================================================================================
mk_legend <- function(size_vals, size_labs, size_title, colour_title, col_breaks, col_labs) {
  n  <- length(size_vals)
  df <- data.frame(x = seq_len(n + 1),
                   y = c(rep(1, n), 1.7),
                   size = c(size_vals, min(size_vals)),
                   score = c(seq(min(col_breaks), max(col_breaks), length.out = n), 0),
                   key = c(rep("sig", n), "ns"), stringsAsFactors = FALSE)
  df_sig <- df[df$key == "sig", ]
  df_ns  <- df[df$key == "ns", ]
  ggplot() +
    geom_point(data = df_sig, aes(x = x, y = y, size = size, colour = score)) +
    scale_colour_gradient(low = RED_LOW, high = RED_HIGH, breaks = col_breaks,
                          labels = col_labs, name = colour_title) +
    ggnewscale::new_scale_colour() +
    geom_point(data = df_ns, aes(x = x, y = y), size = 1.4, colour = COL_NS) +
    scale_size_continuous(range = c(0.10, 3.00), breaks = size_vals, name = size_title,
                          labels = function(b) {
                            i <- match(b, size_vals)
                            ifelse(is.na(i), format(b, trim = TRUE), size_labs[i])
                          }) +
    theme_void(base_size = BASE_PT) +
    legend_bottom +
    theme(legend.key.size = unit(3, "mm"))
}
mk_legend_binary <- function(size_vals, size_labs, size_title) {
  n  <- length(size_vals)
  df <- data.frame(x = c(seq_len(n), 1, 1.6), y = 1,
                   size = c(size_vals, min(size_vals), min(size_vals)),
                   sig  = c(rep("sig", n), "ns", "sig"), stringsAsFactors = FALSE)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(size = size, colour = sig)) +
    scale_colour_manual(values = c(sig = COL_SIG, ns = COL_NS), name = NULL,
                        breaks = c("sig", "ns"), labels = c("P < 0.05", "n.s.")) +
    scale_size_continuous(range = c(0.10, 3.00), breaks = size_vals, name = size_title,
                          labels = function(b) {
                            i <- match(b, size_vals)
                            ifelse(is.na(i), format(b, trim = TRUE), size_labs[i])
                          }) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 1.4)),
           size   = guide_legend(order = 2)) +
    theme_void(base_size = BASE_PT) +
    legend_bottom +
    theme(legend.key.size = unit(3, "mm"))
}

ggsave(file.path(outdir, "fig2I_legend.pdf"),
       mk_legend(c(0.001, 0.002, 0.003), c("0.001", "0.002", "0.003"),
                 expression(italic(Delta)*pi), expression(-log[10]~"BH-FDR adjusted "~italic(P)),
                 c(1.3, 2.5, 4.0), c("1.3", "2.5", "4.0")),
       width = 62, height = 19, units = "mm", device = grDevices::cairo_pdf)
ggsave(file.path(outdir, "fig2J_legend.pdf"),
       mk_legend(c(0.1, 0.3, 0.5, 0.8), c("0.1", "0.3", "0.5", "0.8"),
                 expression(italic(F)[ST]), expression(-log[10]~italic(P)),
                 c(20, 100, 200), c("20", "100", "200")),
       width = 62, height = 19, units = "mm", device = grDevices::cairo_pdf)

message("done")
