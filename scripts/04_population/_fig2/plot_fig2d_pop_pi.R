#!/usr/bin/env Rscript
# Panel d of Supplementary Fig. S2 — nucleotide diversity (pi) per ADMIXTURE population.
#
# Same visual language as panel c (fig2H): a boxplot of per-window pi, faceted by species,
# with the species backdrop tint and a per-stratum box fill. The only difference from panel
# c is the grouping: here the x-axis carries the ADMIXTURE populations (all / pop1 / pop2 /
# pop3) — the same 11 groups as panel a (the Fst matrix, main Fig. 2G) — instead of the
# isolation-source strata (all / env / clinical).
#
# Run from 04.population/_fig2:
#   Rscript plot_fig2d_pop_pi.R
# Output: panels/fig2d_pop_pi.{pdf,png}
suppressMessages({library(data.table); library(ggplot2); library(ggnewscale)})

COL_INK <- "#1A1A1A"
SP_BG <- c(solani = "#D6E4F0", falciformis = "#FFF3CC", keratoplastica = "#E6DAF2")
SP_ORDER <- c("solani", "falciformis", "keratoplastica")   # same facet order as panel c
sp_plot <- function(x) sprintf('italic("%s")', x)
STRATUM <- c(all = "all", pop1 = "pop1", pop2 = "pop2", pop3 = "pop3")
# all keeps panel c's blue; pop1/pop2/pop3 carry the Fig. 2D-F population colours
POP_BREW <- c(all = "#4C72B0", pop1 = "#DB5243", pop2 = "#4BB2C8", pop3 = "#2D9D89")

w <- fread("data/panelJ_pop_pi_windows.tsv")
w[, species := factor(species, levels = SP_ORDER)]
w[, stratum := factor(sub("^[a-z]+-", "", group), levels = names(STRATUM))]
w[, sp_plot := factor(sp_plot(as.character(species)), levels = sp_plot(SP_ORDER))]

sp_keys <- sp_plot(SP_ORDER)
bg <- data.table(sp_plot = factor(sp_keys, levels = sp_keys))

y_max <- 0.025
p <- ggplot(w, aes(x = stratum, y = pi)) +
  geom_rect(data = bg, aes(fill = sp_plot), xmin = -Inf, xmax = Inf,
            ymin = -Inf, ymax = Inf, alpha = 0.55, inherit.aes = FALSE) +
  scale_fill_manual(values = setNames(SP_BG, sp_keys), guide = "none") +
  new_scale_fill() +
  geom_boxplot(aes(fill = stratum), outlier.shape = 16, outlier.size = 0.30,
               outlier.alpha = 0.45, outlier.stroke = 0,
               linewidth = 0.22, width = 0.80, colour = COL_INK) +
  scale_fill_manual(values = POP_BREW, name = NULL) +
  scale_x_discrete(labels = STRATUM, drop = TRUE) +
  scale_y_continuous(labels = function(v) format(v, scientific = FALSE),
                     expand = expansion(mult = c(0.02, 0.03))) +
  labs(x = NULL, y = expression("Nucleotide diversity ("*pi*")")) +
  coord_cartesian(ylim = c(0, y_max)) +
  facet_grid(. ~ sp_plot, scales = "free_x", space = "free_x",
             labeller = labeller(sp_plot = label_parsed)) +
  theme_minimal(base_size = 5.5, base_family = "sans") +
  theme(
    text = element_text(colour = COL_INK, size = 5.5),
    axis.text = element_text(colour = COL_INK, size = 5.2),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.2),
    axis.title = element_text(size = 5.5),
    axis.ticks = element_line(colour = COL_INK, linewidth = 0.25),
    axis.ticks.length = unit(0.6, "mm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = COL_INK, fill = NA, linewidth = 0.3),
    panel.spacing = unit(1.1, "mm"),
    strip.text = element_text(size = 5.9, colour = COL_INK),
    strip.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"),
    legend.position = "bottom", legend.direction = "horizontal",
    legend.key.size = unit(2.8, "mm"),
    legend.text = element_text(size = 5.8),
    legend.title = element_text(size = 5.8),
    legend.box.spacing = unit(0.8, "mm")
  )

ggsave("panels/fig2d_pop_pi.pdf", p, width = 61.7, height = 54, units = "mm",
       device = grDevices::cairo_pdf)
ggsave("panels/fig2d_pop_pi.png", p, width = 61.7, height = 54, units = "mm", dpi = 600)
cat("wrote panels/fig2d_pop_pi.{pdf,png}\n")
