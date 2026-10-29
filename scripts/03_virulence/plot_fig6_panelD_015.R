#!/usr/bin/env Rscript
# Figure 6 panel D — per-strain temperature preference, |Δ mortality| > 15 pp threshold.
#
# Replaces the mean ± 1 s.d. band of the collaborator-supplied composite panel with an
# explicit, transportable threshold so that the panel matches the main text, which names
# 5 x 28 °C-biased and 2 x 37 °C-biased isolates.
#
#   Δ mortality = mortality(37 °C) − mortality(28 °C), averaged over the three inoculum doses
#               = −temp_effect  (phenotype_metrics.tsv stores the opposite sign convention)
#
# backend = r (nature-figure): cairo_pdf + Arial/NimbusSans, Arial-safe glyphs only.

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

ROOT <- "/share/pasteur/gaosh/Fusarium/phase2/solaniPathogenicityPaper"
OUT  <- file.path(ROOT, "00.manuscript/figures_extracted")
THR  <- 15            # percentage points

C37 <- "#C92626"; C28 <- "#4D8AC1"; CNEU <- "#BEBEBE"
CGRID <- "#EAEAEA"

pheno <- read.delim(file.path(ROOT, "07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv"),
                    check.names = FALSE, row.names = 1)
pheno$Isolate <- rownames(pheno)
pheno$delta <- (pheno$mort_37_mean - pheno$mort_28_mean) * 100

pheno$pref <- ifelse(pheno$delta >  THR, "37",
              ifelse(pheno$delta < -THR, "28", "flat"))
pheno <- pheno[order(pheno$delta), ]
pheno$Isolate <- factor(pheno$Isolate, levels = pheno$Isolate)

stopifnot(nrow(pheno) == 21)
cat("threshold ±", THR, "pp\n", sep = "")
print(table(pheno$pref)[c("37", "28", "flat")])
cat("\nclassified strains:\n")
print(pheno[pheno$pref != "flat", c("Isolate", "delta", "pref")], row.names = FALSE)
cat("\nmargin to the cut, in pp (inside = classified, outside = nearest unclassified):\n")
d <- pheno$delta
inside_up   <- d[d >  THR]; outside_up   <- d[d <= THR]
inside_dn   <- d[d < -THR]; outside_dn   <- d[d >= -THR]
cat("  37 °C side: closest inside", min(inside_up), " | closest outside", max(outside_up), "\n")
cat("  28 °C side: closest inside", max(inside_dn), " | closest outside", min(outside_dn), "\n")

write.table(pheno[, c("Isolate", "delta", "pref", "mort_28_mean", "mort_37_mean", "temp_effect")],
            file.path(OUT, "figure6-panelD-015-classification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lab_dn <- pheno[pheno$pref == "28", ]
lab_up <- pheno[pheno$pref == "37", ]

p <- ggplot(pheno, aes(x = Isolate, y = delta)) +
  geom_hline(yintercept =  THR, colour = "#FF8989", linetype = "dotted", linewidth = 0.45) +
  geom_hline(yintercept = -THR, colour = "#8A8AFF", linetype = "dotted", linewidth = 0.45) +
  geom_col(aes(fill = pref), width = 0.8, colour = "#333333", linewidth = 0.28) +
  # 28 °C-biased: label right of the bar end, repelled vertically (five bars crowd -20…-15 pp)
  geom_text_repel(data = lab_dn, aes(label = Isolate),
                  size = 1.22, family = "Arial", colour = "black",
                  nudge_x = 0.5, hjust = 0, direction = "y",
                  segment.size = 0.16, segment.colour = "grey40",
                  min.segment.length = 0, box.padding = 0.05,
                  point.padding = 0, max.overlaps = Inf, seed = 11) +
  # 37 °C-biased: label above the bar top, repelled horizontally (they sit at the right edge)
  geom_text_repel(data = lab_up, aes(label = Isolate),
                  size = 1.22, family = "Arial", colour = "black",
                  nudge_y = 1.6, vjust = 0, direction = "x",
                  segment.size = 0.16, segment.colour = "grey40",
                  min.segment.length = 0, box.padding = 0.05,
                  point.padding = 0, max.overlaps = Inf, seed = 11) +
  scale_fill_manual(values = c("37" = C37, "flat" = CNEU, "28" = C28),
                    breaks = c("37", "flat", "28"),
                    labels = c("37" = "37°C-biased", "flat" = "Neutral", "28" = "28°C-biased"),
                    name = "Temperature preference") +
  scale_y_continuous(limits = c(-50, 50), breaks = seq(-50, 50, 25),
                     labels = function(b) paste0(b, "%"),
                     expand = expansion(mult = 0.05)) +
  scale_x_discrete(expand = expansion(add = 0.6)) +
  labs(x = "Strain", y = "Δ mortality (37°C − 28°C)") +
  theme_bw(base_size = 5, base_family = "Arial") +
  theme(
    panel.grid.major   = element_line(colour = CGRID, linewidth = 0.28),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(colour = "black", linewidth = 0.35, fill = NA),
    axis.ticks         = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length  = unit(0.55, "mm"),
    axis.text          = element_text(size = 4.5, colour = "black"),
    axis.text.x        = element_blank(),
    axis.title         = element_text(size = 5.2, colour = "black"),
    axis.title.x       = element_text(margin = margin(t = 1.1)),
    axis.title.y       = element_text(margin = margin(r = 1.1)),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 4.6, margin = margin(r = 1.5)),
    legend.text        = element_text(size = 4.4),
    legend.key         = element_blank(),
    legend.key.size    = unit(2.6, "mm"),
    legend.spacing.x   = unit(0.9, "mm"),
    legend.background  = element_blank(),
    legend.margin      = margin(t = 1),
    legend.box.spacing = unit(1.2, "mm"),
    plot.margin        = margin(1.4, 1.8, 2.0, 1.4, unit = "mm")
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE,
                             override.aes = list(linewidth = 0.28, colour = "#333333")))

W <- 74; H <- 44   # mm, plot box + axis titles + legend, panel letter excluded

cairo_pdf(file.path(OUT, "figure6-panelD-015.pdf"), width = W / 25.4, height = H / 25.4,
          family = "Arial")
print(p)
invisible(dev.off())

png(file.path(OUT, "figure6-panelD-015.png"), width = W / 25.4, height = H / 25.4,
    units = "in", res = 600, type = "cairo")
print(p)
invisible(dev.off())

cat("\nwrote figure6-panelD-015.{pdf,png} (", W, "x", H, " mm )\n", sep = "")
