#!/usr/bin/env Rscript
# Figure 6 panel D — per-strain temperature preference.
#
# Classification is the per-strain Cochran-Mantel-Haenszel test stratified by inoculum
# dose (scripts/cmh_temperature_preference.py, P < 0.05), NOT a cutoff on the continuous
# difference in mean mortality. A test criterion has no cut point to draw, so the panel
# carries no threshold lines: only the four significant strains are coloured and labelled.
#
#   Δ mortality = mortality(37 °C) − mortality(28 °C), mean over the three inoculum doses
#               = −temp_effect  (phenotype_metrics.tsv stores the opposite sign convention)
#
# backend = r (nature-figure): cairo_pdf + Arial/NimbusSans.

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

ROOT <- "/share/pasteur/gaosh/Fusarium/phase2/solaniPathogenicityPaper"
OUT  <- file.path(ROOT, "00.manuscript/figures_extracted")

C37 <- "#C92626"; C28 <- "#4D8AC1"; CNEU <- "#BEBEBE"
CGRID <- "#EAEAEA"

pheno <- read.delim(file.path(ROOT, "07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv"),
                    check.names = FALSE, row.names = 1)
bias <- read.delim(file.path(ROOT, "07.lethal_groups/phenotype_analysis/temperature_bias_cmh.tsv"),
                   stringsAsFactors = FALSE)

pheno$Isolate <- rownames(pheno)
pheno$delta <- (pheno$mort_37_mean - pheno$mort_28_mean) * 100
pheno$pref <- bias$class[match(pheno$Isolate, bias$strain)]
pheno$P <- bias$P[match(pheno$Isolate, bias$strain)]
pheno$q <- bias$q_BH[match(pheno$Isolate, bias$strain)]
pheno$pref <- factor(ifelse(pheno$pref == "37C_biased", "37",
                     ifelse(pheno$pref == "28C_biased", "28", "flat")),
                     levels = c("37", "flat", "28"))
pheno <- pheno[order(pheno$delta), ]
pheno$Isolate <- factor(pheno$Isolate, levels = pheno$Isolate)

stopifnot(nrow(pheno) == 21)
cat("classification, CMH P < 0.05:\n")
print(table(pheno$pref)[c("37", "28", "flat")])
cat("\nclassified strains:\n")
print(pheno[pheno$pref != "flat", c("Isolate", "delta", "P", "q")], row.names = FALSE)
cat("\nall P values:\n")
print(pheno[order(pheno$P), c("Isolate", "delta", "P", "q")], row.names = FALSE)

write.table(pheno[, c("Isolate", "delta", "pref", "P", "q", "mort_28_mean", "mort_37_mean")],
            file.path(OUT, "figure6-panelD-cmh-classification.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lab_dn <- pheno[pheno$pref == "28", ]
lab_up <- pheno[pheno$pref == "37", ]

p <- ggplot(pheno, aes(x = Isolate, y = delta)) +
  geom_col(aes(fill = pref), width = 0.8, colour = "#333333", linewidth = 0.28) +
  # 28 °C-biased: label right of the bar end
  geom_text_repel(data = lab_dn, aes(label = Isolate),
                  size = 1.22, family = "Arial", colour = "black",
                  nudge_x = 0.5, hjust = 0, direction = "y",
                  segment.size = 0.16, segment.colour = "grey40",
                  min.segment.length = 0, box.padding = 0.05,
                  point.padding = 0, max.overlaps = Inf, seed = 11) +
  # 37 °C-biased: label above the bar tops, which sit at the right edge
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

cairo_pdf(file.path(OUT, "figure6-panelD-cmh.pdf"), width = W / 25.4, height = H / 25.4,
          family = "Arial")
print(p)
invisible(dev.off())

png(file.path(OUT, "figure6-panelD-cmh.png"), width = W / 25.4, height = H / 25.4,
    units = "in", res = 600, type = "cairo")
print(p)
invisible(dev.off())

cat("\nwrote figure6-panelD-cmh.{pdf,png} (", W, "x", H, " mm )\n", sep = "")
