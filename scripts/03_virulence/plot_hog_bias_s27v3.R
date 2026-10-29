#!/usr/bin/env Rscript
# S27 v3 (R, ggplot2 + patchwork) — temperature-bias strain-specific HOGs.
# Two stacked heatmaps (37C toolkit warm / 28C toolkit cool) + strain-group
# strip on top + two-line large row labels + per-cell copy numbers.

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

theme_set(theme_minimal(base_size = 13))
B37 <- c("BMU01287", "BMU02031")
B28 <- c("BMU03976", "BMU12685", "BMU02028", "BMU01936")
MID <- c("BMU12684", "BMU02407", "BMU02479", "BMU01973", "BMU01935",
         "BMU12686", "BMU12683", "BMU10147", "BMU02745")
STRAINS <- c(B37, MID, B28)
strain_group <- data.frame(
  strain = STRAINS,
  group  = c(rep("37C", 2), rep("mid", 9), rep("28C", 4)),
  stringsAsFactors = FALSE
)
grp_col <- c("37C" = "#C44E52", "mid" = "#BBBBBB", "28C" = "#3B7DD8")

d <- read.delim("tmp/s27_long.tsv", stringsAsFactors = FALSE)
m <- read.delim("tmp/s27_meta.tsv", stringsAsFactors = FALSE)
d <- d %>% left_join(m, by = c("HOG", "block"))

d$strain <- factor(d$strain, levels = STRAINS)
d$label <- ifelse(d$chit, paste0("\u25cf ", d$HOG, "\n", d$tag),
                  paste0(d$HOG, "\n", d$tag))
d <- d %>% arrange(block, chit, inf) %>%
  mutate(HOG = factor(HOG, levels = unique(HOG)))
d <- d %>% group_by(HOG) %>% mutate(ypos = cur_group_id()) %>% ungroup()
d <- d %>% mutate(ypos = ifelse(block == "37C", ypos, ypos + 46))

# per-cell colour by block cmap
cmap37 <- colorRampPalette(c("#FFF7EC", "#FDBE85", "#E6550D", "#A63603"))(6)
cmap28 <- colorRampPalette(c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"))(6)
d$fill <- ifelse(d$block == "37C",
                 cmap37[pmin(d$copy + 1, 6)],
                 cmap28[pmin(d$copy + 1, 6)])

# strain group strip
strip <- strain_group %>%
  mutate(x = seq_along(strain), y = 1) %>%
  ggplot(aes(x, y, fill = group)) +
  geom_tile(width = 1, height = 1) +
  scale_fill_manual(values = grp_col, name = NULL) +
  scale_x_continuous(breaks = seq_along(STRAINS), labels = STRAINS,
                     limits = c(0.5, length(STRAINS) + 0.5), expand = c(0, 0)) +
  theme_void(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 11, face = "bold",
                                   colour = grp_col[strain_group$group]),
        legend.position = "bottom",
        plot.margin = margin(0, 108, -4, 0))

# heatmap
p_heat <- ggplot(d, aes(strain, ypos)) +
  geom_tile(aes(fill = fill), width = 0.96, height = 0.96) +
  scale_fill_identity() +
  geom_text(aes(label = ifelse(copy > 0, copy, "")), size = 3.8,
            colour = "#1a1a1a", fontface = "bold") +
  geom_hline(yintercept = 46.5, linewidth = 1.1, colour = "black") +
  scale_y_continuous(
    breaks = d$ypos, labels = d$label, position = "right",
    limits = c(-0.5, 56.8), expand = c(0, 0)) +
  scale_x_discrete(breaks = STRAINS, labels = rep("", length(STRAINS)),
                   expand = c(0, 0)) +
  theme(axis.text.y = element_text(size = 8, hjust = 0, vjust = 0.5,
                                   lineheight = 0.85, face = "plain"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) +
  annotate("text", x = 8, y = 56.4, label = "37 \u00b0C-biased toolkit (n = 46)",
           size = 5.8, hjust = 0.5, colour = "#A63603", fontface = "bold") +
  annotate("text", x = 8, y = 49, label = "28 \u00b0C-biased toolkit (n = 10)",
           size = 5.8, hjust = 0.5, colour = "#08306B", fontface = "bold")

# copy-number legend (single tile strip)
leg_df <- data.frame(copy = 0:5)
p_leg <- ggplot(leg_df, aes(copy, 1)) +
  geom_tile(aes(fill = copy), width = 1, height = 1) +
  scale_fill_gradientn(colours = cmap37, limits = c(0, 5)) +
  scale_x_continuous(breaks = 0:5, labels = 0:5) +
  theme_void(base_size = 11) +
  theme(axis.text.x = element_text(size = 10),
        legend.position = "none",
        plot.margin = margin(2, 0, 0, 0)) +
  coord_cartesian(xlim = c(-0.5, 5.5))

layout <- strip / p_heat / p_leg + plot_layout(heights = c(0.35, 4, 0.12))
ggsave("07.lethal_groups/genomewide_assoc/figure_hog_bias_S27v3.pdf",
       layout, width = 9.5, height = 14, dpi = 300)
ggsave("07.lethal_groups/genomewide_assoc/figure_hog_bias_S27v3.png",
       layout, width = 9.5, height = 14, dpi = 300)
cat("S27 v3 saved\n")
