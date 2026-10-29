#!/usr/bin/env Rscript
# s20ef_panel.R — S20 panels e and f as standalone panels.
#
#   e = BMU02483 shared k-mers along the clinical accessory chromosomes
#       (sliding 10 kb windows, x = chromosome position 0-100 %)
#   f = whole-AC k-mer sharing per chromosome, coloured by the k-mer class
#
# Both panels carry the same chromosome order and the same y axis, so they can be
# placed side by side and line up. Panel e is drawn 2:1, panel f square.
#
# Sources: 06.accessory-chrs/BMU02483_gene_content/sliding/sliding_2483_prop.tsv
#          06.accessory-chrs/kmer_analysis/kmer_matrix.ACx211.tsv
#          06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv

suppressMessages({library(ggplot2)})

S20_RED <- "#C44E52"; S20_BLUE <- "#3B7DD8"

s20e_data <- function(path) {
  sl <- read.delim(path, stringsAsFactors = FALSE)
  ac_order <- sort(unique(sl$AC), decreasing = TRUE)
  sl$AC <- factor(sl$AC, levels = ac_order)
  ac_len <- tapply(sl$end, sl$AC, max)
  sl$xmin <- sl$start / ac_len[as.character(sl$AC)]
  sl$xmax <- sl$end / ac_len[as.character(sl$AC)]
  sl
}

s20f_data <- function(matrix_path, v3_path, sl) {
  mat <- read.delim(matrix_path, row.names = 1, check.names = FALSE)
  v3 <- read.delim(v3_path, stringsAsFactors = FALSE)
  ac_order <- levels(sl$AC)
  b <- data.frame(AC = rownames(mat), share = mat[, "BMU02483"])
  b <- b[b$AC %in% ac_order, ]
  b$AC <- factor(b$AC, levels = ac_order)
  b$type <- ifelse(b$AC %in% v3$AC[v3$kmer_type == "type1_stable"],
                   "Block", "Mosaic")
  b$type <- factor(b$type, levels = c("Block", "Mosaic"))
  b
}

s20e_panel <- function(sl, family = "Arial", fontsize = 9,
                       title = "BMU02483 shared k-mers along clinical ACs") {
  ggplot(sl, aes(xmin = xmin, xmax = xmax, y = AC, fill = prop)) +
    geom_rect(height = 0.85) +
    scale_fill_gradientn(colours = c("#4575B4", "#E0F3F8", "#FFFFBF",
                                     "#FDAE61", "#D73027"),
                         values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                         limits = c(0, 1), name = "Shared\nk-mer") +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0", "25", "50", "75", "100%"),
                       expand = c(0, 0)) +
    labs(x = "Chromosome position (%)", y = NULL, title = title) +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(size = fontsize, colour = "black"),
      axis.text.y = element_text(size = fontsize * 0.8, colour = "black"),
      axis.title.x = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold", hjust = 0),
      legend.position = "right",
      legend.key.height = unit(14, "mm"), legend.key.width = unit(3.5, "mm"),
      legend.title = element_text(size = fontsize - 1),
      legend.text = element_text(size = fontsize - 1.5),
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 4, 4)
    )
}

s20f_panel <- function(b, family = "Arial", fontsize = 9,
                       title = "Overall",
                       ylabels = TRUE, legend_position = "bottom") {
  p <- ggplot(b, aes(x = share, y = AC, fill = type)) +
    geom_col(width = 0.72, alpha = 0.95) +
    scale_fill_manual(values = c("Block" = S20_BLUE, "Mosaic" = S20_RED), name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Whole-AC k-mer sharing", y = NULL, title = title) +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(size = fontsize, colour = "black"),
      axis.title.x = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold", hjust = 0),
      legend.position = legend_position,
      legend.text = element_text(size = fontsize),
      legend.key.size = unit(4, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 3, 4)
    )
  if (!ylabels) p <- p + theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank())
  p
}
