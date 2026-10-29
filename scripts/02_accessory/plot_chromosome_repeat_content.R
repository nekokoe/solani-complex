#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

option_list <- list(
  make_option(c("-r", "--repeat_file"), type="character", default=NULL,
              help="Repeat content file (tsv)", metavar="character"),
  make_option(c("--contig_order"), type="character", default=NULL,
              help="File with list of contigs in desired order (tsv, single column)", metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default="chromosome_repeat_content",
              help="Output file prefix", metavar="character"),
  make_option(c("--width"), type="numeric", default=6,
              help="Plot width in inches", metavar="numeric"),
  make_option(c("--height"), type="numeric", default=3,
              help="Plot height in inches", metavar="numeric"),
  make_option(c("--title"), type="character", default="Repeat Content Distribution by Chromosome",
              help="Plot title", metavar="character"),
  make_option(c("--font_title"), type="numeric", default=10,
              help="Title font size", metavar="numeric"),
  make_option(c("--font_axis"), type="numeric", default=7,
              help="Axis text font size", metavar="numeric"),
  make_option(c("--font_legend_title"), type="numeric", default=8,
              help="Legend title font size", metavar="numeric"),
  make_option(c("--font_legend_text"), type="numeric", default=7,
              help="Legend text font size", metavar="numeric"),
  make_option(c("--font_label"), type="numeric", default=1.5,
              help="Bar label font size", metavar="numeric"),
  make_option(c("--bar_width"), type="numeric", default=0.8,
              help="Bar width (default: 0.8)", metavar="numeric"),
  make_option(c("--min_percentage"), type="numeric", default=0.05,
              help="Minimum percentage to keep (default: 0.05)", metavar="numeric")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

if (is.null(args$repeat_file)) {
  print_help(parser)
  stop("Repeat file must be provided", call.=FALSE)
}

message("Reading repeat content file...")
first_line <- readLines(args$repeat_file, n = 1)
has_header <- grepl("^ContigID", first_line)

if (has_header) {
  repeat_data <- read.table(args$repeat_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  repeat_data <- read.table(args$repeat_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(repeat_data) <- c("ContigID", "RepeatType", "Count", "Percentage")
}

message(paste("Total rows:", nrow(repeat_data)))
message(paste("Unique contigs:", length(unique(repeat_data$ContigID))))
message(paste("Unique repeat types:", length(unique(repeat_data$RepeatType))))

filtered_before <- nrow(repeat_data)
repeat_data <- repeat_data %>% filter(Percentage >= args$min_percentage)
filtered_after <- nrow(repeat_data)
message(paste("Filtered", filtered_before - filtered_after, "rows with Percentage <", args$min_percentage))

repeat_types <- sort(unique(repeat_data$RepeatType))
n_types <- length(repeat_types)
spectral_colors <- colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(n_types)
names(spectral_colors) <- repeat_types
message(paste("Using", n_types, "repeat types with gradient coloring"))

contig_order <- unique(repeat_data$ContigID)
if (!is.null(args$contig_order)) {
  contig_order_df <- read.table(args$contig_order, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  contig_order <- contig_order_df[[1]]
  repeat_data <- repeat_data %>% filter(ContigID %in% contig_order)
  repeat_data$ContigID <- factor(repeat_data$ContigID, levels = contig_order)
  message(paste("Using custom contig order from:", args$contig_order))
} else {
  repeat_data$ContigID <- factor(repeat_data$ContigID, levels = contig_order)
}

p <- ggplot(repeat_data, aes(x = ContigID, y = Percentage, fill = RepeatType)) +
  geom_bar(stat = "identity", position = "stack", width = args$bar_width) +
  geom_text(aes(label = ifelse(Percentage >= 5, sprintf("%.1f%%", Percentage), "")),
            position = position_stack(vjust = 0.5),
            size = args$font_label * 1.5, color = "black", fontface = "bold", angle = 90, hjust = 0.5) +
  scale_fill_manual(values = spectral_colors, name = "Repeat Type") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = c(0, 0), limits = c(0, NA)) +
  labs(title = args$title, x = "Chromosome (ContigID)", y = "Repeat Percentage (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = args$font_title, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = args$font_axis),
    axis.text.y = element_text(size = args$font_axis),
    legend.title = element_text(size = args$font_legend_title),
    legend.text = element_text(size = args$font_legend_text),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")
  )

message(paste("Plotting to", args$output_prefix))

png_file <- paste0(args$output_prefix, ".png")
pdf_file <- paste0(args$output_prefix, ".pdf")

png(png_file, width = args$width, height = args$height, units = "in", res = 300)
print(p)
dev.off()

pdf(pdf_file, width = args$width, height = args$height)
print(p)
dev.off()

message(paste("Output saved to:", png_file, "and", pdf_file))
message("Done!")
