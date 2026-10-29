#!/usr/bin/env Rscript

library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

option_list <- list(
  make_option(c("-g", "--gene_chrom_file"), type="character", default=NULL,
              help="Gene to chromosome affiliation file (tsv)", metavar="character"),
  make_option(c("-m", "--gene_hog_file"), type="character", default=NULL,
              help="Gene to HOG mapping file (tsv)", metavar="character"),
  make_option(c("-p", "--hog_interval_file"), type="character", default=NULL,
              help="HOG to presence rate interval file (tsv)", metavar="character"),
  make_option(c("--selected_contig"), type="character", default=NULL,
              help="File with list of contigs to plot (tsv, single column)", metavar="character"),
  make_option(c("-o", "--output_prefix"), type="character", default="chromosome_presence_rate",
              help="Output file prefix", metavar="character"),
  make_option(c("--width"), type="numeric", default=12,
              help="Plot width in inches", metavar="numeric"),
  make_option(c("--height"), type="numeric", default=6,
              help="Plot height in inches", metavar="numeric"),
  make_option(c("--title"), type="character", default="Gene Presence Rate Distribution by Chromosome",
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
              help="Bar width (default: 0.8)", metavar="numeric")
)

parser <- OptionParser(option_list=option_list)
args <- parse_args(parser)

if (is.null(args$gene_chrom_file) || is.null(args$gene_hog_file) || is.null(args$hog_interval_file)) {
  print_help(parser)
  stop("All three input files must be provided", call.=FALSE)
}

read_tsv_auto_header <- function(file_path, expected_cols, col_names) {
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  first_line <- readLines(file_path, n = 1)
  has_header <- grepl(paste0("^", expected_cols[1]), first_line)
  
  if (has_header) {
    read.table(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  } else {
    df <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    colnames(df) <- col_names
    df
  }
}

message("Reading gene-chromosome affiliation file...")
gene_chrom <- read_tsv_auto_header(args$gene_chrom_file,
                                   expected_cols = c("TranscriptID", "BMU"),
                                   col_names = c("TranscriptID", "GeneID", "ContigID"))

message("Reading gene-HOG mapping file...")
gene_hog <- read_tsv_auto_header(args$gene_hog_file,
                                 expected_cols = c("GeneID", "BMU"),
                                 col_names = c("TranscriptID", "SampleID", "HOG_ID", "OG_ID"))
colnames(gene_hog)[1] <- "TranscriptID"

message("Reading HOG-presence rate interval file...")
hog_interval <- read_tsv_auto_header(args$hog_interval_file,
                                     expected_cols = c("HOG", "N0"),
                                     col_names = c("HOG", "Presence_Rate", "Presence_Rate_Interval"))

merged_data <- gene_chrom %>%
  select(TranscriptID, ContigID) %>%
  inner_join(gene_hog %>% select(TranscriptID, HOG_ID), by = "TranscriptID") %>%
  inner_join(hog_interval %>% select(HOG_ID = HOG, Presence_Rate_Interval), by = "HOG_ID") %>%
  distinct(TranscriptID, ContigID, HOG_ID, Presence_Rate_Interval)

if (!is.null(args$selected_contig)) {
  selected_contig_df <- read.table(args$selected_contig, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  selected_contigs <- selected_contig_df[[1]]
  merged_data <- merged_data %>% filter(ContigID %in% selected_contigs)
  message(paste("Filtered to", length(selected_contigs), "selected contigs"))
}

message(paste("Total transcripts with mapping:", nrow(merged_data)))
message(paste("Unique contigs:", length(unique(merged_data$ContigID))))

interval_order <- c("0%-20%", "20%-40%", "40%-60%", "60%-80%", "80%-100%")
merged_data$Presence_Rate_Interval <- factor(merged_data$Presence_Rate_Interval, levels = interval_order)

contig_gene_counts <- merged_data %>%
  group_by(ContigID) %>%
  summarise(Total_Genes = n())

interval_counts <- merged_data %>%
  group_by(ContigID, Presence_Rate_Interval) %>%
  summarise(Count = n(), .groups = "drop")

interval_counts <- interval_counts %>%
  left_join(contig_gene_counts, by = "ContigID") %>%
  mutate(Percentage = Count / Total_Genes * 100)

contig_order <- interval_counts %>%
  filter(Presence_Rate_Interval == "80%-100%") %>%
  arrange(desc(Percentage)) %>%
  pull(ContigID)

interval_counts$ContigID <- factor(interval_counts$ContigID, levels = contig_order)

n_contigs <- length(contig_order)

contig_order_file <- paste0(args$output_prefix, "_contig_order.txt")
write.table(data.frame(ContigID = contig_order), contig_order_file,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message(paste("Contig order saved to:", contig_order_file))
spectral_colors <- rev(brewer.pal(n = 5, name = "Spectral"))
names(spectral_colors) <- interval_order

p <- ggplot(interval_counts, aes(x = ContigID, y = Percentage, fill = Presence_Rate_Interval)) +
  geom_bar(stat = "identity", position = position_fill(vjust = 0.5), width = args$bar_width) +
  geom_text(aes(label = ifelse(Percentage >= 5, sprintf("%.1f%%", Percentage), "")),
            position = position_fill(vjust = 0.5),
            size = args$font_label * 1.5, color = "black", fontface = "bold", angle = 90, hjust = 0.5) +
  scale_fill_manual(values = spectral_colors, name = "Presence Rate\nInterval") +
  scale_y_continuous(labels = scales::percent_format(scale = 100), expand = c(0, 0)) +
  labs(title = args$title, x = "Chromosome (ContigID)", y = "Gene Percentage (%)") +
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
