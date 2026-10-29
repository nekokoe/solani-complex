#!/usr/bin/env Rscript

options(error = function() traceback(sys.calls(), max.lines = 200))
options(warn = 1)

library(optparse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(dplyr)

option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL,
              help="Path to input data file (contig x contig matrix)", metavar="character"),
  make_option(c("--output_prefix"), type="character", default=NULL,
              help="Path to output heatmap file prefix (will output both png and pdf)", metavar="character"),
  make_option(c("--width"), type="numeric", default=10,
              help="Heatmap width in inches (default: 10)", metavar="numeric"),
  make_option(c("--height"), type="numeric", default=10,
              help="Heatmap height in inches (default: 10)", metavar="numeric"),
  make_option(c("--long_format"), type="character", default=NULL,
              help="Output file for long format matrix (Query, Ref, Percent)", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

main <- function() {
  tryCatch({
    cat(paste("Reading input file:", opt$input_file, "\n"))
    d <- read.table(opt$input_file, sep = ",", row.names = 1, header = FALSE)

    header_file <- paste(opt$input_file, ".header", sep = "")
    cat(paste("Reading header file:", header_file, "\n"))
    header <- read.table(header_file, sep = ",")
    header <- as.matrix(header)
    colnames(d) <- header[1,]

    d <- d[, -dim(d)[2]]
    cat(paste("Data dimensions:", dim(d)[1], "x", dim(d)[2], "\n"))

    cat("Transposing matrix...\n")
    d <- t(d)
    cat(paste("Transposed data dimensions:", dim(d)[1], "x", dim(d)[2], "\n"))

    if (!is.null(opt$long_format)) {
      cat(paste("Writing long format matrix to:", opt$long_format, "\n"))
      long_df <- melt(t(d), varnames = c("Query", "Ref"), value.name = "Percent")
      long_df <- long_df %>%
        mutate(Query = as.character(Query), Ref = as.character(Ref)) %>%
        filter(Query != Ref, Percent >= 0.20)
      write.table(long_df, opt$long_format, sep = "\t", quote = FALSE, row.names = FALSE)
      cat(paste("Long format written:", nrow(long_df), "rows\n"))
    }

    breaks <- seq(0, 1, by = 0.001)

    output_prefix <- ifelse(is.null(opt$output_prefix), opt$input_file, opt$output_prefix)
    output_pdf <- paste0(output_prefix, ".pdf")
    output_png <- paste0(output_prefix, ".png")
    cat(paste("Output files:", output_pdf, "and", output_png, "\n"))

    cat("Generating heatmap...\n")
    pheatmap(d,
             scale = "none",
             cluster_rows = F,
             cluster_cols = F,
             clustering_method = "ward.D2",
             color = colorRampPalette(rev(brewer.pal(11, "Spectral")))(length(breaks)),
             legend_breaks = seq(0, 1, by = 0.2),
             breaks = breaks,
             show_rownames = TRUE,
             show_colnames = TRUE,
             filename = output_pdf,
             width = opt$width,
             height = opt$height)

    pheatmap(d,
             scale = "none",
             cluster_rows = F,
             cluster_cols = F,
             clustering_method = "ward.D2",
             color = colorRampPalette(rev(brewer.pal(11, "Spectral")))(length(breaks)),
             legend_breaks = seq(0, 1, by = 0.2),
             breaks = breaks,
             show_rownames = TRUE,
             show_colnames = TRUE,
             filename = output_png,
             width = opt$width,
             height = opt$height)

    cat(paste("Heatmaps saved to:", output_pdf, "and", output_png, "\n"))

  }, error = function(e) {
    cat("Error occurred:", conditionMessage(e), "\n")
    traceback()
    stop(e)
  })
}

main()
