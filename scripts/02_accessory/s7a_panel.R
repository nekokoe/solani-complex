#!/usr/bin/env Rscript
# s7a_panel.R — S7 panel a as a standalone square panel: repeat content by
# compartment (whole core vs core terminal 100 kb vs accessory).
#
# Same conventions as s13a_panel.R / s13e_panel.R: square canvas, legible font
# sizes driven by one number, NimbusSans (Arial resolves to it on this host), and
# the colours the S7 assembly already uses.
#
# The message the panel carries: the terminal 100 kb of a core chromosome is as
# repeat-rich as an accessory chromosome, while the core chromosome as a whole is
# not (medians 42.1 % / 44.8 % / 8.0 %).
#
# Source: 06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_content.tsv
# (166 core chromosomes, 78 accessory chromosomes, 13 TGS assemblies)

suppressMessages({library(ggplot2); library(ggsignif)})

S7_GREY <- "#BBBBBB"; S7_GOLD <- "#E0B74E"; S7_RED <- "#C44E52"

read_repeat_content <- function(path) {
  cc <- read.delim(path, stringsAsFactors = FALSE)
  need <- c("Type", "Repeat_pct", "Terminal100kb_pct")
  if (!all(need %in% names(cc))) stop("missing columns: ",
                                      paste(setdiff(need, names(cc)), collapse = ", "))
  cc
}

s7a_data <- function(cc) {
  core <- cc$Type == "core"; acc <- cc$Type == "accessory"
  lv <- c("Core\n(whole)", "Core terminal\n(100 kb/end)", "Accessory")
  df <- data.frame(
    compartment = c(rep(lv[1], sum(core)), rep(lv[2], sum(core)),
                    rep(lv[3], sum(acc))),
    repeat_pct = c(cc$Repeat_pct[core], cc$Terminal100kb_pct[core],
                   cc$Repeat_pct[acc]))
  df$compartment <- factor(df$compartment, levels = lv)
  df
}

# significance labels come from the existing test table rather than being typed in
s7a_stars <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                            ifelse(p < 0.05, "*", "ns")))

# read the pooled Mann-Whitney comparisons and map them onto the x levels
s7a_signif <- function(path, levels) {
  tt <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  tt <- tt[grepl("pooled", tt$Comparison), ]
  key <- c(Core_whole = levels[1],
           Core_terminal100kb = levels[2],
           Accessory = levels[3])
  out <- list(comparisons = list(), annotations = character())
  for (i in seq_len(nrow(tt))) {
    parts <- strsplit(sub("\\s*\\(pooled.*$", "", tt$Comparison[i]), "_vs_")[[1]]
    lab <- key[c(parts[1], parts[2])]
    if (any(is.na(lab))) next
    out$comparisons[[length(out$comparisons) + 1]] <- unname(lab)
    out$annotations <- c(out$annotations, s7a_stars(tt$MannWhitney_p[i]))
  }
  out
}

s7a_panel <- function(df, family = "Arial", fontsize = 9,
                      signif = NULL, title = "Repeat content by compartment") {
  ggplot(df, aes(x = compartment, y = repeat_pct, fill = compartment)) +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.4, width = 0.58) +
    scale_fill_manual(values = c("Core\n(whole)" = S7_GREY,
                                 "Core terminal\n(100 kb/end)" = S7_GOLD,
                                 "Accessory" = S7_RED)) +
    labs(x = NULL, y = "Repeat content (%)", title = title) +
    (if (!is.null(signif)) {
       lv <- levels(df$compartment)
       span <- vapply(signif$comparisons, function(cc)
         abs(match(cc[1], lv) - match(cc[2], lv)) > 1, logical(1))
       hi <- max(df$repeat_pct) * 1.24; lo <- max(df$repeat_pct) * 1.03
       geom_signif(
         comparisons = signif$comparisons, annotations = signif$annotations,
         # a span-2 bracket goes on top, the adjacent ones sit side by side below it
         y_position = ifelse(span, hi, lo),
         tip_length = 0.012, textsize = fontsize / 3.2, vjust = 0.45,
         colour = "black", size = 0.35)
     } else NULL) +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(size = fontsize, colour = "black", lineheight = 0.85),
      axis.text.y = element_text(size = fontsize, colour = "black"),
      axis.title.y = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold", hjust = 0),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 4, 4)
    )
}
