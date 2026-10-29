#!/usr/bin/env Rscript
# s13a_panel.R — S13 panel a: AC-related gene transfer, split by species.
#
# Sourced by plot_s13a_transfer_by_species.R (standalone square panel) and by
# ac_plot_transfer_mechanism_combined_v2.R (the a-e assembly), so the two cannot
# drift apart.
#
# Attribution matters here. The counted genes are the 13 chromosome-level (TGS)
# assemblies' genes; the 211 genomes only supply the mmseqs clusters those genes
# sit in. See scripts/stats_gene_transfer_211.py: ACC_TYPE is built from the 13 TGS
# isolates alone, and only members whose isolate is in tgs_set are recorded.
#
# BOTH classes are normalised over the genes the group carries on accessory
# chromosomes, and the BARS therefore carry that percentage -- the y axis has a single
# scale and there is no secondary axis. The gene COUNT is annotation only: it rides on
# the label's second line and never sets a bar height. That share is what separates
# F. solani (C-A 8.8 %) from F. falciformis (C-A 2.4 %); comparing raw counts instead
# is misleading because the groups differ about eightfold in AC gene content.

suppressMessages({library(ggplot2)})

SPECIES_OF <- c(BMU01933 = "falciformis", BMU02028 = "solani", BMU02474 = "solani",
                BMU03267 = "falciformis", BMU06673 = "keratoplastica",
                BMU12683 = "solani", BMU12684 = "solani", BMU12685 = "solani",
                BMU12686 = "solani", Fu3.1 = "falciformis", Fu6.1 = "keratoplastica",
                SB1 = "solani", SX13 = "solani")

S13_RED <- "#C44E52"; S13_BLUE <- "#3B7DD8"

dual_counts <- read.delim(file.path("06.accessory-chrs", "gene_transfer_211_99",
                                    "suppl_table.transfer_dual_counts.tsv"))

s13a_table <- function(table, with_all = TRUE) {
  df <- read.delim(table, stringsAsFactors = FALSE)
  df$species <- SPECIES_OF[df$Isolate]
  if (any(is.na(df$species)))
    stop("unmapped isolate(s): ", paste(unique(df$Isolate[is.na(df$species)]), collapse = ", "))

  sp_order <- c("solani", "falciformis", "keratoplastica")

  # DENOMINATOR for both classes = the genes the group carries ON ACCESSORY
  # CHROMOSOMES. Only genes sitting on an accessory contig are counted, so the
  # quantity is "of the genes on accessory chromosomes, what fraction is involved
  # in each kind of exchange". Normalising C-A over every C-A gene instead (which
  # includes the ones that have already moved onto a core chromosome) mixes two
  # different denominators and inflates that bar roughly threefold.
  acc <- df[df$Contig_type == "accessory", ]

  one <- function(grp) {
    sub <- if (grp == "all13") acc else acc[acc$species == grp, ]
    den <- nrow(sub)
    # The two numerators are computed INDEPENDENTLY over the same denominator:
    # A-A and C-A are transfer EVENTS, not exclusive classes, so a gene whose
    # cluster contains both accessory-accessory pairs and core-accessory pairs
    # is counted in both bars. The table labels each gene with one Transfer_type
    # (C-A whenever its cluster mixes accessory and core chromosomes), which
    # hides the accessory-accessory pairs of mixed clusters; the dual count
    # recovered from the raw mmseqs clusters is therefore added back to the A-A
    # numerator (suppl_table.transfer_dual_counts.tsv).
    dual_n <- dual_counts$dual[dual_counts$grp == grp]
    if (length(dual_n) == 0) dual_n <- 0
    aa <- sum(sub$Transfer_type == "A-A") + dual_n
    ca <- sum(sub$Transfer_type == "C-A")
    data.frame(grp = grp, type = c("A-A", "C-A"), n = c(aa, ca),
               pct = 100 * c(aa, ca) / den, acc_genes = den, dual = dual_n)
  }
  rows <- do.call(rbind, lapply(if (with_all) c(sp_order, "all13") else sp_order, one))
  if (any(rows$dual != 0))
    warning("genes participating in both A-A and C-A events: ",
            sum(rows$dual[rows$type == "A-A"]), " (they are counted in both bars)")
  rows$grp <- factor(rows$grp, levels = if (with_all) c(sp_order, "all13") else sp_order)
  rows$type <- factor(rows$type, levels = c("A-A", "C-A"))
  rows
}

s13a_panel <- function(rows, family = "Arial", fontsize = 11,
                       show_subtitle = FALSE, show_pct = FALSE) {
  labs_sp <- c(solani = "F. solani", falciformis = "F. falciformis",
               keratoplastica = "F. keratoplastica")

  # Bars carry the PERCENTAGE; the y axis therefore has one scale only. The gene
  # COUNT is annotation, not geometry: it rides on the label's second line and is
  # never plotted, so no bar height and no second axis depend on it.
  p <- ggplot(rows, aes(x = grp, y = pct)) +
    geom_col(aes(fill = type), position = position_dodge(width = 0.72), width = 0.66) +
    geom_text(aes(label = sprintf("%.1f%%\n%s", pct, format(n, big.mark = ",")),
                  group = type),
              position = position_dodge(width = 0.72),
              vjust = -0.35, lineheight = 0.92,
              size = (fontsize - 1.1) / 2.845) +
    scale_fill_manual(values = c("A-A" = S13_RED, "C-A" = S13_BLUE),
                      labels = c("A-A (accessory\u2013accessory)",
                                 "C-A (core\u2013accessory)"), name = NULL) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
    scale_x_discrete(labels = function(x) {
      parse(text = ifelse(x == "all13", "plain('all 13')",
                          paste0("italic('", labs_sp[x], "')")))
    }) +
    scale_y_continuous(
      limits = c(0, max(rows$pct) * 1.34),
      expand = expansion(mult = c(0, 0.02)),
      labels = function(v) paste0(v, "%")) +
    labs(x = NULL, y = "Share of accessory-chromosome genes",
         title = "AC-related gene transfer",
         subtitle = if (show_subtitle)
           "genes on the 13 chromosome-level assemblies\nclusters from 211 genomes (100% identity)" else NULL,
         # The footer states the PERCENTAGE's definition, which is the panel's one
         # easily-misread property: the two classes are counted independently over the
         # same denominator, so a gene involved in both kinds of exchange is counted in
         # both bars and the two therefore need not sum to 100%.
         caption = "percentages are of the genes carried on accessory chromosomes;\nA-A and C-A are counted independently, so a gene in both is counted twice") +
    theme_classic(base_size = fontsize, base_family = family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(size = fontsize, colour = "black",
                                 angle = 60, hjust = 1),
      axis.text.y = element_text(size = fontsize, colour = "black"),
      axis.title.y = element_text(size = fontsize),
      plot.title = element_text(size = fontsize + 1, face = "bold"),
      plot.subtitle = element_text(size = fontsize - 2, colour = "grey25", lineheight = 1.05),
      plot.caption = element_text(size = fontsize - 2.5, colour = "grey35", hjust = 0),
      legend.position = "bottom",
      legend.box.spacing = unit(1.2, "mm"),
      legend.spacing.x = unit(2.2, "mm"),
      legend.text = element_text(size = fontsize),
      legend.key.size = unit(4.4, "mm"),
      legend.margin = margin(t = -1),
      panel.grid = element_blank(),
      plot.margin = margin(4.2, 4, 4.4, 4)
    )
  p
}
