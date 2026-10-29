#!/usr/bin/env Rscript
# CNV (core-HOG copy number) PERMANOVA — control alongside the PAV PERMANOVA
# (scripts/stats_permanova.R).
# Input: the same sample x core-HOG copy-number matrix that produced the
# manuscript's core PCoA:
#   05.pan-genome/accessory_genes/N0.HOG.hog_count.allSC_core.pcoa/subset_feature_matrix.tsv
#   (SampleID x 6227 core HOGs, values = gene copy number)
# Metadata merged exactly as the PAV script: med_env + species,
# source=='unknown' -> NA -> dropped.
# Distance: euclidean on the copy-number matrix (natural for CNV). Also
# computes bray-curtis for cross-check (the core PCoA used braycurtis).
# Output: 05.stats_pan_by_species_shorter/permanova_results_cnv.tsv
suppressMessages({ library(vegan) })

BASE <- "05.stats_pan_by_species_shorter"
MATRIX <- "05.pan-genome/accessory_genes/N0.HOG.hog_count.allSC_core.pcoa/subset_feature_matrix.tsv"

mat <- read.delim(MATRIX, row.names = 1, check.names = FALSE)
samples <- rownames(mat)

med <- read.delim(file.path(BASE, "bmu_ref_groups.med_env.txt"), header = FALSE,
                  col.names = c("sample", "source"), stringsAsFactors = FALSE)
sp  <- read.delim(file.path(BASE, "bmu_ref_groups.species.txt"), header = FALSE,
                  col.names = c("sample", "species"), stringsAsFactors = FALSE)
meta <- merge(med, sp, by = "sample")
meta <- meta[meta$sample %in% samples, ]
meta$source[meta$source == "unknown"] <- NA
meta <- meta[!is.na(meta$source), ]
mat <- mat[meta$sample, ]

cat("matrix dims (after metadata merge):", nrow(mat), "x", ncol(mat), "\n")

run <- function(method) {
  dist_full <- vegdist(mat, method = method)
  res <- data.frame(analysis = character(), n = integer(), term = character(),
                    Df = integer(), R2 = numeric(), F = numeric(), p = numeric(),
                    stringsAsFactors = FALSE)
  add <- function(analysis, fit, n, term) {
    tab <- as.data.frame(fit)
    rn <- rownames(tab)
    pick <- if (term %in% rn) term else "Model"
    row <- tab[pick, ]
    res <<- rbind(res, data.frame(analysis = analysis, n = n,
                                  term = term, Df = row$Df, R2 = row$R2,
                                  F = row$F, p = row$`Pr(>F)`))
  }
  set.seed(1)
  fitA <- adonis2(dist_full ~ source, data = meta, permutations = 999,
                  strata = meta$species, by = "terms")
  add("FSSC_all_strata_species", fitA, nrow(meta), "source")
  for (spn in c("falciformis", "solani")) {
    idx <- which(meta$species == spn)
    if (length(idx) >= 10) {
      set.seed(1)
      fit <- adonis2(vegdist(mat[idx, , drop = FALSE], method = method) ~ source,
                     data = meta[idx, ], permutations = 999)
      add(paste0(spn, "_within"), fit, length(idx), "source")
    }
  }
  res
}

res_euclidean <- run("euclidean")
res_bray      <- run("bray")

cat("=== euclidean (primary) ===\n")
print(res_euclidean)
cat("=== bray (cross-check) ===\n")
print(res_bray)

write.table(res_euclidean, file.path(BASE, "permanova_results_cnv.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", file.path(BASE, "permanova_results_cnv.tsv"), "\n")
