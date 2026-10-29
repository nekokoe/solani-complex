#!/usr/bin/env Rscript
# T2 审稿回应 #505/#61: PERMANOVA (adonis2) 支持 PCoA 的临床/环境分离结论。
# 输入: 预处理的 PAV 矩阵 (tmp/permanova_pav_matrix.tsv, 由 prepare_permanova.py 生成)
#       + med_env 分组 + species 分组
# 分析:
#   A. 全 FSSC (211 株): adonis2(vegdist(pav, method="bray") ~ Source, strata=Species)
#   B. F. falciformis 内: adonis2(~ Source)
#   C. F. solani 内:     adonis2(~ Source)
# 输出: 05.stats_pan_by_species_shorter/permanova_results.tsv
suppressMessages({
  library(vegan)
})

BASE <- "05.stats_pan_by_species_shorter"
mat <- read.delim("tmp/permanova_pav_matrix.tsv", row.names = 1, check.names = FALSE)
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

dist_full <- vegdist(mat, method = "bray")
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

# A. full FSSC, strata = species
fitA <- adonis2(dist_full ~ source, data = meta, permutations = 999,
                strata = meta$species, by = "terms")
add("FSSC_all_strata_species", fitA, nrow(meta), "source")

# B/C. per species
for (spn in c("falciformis", "solani")) {
  idx <- which(meta$species == spn)
  if (length(idx) >= 10) {
    fit <- adonis2(vegdist(mat[idx, , drop = FALSE], method = "bray") ~ source,
                   data = meta[idx, ], permutations = 999)
    add(paste0(spn, "_within"), fit, length(idx), "source")
  }
}
write.table(res, file.path(BASE, "permanova_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
