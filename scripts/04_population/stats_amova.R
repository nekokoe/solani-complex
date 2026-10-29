#!/usr/bin/env Rscript
# T2 审稿回应 #61 — AMOVA (pegas) 量化临床/环境分组的遗传分化。
# 输入: tmp/amova/<sp>_dist.npy (numpy 距离矩阵) + <sp>_samples.tsv + med_env 分组
# 输出: 05.stats_pan_by_species_shorter/amova_results.tsv
suppressMessages({
  library(pegas)
})

dist2 <- function(path) {
  m <- as.matrix(read.table(path))
  d <- as.dist(m)
  d
}

read_npy_dummy <- function(sp) {
  m <- as.matrix(read.delim(sprintf("tmp/amova/%s_dist_matrix.tsv", sp),
                            row.names = 1, check.names = FALSE))
  m
}

med <- read.delim("05.stats_pan_by_species_shorter/bmu_ref_groups.med_env.txt",
                  header = FALSE, col.names = c("sample", "source"), stringsAsFactors = FALSE)

res <- data.frame(species = character(), n = integer(), n_med = integer(),
                  n_env = integer(), phi_ST = numeric(),
                  sigma2_source = numeric(), sigma2_error = numeric(),
                  p = numeric(), stringsAsFactors = FALSE)

for (sp in c("solani", "falciformis", "keratoplastica")) {
  dd <- read_npy_dummy(sp)
  samples <- rownames(dd)
  grp <- data.frame(sample = samples, stringsAsFactors = FALSE)
  grp <- merge(grp, med, by = "sample", all.x = TRUE)
  grp$source[is.na(grp$source)] <- "unknown"
  # keep only med/env for the test (unknown dropped from AMOVA)
  keep <- grp$source %in% c("med", "env")
  if (sum(keep) < 6) { cat(sp, ": too few grouped samples, skip\n"); next }
  dsub <- as.dist(dd[keep, keep])
  gsub <- factor(grp$source[keep])
  D <- dsub
  fit <- amova(D ~ gsub, nperm = 999)
  vc <- fit$varcomp
  s2_g <- vc["gsub", "sigma2"]
  s2_e <- vc["Error", "sigma2"]
  phi_st <- s2_g / (s2_g + s2_e)
  res <- rbind(res, data.frame(species = sp, n = sum(keep),
                               n_med = sum(gsub == "med"), n_env = sum(gsub == "env"),
                               phi_ST = phi_st,
                               sigma2_source = s2_g, sigma2_error = s2_e,
                               p = vc["gsub", "P.value"],
                               stringsAsFactors = FALSE))
}
write.table(res, "05.stats_pan_by_species_shorter/amova_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
print(res)
