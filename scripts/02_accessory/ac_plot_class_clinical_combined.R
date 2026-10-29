#!/usr/bin/env Rscript
# 组3: AC 分类与临床组合图 (原 S20 + S21 + S22 -> 新 S20)
# nature-figure, backend = R
# 核心结论: FSSC AC 按 kmer 群体分布分为 block(种特异) / mosaic(跨种镶嵌) 两类;
#   mosaic AC 临床富集, 物种校正后存在真正临床 AC; BMU02483 携带多个近完整临床 mosaic AC。
# 面板:
#   a 分类景观: 78 AC 的种内共享广度 x 跨种共享广度 (block 蓝 / mosaic 红, 点大小=携带株数)
#   b mosaic 亚型: 3 亚型的跨种共享谱 (subtype1 深SC-wide / subtype2 broad / subtype3 pair-restricted)
#   c 临床 AC 鉴定: 观察 vs 物种校正期望临床比例 (对角线上方 = 临床富集, 22 个显著)
#   d RF 临床 HOG 分布: clinical-biased vs non-clinical-biased 在 mosaic/block/core 的比例
#   e BMU02483 图谱: 78 AC 上 kmer 共享谱 (临床 mosaic 红 / 非临床 mosaic 粉 / block 蓝)
suppressMessages({library(ggplot2); library(patchwork); library(scales)})

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 5.8),
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid = element_blank()))

RED <- "#C44E52"; BLUE <- "#3B7DD8"; PINK <- "#E8A2A6"; GREY <- "#BBBBBB"; DARK <- "#444444"

# ---- a: classification landscape ----
dd <- read.delim("06.accessory-chrs/kmer_analysis/kmer_type.dual_domain.tsv")
v3 <- read.delim("06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv")
mat <- read.delim("06.accessory-chrs/kmer_analysis/kmer_matrix.ACx211.tsv", row.names = 1, check.names = FALSE)
n_carriers <- apply(mat >= 0.4, 1, sum)
a_df <- merge(dd, v3[, c("AC", "kmer_type")], by = "AC")
a_df$n_carriers <- n_carriers[a_df$AC]
a_df$type_lab <- ifelse(a_df$kmer_type == "type1_stable", "Block (species-specific)", "Mosaic (cross-species)")
a_df$type_lab <- factor(a_df$type_lab, levels = c("Block (species-specific)", "Mosaic (cross-species)"))
# add transfer-concentration axis (genes per partner) and bulk/scatter shape
tf <- read.delim("06.accessory-chrs/gene_transfer_211/suppl_table.ac_transfer_profile.tsv")
a_df <- merge(a_df, tf[, c("AC", "genes_per_partner", "AC_class")], by = "AC", all.x = TRUE)
a_df$genes_per_partner[is.na(a_df$genes_per_partner)] <- 0.1
a_df$transfer_lab <- ifelse(a_df$AC_class == "Type1_block", "Block transfer",
                            "Scattered transfer")
p_a <- ggplot(a_df, aes(x = intra_breadth_0.4, y = inter_breadth_0.4, colour = type_lab)) +
  geom_point(aes(size = n_carriers, shape = transfer_lab), alpha = 0.8) +
  scale_colour_manual(values = c("Block (species-specific)" = BLUE, "Mosaic (cross-species)" = RED)) +
  scale_shape_manual(values = c("Block transfer" = 17, "Scattered transfer" = 16)) +
  scale_size_continuous(range = c(1.2, 4.5), name = "Carriers\n(k-mer \u2265 0.4)") +
  labs(x = "Intra-species sharing breadth", y = "Cross-species sharing breadth",
       title = "K-mer classification landscape (78 ACs)") +
  theme(legend.position = "right", legend.key.size = unit(0.3, "cm"))

# ---- b: mosaic subtypes ----
sub <- read.delim("06.accessory-chrs/kmer_analysis/kmer_type2.subtypes.tsv")
sub_b <- sub[, c("AC", "br_solani", "br_falciformis", "br_keratoplastica", "br_pisi", "br_other", "subtype")]
sub_b <- reshape(sub_b, direction = "long",
                 varying = list(c("br_solani", "br_falciformis", "br_keratoplastica", "br_pisi", "br_other")),
                 v.names = "breadth", timevar = "sp", times = c("solani", "falciformis", "keratoplastica", "pisi", "other"))
sub_b$subtype <- factor(sub_b$subtype, levels = c("subtype1", "subtype2", "subtype3"))
sub_b$sp <- factor(sub_b$sp, levels = c("solani", "falciformis", "keratoplastica", "pisi", "other"))
p_b <- ggplot(sub_b, aes(x = sp, y = AC, fill = breadth)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = RED, name = "Sharing\nbreadth") +
  facet_grid(subtype ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "Mosaic subtypes: cross-species sharing") +
  theme(strip.text.y = element_text(angle = 0, size = 6, face = "bold",
                                    margin = margin(l = 1, r = 1)),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
        legend.position = "right", legend.key.height = unit(0.3, "cm"))

# ---- c: clinical AC identification (x/y same scale) ----
ca <- read.delim("06.accessory-chrs/kmer_analysis/kmer_clinical_ac.burden_controlled.tsv")
ca$type <- ifelse(ca$AC %in% v3$AC[v3$kmer_type == "type1_stable"], "Block", "Mosaic")
ca$excess <- ca$obs_minus_exp
p_c <- ggplot(ca, aes(x = expected_med_frac_burden, y = med_frac, colour = type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = GREY) +
  geom_point(aes(size = n_carriers), alpha = 0.8) +
  scale_colour_manual(values = c("Block" = BLUE, "Mosaic" = RED)) +
  scale_size_continuous(range = c(1.2, 4), name = "Carriers") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_equal() +
  labs(x = "Species-mix-adjusted expected clinical fraction",
       y = "Observed clinical fraction",
       title = "Clinical AC identification (76 ACs)") +
  theme(legend.position = "right")

# ---- d: RF HOG distribution ----
# 数据读自 rf_candidate_location_by_ac_class.tsv, 由 scripts/plot_kmer_analysis_figures.py 的
# plot_rf_candidate_location()
# 用 PAV 候选集算出后写出 (join 表 = 06.accessory-chrs/hog_to_ac.tsv)。
# 不要在这里写死百分比: 上一版把 133/296、后来的 177/438 硬编码进来, 候选集一换
# 图就和图注/正文脱节, 而且不会报错。
loc_tab <- read.delim("06.accessory-chrs/kmer_analysis/rf_candidate_location_by_ac_class.tsv")
loc_lab <- c(type2_AC = "Mosaic AC", type1_AC = "Block AC",
             core_or_no13 = "Core chromosome")
loc_tab$cls <- factor(loc_lab[loc_tab$location],
                  levels = c("Core chromosome", "Block AC", "Mosaic AC"))
clin_lab <- sprintf("Clinical-biased (n = %d)", loc_tab$clinical_total[1])
non_lab <- sprintf("Non-clinical-biased (n = %d)", loc_tab$nonclinical_total[1])
rf_l <- rbind(
  data.frame(location = loc_tab$cls, bias = clin_lab, pct = loc_tab$clinical_pct),
  data.frame(location = loc_tab$cls, bias = non_lab, pct = loc_tab$nonclinical_pct))
rf_l$bias <- factor(rf_l$bias, levels = c(clin_lab, non_lab))
p_d <- ggplot(rf_l, aes(x = location, y = pct, fill = bias)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = setNames(c(RED, BLUE), c(clin_lab, non_lab))) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 2.0) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "RF clinical-associated HOGs (%)",
       title = "RF clinical HOGs across chromosome classes") +
  theme(axis.text.x = element_text(size = 6), legend.position = "bottom",
        legend.title = element_blank())

# ---- e: BMU02483 sliding-window heatmap + overall profile ----
# 左侧: 临床 AC 沿染色体位置 (0-100%) 与 BMU02483 共有 k-mer 比例热图 (RdYlBu)
#   每窗口用 geom_rect 按其真实 [start/len, end/len] 区间绘制, 铺满 0-100%
# 右侧: 各 AC 整体共有比例柱状图, 按 block/mosaic 着色 (与热图 Y 轴顺序一致)
sl <- read.delim("06.accessory-chrs/BMU02483_gene_content/sliding/sliding_2483_prop.tsv")
ac_order <- sort(unique(sl$AC), decreasing = TRUE)
sl$AC <- factor(sl$AC, levels = ac_order)
# 长度归一化区间 [xmin, xmax]
ac_len <- tapply(sl$end, sl$AC, max)
sl$xmin <- sl$start / ac_len[as.character(sl$AC)]
sl$xmax <- sl$end / ac_len[as.character(sl$AC)]
# 整体比例 (整 AC kmer 共享), 同一 AC 顺序 + block/mosaic 分类
b_prof <- data.frame(AC = rownames(mat), share = mat[, "BMU02483"])
b_prof <- b_prof[b_prof$AC %in% ac_order, ]
b_prof$AC <- factor(b_prof$AC, levels = ac_order)
b_prof$type <- ifelse(b_prof$AC %in% v3$AC[v3$kmer_type == "type1_stable"],
                      "Block", "Mosaic")
b_prof$type <- factor(b_prof$type, levels = c("Block", "Mosaic"))

p_e_left <- ggplot(sl, aes(xmin = xmin, xmax = xmax, y = AC, fill = prop)) +
  geom_rect(height = 0.85) +
  scale_fill_gradientn(colours = c("#4575B4", "#E0F3F8", "#FFFFBF",
                                   "#FDAE61", "#D73027"),
                       values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
                       limits = c(0, 1), name = "Shared\nk-mer") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "25", "50", "75", "100%")) +
  labs(x = "Chromosome position (%)", y = NULL,
       title = "BMU02483 shared k-mers along clinical ACs") +
  theme(axis.text.y = element_text(size = 6.5),
        legend.position = "right", legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(size = 5.5))

p_e_right <- ggplot(b_prof, aes(x = share, y = AC, fill = type)) +
  geom_col(width = 0.7, alpha = 0.9) +
  scale_fill_manual(values = c("Block" = BLUE, "Mosaic" = RED)) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Whole-AC k-mer\nsharing", y = NULL, fill = NULL,
       title = "Overall") +
  theme(axis.text.y = element_text(size = 6.5),
        axis.text.x = element_text(size = 6),
        legend.position = "bottom", legend.key.size = unit(0.25, "cm"))

p_e <- p_e_left + p_e_right + plot_layout(widths = c(3, 1))

# ---- combine: 2x3 layout (e = sliding heatmap + overall, full bottom row) ----
fig <- (p_a | p_b) / (p_c | p_d) / p_e +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))

out <- "06.accessory-chrs/kmer_analysis/figure_combined_ac_class_clinical"
ggsave(paste0(out, ".pdf"), fig, width = 11, height = 12, device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 11, height = 12, dpi = 300)
cat("saved ->", out, ".pdf/.png\n")
