#!/usr/bin/env Rscript
# S7 重绘 (6 面板 a-f, 与原 S7 内容一致, R nature-figure 统一风格)
# 面板:
#   a 染色体 repeat 含量: core whole vs core terminal 100kb vs accessory (箱线)
#   b 分家族分布: terminal 100kb vs accessory (箱线 + 星号显著性)
#   c 按株配对 (13 isolates)
#   d 固定窗口扫描 (40-1000 kb/end, log): terminal vs accessory
#   e 按染色体大小分层: 百分比(3%) vs 固定(100kb)窗口比较
#   f 分家族窗口扫描: 各 family terminal-accessory 差
suppressMessages({library(ggplot2); library(patchwork)})

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6),
        plot.title = element_text(size = 7.5, face = "bold"),
        panel.grid = element_blank()))

RED <- "#C44E52"; BLUE <- "#3B7DD8"; GOLD <- "#E0B74E"; GREY <- "#BBBBBB"; DARK <- "#444444"

cc <- read.delim("06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_content.tsv")
fs <- read.delim("06.accessory-chrs/repeat_terminal/suppl_table.terminal_fixed_scan.tsv")
fw <- read.delim("06.accessory-chrs/repeat_terminal/suppl_table.repeat_family_window_scan.tsv")
fam <- read.delim("06.accessory-chrs/repeat_terminal/suppl_table.repeat_terminal_family_100kb.tsv")
strat <- read.delim("06.accessory-chrs/repeat_terminal/suppl_table.terminal_window_stratified.tsv")

# ---- a: overall repeat content by compartment ----
core_idx <- cc$Type == "core"
acc_idx <- cc$Type == "accessory"
a_df <- data.frame(
  compartment = c(rep("Core\n(whole)", sum(core_idx)),
                  rep("Core terminal\n(100 kb/end)", sum(core_idx)),
                  rep("Accessory", sum(acc_idx))),
  repeat_pct = c(cc$Repeat_pct[core_idx],
                 cc$Terminal100kb_pct[core_idx],
                 cc$Repeat_pct[acc_idx]))
a_df$compartment <- factor(a_df$compartment,
                           levels = c("Core\n(whole)", "Core terminal\n(100 kb/end)",
                                      "Accessory"))
p_a <- ggplot(a_df, aes(x = compartment, y = repeat_pct, fill = compartment)) +
  geom_boxplot(outlier.size = 0.6, width = 0.55) +
  scale_fill_manual(values = c("Core\n(whole)" = GREY,
                               "Core terminal\n(100 kb/end)" = GOLD,
                               "Accessory" = RED)) +
  labs(x = NULL, y = "Repeat content (%)",
       title = "Repeat content by compartment") +
  theme(legend.position = "none", axis.text.x = element_text(size = 6, lineheight = 0.8))

# ---- b: per-family distribution (terminal 100kb vs accessory) ----
fam$lab <- paste0(fam$Family, ifelse(fam$MW_p < 0.001, " ***",
                              ifelse(fam$MW_p < 0.01, " **",
                              ifelse(fam$MW_p < 0.05, " *", ""))))
fam$Family <- factor(fam$Family, levels = fam$Family[order(fam$terminal_mean)])
fam_l <- reshape(fam[, c("Family", "terminal_mean", "accessory_mean", "lab")],
                 direction = "long",
                 varying = list(c("terminal_mean", "accessory_mean")),
                 v.names = "mean_pct", timevar = "comp",
                 times = c("Core terminal (100 kb)", "Accessory"))
p_b <- ggplot(fam_l, aes(x = Family, y = mean_pct, fill = comp)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("Core terminal (100 kb)" = GOLD, "Accessory" = RED)) +
  geom_text(data = fam, aes(x = Family, y = pmax(terminal_mean, accessory_mean),
                            label = lab),
            vjust = -0.3, size = 2.2, inherit.aes = FALSE) +
  labs(x = NULL, y = "Mean repeat content (%)", fill = NULL,
       title = "Per-family repeat content") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom", legend.key.size = unit(0.25, "cm"))

# ---- c: per-isolate paired comparison ----
term_mean <- aggregate(Terminal100kb_pct ~ Isolate, data = cc[core_idx, ], FUN = mean)
acc_mean <- aggregate(Repeat_pct ~ Isolate, data = cc[acc_idx, ], FUN = mean)
c_df <- merge(term_mean, acc_mean, by = "Isolate")
colnames(c_df) <- c("Isolate", "terminal", "accessory")
p_c <- ggplot(c_df, aes(x = terminal, y = accessory)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = GREY) +
  geom_point(size = 2.6, colour = DARK) +
  ggrepel::geom_text_repel(aes(label = Isolate), size = 2.3,
                           max.overlaps = 40, box.padding = 0.7,
                           point.padding = 0.6, min.segment.length = 0,
                           segment.size = 0.25, seed = 42, force = 3,
                           force_pull = 0.3) +
  scale_x_continuous(expand = expansion(mult = 0.22)) +
  scale_y_continuous(expand = expansion(mult = 0.22)) +
  labs(x = "Core terminal repeat (%)", y = "Accessory repeat (%)",
       title = "Per-isolate paired comparison") +
  coord_equal()

# ---- d: fixed-window scan (log) ----
fs$Window_kb_per_end <- as.numeric(fs$Window_kb_per_end)
p_d <- ggplot(fs, aes(x = Window_kb_per_end)) +
  geom_line(aes(y = terminal_mean, colour = "Core terminal"), linewidth = 0.7) +
  geom_line(aes(y = accessory_mean, colour = "Accessory"), linewidth = 0.7) +
  geom_point(aes(y = terminal_mean, colour = "Core terminal"), size = 1.5) +
  geom_point(aes(y = accessory_mean, colour = "Accessory"), size = 1.5) +
  scale_colour_manual(values = c("Core terminal" = GOLD, "Accessory" = RED)) +
  scale_x_continuous(breaks = c(40, 60, 100, 200, 400, 1000), trans = "log10") +
  labs(x = "Terminal window (kb/end, log)", y = "Mean repeat content (%)",
       colour = NULL, title = "Fixed-window scan") +
  theme(legend.position = "bottom")

# ---- e: fixed vs percentage window stratified by size ----
strat$Window_lab <- ifelse(strat$Window == "fixed100kb_terminal_pct", "Fixed 100 kb",
                           "Percentage (3%)")
strat$Size_class <- factor(strat$Size_class,
                           levels = c("small (<2 Mb)", "mid (2-4 Mb)", "large (>4 Mb)"))
p_e <- ggplot(strat, aes(x = Size_class, y = terminal_mean, fill = Window_lab)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 41.45, linetype = "dashed", colour = DARK) +
  scale_fill_manual(values = c("Fixed 100 kb" = GOLD, "Percentage (3%)" = BLUE)) +
  geom_text(aes(label = sprintf("%.1f%%", terminal_mean),
                y = terminal_mean + 2.2),
            position = position_dodge(width = 0.7), vjust = 0, size = 2.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Terminal repeat content (%)", fill = NULL,
       title = "Fixed vs percentage window by chromosome size") +
  theme(axis.text.x = element_text(size = 5.5),
        legend.position = "bottom", legend.key.size = unit(0.25, "cm"))

# ---- f: family window scan (log) ----
fw$Window_kb <- as.numeric(fw$Window_kb)
fw$sig <- fw$MW_p < 0.05
p_f <- ggplot(fw, aes(x = Window_kb, y = diff, colour = Family)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = GREY) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(shape = sig), size = 1.3) +
  scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1), guide = "none") +
  scale_x_continuous(breaks = c(40, 60, 100, 150, 200), trans = "log10") +
  labs(x = "Terminal window (kb/end, log)", y = "Terminal \u2212 accessory repeat (%)",
       colour = "Family", title = "Family-specific window scan") +
  theme(legend.position = "bottom")

# ---- combine: 2x3 ----
fig <- (p_a | p_b | p_c) / (p_d | p_e | p_f) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))

out <- "06.accessory-chrs/repeat_terminal/figure_repeat_panels_r"
ggsave(paste0(out, ".pdf"), fig, width = 12, height = 8, device = cairo_pdf)
ggsave(paste0(out, ".png"), fig, width = 12, height = 8, dpi = 300)
cat("saved ->", out, ".pdf/.png (6 panels)\n")
