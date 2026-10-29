#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Combine all repeat-analysis figures into ONE multi-panel supplementary
figure (S7, panels a-f):
  a = overall (whole core vs terminal 100 kb vs accessory)
  b = per-family (terminal vs accessory boxplots)
  c = per-isolate (terminal vs accessory boxplots)
  d = fixed-length window scan
  e = chromosome-size stratification (fixed vs percentage windows)
  f = family x window scan
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

C_TERM = "#E69F00"
C_ACC = "#C44E52"
C_CORE = "#6B6B6B"

cont = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_content.tsv", sep="\t")
byf = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_by_family.tsv", sep="\t")
accf = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.accessory_repeat_by_family.tsv", sep="\t")
mwp = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.repeat_terminal_family_100kb.tsv", sep="\t")
fixed = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.terminal_fixed_scan.tsv", sep="\t")
strat = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.terminal_window_stratified.tsv", sep="\t")
famwin = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.repeat_family_window_scan.tsv", sep="\t")

fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.8), dpi=300)
(axa, axb, axc), (axd, axe, axf) = axes

def style_box(bp, cols):
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.75)
        patch.set_edgecolor("black"); patch.set_linewidth(0.6)
    for w, c in zip(bp["whiskers"], cols + cols):
        w.set_color(c); w.set_linewidth(0.7)
    for cap, c in zip(bp["caps"], cols + cols):
        cap.set_color(c); cap.set_linewidth(0.7)

def clean(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)

# ---------------- a: overall ----------------
core_whole = cont[cont.Type == "core"].Repeat_pct
core_term = cont[cont.Type == "core"].Terminal100kb_pct
acc_all = cont[cont.Type == "accessory"].Repeat_pct
bp = axa.boxplot([core_whole, core_term, acc_all], positions=[0, 1, 2], widths=0.6,
                 patch_artist=True, showfliers=False,
                 medianprops=dict(color="black", lw=0.9))
style_box(bp, [C_CORE, C_TERM, C_ACC])
axa.set_xticks([0, 1, 2])
axa.set_xticklabels(["Core\nwhole", "Core\nterminal 100 kb", "Accessory"], fontsize=6.5)
axa.set_ylabel("Repeat content (%)", fontsize=7)
axa.tick_params(labelsize=6)
axa.text(0.5, max(core_whole.max(), core_term.max(), acc_all.max()) * 0.97,
         "***", ha="center", fontsize=7)
axa.set_title("a  Repeat content: whole core vs terminal 100 kb vs accessory",
              fontsize=7.5, loc="left")
clean(axa)

# ---------------- b: per-family ----------------
fams = ["DNA", "LTR", "LINE", "RC", "Simple_repeat", "Low_complexity", "Unknown", "Other"]
p_map = dict(zip(mwp.Family, mwp.MW_p))
xs, labels, max_whisk = [], [], 0
for i, fam in enumerate(fams):
    core = byf[(byf.Family == fam) & (byf.Type == "core")].Terminal100kb_pct
    acc = accf[accf.Family == fam].Accessory_pct
    if len(core) == 0 or len(acc) == 0:
        continue
    x0 = i * 2
    bp = axb.boxplot([core, acc], positions=[x0, x0 + 0.75], widths=0.62,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", lw=0.9))
    style_box(bp, [C_TERM, C_ACC])
    whisk_top = max(w.get_ydata()[1] for w in bp["whiskers"])
    max_whisk = max(max_whisk, whisk_top)
    p = p_map.get(fam, 1)
    star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "n.s."))
    axb.annotate("", xy=(x0 + 0.75, whisk_top + 1.2), xytext=(x0, whisk_top + 1.2),
                 arrowprops=dict(arrowstyle="-", lw=0.5, color="grey"))
    axb.text(x0 + 0.375, whisk_top + 1.8, star, ha="center", fontsize=6)
    xs.append(x0 + 0.375); labels.append(fam)
axb.set_xticks(xs)
axb.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.2)
axb.set_xlim(-0.5, xs[-1] + 1.6)
axb.set_ylim(0, max_whisk * 1.22)
axb.set_ylabel("Repeat content (%)", fontsize=7)
axb.tick_params(labelsize=6)
axb.set_title("b  Per-family repeat content", fontsize=7.5, loc="left")
clean(axb)

# ---------------- c: per-isolate ----------------
ISOS = ["BMU01933","BMU02028","BMU02474","BMU03267","BMU06673","BMU12683",
        "BMU12684","BMU12685","BMU12686","Fu3.1","Fu6.1","SB1","SX13"]
sp = {}
for line in open("06.accessory-chrs/tgs_assemblies/species.211.map.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2: sp[p[0]] = p[1]
SO = {"falciformis": 0, "solani": 1, "keratoplastica": 2}
order = sorted(ISOS, key=lambda i: (SO.get(sp.get(i, "solani"), 9), i))
xs, labels = [], []
for i, iso in enumerate(order):
    sub = cont[cont.Isolate == iso]
    core = sub[sub.Type == "core"].Terminal100kb_pct
    acc = sub[sub.Type == "accessory"].Repeat_pct
    if len(core) == 0 or len(acc) == 0:
        continue
    x0 = i * 2
    bp = axc.boxplot([core, acc], positions=[x0, x0 + 0.75], widths=0.62,
                     patch_artist=True, showfliers=False,
                     medianprops=dict(color="black", lw=0.9))
    style_box(bp, [C_TERM, C_ACC])
    xs.append(x0 + 0.375); labels.append(iso.replace("BMU", ""))
axc.set_xticks(xs)
axc.set_xticklabels(labels, rotation=60, ha="right", fontsize=5)
axc.set_xlim(-0.5, xs[-1] + 1.8)
axc.set_ylabel("Repeat content (%)", fontsize=7)
axc.tick_params(labelsize=6)
axc.set_title("c  Per-isolate: terminal 100 kb vs accessory", fontsize=7.5, loc="left")
clean(axc)

# ---------------- d: fixed-length window scan ----------------
acc_ref = fixed.accessory_mean.iloc[0]
axd.plot(fixed.Window_kb_per_end, fixed.terminal_mean, "o-", color=C_TERM,
         ms=3.5, lw=1.2, label="Core terminal")
axd.axhline(acc_ref, color=C_ACC, lw=1.2, ls="--", label="Accessory (41.5%)")
for _, r in fixed.iterrows():
    if r.MW_p > 0.05:
        axd.plot(r.Window_kb_per_end, r.terminal_mean, "o", color="white",
                 mec=C_TERM, ms=5.5, zorder=3)
axd.axvspan(80, 120, color="grey", alpha=0.12)
axd.text(165, 41.2, "P > 0.05\n(80-120 kb)", ha="left",
         fontsize=5.8, va="bottom")
axd.set_xlabel("Terminal window per end (kb)", fontsize=6.5)
axd.set_ylabel("Repeat content (%)", fontsize=7)
axd.tick_params(labelsize=6)
axd.set_title("d  Fixed-length window scan", fontsize=7.5, loc="left")
axd.legend(fontsize=5.8, frameon=False)
clean(axd)

# ---------------- e: stratification ----------------
sizes = ["small (<2 Mb)", "mid (2-4 Mb)", "large (>4 Mb)"]
slabels = ["small\n(<2 Mb)", "mid\n(2-4 Mb)", "large\n(>4 Mb)"]
xs = []
for i, sz in enumerate(sizes):
    for j, win in enumerate(["fixed100kb_terminal_pct", "pct3_terminal_pct"]):
        r = strat[(strat.Size_class == sz) & (strat.Window == win)]
        if len(r) == 0: continue
        x = i * 2 + j * 0.7
        ax_ = axe
        col = C_TERM if j == 0 else "#A9A9A9"
        axe.bar(x, r.terminal_mean.iloc[0], 0.6, color=col,
                edgecolor=("black" if j == 0 else "#666666"), lw=0.5,
                alpha=1.0 if j == 0 else 0.45)
        p = r.MW_p.iloc[0]
        star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
        if star:
            axe.text(x, r.terminal_mean.iloc[0] + 0.8, star, ha="center", fontsize=6)
        xs.append((x, sz, win))
axe.axhline(acc_ref, color=C_ACC, lw=1.4, ls="--", zorder=5)
axe.set_xticks([i * 2 + 0.35 for i in range(3)])
axe.set_xticklabels(slabels, fontsize=5.8)
axe.set_ylabel("Repeat content (%)", fontsize=7)
axe.tick_params(labelsize=6)
axe.set_title("e  Window behaviour by chromosome size", fontsize=7.5, loc="left")
from matplotlib.patches import Patch
axe.legend(handles=[Patch(color=C_TERM, label="Fixed 100 kb"),
                    Patch(facecolor="white", edgecolor="#666666", label="3% per end"),
                    Line2D([0], [0], color=C_ACC, ls="--", lw=1.2, label="Accessory")],
           fontsize=5.8, frameon=False)
clean(axe)

# ---------------- f: family x window ----------------
for fam, col, mark in [("DNA", C_TERM, "o"), ("LTR", "#377EB8", "s"),
                       ("LINE", "#4DAF4A", "^"), ("RC", "#984EA3", "v"),
                       ("Unknown", "#777777", "d")]:
    sub = famwin[famwin.Family == fam]
    axf.plot(sub.Window_kb, sub.terminal_mean, "-", color=col, lw=1.1, marker=mark,
             ms=2.8, label=fam)
    axf.axhline(sub.accessory_mean.iloc[0], color=col, lw=0.9, ls=":")
axf.set_xlabel("Terminal window per end (kb)", fontsize=6.5)
axf.set_ylabel("Repeat content (%)", fontsize=7)
axf.tick_params(labelsize=6)
axf.set_title("f  Per-family window scan (terminal vs accessory)", fontsize=7.5, loc="left")
axf.set_xlim(35, 245)
from matplotlib.lines import Line2D as L2
handles_f, labels_f = axf.get_legend_handles_labels()
handles_f.append(L2([0], [0], color="grey", lw=0.9, ls=":"))
labels_f.append("Accessory level (dotted)")
axf.legend(handles_f, labels_f, fontsize=5.2, frameon=False, ncol=1, loc="upper right",
           bbox_to_anchor=(0.995, 1.0))
clean(axf)

leg = [Line2D([0], [0], marker="s", color="w", markerfacecolor=C_TERM, markersize=8, label="Core terminal 100 kb"),
       Line2D([0], [0], marker="s", color="w", markerfacecolor=C_ACC, markersize=8, label="Accessory chromosomes")]
fig.legend(handles=leg, loc="lower center", ncol=2, fontsize=7, frameon=False,
           bbox_to_anchor=(0.5, -0.015))
fig.tight_layout(rect=[0, 0.02, 1, 1])
fig.savefig("06.accessory-chrs/repeat_terminal/figure_repeat_panels.pdf")
fig.savefig("06.accessory-chrs/repeat_terminal/figure_repeat_panels.png")
print("combined repeat panel figure saved")
