#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Redesign S5/S6 as boxplot comparisons and S15 as a flanking-TE-pair
heatmap (user request 2026-08-12)."""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

C_TERM = "#E69F00"   # orange = core terminal 100 kb
C_ACC = "#C44E52"    # red = accessory
ISOS = ["BMU01933","BMU02028","BMU02474","BMU03267","BMU06673","BMU12683",
        "BMU12684","BMU12685","BMU12686","Fu3.1","Fu6.1","SB1","SX13"]
iso_set = set(ISOS)
SPECIES_ORDER = {"falciformis": 0, "solani": 1, "keratoplastica": 2}

# ================= S5: per-family boxplots =================
byf = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_by_family.tsv", sep="\t")
accf = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.accessory_repeat_by_family.tsv", sep="\t")
mwp = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.repeat_terminal_family_100kb.tsv", sep="\t")
p_map = dict(zip(mwp.Family, mwp.MW_p))

fams = ["DNA", "LTR", "LINE", "RC", "Simple_repeat", "Low_complexity", "Unknown", "Other"]
fig, ax = plt.subplots(figsize=(6.4, 2.9), dpi=300)
xs, labels, max_whisk = [], [], 0
for i, fam in enumerate(fams):
    core = byf[(byf.Family == fam) & (byf.Type == "core")].Terminal100kb_pct
    acc = accf[accf.Family == fam].Accessory_pct
    if len(core) == 0 or len(acc) == 0:
        continue
    x0 = i * 2
    bp1 = ax.boxplot([core, acc], positions=[x0, x0 + 0.75], widths=0.62,
                     patch_artist=True, showfliers=False, medianprops=dict(color="black", lw=0.9))
    for patch, c in zip(bp1["boxes"], [C_TERM, C_ACC]):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("black"); patch.set_linewidth(0.6)
    for w, c in zip(bp1["whiskers"], [C_TERM, C_ACC, C_TERM, C_ACC]):
        w.set_color(c); w.set_linewidth(0.7)
    for cap, c in zip(bp1["caps"], [C_TERM, C_ACC, C_TERM, C_ACC]):
        cap.set_color(c); cap.set_linewidth(0.7)
    # whisker top = max of the two boxes' upper whisker (exclude outliers already)
    whisk_top = max(w.get_ydata()[1] for w in bp1["whiskers"])
    max_whisk = max(max_whisk, whisk_top)
    p = p_map.get(fam, 1)
    if p < 0.001:
        star = "***"
    elif p < 0.01:
        star = "**"
    elif p < 0.05:
        star = "*"
    else:
        star = "n.s."
    ax.annotate("", xy=(x0 + 0.75, whisk_top + 1.2), xytext=(x0, whisk_top + 1.2),
                arrowprops=dict(arrowstyle="-", lw=0.5, color="grey"))
    ax.text(x0 + 0.375, whisk_top + 1.8, star, ha="center", fontsize=6.5)
    xs.append(x0 + 0.375); labels.append(fam)
ax.set_xticks(xs)
ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=6.5)
ax.set_ylim(0, max_whisk * 1.22)
ax.set_ylabel("Repeat content (%)", fontsize=7.5)
ax.tick_params(labelsize=6.5)
from matplotlib.lines import Line2D
leg = [Line2D([0], [0], marker="s", color="w", markerfacecolor=C_TERM, markersize=7, label="Core terminal 100 kb"),
       Line2D([0], [0], marker="s", color="w", markerfacecolor=C_ACC, markersize=7, label="Accessory chromosomes")]
ax.legend(handles=leg, fontsize=6, frameon=False, loc="upper right")
ax.set_title("Per-family repeat content: core terminal 100 kb vs accessory", fontsize=7.8)
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("06.accessory-chrs/repeat_terminal/figure_repeat_terminal.families_100kb.pdf")
print("S5 boxplot saved")

# ================= S6: per-isolate boxplots =================
cont = pd.read_csv("06.accessory-chrs/repeat_terminal/suppl_table.chromosome_repeat_content.tsv", sep="\t")
sp = {}
for line in open("06.accessory-chrs/tgs_assemblies/species.211.map.tsv"):
    p = line.rstrip("\n").split("\t")
    if len(p) >= 2: sp[p[0]] = p[1]
order = sorted(ISOS, key=lambda i: (SPECIES_ORDER.get(sp.get(i, "solani"), 9), i))

fig, ax = plt.subplots(figsize=(6.8, 3.0), dpi=300)
xs, labels = [], []
for i, iso in enumerate(order):
    sub = cont[cont.Isolate == iso]
    core = sub[(sub.Type == "core")].Terminal100kb_pct
    acc = sub[(sub.Type == "accessory")].Repeat_pct
    if len(core) == 0 or len(acc) == 0:
        continue
    x0 = i * 2
    bp1 = ax.boxplot([core, acc], positions=[x0, x0 + 0.75], widths=0.62,
                     patch_artist=True, showfliers=False, medianprops=dict(color="black", lw=0.9))
    for patch, c in zip(bp1["boxes"], [C_TERM, C_ACC]):
        patch.set_facecolor(c); patch.set_alpha(0.75); patch.set_edgecolor("black"); patch.set_linewidth(0.6)
    for w, c in zip(bp1["whiskers"], [C_TERM, C_ACC, C_TERM, C_ACC]):
        w.set_color(c); w.set_linewidth(0.7)
    for cap, c in zip(bp1["caps"], [C_TERM, C_ACC, C_TERM, C_ACC]):
        cap.set_color(c); cap.set_linewidth(0.7)
    xs.append(x0 + 0.375); labels.append(iso.replace("BMU", ""))
ax.set_xticks(xs)
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.8)
ax.set_ylabel("Repeat content (%)", fontsize=7.5)
ax.tick_params(labelsize=6.5)
ax.set_title("Per-isolate repeat content: core terminal 100 kb vs accessory", fontsize=7.8)
ax.legend(handles=leg, fontsize=6, frameon=False, loc="upper right")
for s in ["top", "right"]:
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig("06.accessory-chrs/repeat_terminal/figure_repeat_terminal.paired.pdf")
print("S6 boxplot saved")

# ================= S15: flanking-TE-pair heatmap =================
aff, gpos = {}, {}
for iso in ISOS:
    with open(f"02.annotation/per_isolates/{iso}/funannotate/predict_results/fusarium.gff3") as f:
        for line in f:
            if line.startswith("#"): continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "mRNA": continue
            tid = None
            for a in p[8].split(";"):
                if a.startswith("ID="): tid = a[3:]; break
            if tid:
                aff[tid] = p[0]
                gpos[f"{iso}|{tid}"] = (int(p[3]), int(p[4]))
ct = {}
for iso in ISOS:
    acc = set(open(f"06.accessory-chrs/tgs_assemblies/{iso}.accessory_chrs.id").read().split())
    for line in open(f"02.annotation/per_isolates/{iso}/repeatmasker/{iso}.assembly.fasta.masked.length"):
        p = line.split()
        if len(p) >= 2:
            c, L = p[0], int(p[1])
            ct[c] = "acc" if c in acc else ("core" if L >= 500000 else "scaf")
reps = {}
for iso in ISOS:
    rep = defaultdict(list)
    with open(f"06.accessory-chrs/tgs_assemblies/{iso}.assembly.fasta.out") as f:
        for i, line in enumerate(f):
            if i < 3: continue
            p = line.split()
            if not p or len(p) < 11: continue
            if len(p) > 15 and p[15] == "*": continue
            rep[p[4]].append((int(p[5]), int(p[6]), p[10]))
    reps[iso] = {c: sorted(v) for c, v in rep.items()}

def flank_pair(iso, contig, gs, ge, win=5000):
    iv = reps[iso].get(contig, [])
    left = [fam for s, e, fam in iv if e < gs and e >= gs - win]
    right = [fam for s, e, fam in iv if s > ge and s <= ge + win]
    return (left[-1] if left else None, right[0] if right else None)

NON_TE = {"Unknown", "Simple_repeat", "Low_complexity"}
rep_members = defaultdict(list)
with open("06.accessory-chrs/tgs_assemblies/all_solaniSC.cluster100_cluster.tsv") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 2: continue
        iso, _, trans = p[1].partition("|")
        if iso in iso_set and trans in aff:
            rep_members[p[0]].append(p[1])

pair_counts, fam_counts = Counter(), Counter()
for rep_id, members in rep_members.items():
    contigs = {aff[m.split("|")[1]] for m in members}
    if len(contigs) < 2: continue
    ctypes = {ct.get(c) for c in contigs if c}
    if "acc" not in ctypes: continue
    for m in members:
        iso = m.split("|")[0]; tid = m.split("|")[1]
        pos = gpos.get(m)
        if pos is None: continue
        lf, rf = flank_pair(iso, aff[tid], pos[0], pos[1])
        if not lf or not rf: continue
        if lf in NON_TE or rf in NON_TE: continue
        pair_counts[(lf, rf)] += 1
        fam_counts[lf] += 1; fam_counts[rf] += 1

top_fams = [f for f, _ in fam_counts.most_common(12)]
mat = np.zeros((12, 12), dtype=int)

for (l, r), c in pair_counts.items():
    if l in top_fams and r in top_fams:
        mat[top_fams.index(l), top_fams.index(r)] = c

import pandas as _pd
_pd.DataFrame(mat, index=top_fams, columns=top_fams).to_csv(
    "06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.tsv", sep="\t")

# Denominator for the percentage version of this panel. It has to be the number of
# TE-BRACKETED genes (both flanks a repeat family other than Unknown/Simple_repeat/
# Low_complexity), NOT sum(mat) -- the matrix only keeps pairs whose two families are
# both inside the top 12, so sum(mat) is a subset of the bracketed genes.
_pd.DataFrame([{
    "te_bracketed_genes": int(sum(pair_counts.values())),
    "in_top12_grid": int(mat.sum()),
    "n_top_families": len(top_fams),
    "n_distinct_pairs": len(pair_counts),
}]).to_csv("06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.total.tsv",
           sep="\t", index=False)
print("te-bracketed genes:", int(sum(pair_counts.values())),
      "| in top-12 grid:", int(mat.sum()), "| distinct pairs:", len(pair_counts))

fig, ax = plt.subplots(figsize=(5.8, 4.6), dpi=300)
im = ax.imshow(mat + 1, norm=LogNorm(vmin=1, vmax=mat.max() + 1), cmap="YlOrRd")
ax.set_xticks(range(12)); ax.set_yticks(range(12))
ax.set_xticklabels([f.replace("DNA/", "").replace("RC/", "").replace("LINE/", "").replace("LTR/", "")
                    for f in top_fams], rotation=45, ha="right", fontsize=6)
ax.set_yticklabels([f.replace("DNA/", "").replace("RC/", "").replace("LINE/", "").replace("LTR/", "")
                    for f in top_fams], fontsize=6)
for i in range(12):
    for j in range(12):
        v = mat[i, j]
        if v >= 3:
            ax.text(j, i, str(v), ha="center", va="center", fontsize=5.5,
                    color="white" if v >= mat.max() * 0.45 else "black")
ax.set_xlabel("Downstream flanking TE family", fontsize=7.5)
ax.set_ylabel("Upstream flanking TE family", fontsize=7.5)
ax.set_title("Genes bracketed by TE pairs (upstream x downstream, 5 kb)",
             fontsize=7.8)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03,
                      ticks=[1, 3, 7, 14, 28])
cbar.ax.set_yticklabels(["1", "3", "7", "14", "28"])
cbar.ax.tick_params(labelsize=5.5)
cbar.set_label("Gene count", fontsize=6.5)
fig.tight_layout()
fig.savefig("06.accessory-chrs/gene_transfer/figure_transfer_cassette.pdf")
print("S15 heatmap saved; matrix diag:", [mat[i, i] for i in range(12)])
