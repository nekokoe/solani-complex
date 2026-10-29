#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S24 v7 — final layout.

a: EVERY trait of every HOG plotted; each HOG appears on BOTH sides
   (28C traits on the left half, 37C traits on the right half), ranked by
   that temperature's association strength (strongest at the outer edges,
   weakest at the centre). x = rank, not genomic position.
   Colour strip below (28C-specific blue / 37C-specific red / shared grey /
   weak white) with its legend in the same row.
b: P<0.001 orthogroups that are 28C-specific (spec > 0.25) — both 28 C rho.
c: P<0.001 orthogroups that are 37C-specific (spec < -0.25) — both 37 C rho.
   Shade = inoculum concentration (light 1e3 to dark 1e5), x labels carry
   a short functional annotation.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(BASE)
OUT = "07.lethal_groups/genomewide_assoc"
C28, C37 = "#3B7DD8", "#C44E52"
SHADES28 = ["#C9E1F7", "#1F4E9C"]   # 1e3 (light) -> 1e4 (dark)
SHADES37 = ["#F7C9D0", "#A61E38"]

# Phenotypes are 1 - AUC at 10^3 and 10^4 plus the two-dose mean; the 10^5 dose is
# excluded because 16 of 21 isolates sit on the mortality ceiling there.
assoc = pd.read_csv(f"{OUT}/hog_phenotype_assoc_auc_no10e5.tsv", sep="\t")
T28 = ["vir_28Ce3", "vir_28Ce4", "vir_28_mean2"]
T37 = ["vir_37Ce3", "vir_37Ce4", "vir_37_mean2"]
C28C = T28[:2]
C37C = T37[:2]

piv_rho = assoc.pivot(index="HOG_ID", columns="trait", values="rho")
piv_p = assoc.pivot(index="HOG_ID", columns="trait", values="p")

def best_rho(r, p, ts):
    sub = r[ts].dropna()
    if sub.empty:
        return np.nan, np.nan
    i = sub.abs().idxmax()
    return sub[i], p[ts][i]

rows = {}
for h in piv_rho.index:
    r28, p28 = best_rho(piv_rho.loc[h], piv_p.loc[h], T28)
    r37, p37 = best_rho(piv_rho.loc[h], piv_p.loc[h], T37)
    rows[h] = (r28, p28, r37, p37)
df = pd.DataFrame.from_dict(rows, orient="index",
                            columns=["rho28", "p28", "rho37", "p37"])
df = df.dropna(subset=["rho28", "rho37"]).fillna({"p28": 1, "p37": 1})
df["spec"] = df.rho28.abs() - df.rho37.abs()

# functional short names for the significant HOGs
annot = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv",
                    sep="\t")
GENE_CLASS = {
    "VPS1": "Dynamin", "MDR1": "ABC transporter", "RSM10": "Ribosomal protein",
    "RPL16A": "Ribosomal protein", "AOX2": "Alternative oxidase",
    "NIT3": "Nitrilase", "SEC62": "ER translocation",
    "CDC54": "Helicase", "MRH4": "Helicase", "PMC1": "Ca2+ ATPase",
    "CDC55": "Protein phosphatase", "EPT1": "Phospholipid metabolism",
    "ARO4": "Aromatic AA biosynthesis", "EGD1": "GMC oxidoreductase",
    "SGS1": "Helicase", "BEM1": "Cell polarity scaffold",
    "GCV3": "Glycine cleavage", "CYM1": "Metallopeptidase",
    "DNM1": "Dynamin", "SOD2": "Mn superoxide dismutase",
    "FCR1": "Transcription factor", "HMT1": "Methyltransferase",
    "CDC4": "F-box/Skp1 complex",
}
CLASSES = [
    ("Ribosomal protein", ["ribosomal", "rpl", "rps", "rsm", "ribo"]),
    ("Helicase", ["helicase", "dead", "deah", "res_iii", "mcm"]),
    ("GTPase/dynamin", ["dynamin", "ged", "arp", "ras_"]),
    ("Transporter", ["mfs", "permease", "transp_", "transport", "cotransp",
                     "sugar_tr", "exchanger", "abc_tran", "abc2"]),
    ("Oxidoreductase", ["adh_", "p450", "oxido", "reductase", "dehydrogenase",
                        "monooxygenase", "nad_binding", "fad_binding", "pyr_redox",
                        "flavin", "ferredoxin", "cytochrom", "gmc_oxred", "aox"]),
    ("Transcription factor", ["fungal_trans", "zn_clus", "bzip", "homeobox"]),
    ("Cell cycle", ["cyclin", "cdc"]),
    ("Protein kinase", ["pkinase", "kinase", "hiska"]),
    ("Peptidase", ["peptidase", "protease", "trypsin"]),
    ("Methyltransferase", ["methyltransf", "methylase"]),
    ("Ankyrin repeat", ["ank"]),
    ("WD40/TPR repeat", ["wd40", "tpr"]),
    ("HET domain", ["^het"]),
    ("Mitochondrial", ["cox", "atp_synt", "mrpl", "mrp"]),
    ("RNA-binding", ["rrm", "rrp", "zf-cchc"]),
]
def domain_class(doms_str):
    doms = [d.strip().lower() for d in str(doms_str).split(";") if d.strip()]
    for name, kws in CLASSES:
        for d in doms:
            for kw in kws:
                if kw == "^het":
                    if d.startswith("het"):
                        return name
                elif kw in d:
                    return name
    return "Other function"

def short_func(h):
    r = annot[annot.HOG_ID == h]
    if len(r) == 0:
        return ("", "")
    gn = str(r.Gene_Name.iloc[0]).split(";")[0].strip()
    if not gn or gn.lower() == "nan" or gn == "-":
        gn = ""
    cls = GENE_CLASS.get(gn, "")
    if not cls:
        cls = domain_class(r.Pfam_Domain.iloc[0])
    return (cls, gn)
piv_rho[T37].abs().mean(axis=1)
mean_rho28 = piv_rho[T28].abs().mean(axis=1)
mean_rho37 = piv_rho[T37].abs().mean(axis=1)
side_L = mean_rho28 >= mean_rho37          # 28C stronger -> left, once only
left = df[side_L].copy()
right = df[~side_L].copy()
left["m"] = mean_rho28[left.index]
right["m"] = mean_rho37[right.index]
order28 = left.sort_values("m", ascending=False).index.tolist()
order37 = right.sort_values("m", ascending=False).index.tolist()
nL, nR = len(order28), len(order37)
N = len(df)
N2 = N / 2
pos28 = {h: i / nL * N2 for i, h in enumerate(order28)}          # 0 .. N/2
pos37 = {h: N2 + (1 - (i + 1) / nR) * N2 for i, h in enumerate(order37)}  # N/2 .. N

# ---- d: HOGs linked to the overall temperature effect (old fig d style) ----
# (removed: a stale TE_HOGS/te_carry block baked from the 15-strain run - its only
#  consumers were gone, and it carried hard-coded rho/P values that no longer
#  matched the data.)
fig = plt.figure(figsize=(11.5, 8.4), dpi=300)
gs = fig.add_gridspec(3, 2, height_ratios=[1.7, 0.28, 1.25], hspace=0.5, wspace=0.35)
ax1 = fig.add_subplot(gs[0, :])
axstrip = fig.add_subplot(gs[1, :])
axb = fig.add_subplot(gs[2, 0])
axc = fig.add_subplot(gs[2, 1])

alphas = [0.25, 0.50, 0.75, 0.9]   # 1e3, 1e4, 1e5, temperature mean
for j, cond in enumerate(T28):
    p = piv_p[cond].reindex(order28).values
    xs = np.array([pos28[h] for h in order28])
    ax1.scatter(xs, -np.log10(np.clip(p, 1e-8, 1)), s=3, c=C28, alpha=alphas[j],
                rasterized=True)
for j, cond in enumerate(T37):
    p = piv_p[cond].reindex(order37).values
    xs = np.array([pos37[h] for h in order37])
    ax1.scatter(xs, -np.log10(np.clip(p, 1e-8, 1)), s=3, c=C37, alpha=alphas[j],
                rasterized=True)
ax1.axhline(-np.log10(0.001), ls=":", color="grey", lw=0.8)
ax1.text(N * 0.985, -np.log10(0.001) + 0.12, "P = 0.001", fontsize=5.5,
         color="grey", ha="right")
ax1.axvline(N2, color="grey", lw=0.8, alpha=0.6)
ax1.set_xlim(0, N)
ax1.set_ylim(0, 5.2)
ax1.tick_params(labelbottom=False)
ax1.set_ylabel("-log10 P per trait")
ax1.set_title("a  Distribution of HOG\u2013temperature associations "
              "(x = rank, no genomic position)", fontsize=8, loc="left")
ax1.text(N * 0.25, 4.85, "28 \u00b0C (strongest at left)", ha="center",
         fontsize=6.5, color="#1F4E9C", bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
ax1.text(N * 0.75, 4.85, "37 \u00b0C (strongest at right)", ha="center",
         fontsize=6.5, color="#A61E38", bbox=dict(facecolor="white", edgecolor="none", pad=0.5))
ax1.text(N2, 4.75, "weakest", ha="center", fontsize=5.5, color="#888888",
         bbox=dict(facecolor="white", edgecolor="none", pad=0.5))

# ---------- strip ----------
for i, h in enumerate(order28):
    x0 = i / nL * N2
    s = df.loc[h, "spec"]
    m = max(abs(df.loc[h, "rho28"]), abs(df.loc[h, "rho37"]))
    c = C28 if s > 0.25 else ("#999999" if m >= 0.55 else "#EEEEEE")
    axstrip.add_patch(Rectangle((x0, 0.55), N2 / nL + 0.001, 0.45, facecolor=c,
                                edgecolor="none"))
for i, h in enumerate(order37):
    x0 = N2 + (1 - (i + 1) / nR) * N2
    s = df.loc[h, "spec"]
    m = max(abs(df.loc[h, "rho28"]), abs(df.loc[h, "rho37"]))
    c = C37 if s < -0.25 else ("#999999" if m >= 0.55 else "#EEEEEE")
    axstrip.add_patch(Rectangle((x0, 0.55), N2 / nR + 0.001, 0.45, facecolor=c,
                                edgecolor="none"))
axstrip.set_xlim(0, N)
axstrip.set_ylim(0, 1)
axstrip.set_yticks([])
axstrip.set_xticks([])
for s_ in ["top", "right", "left", "bottom"]:
    axstrip.spines[s_].set_visible(False)
l2 = [Line2D([0],[0], marker="s", color="w", markerfacecolor=C28, markersize=8,
             label="28 \u00b0C-specific"),
      Line2D([0],[0], marker="s", color="w", markerfacecolor=C37, markersize=8,
             label="37 \u00b0C-specific"),
      Line2D([0],[0], marker="s", color="w", markerfacecolor="#999999", markersize=8,
             label="shared association"),
      Line2D([0],[0], marker="s", color="w", markerfacecolor="#EEEEEE", markersize=8,
             label="weak")]
axstrip.legend(handles=l2, fontsize=6, ncol=4, loc="lower center",
               bbox_to_anchor=(0.5, 0.02))

# ---------- b / c: significant HOGs split by specificity ----------
sig = df[(df.p28 < 0.001) | (df.p37 < 0.001)].copy()
sig["diff"] = sig.rho28 - sig.rho37
b28 = sig[sig.spec >= 0].sort_values("diff", key=lambda c: c.abs(), ascending=False)
b37 = sig[sig.spec < 0].sort_values("diff", key=lambda c: c.abs(), ascending=False)

def plot_hog_panel(ax, panel_df, title, rot=0):
    sx = np.arange(len(panel_df))
    ax.plot([0, len(panel_df)], [0, 0], ls="--", color="grey", lw=0.8)
    for j, cond in enumerate(C28C):
        r = piv_rho[cond].reindex(panel_df.index).values
        ax.scatter(sx, r, s=24, c=SHADES28[j], marker="o", zorder=3)
    for j, cond in enumerate(C37C):
        r = piv_rho[cond].reindex(panel_df.index).values
        ax.scatter(sx, r, s=24, c=SHADES37[j], marker="^", zorder=3)
    labels = []
    for h in panel_df.index:
        cls, gn = short_func(h)
        ann = f"{cls} ({gn})" if gn else cls
        labels.append(f"{h.replace('N0.HOG','')}\n{ann}")
    ax.set_xticks(sx)
    ax.set_xticklabels(labels, rotation=rot, fontsize=5.0, ha="center" if rot == 0 else "right")
    ax.set_ylim(-1.1, 1.1)
    ax.set_ylabel("Spearman \u03c1 (copy vs trait)")
    ax.set_title(title, fontsize=8, loc="left")
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)

plot_hog_panel(axb, b28,
               f"b  P < 0.001 orthogroups, 28 \u00b0C-specific (n = {len(b28)})",
               rot=45)
plot_hog_panel(axc, b37,
               f"c  P < 0.001 orthogroups, 37 \u00b0C-specific (n = {len(b37)})",
               rot=60)

l3 = [Line2D([0],[0], marker="o", color="w", markerfacecolor="#6BA3D6", markersize=6,
             label="28 \u00b0C (o)"),
      Line2D([0],[0], marker="^", color="w", markerfacecolor="#E07B8F", markersize=6,
             label="37 \u00b0C (\u25b3)"),
      Line2D([0],[0], marker="o", color="w", markerfacecolor=SHADES28[0], markersize=5,
             label="1\u00d710\u00b3"),
      Line2D([0],[0], marker="o", color="w", markerfacecolor=SHADES28[1], markersize=5,
             label="1\u00d710\u2074")]
fig.legend(handles=l3, fontsize=6, ncol=4, loc="lower center",
            bbox_to_anchor=(0.5, 0.006))

fig.subplots_adjust(bottom=0.22, top=0.97)
fig.savefig(f"{OUT}/figure_gwa_s24.pdf")
fig.savefig(f"{OUT}/figure_gwa_s24.png", dpi=300)
print(f"S24 v7 saved: {N} HOGs; sig 28C={len(b28)} 37C={len(b37)}")
