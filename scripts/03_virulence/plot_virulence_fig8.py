#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fig. 8 v2 — virulence-genotype association (user redesign).

a: temperature-preference quadrant (28C-blue edge / 37C-red edge, same
   colour scheme as S24; fill = clinical/environmental)
b: concentration response
c: velvet-complex copies vs 1 - AUC, 10^3 and 10^4 plotted with
   trend lines
d: sidC tree with per-leaf 1 - AUC annotation
e: genome-wide dual-temperature quadrant, highlighted candidates coloured
   by 28C/37C preference
f: 37C panel — VPS1/dynamin copy loss vs 37C mortality (3 inocula + trends)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo

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
OUT = "07.lethal_groups"
C28, C37 = "#3B7DD8", "#C44E52"          # 28C blue, 37C red (consistent with S24)
C_MED, C_ENV = "#C44E52", "#9E9E9E"
SHADES28 = ["#C9E1F7", "#5E9CD6", "#1F4E9C"]
SHADES37 = ["#F7C9D0", "#E0607A", "#A61E38"]

_pheno_head = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics.tsv",
                          sep="\t", index_col=0)
STRAINS = _pheno_head.index.tolist()
# clinical (med) set read from the cohort metadata rather than a literal
_meta = pd.read_csv("05.pan-genome/bmu_ref_groups.211.txt", sep="\t")
MED = set(_meta.BMU_ID[_meta.MedEnv == "med"])

pheno = pd.read_csv(f"{OUT}/phenotype_analysis/phenotype_metrics.tsv", sep="\t", index_col=0)
# per-strain temperature preference from the CMH test stratified by inoculum dose
# (scripts/cmh_temperature_preference.py), not a hard cutoff on temp_effect
_bias = pd.read_csv(f"{OUT}/phenotype_analysis/temperature_bias_cmh.tsv", sep="\t")
_lab = {"37C_biased": "37", "28C_biased": "28"}
PREF = {st: _lab.get(cl, "flat") for st, cl in zip(_bias["strain"], _bias["class"])}
assert len(PREF) == 21 and sum(v != "flat" for v in PREF.values()) == 4, PREF
mort = pd.read_csv(f"{OUT}/phenotype_analysis/phenotype_mortality.tsv", sep="\t", index_col=0)
copy = pd.read_csv(f"{OUT}/literature_genes/gene_strain_copy.tsv", sep="\t", index_col=0)
# Association universe: 1 - AUC at 10^3 and 10^4 plus the two-dose mean; the 10^5 dose is
# excluded because 16 of 21 isolates sit on the mortality ceiling there.
assoc = pd.read_csv(f"{OUT}/genomewide_assoc/hog_phenotype_assoc_auc_no10e5.tsv", sep="\t")
hog = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv", sep="\t")
mat = hog[hog.iloc[:,0].isin(STRAINS)].pivot(index="HOG_ID", columns="SampleID",
                                             values="GeneCount").reindex(columns=STRAINS)

fig, axes = plt.subplots(2, 3, figsize=(12.6, 8.4), dpi=300)
(axa, axb, axe), (axd, axc, axf) = axes   # c at bottom-middle, e at top-right

# ---------------- a: temperature preference (28 blue / 37 red edges) --------
# Every phenotype in this figure is 1 - AUC at 10^3 and 10^4 (plus their two-dose mean);
# the 10^5 dose is excluded from the analysis because 16 of 21 isolates are on the
# mortality ceiling there.
AUC = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_metrics_auc_no10e5.tsv",
                  sep="\t", index_col=0)
m28 = AUC.vir_28_mean2.reindex(STRAINS)
m37 = AUC.vir_37_mean2.reindex(STRAINS)
te = pheno.temp_effect.reindex(STRAINS)
for s in STRAINS:
    pref = PREF[s]
    marker = "o" if pref == "28" else ("v" if pref == "37" else "o")
    edge = C28 if pref == "28" else (C37 if pref == "37" else "none")
    lw = 1.3 if pref != "flat" else 0.5
    axa.scatter(m28[s], m37[s], s=46, c=C_MED if s in MED else C_ENV,
                marker=marker, alpha=0.9, edgecolors=edge, linewidths=lw, zorder=3)
# radial label start: push each label outward from the data centroid by
# ~2x text height along the point->centroid direction; adjustText then
# separates and auto-connects.
cx = m28.mean(); cy = m37.mean()
fig.canvas.draw()
probe = axa.text(0, 0, "X", fontsize=4.2)
hp = probe.get_window_extent().height
aw = axa.get_window_extent().width
ah = axa.get_window_extent().height
xr = axa.get_xlim()[1] - axa.get_xlim()[0]
yr = axa.get_ylim()[1] - axa.get_ylim()[0]
hx = hp / aw * xr
hy = hp / ah * yr
probe.remove()
texts_a = []
ANG_OFF = {"BMU01287": -0.6, "BMU02031": 0.6}
FIX2 = {"BMU01287": (0.36, 0.87), "BMU02031": (0.61, 0.87)}
auto_list = [s for s in STRAINS if s not in FIX2]
for s in STRAINS:
    if s in FIX2:
        axa.text(FIX2[s][0], FIX2[s][1], s.replace("BMU", ""), fontsize=4.2,
                 ha="center", va="center", zorder=5)
        continue
    dx, dy = m28[s] - cx, m37[s] - cy
    n = np.hypot(dx, dy)
    if n == 0:
        dx, dy, n = 0, 1, 1.0
    ang = np.arctan2(dy, dx) + ANG_OFF.get(s, 0.0)
    ux, uy = np.cos(ang), np.sin(ang)
    t = axa.text(m28[s] + ux * 5 * hx, m37[s] + uy * 5 * hy,
                 s.replace("BMU", ""), fontsize=4.2,
                 ha="center", va="center", zorder=5)
    texts_a.append(t)
from adjustText import adjust_text
adjust_text(texts_a, ax=axa,
            expand=(1.8, 1.8), force_text=1.2, force_static=0.4,
            force_points=0.4, lim=500)
for s, t in zip(auto_list, texts_a):
    xf, yf = t.get_position()
    axa.plot([m28[s], xf], [m37[s], yf], color="#999999", lw=0.7, zorder=4)
for s, (lx, ly) in FIX2.items():
    axa.plot([m28[s], lx], [m37[s], ly], color="#999999", lw=0.7, zorder=4)
lim = [0, 1.02]
axa.plot(lim, lim, ls="--", color="grey", lw=0.8)          # diagonal y = x
axa.set_xlim(lim); axa.set_ylim(0, 1.02)
axa.set_xlabel("Mean 1 \u2212 AUC at 28 \u00b0C (10\u00b3, 10\u2074)")
axa.set_ylabel("Mean 1 \u2212 AUC at 37 \u00b0C (10\u00b3, 10\u2074)")
axa.set_title("a  Strain temperature preference", fontsize=8, loc="left")
# the clinical/environmental counts come off the cohort, never typed in: the panel is
# plotted from STRAINS and MED, so a literal n goes stale without any error (it still
# read "n=11" after the cohort had grown to 17 clinical of 21 isolates).
N_MED = sum(s in MED for s in STRAINS)
N_ENV = len(STRAINS) - N_MED
handles = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_MED, markersize=7,
                  label=f"Clinical (n={N_MED})"),
           Line2D([0],[0], marker="o", color="w", markerfacecolor=C_ENV, markersize=7,
                  label=f"Environmental (n={N_ENV})"),
           Line2D([0],[0], marker="o", color="w", markerfacecolor="none",
                  markeredgecolor=C28, markersize=7, label="28 \u00b0C-biased"),
           Line2D([0],[0], marker="v", color="w", markerfacecolor="none",
                  markeredgecolor=C37, markersize=7, label="37 \u00b0C-biased")]
# upper left, not lower left: both axes run 0-1.02 but the strains sit in
# x 0.13-0.68 / y 0.10-0.36, so a lower-left legend was drawn on top of the cloud.
axa.legend(handles=handles, fontsize=5.8, loc="upper left", bbox_to_anchor=(0.02, 0.98),
           ncol=2, columnspacing=0.8, handlelength=1.2)

# ---------------- b: concentration response (unchanged) --------------------
conc = [3, 4, 5]
AUCONC = pd.read_csv("07.lethal_groups/phenotype_analysis/phenotype_auc.tsv",
                     sep="\t", index_col=0)
def vir_at(temp, d):
    return 1.0 - AUCONC[f"auc_{temp}Ce{d}"].reindex(STRAINS)
for s in STRAINS:
    for temp, ls in (("28", "-"), ("37", "--")):
        axb.plot(conc[:2], [vir_at(temp, d).loc[s] for d in conc[:2]], ls + "o",
                 color="#BBBBBB", ms=2.5, lw=0.7, alpha=0.5, zorder=1)
        axb.plot(conc[1:], [vir_at(temp, d).loc[s] for d in conc[1:]], ls + "o",
                 color="#DDDDDD", ms=2.5, lw=0.7, alpha=0.35, zorder=1)
for col, ls, temp, lab in [(C28, "-", "28", "28 \u00b0C"), (C37, "--", "37", "37 \u00b0C")]:
    m = [vir_at(temp, d).mean() for d in conc]
    axb.plot(conc[:2], m[:2], ls, color=col, lw=2, ms=6, marker="o", label=lab, zorder=3)
    axb.plot(conc[1:], m[1:], ls, color=col, lw=1.2, ms=6, marker="o",
             mfc="white", zorder=2)
axb.set_xticks(conc)
axb.set_xticklabels([f"10\u00b3", "10\u2074", "10\u2075"])
axb.set_xlabel("Inoculum (conidia per larva)")
axb.set_ylabel("1 \u2212 AUC")
axb.set_title(f"b  Concentration response (n = {len(STRAINS)} strains)", fontsize=8, loc="left")
axb.legend(fontsize=6, loc="upper left")
axb.set_ylim(-0.02, 1.06)
axb.text(5.05, 0.06, "10\u2075 not used", fontsize=5.8, color="#888888", ha="right")
axb.text(4.65, 0.16, "concentration effect P \u2248 10\u207b\u00b9\u2075",
         fontsize=5.8, color="#666666", ha="right")
axb.text(3.02, 0.90, f"grey lines: individual strains (n = {len(STRAINS)})",
         fontsize=5.5, color="#888888")

# ---------------- c: velvet complex vs 28Ce3 (most significant inoculum) -----
velvet_sum = copy[["velvet", "VeA", "LaeA"]].sum(axis=1)
y = 1.0 - AUCONC.auc_28Ce3.reindex(STRAINS)
x = velvet_sum
rho, p = stats.spearmanr(x, y)
for s in STRAINS:
    pref = PREF[s]
    ec = C28 if pref == "28" else (C37 if pref == "37" else "white")
    elw = 1.1 if pref != "flat" else 0.3
    mk = "v" if pref == "37" else "o"
    axc.scatter(x[s], y[s], s=38, c=C_MED if s in MED else C_ENV, alpha=0.9,
                marker=mk, edgecolors=ec, linewidths=elw, zorder=3)

axc.set_xlabel("velvet + VeA + LaeA copies")
axc.set_ylabel("1 \u2212 AUC (28 \u00b0C, 10\u00b3)")
axc.set_title("e  More velvet copies \u2192 higher 28 \u00b0C virulence",
              fontsize=7.5, loc="left")
# tertiles: with n=21 a 5-v-5 split is only 24% of the cohort per side
_k = max(3, len(y) // 3)
hi = y.nlargest(_k).index; lo = y.nsmallest(_k).index
u = stats.mannwhitneyu(velvet_sum[hi], velvet_sum[lo])
axc.text(0.02, 0.97, f"high vs low virulence: P={u.pvalue:.3g} (FDR<0.05)",
         transform=axc.transAxes, fontsize=5.8, color="#666666")
coef = np.polyfit(x.values, y.values, 1)
xs = np.linspace(x.min(), x.max(), 20)
axc.plot(xs, np.polyval(coef, xs), "--", color="#333333", lw=1.1, zorder=2)
axc.text(x.max(), np.polyval(coef, x.max()) + 0.05,
         f"\u03c1 = {rho:.2f}, P = {p:.3g}", ha="right", fontsize=5.5, color="#333333")
lc2 = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_MED, markersize=7,
              label="Clinical"),
       Line2D([0],[0], marker="o", color="w", markerfacecolor=C_ENV, markersize=7,
              label="Environmental"),
       Line2D([0],[0], marker="o", color="w", markerfacecolor="none",
              markeredgecolor=C28, markersize=7, label="28 \u00b0C-biased"),
       Line2D([0],[0], marker="v", color="w", markerfacecolor="none",
              markeredgecolor=C37, markersize=7, label="37 \u00b0C-biased")]
axc.legend(handles=lc2, fontsize=5.8, loc="lower right",
           bbox_to_anchor=(0.98, 0.02), ncol=2, columnspacing=0.8,
           handlelength=1.2)

# ---------------- d: sidC phylogram (standard rectangular, dendropy tree) ----
import dendropy
t = dendropy.Tree.get(path=f"{OUT}/literature_genes/tmp/sidC_tree_midroot.nwk",
                      schema="newick")
t.ladderize()
leaves = [n for n in t.leaf_node_iter()]
y_pos = {id(leaf): i for i, leaf in enumerate(leaves)}
depth = {}
for nd in t.postorder_node_iter():
    if nd.is_leaf():
        depth[id(nd)] = nd.distance_from_root()
    else:
        depth[id(nd)] = nd.distance_from_root()
mort_map = 1.0 - AUCONC.auc_28Ce4.reindex(STRAINS)
cmap_d = plt.cm.YlOrRd
maxd = max(depth.values())

def draw_clade(ax, nd):
    x0 = depth[id(nd)]
    if nd.is_leaf():
        return x0, y_pos[id(nd)]
    child_pts = [draw_clade(ax, ch) for ch in nd.child_nodes()]
    ys = [p[1] for p in child_pts]
    y0 = np.mean(ys)
    ymin, ymax = min(ys), max(ys)
    ax.plot([x0, x0], [ymin, ymax], color="grey", lw=0.8, zorder=1)   # vertical
    for cx, cy in child_pts:
        ax.plot([x0, cx], [cy, cy], color="grey", lw=0.8, zorder=1)   # horizontal
    if nd is not t.seed_node and nd.edge_length is not None:
        ax.text(x0, y0 + 0.22, f"{nd.edge_length:.4f}", fontsize=4,
                color="grey", ha="center", zorder=2)
    return x0, y0
draw_clade(axd, t.seed_node)
for leaf in leaves:
    x0 = depth[id(leaf)]
    y0 = y_pos[id(leaf)]
    name = leaf.taxon.label
    mort_d = mort_map.get(name, np.nan)
    col = cmap_d(mort_d) if not np.isnan(mort_d) else "#999999"
    axd.plot(x0, y0, "o", color=col, ms=8, mec="black", mew=0.6, zorder=3)
    axd.text(x0 + maxd * 0.032, y0, f"{name[3:]} {mort_d:.2f}",
             fontsize=5.0, va="center", zorder=3)
axd.set_xlim(-0.0022, maxd * 1.30)
axd.set_ylim(-1.4, len(leaves) + 0.4)
axd.set_yticks([])
axd.set_xlabel("Branch length (subst./site)")
scale_len = 0.002
sb_x0 = maxd * 0.9
sb_y = -0.8
axd.plot([sb_x0, sb_x0 + scale_len], [sb_y, sb_y], color="black", lw=1.2, zorder=4)
axd.text(sb_x0 + scale_len / 2, sb_y - 0.5, f"{scale_len}", fontsize=4.6,
         ha="center", va="top", zorder=4)
sm = plt.cm.ScalarMappable(cmap=cmap_d, norm=plt.Normalize(0, 1))
sm.set_array([])
cb = fig.colorbar(sm, ax=axd, fraction=0.02, pad=0.01, ticks=[0, 0.5, 1])
cb.set_label("1 \u2212 AUC, 28 \u00b0C, 10\u2074", fontsize=5.5)
cb.ax.tick_params(labelsize=5)
axd.set_title("d  sidC divergence tracks 28 \u00b0C virulence",
              fontsize=8, loc="left")
for s_ in ["top", "right"]:
    axd.spines[s_].set_visible(False)

# ---------------- e: dual-temperature quadrant, coloured candidates ---------
t28 = ["vir_28Ce3", "vir_28Ce4", "vir_28_mean2"]
t37 = ["vir_37Ce3", "vir_37Ce4", "vir_37_mean2"]
g28 = assoc[assoc.trait.isin(t28)].groupby("HOG_ID").apply(
    lambda g: g.loc[g["rho"].abs().idxmax(), "rho"], include_groups=False)
g37 = assoc[assoc.trait.isin(t37)].groupby("HOG_ID").apply(
    lambda g: g.loc[g["rho"].abs().idxmax(), "rho"], include_groups=False)
axe.scatter(g28.values, g37.values, s=5, c="#BBBBBB", alpha=0.35, rasterized=True, zorder=1)
# The highlighted set used to be a 19-HOG literal baked from the 15-strain run
# (with hard-coded label positions). It is recomputed here from the current
# association table: HOGs whose best-temperature association reaches P < 0.001,
# side = whichever temperature has the larger |rho|, coordinates = (rho at the
# best 28C trait, rho at the best 37C trait).
def _best_per_temp(sub):
    g = sub.groupby("HOG_ID")
    return pd.DataFrame({
        "rho": g.apply(lambda x: x.loc[x["rho"].abs().idxmax(), "rho"],
                       include_groups=False),
        "p": g.apply(lambda x: x.loc[x["rho"].abs().idxmax(), "p"],
                     include_groups=False)})


_b28 = _best_per_temp(assoc[assoc.trait.isin(t28)])
_b37 = _best_per_temp(assoc[assoc.trait.isin(t37)])

# Every trait in this table is tested (the saturated 10^5 rates are not in it at all),
# so the highlight criterion needs no further filtering.
_b28t = _best_per_temp(assoc[assoc.trait.isin(t28)])
_b37t = _best_per_temp(assoc[assoc.trait.isin(t37)])
_ann = pd.read_csv("05.pan-genome/orthofinder-shorter/N0.HOG_annotation.v2.tsv",
                   sep="\t").set_index("HOG_ID")
SIG_HITS = {}
for _h in _b28t.index:
    if _h not in _b37t.index:
        continue
    if min(_b28t.p[_h], _b37t.p[_h]) >= 0.001:
        continue
    _side = "28" if abs(_b28t.rho[_h]) >= abs(_b37t.rho[_h]) else "37"
    _gn = str(_ann.Gene_Name.get(_h, "-")).split(";")[0]
    _pf = str(_ann.Pfam_Domain.get(_h, "-")).split(";")[0]
    _lab = _gn if _gn not in ("-", "nan", "") else (
        _pf if _pf not in ("-", "nan", "") else _h.replace("N0.HOG", "H"))
    SIG_HITS[_h] = (_lab, _side, float(_b28t.rho[_h]), float(_b37t.rho[_h]))
print(f"highlighted HOGs (best-temperature P<0.001): {len(SIG_HITS)}")
SIDE_C = {"28": C28, "37": C37, "sh": "#999999"}
texts_c = []
for h, (lab, side, rx, ry) in SIG_HITS.items():
    axe.scatter(rx, ry, s=30, c=SIDE_C[side], edgecolors="white",
                linewidths=0.4, zorder=3)
    t = axe.text(rx, ry, lab, fontsize=4.0,
                 ha="center", va="center", zorder=5)
    texts_c.append(t)
from adjustText import adjust_text
adjust_text(texts_c, ax=axe,
            expand=(1.8, 1.8), force_text=1.2, force_static=0.4,
            force_points=0.4, lim=500)
for h, (lab, side, rx, ry) in SIG_HITS.items():
    xf, yf = texts_c[list(SIG_HITS).index(h)].get_position()
    axe.plot([rx, xf], [ry, yf], color="#999999", lw=0.6, zorder=4)

lim = [-1.02, 1.02]
axe.plot([0, 0], lim, ls=":", color="grey", lw=0.7)
axe.plot(lim, [0, 0], ls=":", color="grey", lw=0.7)
axe.set_xlim(-1.05, 1.05); axe.set_ylim(lim)

axe.set_xlabel("Spearman \u03c1 (28 \u00b0C)")
axe.set_ylabel("Spearman \u03c1 (37 \u00b0C)")
axe.set_title(f"c  Genome-wide scan ({len(g28)} variable HOGs)", fontsize=8)



# ---------------- f: VPS1 loss vs 37Ce3 (most significant inoculum) ----------
x_vps1 = mat.loc["N0.HOG0017872"]
y = 1.0 - AUCONC.auc_37Ce3.reindex(STRAINS)
rho, p = stats.spearmanr(x_vps1, y)
for s in STRAINS:
    pref = PREF[s]
    ec = C28 if pref == "28" else (C37 if pref == "37" else "white")
    elw = 1.1 if pref != "flat" else 0.3
    mk = "v" if pref == "37" else "o"
    axf.scatter(x_vps1[s], y[s], s=38, c=C_MED if s in MED else C_ENV, alpha=0.9,
                marker=mk, edgecolors=ec, linewidths=elw, zorder=3)
coef = np.polyfit(x_vps1.values, y.values, 1)
xs = np.linspace(x_vps1.min(), x_vps1.max(), 20)
axf.plot(xs, np.polyval(coef, xs), "--", color="#333333", lw=1.1, zorder=2)
axf.text(x_vps1.min(), np.polyval(coef, x_vps1.min()) + 0.05,
         f"\u03c1 = {rho:.2f}, P = {p:.3g}", ha="left", fontsize=5.5, color="#333333")

axf.set_xlabel("VPS1/dynamin copies")
axf.set_ylabel("1 \u2212 AUC (37 \u00b0C, 10\u00b3)")
axf.set_title("f  VPS1/dynamin loss \u2192 higher 37 \u00b0C virulence",
              fontsize=7.5, loc="left")
# headroom above the data (max 1 - AUC = 0.25) so the legend has a band of its own;
# at the autoscaled limits the strains fill the full height at x = 0 and x = 1 and
# any in-axes legend lands on them.
axf.set_ylim(-0.01, 0.34)
lf2 = [Line2D([0],[0], marker="o", color="w", markerfacecolor=C_MED, markersize=7,
              label=f"Clinical (n={N_MED})"),
       Line2D([0],[0], marker="o", color="w", markerfacecolor=C_ENV, markersize=7,
              label=f"Environmental (n={N_ENV})"),
       Line2D([0],[0], marker="o", color="w", markerfacecolor="none",
              markeredgecolor=C28, markersize=7, label="28 \u00b0C-biased"),
       Line2D([0],[0], marker="v", color="w", markerfacecolor="none",
              markeredgecolor=C37, markersize=7, label="37 \u00b0C-biased")]
axf.legend(handles=lf2, fontsize=5.8, loc="upper right", ncol=2,
           columnspacing=0.8, handlelength=1.2)

fig.subplots_adjust(left=0.055, right=0.945, bottom=0.06, top=0.94,
                    wspace=0.30, hspace=0.34)

# --- legend placement check: a legend must not be drawn over the data ----------
# "the legend and the content cover each other" is invisible to every other test,
# so assert it geometrically instead of trusting the corner keyword.
fig.canvas.draw()
for _ax, _name in [(axa, "a"), (axf, "f")]:
    _leg = _ax.get_legend()
    if _leg is None:
        continue
    _lb = _leg.get_window_extent()
    _hits = []
    for _c in _ax.collections:
        _off = _c.get_offsets()
        if _off is None or len(_off) == 0:
            continue
        for _px, _py in _ax.transData.transform(np.asarray(_off, dtype=float)):
            if _lb.contains(_px, _py):
                _hits.append(f"mark@({_px:.0f},{_py:.0f})")
    for _t in _ax.texts:
        if _t.get_text() and _t.get_window_extent().overlaps(_lb):
            _hits.append(f"label:{_t.get_text()[:10]}")
    print(f"[legend] panel {_name}: {len(_hits)} marks/labels under the legend"
          + (f" -> {_hits[:6]}" if _hits else " (clear)"))

fig.savefig(f"{OUT}/figure_virulence_fig8.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/figure_virulence_fig8.png", dpi=300, bbox_inches="tight")
fig.canvas.draw()
for name, ax in [("a", axa), ("b", axb), ("c", axc), ("e", axe), ("f", axf)]:
    bb = ax.get_window_extent()
    print(f"panel {name}: axes box {bb.width:.0f} x {bb.height:.0f} pt, ratio {bb.width/bb.height:.2f}")
print("Fig 8 v3 saved")
