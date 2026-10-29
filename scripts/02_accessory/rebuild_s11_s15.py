#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rebuild S11 (method abstract, cleaned text) and S15 (cassette, cluster100
only) after the user decision to drop the 90% analysis from the paper."""

import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ISOS = ["BMU01933","BMU02028","BMU02474","BMU03267","BMU06673","BMU12683",
        "BMU12684","BMU12685","BMU12686","Fu3.1","Fu6.1","SB1","SX13"]
iso_set = set(ISOS)

# ---------- shared data loading ----------
aff = {}
gpos = {}
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

# ---------- S15: cassette (cluster100 only) ----------
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

NON_TE = {"Unknown", "Simple_repeat", "Low_complexity"}

def flank_pair(iso, contig, gs, ge, win=1000):
    iv = reps[iso].get(contig, [])
    left = [fam for s, e, fam in iv if e < gs and e >= gs - win]
    right = [fam for s, e, fam in iv if s > ge and s <= ge + win]
    lf = left[-1] if left else None
    rf = right[0] if right else None
    return lf, rf

# transfer clusters (cluster100): 13-TGS members, >=2 contigs, accessory involved
rep_members = defaultdict(list)
with open("06.accessory-chrs/tgs_assemblies/all_solaniSC.cluster100_cluster.tsv") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 2: continue
        iso, _, trans = p[1].partition("|")
        if iso in iso_set and trans in aff:
            rep_members[p[0]].append(p[1])

pair_counter = Counter()
n_all, n_all_shared, n_spec, n_spec_shared = 0, 0, 0, 0
for rep_id, members in rep_members.items():
    contigs = {aff[m.split("|")[1]] for m in members}
    if len(contigs) < 2: continue
    ctypes = {ct.get(c) for c in contigs if c}
    if "acc" not in ctypes: continue
    n_all += 1
    pairs = []
    for m in members:
        iso = m.split("|")[0]; tid = m.split("|")[1]
        pos = gpos.get(m)
        if pos is None: continue
        lf, rf = flank_pair(iso, aff[tid], pos[0], pos[1])
        if lf and rf:
            pairs.append((lf, rf))
    if len(pairs) >= 2:
        cnt = Counter(pairs)
        if max(cnt.values()) >= 2:
            n_all_shared += 1
            spec_pairs = {pr: c for pr, c in cnt.items()
                          if pr[0] not in NON_TE and pr[1] not in NON_TE}
            if spec_pairs:
                n_spec_shared += 1
                best = max(spec_pairs, key=spec_pairs.get)
                pair_counter[best] += 1
    n_spec += 1

# left-panel values from suppl_table.cassette_test.tsv (cluster100 rows)
ct_test = pd.read_csv("06.accessory-chrs/gene_transfer/suppl_table.cassette_test.tsv", sep="\t")
c100 = ct_test[ct_test.Cluster == "cluster100"]
obs_all = c100[c100.Specific_TE_only == False].observed_frac.iloc[0]
null_all = c100[c100.Specific_TE_only == False].null_mean_frac.iloc[0]
obs_spec = c100[c100.Specific_TE_only == True].observed_frac.iloc[0]
null_spec = c100[c100.Specific_TE_only == True].null_mean_frac.iloc[0]

top_pairs = pair_counter.most_common(8)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(6.4, 2.7), dpi=300,
                               gridspec_kw={"width_ratios": [1, 1.25]})
# left: observed vs permutation
x = np.arange(2)
w = 0.32
axL.bar(x - w/2, [obs_all*100, obs_spec*100], w, color="#C44E52", label="Observed")
axL.bar(x + w/2, [null_all*100, null_spec*100], w, color="#9E9E9E", label="Permutation expectation")
axL.set_xticks(x)
axL.set_xticklabels(["cluster100\n(all families)", "cluster100\n(specific TEs)"], fontsize=6.5)
axL.set_ylabel("Clusters with ≥2 copies sharing\nan identical TE bracket (%)", fontsize=6.5)
axL.set_ylim(0, 70)
axL.tick_params(labelsize=6.5)
axL.legend(fontsize=6, frameon=False)
for xi, o in zip(x - w/2, [obs_all*100, obs_spec*100]):
    axL.text(xi, o + 1, f"{o:.0f}%", ha="center", fontsize=6)
axL.set_title("TE-bracketed cassette sharing", fontsize=7.5)
for s in ["top", "right"]:
    axL.spines[s].set_visible(False)

# right: top specific TE bracket pairs (cluster100)
labels = [f"{a}\n{b}" for a, b in top_pairs]
vals = [v for _, v in top_pairs]
axR.barh(range(len(top_pairs)), vals, color="#3B7DD8", height=0.6)
axR.set_yticks(range(len(top_pairs)))
axR.set_yticklabels(labels, fontsize=5.2)
axR.invert_yaxis()
axR.set_xlabel("Clusters sharing the specific TE bracket", fontsize=6.5)
axR.set_title("Most shared specific TE brackets (cluster100)", fontsize=7.5)
axR.tick_params(labelsize=6.5)
for i, v in enumerate(vals):
    axR.text(v + 0.2, i, str(v), va="center", fontsize=6)
for s in ["top", "right"]:
    axR.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig("06.accessory-chrs/gene_transfer/figure_transfer_cassette.pdf")
print("S15 rebuilt (cluster100 only):", len(top_pairs), "top pairs")

# ---------- S11: method abstract (cleaned) ----------
fig, ax = plt.subplots(figsize=(6.2, 4.4), dpi=300)
ax.set_xlim(0, 12); ax.set_ylim(0, 12); ax.axis("off")

def box(x, y, w, h, text, fc="#F4F4F4", fs=6.2, bold=False, align="center"):
    ax.add_patch(plt.Rectangle((x, y), w, h, fc=fc, ec="#333333", lw=0.7, zorder=1))
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs,
            fontweight="bold" if bold else "normal", zorder=2)

def arr(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#333333"))

ax.text(6, 11.6, "Method: are transferred genes mediated by specific repeat families?",
        ha="center", fontsize=8, fontweight="bold")

box(0.6, 10.3, 5.0, 0.95, "All 211 FSSC genomes\nmmseqs protein clusters, 100% identity")
box(6.4, 10.3, 5.0, 0.95, "13 TGS assemblies\nGFF3 genes + RepeatMasker (non-nested)\n+ accessory-chromosome calls")
arr(3.1, 10.8, 6.4, 10.8); arr(6.15, 10.3, 6.15, 9.78)

box(2.0, 8.6, 8.0, 1.1,
    "1. Transfer detection (211-genome context)\ncluster spans ≥2 chromosomes AND involves ≥1 accessory chromosome\n→ AC-related transfer candidates (A-A: accessory-accessory; C-A: core-accessory)", fs=5.6)

box(2.0, 6.6, 8.0, 1.15,
    "2. Filter candidates on the 13 TGS assemblies\n11,895 transfer-involved genes (A-A 8,988 + C-A 2,907)\nvs whole-genome background: 199,994 other genes", fs=5.6)

box(2.0, 4.6, 8.0, 1.15,
    "3. Flanking test (1 kb and 5 kb)\nfor each gene: repeat families overlapping or flanking the gene\nper-family Fisher exact test: transfer candidates vs genome background", fs=5.6)

box(2.0, 2.6, 8.0, 1.35,
    "4. Conclusion\ntransferred genes are significantly TE-flanked: LINE (×17), RC (×15),\nLTR & DNA (×13; hAT, IS3EU, PiggyBac, TcMar-Fot1, Tad1)\nSimple repeats / low-complexity depleted (×0.5)", fs=5.6)

box(2.0, 0.5, 8.0, 1.35,
    "Robustness: core-resident matched control\naccessory-derived genes on core chromosomes vs other core genes:\nDNA ×7.0, LINE ×7.1, LTR ×5.1 (P < 0.001)", fs=5.6)

# flow arrows along the LEFT margin, outside the text boxes
arr(1.2, 8.6, 1.2, 7.8)
arr(1.2, 6.6, 1.2, 5.8)
arr(1.2, 4.6, 1.2, 4.0)
arr(1.2, 2.6, 1.2, 1.9)

fig.savefig("06.accessory-chrs/gene_transfer_211/figure_method_abstract.pdf")
print("S11 rebuilt (cleaned text)")
