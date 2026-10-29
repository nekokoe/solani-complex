#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Statistical tests for the two routes by which BMU02483 acquired its clinical-only genes.

Population (same as bmu2483_clinical_only_description.py): the 1,138 BMU02483 genes whose
donors (>=99% identity over 211 isolates) include clinical isolates and none of the 17
environmental isolates closest to it in the core-genome tree order.

AC route, two questions
  A1  is the accessory-chromosome contribution larger than the genome's compartment
      composition would give?  229 AC-only vs 164 CC-only against 13,825 AC-resident vs
      198,064 core-resident genes (Fisher exact).
  A2  within the AC side, is the mosaic class over-represented relative to its share of
      accessory-chromosome gene content?  (pooled Fisher + per-AC Mann-Whitney)

CC route
  B1  are the core-side donor copies closer to a chromosome end than core genes at large?
      Fisher exact at 10 / 50 / 100 kb, plus Mann-Whitney on the distance distributions.
"""

import collections
import csv
import statistics

from scipy import stats as st

TGS = "06.accessory-chrs/tgs_assemblies"
ANN = "02.annotation/per_isolates"
CLASS = "06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv"
ISO = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
       "BMU12683", "BMU12684", "BMU12685", "BMU12686",
       "Fu3.1", "Fu6.1", "SB1", "SX13"]

aclass = {r["AC"]: r["kmer_type"] for r in csv.DictReader(open(CLASS), delimiter="\t")}

gene = {}
clen = {}
ac_genes = collections.Counter()
for iso in ISO:
    acc = set(open(f"{TGS}/{iso}.accessory_chrs.id").read().split())
    for line in open(f"{ANN}/{iso}/repeatmasker/{iso}.assembly.fasta.masked.length"):
        q = line.split()
        if len(q) >= 2:
            clen[q[0]] = int(q[1])
    for line in open(f"{ANN}/{iso}/funannotate/predict_results/fusarium.gff3"):
        if line.startswith("#"):
            continue
        f = line.rstrip().split("\t")
        if len(f) < 9 or f[2] != "mRNA":
            continue
        for at in f[8].split(";"):
            if at.startswith("ID="):
                tid = at[3:]
                break
        else:
            continue
        L = clen.get(f[0], 0)
        s, e = int(f[3]), int(f[4])
        gene[tid] = (f[0], f[0] if f[0] in acc else None, min(s, max(L - e, 0)))
        if f[0] in acc:
            ac_genes[f[0]] += 1

grp = {r["BMU_ID"]: r for r in csv.DictReader(
    open("05.pan-genome/bmu_ref_groups.211.txt"), delimiter="\t")}
names = sorted([l.strip() for l in open("05.pan-genome/sample.solaniSC.211.list") if l.strip()],
               key=len, reverse=True)


def iso_of(g):
    for n in names:
        if g.startswith(n + "_"):
            return n
    return None


order = [l.strip() for l in open(f"{TGS}/tree_order_allSC.txt") if l.strip()]
b = order.index("BMU02483")
near = order[max(0, b - 10):b + 11]
NEAR = {s for s in near if grp.get(s, {}).get("Species") == "falciformis"
        and grp[s]["MedEnv"] == "env"}
CLIN = {k for k, v in grp.items() if v["MedEnv"] == "med"}

clusters = collections.defaultdict(list)
for line in open("tmp/mmseqs99/out99_cluster.tsv"):
    p = line.rstrip().split("\t")
    if len(p) > 1:
        clusters[p[0]].append(p[1])

comp = collections.Counter()
attr = collections.Counter()
cc_dist = []
for rep, mem in clusters.items():
    mine = [m for m in mem if m.startswith("BMU02483_")]
    if not mine:
        continue
    others = [m for m in mem if not m.startswith("BMU02483_")]
    isos = {iso_of(m) for m in others}
    isos.discard(None)
    if (isos & NEAR) or not (isos & CLIN):
        continue
    acm = [m for m in others if m in gene and gene[m][1]]
    ccm = [m for m in others if m in gene and not gene[m][1]]
    if acm and ccm:
        comp["both"] += len(mine)
    elif acm:
        comp["AC only"] += len(mine)
    elif ccm:
        comp["CC only"] += len(mine)
    else:
        comp["unknown"] += len(mine)
    for m in acm:
        attr[gene[m][1]] += len(mine)
    for m in ccm:
        cc_dist.append(gene[m][2])

n_res = comp["AC only"] + comp["CC only"] + comp["both"]
print(f"clinical-only genes {sum(comp.values())}: " +
      " · ".join(f"{k} {v}" for k, v in comp.most_common()))

# ---- A1: AC contribution vs the genome's compartment composition -------------------
ACRES = sum(ac_genes.values())                       # genes on the 78 ACs
CCRES = sum(1 for t in gene if not gene[t][1])       # genes on core chromosomes
obs_ac = comp["AC only"] + comp["both"]
obs_cc = comp["CC only"]
orr, p = st.fisher_exact([[obs_ac, obs_cc], [ACRES, CCRES]])
print(f"\nA1  AC vs CC contribution")
print(f"    observed  AC {obs_ac} / CC {obs_cc}   ({100*obs_ac/(obs_ac+obs_cc):.1f}% AC)")
print(f"    available AC {ACRES} / CC {CCRES}       ({100*ACRES/(ACRES+CCRES):.2f}% AC)")
print(f"    OR = {orr:.1f}   Fisher P = {p:.3e}")

# ---- A2: mosaic over-representation within the AC side ---------------------------
mos_attr = sum(v for k, v in attr.items() if aclass.get(k) == "type2_mosaic")
blk_attr = sum(v for k, v in attr.items() if aclass.get(k) == "type1_stable")
mos_g = sum(v for k, v in ac_genes.items() if aclass.get(k) == "type2_mosaic")
blk_g = sum(v for k, v in ac_genes.items() if aclass.get(k) == "type1_stable")
orr2, p2 = st.fisher_exact([[mos_attr, blk_attr], [mos_g, blk_g]])
rate = lambda a, n: 100 * a / n
print(f"\nA2  mosaic vs block within the AC side")
print(f"    attributed genes  mosaic {mos_attr} / block {blk_attr}   ({100*mos_attr/(mos_attr+blk_attr):.1f}% mosaic)")
print(f"    AC gene content   mosaic {mos_g} / block {blk_g}   ({100*mos_g/(mos_g+blk_g):.1f}% mosaic)")
print(f"    per-AC rate       mosaic {rate(mos_attr, mos_g):.2f}%   block {rate(blk_attr, blk_g):.2f}%")
print(f"    pooled OR = {orr2:.2f}   Fisher P = {p2:.3e}")
mr = [attr.get(a, 0) / n for a, n in ac_genes.items() if aclass.get(a) == "type2_mosaic"]
br = [attr.get(a, 0) / n for a, n in ac_genes.items() if aclass.get(a) == "type1_stable"]
u, pu = st.mannwhitneyu(mr, br, alternative="greater")
print(f"    per-AC Mann-Whitney (n={len(mr)} mosaic, {len(br)} block): "
      f"median {100*statistics.median(mr):.3f}% vs {100*statistics.median(br):.3f}%  P = {pu:.3e}")

# ---- B1: are the CC donors terminal? ---------------------------------------------
bg = [gene[t][2] for t in gene if not gene[t][1]]
print(f"\nB1  CC route: distance to the nearest chromosome end")
print(f"    donor copies n={len(cc_dist)}, median {statistics.median(cc_dist)/1000:.1f} kb"
      f"   background n={len(bg)}, median {statistics.median(bg)/1000:.1f} kb")
for t in (10_000, 50_000, 100_000):
    a1 = sum(1 for d in cc_dist if d <= t)
    a0 = sum(1 for d in bg if d <= t)
    o, pp = st.fisher_exact([[a1, len(cc_dist) - a1], [a0, len(bg) - a0]])
    print(f"    within {t//1000:3d} kb: {a1:3d}/{len(cc_dist)} = {100*a1/len(cc_dist):4.1f}%"
          f"  vs {a0}/{len(bg)} = {100*a0/len(bg):4.1f}%"
          f"   OR = {o:5.2f}   P = {pp:.2e}")
u2, p2b = st.mannwhitneyu(cc_dist, bg, alternative="less")
print(f"    Mann-Whitney (donors closer to an end): U = {u2:.0f}   P = {p2b:.2e}")
