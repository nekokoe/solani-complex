#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Describe, without any hypothesis test, where BMU02483's clinical-only genes come from.

Population: the 1,138 BMU02483 genes whose donors (>=99% identity, all 211 isolates)
include clinical isolates and include NONE of the 17 environmental isolates closest to
BMU02483 in the core-genome tree order -- i.e. the genes its own clade lacks.

Three descriptions
  1. donor compartment: accessory chromosome vs core chromosome vs unresolvable
  2. the AC donors: block (type1_stable) or mosaic (type2_mosaic) accessory chromosome
  3. the CC donors: distance of the donor copy to the nearest end of its chromosome,
     against the same measurement for every core-chromosome gene of the 13 assemblies
"""

import collections
import csv
import os
import statistics

TGS = "06.accessory-chrs/tgs_assemblies"
ANN = "02.annotation/per_isolates"
CLASS = "06.accessory-chrs/kmer_analysis/kmer_type.classification.v3.tsv"
ISO = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
       "BMU12683", "BMU12684", "BMU12685", "BMU12686",
       "Fu3.1", "Fu6.1", "SB1", "SX13"]

aclass = {r["AC"]: r["kmer_type"] for r in
          csv.DictReader(open(CLASS), delimiter="\t")}

# gene -> (contig, ac_id or None, distance to the nearest contig end)
gene = {}
clen = {}
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
ac_hit = collections.Counter()            # AC id -> clinical-only genes
cc_end = []                               # distances of CC donor copies
n_genes = 0
CC_ENDS = {100_000: 0, 50_000: 0, 10_000: 0}

for rep, mem in clusters.items():
    mine = [m for m in mem if m.startswith("BMU02483_")]
    if not mine:
        continue
    others = [m for m in mem if not m.startswith("BMU02483_")]
    isos = {iso_of(m) for m in others}
    isos.discard(None)
    if (isos & NEAR) or not (isos & CLIN):
        continue
    n_genes += len(mine)
    has_ac = [m for m in others if m in gene and gene[m][1]]
    has_cc = [m for m in others if m in gene and not gene[m][1]]
    if has_ac and has_cc:
        comp["AC + CC"] += len(mine)
    elif has_ac:
        comp["AC only"] += len(mine)
    elif has_cc:
        comp["CC only"] += len(mine)
    else:
        comp["unknown (draft donors)"] += len(mine)
    for m in has_ac:
        ac_hit[gene[m][1]] += len(mine)
    for m in has_cc:
        d = gene[m][2]
        cc_end.append(d)
        for t in CC_ENDS:
            if d <= t:
                CC_ENDS[t] += 1

print(f"BMU02483 clinical-only genes: {n_genes}\n")
print("1) donor compartment")
for k, v in comp.most_common():
    print(f"   {k:26s} {v:5d}  {100*v/n_genes:5.1f}%")

print("\n2) the AC donors: block vs mosaic")
bt = collections.Counter()
for ac, n in ac_hit.items():
    bt[aclass.get(ac, "unclassified")] += n
for k, v in bt.most_common():
    nacs = sum(1 for ac in ac_hit if aclass.get(ac, "unclassified") == k)
    print(f"   {k:18s} {v:5d} genes across {nacs:2d} ACs  ({100*v/sum(bt.values()):4.1f}% of AC-attributed genes)")
print("   per-AC detail (top 10):")
for ac, n in ac_hit.most_common(10):
    print(f"     {ac:16s} {aclass.get(ac,'?'):14s} {n:4d}")

print("\n3) the CC donors: distance to the nearest chromosome end")
bg = [gene[t][2] for t in gene if not gene[t][1]]
print(f"   CC donor copies of clinical-only genes: n={len(cc_end)}  median {statistics.median(cc_end)/1000:.1f} kb")
print(f"   all core-chromosome genes (13 assemblies): n={len(bg)}  median {statistics.median(bg)/1000:.1f} kb")
for t in sorted(CC_ENDS, reverse=True):
    print(f"   within {t//1000:3d} kb of an end: donors {100*CC_ENDS[t]/len(cc_end):5.1f}%"
          f"   background {100*sum(1 for d in bg if d <= t)/len(bg):5.1f}%")
