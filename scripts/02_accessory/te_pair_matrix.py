#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TE-pair matrix for transferred gene clusters (the S13e / S15 heatmap input).

For every transmitted cluster (>=2 contigs, at least one accessory), each member gene
contributes its flanking TE PAIR (nearest upstream family, nearest downstream family
within --window bp). Pairs with a non-TE flank (Unknown / Simple_repeat /
Low_complexity) are dropped. Families are ranked by how many pairs they take part in;
the top --top families form a square matrix of pair counts, rows = upstream,
columns = downstream.

Provenance
    Previously this was computed inside plot_s5_s6_boxplots_s15_heatmap.py, which also
    re-renders the unrelated S5/S6 repeat-content panels and reads one hard-coded
    cluster file. It is extracted here so the matrix can be recomputed for a different
    identity level without touching S5/S6, and so that BOTH id forms are accepted
    ("ISOLATE|transcript" from the archived clustering, bare "ISOLATE_transcript" from
    a fresh mmseqs run).

    Validated by reproducing 06.accessory-chrs/gene_transfer/suppl_table.te_pair_matrix_12.tsv
    exactly at the 100% level.

Usage
    python3 scripts/te_pair_matrix.py --clusters <clusters.tsv> --out-dir <dir> \
        [--window 5000] [--top 12]
"""

import argparse
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

ISOS = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
        "BMU12683", "BMU12684", "BMU12685", "BMU12686",
        "Fu3.1", "Fu6.1", "SB1", "SX13"]
NON_TE = {"Unknown", "Simple_repeat", "Low_complexity"}
ANN = "02.annotation/per_isolates"
TGS = "06.accessory-chrs/tgs_assemblies"


def load_genome():
    aff, gpos = {}, {}
    for iso in ISOS:
        with open(f"{ANN}/{iso}/funannotate/predict_results/fusarium.gff3") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 9 or p[2] != "mRNA":
                    continue
                tid = None
                for a in p[8].split(";"):
                    if a.startswith("ID="):
                        tid = a[3:]
                        break
                if tid:
                    aff[f"{iso}|{tid}"] = p[0]
                    gpos[f"{iso}|{tid}"] = (int(p[3]), int(p[4]))
    ct = {}
    for iso in ISOS:
        acc = set(open(f"{TGS}/{iso}.accessory_chrs.id").read().split())
        for line in open(f"{ANN}/{iso}/repeatmasker/{iso}.assembly.fasta.masked.length"):
            q = line.split()
            if len(q) >= 2:
                c, L = q[0], int(q[1])
                ct[c] = "acc" if c in acc else ("core" if L >= 500000 else "scaf")
    reps = {}
    for iso in ISOS:
        rep = defaultdict(list)
        with open(f"{TGS}/{iso}.assembly.fasta.out") as f:
            for i, line in enumerate(f):
                if i < 3:
                    continue
                q = line.split()
                if not q or len(q) < 11:
                    continue
                if len(q) > 15 and q[15] == "*":
                    continue
                rep[q[4]].append((int(q[5]), int(q[6]), q[10]))
        reps[iso] = {c: sorted(v) for c, v in rep.items()}
    return aff, gpos, ct, reps


def flank_pair(reps, iso, contig, gs, ge, win):
    iv = reps[iso].get(contig, [])
    left = [fam for s, e, fam in iv if e < gs and e >= gs - win]
    right = [fam for s, e, fam in iv if s > ge and s <= ge + win]
    return (left[-1] if left else None, right[0] if right else None)


def split_member(member, names):
    if "|" in member:
        return member
    for n in names:
        if member.startswith(n + "_"):
            return f"{n}|{member}"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--window", type=int, default=5000)
    ap.add_argument("--top", type=int, default=12)
    a = ap.parse_args()

    aff, gpos, ct, reps = load_genome()
    names = sorted(ISOS, key=len, reverse=True)

    members = defaultdict(list)
    with open(a.clusters) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            g = split_member(p[1], names)
            if g and g in aff:
                members[p[0]].append(g)

    pair_counts, fam_counts = Counter(), Counter()
    for rep_id, mem in members.items():
        contigs = {aff[m] for m in mem}
        if len(contigs) < 2:
            continue
        if "acc" not in {ct.get(c) for c in contigs}:
            continue
        for m in mem:
            pos = gpos.get(m)
            if pos is None:
                continue
            iso = m.split("|")[0]
            lf, rf = flank_pair(reps, iso, aff[m], pos[0], pos[1], a.window)
            if not lf or not rf:
                continue
            if lf in NON_TE or rf in NON_TE:
                continue
            pair_counts[(lf, rf)] += 1
            fam_counts[lf] += 1
            fam_counts[rf] += 1

    top_fams = [f for f, _ in fam_counts.most_common(a.top)]
    mat = np.zeros((len(top_fams), len(top_fams)), dtype=int)
    for (l, r), c in pair_counts.items():
        if l in top_fams and r in top_fams:
            mat[top_fams.index(l), top_fams.index(r)] = c

    os.makedirs(a.out_dir, exist_ok=True)
    pd.DataFrame(mat, index=top_fams, columns=top_fams).to_csv(
        os.path.join(a.out_dir, "suppl_table.te_pair_matrix_12.tsv"), sep="\t")
    pd.DataFrame([{
        "te_bracketed_genes": int(sum(pair_counts.values())),
        "in_top12_grid": int(mat.sum()),
        "n_top_families": len(top_fams),
        "n_distinct_pairs": len(pair_counts),
    }]).to_csv(os.path.join(a.out_dir, "suppl_table.te_pair_matrix_12.total.tsv"),
               sep="\t", index=False)
    print("te-bracketed genes:", int(sum(pair_counts.values())),
          "| in top-%d grid:" % a.top, int(mat.sum()),
          "| distinct pairs:", len(pair_counts))
    print("top families:", ", ".join(top_fams))


if __name__ == "__main__":
    main()
