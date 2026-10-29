#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TE-bracketed cassette test for transferred gene clusters.

Question
    Transferred gene families should share the SAME transposable-element bracket
    across their chromosomal copies if genes travel as TE-bracketed cassettes,
    rather than merely sitting in repeat-rich regions.

Definitions (as used for the published 90%/100% table)
    unit          a TRANSFERRED cluster: its member copies occupy >=2 contigs and at
                  least one of those is an accessory chromosome (A-A or C-A)
    bracket       per gene copy, the pair (nearest upstream family, nearest downstream
                  family) within --window bp on either side; families from RepeatMasker
                  (non-nested hits). This is the same `flank_pair` used for the S13e/S15
                  TE-pair matrix, so the two artefacts cannot drift apart.
    shared        a cluster counts as "sharing a bracket" when at least two of its
                  copies carry the identical (upstream, downstream) pair
    observed_frac fraction of eligible clusters that share a bracket
    specific      TRUE  -> only pairs whose two flanks are both specific TE families
                           (DNA / LTR / LINE / RC), excluding the ubiquitous
                           Unknown / Simple_repeat / Low_complexity
                  FALSE -> any pair whose two flanks are present
    null          permutation: the qualifying pairs are shuffled across the pooled genes
                  and handed back to clusters with their original sizes; the fraction
                  sharing a bracket is recomputed each time

Provenance
    The original producer of suppl_table.cassette_test.tsv was missing from the
    repository, so this is a re-implementation. It is validated by reproducing the
    archived cluster100 row exactly (specific 1,181 clusters / 0.5588; all-families
    2,170 clusters / 0.6041) before it is used for any other identity level.

Usage
    python3 scripts/te_cassette_test.py --clusters <clusters.tsv> --label cluster99 \
        --out <dir> [--window 10000] [--nperm 1000]
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
                if i < 3:                       # RepeatMasker header
                    continue
                q = line.split()
                if not q or len(q) < 11:
                    continue
                if len(q) > 15 and q[15] == "*":  # nested hit
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
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--window", type=int, default=10000)
    ap.add_argument("--nperm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261019)
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

    # pairs[cluster] = list of (up, down); only for TRANSFERRED clusters
    pairs = {}
    n_transferred = 0
    for rep, mem in members.items():
        contigs = {aff[m] for m in mem}
        if len(contigs) < 2:
            continue
        if "acc" not in {ct.get(c) for c in contigs}:
            continue
        n_transferred += 1
        ps = []
        for m in mem:
            pos = gpos.get(m)
            if pos is None:
                continue
            iso = m.split("|")[0]
            lf, rf = flank_pair(reps, iso, aff[m], pos[0], pos[1], a.window)
            if lf and rf:
                ps.append((lf, rf))
        if ps:
            pairs[rep] = ps
    print(f"transferred clusters: {n_transferred}; with >=1 complete bracket: {len(pairs)}")

    rng = np.random.default_rng(a.seed)
    rows = []
    for specific in (False, True):
        per_cluster = {}
        for rep, ps in pairs.items():
            keep = [p for p in ps if (not specific) or
                    (p[0] not in NON_TE and p[1] not in NON_TE)]
            # A cluster with a single bracketed copy can never "share" a bracket, so it
            # cannot dilute the denominator: eligibility needs TWO bracketed copies.
            if len(keep) >= 2:
                per_cluster[rep] = keep
        sizes = np.array([len(v) for v in per_cluster.values()])
        flat = [p for v in per_cluster.values() for p in v]

        def shared(values, sizes):
            i = 0
            hit = 0
            for k in sizes:
                chunk = values[i:i + k]
                i += k
                if len(set(chunk)) < len(chunk):
                    hit += 1
            return hit / len(sizes)

        obs = shared(flat, sizes)

        null = np.empty(a.nperm)
        idx = np.arange(len(flat))
        for t in range(a.nperm):
            rng.shuffle(idx)
            null[t] = shared([flat[j] for j in idx], sizes)
        p_emp = (1 + int((null >= obs).sum())) / (1 + a.nperm)

        rows.append({
            "Cluster": a.label, "Specific_TE_only": specific, "n_clusters": len(per_cluster),
            "observed_frac": obs, "null_mean_frac": float(null.mean()),
            "null_sd_frac": float(null.std(ddof=1)), "empirical_p": p_emp,
        })
        print(f"  specific={specific!s:5s} n={len(per_cluster):5d}  observed {obs:.4f}"
              f"  null {null.mean():.4f} ± {null.std(ddof=1):.4f}  p={p_emp:.4f}")

    os.makedirs(a.out, exist_ok=True)
    out = os.path.join(a.out, "suppl_table.cassette_test.tsv")
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
