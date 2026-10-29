#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cutoff sensitivity for the BMU02483 gene-origin analysis.

The origin Sankey infers a donor by asking which other isolates carry a copy of the
gene. That answer depends entirely on how "a copy" is defined, so this script runs the
SAME classification under four definitions and prints how the flows move:

    identity 90 / 95 / 100 %   mmseqs identity clusters over all 211 isolates
                               (all_solaniSC.cluster{90,95,100}_cluster.tsv)
    HOG                        OrthoFinder orthogroups
                               (N0.gene_vs_HOG.tsv + N0.HOG.hog_count.tsv)

The identity clusters ask "does an IDENTICAL / nearly identical protein exist in another
isolate"; the HOG asks "does an ORTHOLOGUE exist there". The first is a strict claim about
gene flow, the second a loose claim about the gene family, so moving from clusters to HOGs
inflates every donor class -- the question is by how much.

Sources and compartments are defined exactly as in bmu2483_gene_origin_sankey.py:
  nearby env only / clinical only / both / no other donor
  ACs (pooled) / CC / no chromosome-level donor

Output: prints the comparison table and writes
  06.accessory-chrs/BMU02483_gene_content/suppl_table.bmu2483_cutoff_sensitivity.tsv

Usage: python3 scripts/bmu2483_cutoff_sensitivity.py
"""

import collections
import csv
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TGS = os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies")
ANN = os.path.join(REPO, "02.annotation", "per_isolates")
OF = os.path.join(REPO, "05.pan-genome", "orthofinder-shorter")
OUT = os.path.join(REPO, "06.accessory-chrs", "BMU02483_gene_content")
FOCUS = "BMU02483"

ISO = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
       "BMU12683", "BMU12684", "BMU12685", "BMU12686",
       "Fu3.1", "Fu6.1", "SB1", "SX13"]
SRC_ORDER = ["nearby env only", "clinical only", "both", "no other donor"]


def iso_of(gid):
    """'BMU01933|BMU01933_000001-T1' -> BMU01933 ; 'BMU01933_000001-T1' -> BMU01933

    The archived identity files store 'isolate|transcript'; a fresh mmseqs easy-cluster
    run writes the bare transcript id. A few isolate names themselves contain an
    underscore (FSSC_6, IMI_366586, NRRL_20438), so the fallback strips the trailing
    _<digits>-T<n> rather than splitting on the first underscore.
    """
    if "|" in gid:
        return gid.split("|", 1)[0]
    return re.sub(r"_\d+-T\d+$", "", gid)
COMP_ORDER = ["ACs", "CC", "no chromosome-level donor"]


def contexts():
    grp = {r["BMU_ID"]: r for r in csv.DictReader(
        open(os.path.join(REPO, "05.pan-genome", "bmu_ref_groups.211.txt")), delimiter="\t")}
    order = [l.strip() for l in open(os.path.join(TGS, "tree_order_allSC.txt")) if l.strip()]
    b = order.index(FOCUS)
    near = order[max(0, b - 10):b + 11]
    f_env = {k for k, v in grp.items()
             if v["Species"] == "falciformis" and v["MedEnv"] == "env"}
    NEAR = {s for s in near if s in f_env}
    CLIN = {k for k, v in grp.items() if v["MedEnv"] == "med"}
    tid2c = {}
    for iso in ISO:
        acc = set(open(os.path.join(TGS, f"{iso}.accessory_chrs.id")).read().split())
        for line in open(os.path.join(ANN, iso, "funannotate", "predict_results",
                                      "fusarium.gff3")):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            tid = None
            for a in f[8].split(";"):
                if a.startswith("ID="):
                    tid = a[3:]
                    break
            if tid:
                # BOTH keys: the identity clusters use "iso|transcript", the OrthoFinder
                # tables use the bare transcript id. Indexing by only one form made every
                # HOG compartment lookup miss (all genes fell into "no chromosome-level
                # donor" -- a silent 100% rather than an error).
                c = "ACs" if f[0] in acc else "CC"
                tid2c[f"{iso}|{tid}"] = c
                tid2c[tid] = c
    return NEAR, CLIN, tid2c


def classify(members_by_focus, NEAR, CLIN, tid2c):
    """members_by_focus: {focus gene -> set of other-isolate gene ids}"""
    src = collections.Counter()
    comp = collections.Counter()
    for gene, others in members_by_focus.items():
        isos = {iso_of(m) for m in others}
        ne = bool(isos & NEAR)
        cli = bool(isos & CLIN)
        s = ("both" if (ne and cli) else "nearby env only" if ne else
             "clinical only" if cli else "no other donor")
        src[s] += 1
        if any(tid2c.get(m) == "ACs" for m in others):
            comp["ACs"] += 1
        elif any(tid2c.get(m) == "CC" for m in others):
            comp["CC"] += 1
        else:
            comp["no chromosome-level donor"] += 1
    return src, comp


def from_identity_clusters(path):
    """two passes: find the clusters that contain a FOCUS member, then read them fully"""
    reps = set()
    with open(path) as fh:
        for line in fh:
            a, _, b = line.rstrip("\n").partition("\t")
            if not b:
                continue
            if iso_of(a) == FOCUS or iso_of(b) == FOCUS:
                reps.add(a)
    print(f"    clusters containing {FOCUS}: {len(reps)}")
    members = collections.defaultdict(set)
    with open(path) as fh:
        for line in fh:
            a, _, b = line.rstrip("\n").partition("\t")
            if a in reps and b:
                members[a].add(b)
    out = {}
    for rep, mem in members.items():
        mine = [m for m in mem if iso_of(m) == FOCUS]
        others = {m for m in mem if iso_of(m) != FOCUS}
        # the representative itself counts as a member of its own cluster
        for m in mine:
            out[m] = others
    return out


def from_hog():
    """focus gene -> HOG -> the samples that carry that HOG (all copied to every focus gene)"""
    g2h = {}
    with open(os.path.join(OF, "N0.gene_vs_HOG.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if row["SampleID"] == FOCUS or row["SampleID"] in ISO:
                g2h[row["GeneID"]] = row["HOG_ID"]
    carriers = collections.defaultdict(set)
    with open(os.path.join(OF, "N0.HOG.hog_count.tsv")) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            if int(row["GeneCount"]) > 0:
                carriers[row["HOG_ID"]].add(row["SampleID"])
    h2genes = collections.defaultdict(list)
    for g, h in g2h.items():
        h2genes[h].append(g)
    out = {}
    for h, genes in h2genes.items():
        mine = [g for g in genes if g.startswith(FOCUS + "_")]
        if not mine:
            continue
        # donor "genes" at HOG level are the iso set; expand to iso|* so the compartment
        # lookup can search the TGS member genes of the same HOG
        others = {m for m in genes if m.partition("|")[0] != FOCUS}
        for g in mine:
            out[g] = others | {f"{s}|*" for s in carriers[h] if s != FOCUS}
    return out


def main():
    NEAR, CLIN, tid2c = contexts()
    print(f"nearby environmental control: {len(NEAR)} isolates | clinical: {len(CLIN)}")
    rows = []
    # ALL identity levels are re-run locally in one protocol:
    #   mmseqs easy-cluster all_proteins.faa out<N> tmp<N> --min-seq-id <N> -c 1.0
    #     --cov-mode 0 --cluster-reassign 1 --threads 48
    # They are NOT the archived all_solaniSC.cluster{90,95,100}_cluster.tsv files. Those
    # cannot be reproduced: the same command yields 977,779 clusters at 100% against the
    # archive's 906,306, with only 34% of the BMU02483 genes getting an identical member
    # set. Mixing a fresh run with the archive would therefore compare two protocols, so
    # the whole series comes from one run and the archive mismatch is reported instead.
    LOCAL = {90: "out90_cluster.tsv", 95: "out95_cluster.tsv",
             99: "out99_cluster.tsv", 100: "out100_cluster.tsv"}
    base = os.path.join(REPO, "tmp", "mmseqs99")
    levels = []
    for lvl in (90, 95, 99, 100):
        fp = os.path.join(base, LOCAL[lvl])
        if os.path.exists(fp):
            levels.append((lvl, fp))
        else:
            print(f"  (missing local {lvl}% run: {fp})")
    for lvl, path in levels:
        print(f"\n=== identity {lvl}% (locally re-run) ===")
        m = from_identity_clusters(path)
        src, comp = classify(m, NEAR, CLIN, tid2c)
        n = sum(src.values())
        rows.append(("identity %d%%" % lvl, n, src, comp))
    print("\n=== HOG (OrthoFinder orthogroups) ===")
    m = from_hog()
    src, comp = classify(m, NEAR, CLIN, tid2c)
    n = sum(src.values())
    rows.append(("HOG", n, src, comp))

    pct = lambda k, n: f"{100*k/n:5.1f}%" if n else "  n/a"
    print(f"\n{'definition':<14}{'genes':>8}   " +
          "".join(f"{s:>19}" for s in SRC_ORDER))
    for name, n, src, _ in rows:
        print(f"{name:<14}{n:>8}   " +
              "".join(f"{str(src[s])+' ('+pct(src[s],n)+')':>19}" for s in SRC_ORDER))
    print(f"\n{'definition':<14}{'genes':>8}   " + "".join(f"{c:>19}" for c in COMP_ORDER))
    for name, n, _, comp in rows:
        print(f"{name:<14}{n:>8}   " +
              "".join(f"{str(comp[c])+' ('+pct(comp[c],n)+')':>19}" for c in COMP_ORDER))

    with open(os.path.join(OUT, "suppl_table.bmu2483_cutoff_sensitivity.tsv"), "w") as fh:
        fh.write("definition\tn_genes\tdimension\tclass\tn\tpct\n")
        for name, n, src, comp in rows:
            for k in SRC_ORDER:
                fh.write(f"{name}\t{n}\tsource\t{k}\t{src[k]}\t{100*src[k]/n:.2f}\n")
            for k in COMP_ORDER:
                fh.write(f"{name}\t{n}\tcompartment\t{k}\t{comp[k]}\t{100*comp[k]/n:.2f}\n")
    print(f"\nwrote suppl_table.bmu2483_cutoff_sensitivity.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
