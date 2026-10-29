#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Did the accessory chromosome act as the VEHICLE for BMU02483's clinical gene uptake?

The question
    Composition alone cannot answer it: "59% of BMU02483's clinical genes sit on an AC" is
    compatible with ACs being the vehicle, but also with AC-resident genes simply being
    more mobile. The discriminating test is an ENTRY RATE, not a share:

        rate(AC)  = of all gene families found ONLY on the accessory chromosomes of
                    CLINICAL isolates, what fraction reached BMU02483?
        rate(CC)  = the same fraction for families found ONLY on core chromosomes

    Under "random" horizontal transfer the compartment should not matter, so rate(AC) >>
    rate(CC) is the signature of the AC acting as a preferential vehicle.

Controls built in
    * "acquired, not inherited": by default only families ABSENT from the 17 environmental
      isolates closest to BMU02483 are counted, so a family it shares with its own clade
      cannot inflate the rate (--no-near-env-filter turns this off).
    * prevalence: a family carried by many clinical isolates has more opportunities to be
      picked up whatever its compartment, so the comparison is ALSO reported stratified by
      the number of clinical carriers (Mantel-Haenszel odds ratio over the strata).
    * only clinical-carried families enter the denominator.

Definition of the compartments (13 chromosome-level assemblies only)
    AC-only   >=1 member on an accessory chromosome, none on any core chromosome
    CC-only   >=1 member on a core chromosome, none on any accessory chromosome
    both      members on both (excluded -- it cannot discriminate)
    unknown   no member in a chromosome-level assembly (its donor is a draft)

Usage
    python3 scripts/bmu2483_ac_vehicle_test.py [--clusters tmp/mmseqs99/out99_cluster.tsv] \
        [--no-near-env-filter]
"""

import argparse
import collections
import csv
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TGS = os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies")
ANN = os.path.join(REPO, "02.annotation", "per_isolates")
FOCUS = "BMU02483"
ISO = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
       "BMU12683", "BMU12684", "BMU12685", "BMU12686",
       "Fu3.1", "Fu6.1", "SB1", "SX13"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", default=os.path.join(REPO, "tmp", "mmseqs99", "out99_cluster.tsv"))
    ap.add_argument("--no-near-env-filter", action="store_true")
    a = ap.parse_args()

    grp = {r["BMU_ID"]: r for r in csv.DictReader(
        open(os.path.join(REPO, "05.pan-genome", "bmu_ref_groups.211.txt")), delimiter="\t")}
    names = sorted([l.strip() for l in
                    open(os.path.join(REPO, "05.pan-genome", "sample.solaniSC.211.list")) if l.strip()],
                   key=len, reverse=True)
    order = [l.strip() for l in open(os.path.join(TGS, "tree_order_allSC.txt")) if l.strip()]
    b = order.index(FOCUS)
    near = order[max(0, b - 10):b + 11]
    NEAR = {s for s in near if grp.get(s, {}).get("Species") == "falciformis"
            and grp[s]["MedEnv"] == "env"}
    CLIN = {k for k, v in grp.items() if v["MedEnv"] == "med"}
    CLIN_TGS = {i for i in ISO if grp[i]["MedEnv"] == "med"}

    # compartment of every gene in the 13 assemblies
    tid2c = {}
    for iso in ISO:
        acc = set(open(os.path.join(TGS, f"{iso}.accessory_chrs.id")).read().split())
        gff = os.path.join(ANN, iso, "funannotate", "predict_results", "fusarium.gff3")
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "mRNA":
                continue
            for at in f[8].split(";"):
                if at.startswith("ID="):
                    tid2c[at[3:]] = (iso, "AC" if f[0] in acc else "CC")
                    break

    def iso_of(g):
        for n in names:
            if g.startswith(n + "_"):
                return n
        return None

    def contig_class(n):
        return tid2c[n][1] if n in tid2c else None

    # The mmseqs cluster file has ONE LINE PER MEMBER (rep<TAB>member), so the members
    # must be aggregated per representative first. Treating each LINE as a cluster
    # (the first version of this script did) inflates the denominators with member
    # pairs and deflates "in BMU02483" -- it reported 25/8294 vs 359/128587, which
    # cannot be squared with the 229 vs 164 genes counted directly.
    clusters = collections.defaultdict(list)
    with open(a.clusters) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                clusters[p[0]].append(p[1])

    stat = collections.defaultdict(lambda: {"n": 0, "in_bmu": 0, "strata": collections.Counter()})
    if True:
        for _rep, members in clusters.items():
            isos = {iso_of(m) for m in members}
            isos.discard(None)
            if not isos:
                continue
            has_focus = any(m.startswith(FOCUS + "_") for m in members)
            classes = {contig_class(m) for m in members}
            classes.discard(None)
            if "AC" in classes and "CC" in classes:
                comp = "both"
            elif "AC" in classes:
                # the AC must belong to a CLINICAL isolate
                owners = {tid2c[m][0] for m in members if m in tid2c and tid2c[m][1] == "AC"}
                comp = "AC-only" if owners & CLIN_TGS else "AC-env"
            elif "CC" in classes:
                comp = "CC-only"
            else:
                comp = "unknown"
            if comp not in ("AC-only", "CC-only"):
                continue
            n_clin = len(isos & CLIN)
            if n_clin < 1:                       # must be carried by clinical isolates
                continue
            if not a.no_near_env_filter and (isos & NEAR):
                continue                          # present in its own clade -> not acquired
            s = stat[comp]
            s["n"] += 1
            s["in_bmu"] += 1 if has_focus else 0
            s["strata"][n_clin] += 1
            if has_focus:
                s["focus" + str(n_clin)] = s.get("focus" + str(n_clin), 0) + 1

    print("filter: " + ("OFF (all clinical-carried families)" if a.no_near_env_filter
                        else "families absent from the 17 nearby environmental isolates"))
    print(f"\n{'compartment':12s} {'families':>9s} {'in BMU02483':>12s} {'entry rate':>11s}")
    for comp in ("AC-only", "CC-only"):
        s = stat[comp]
        r = s["in_bmu"] / s["n"] if s["n"] else float("nan")
        print(f"{comp:12s} {s['n']:9d} {s['in_bmu']:12d} {100*r:10.2f}%")

    ac, cc = stat["AC-only"], stat["CC-only"]
    from scipy import stats as st
    orr, p = st.fisher_exact([[ac["in_bmu"], ac["n"] - ac["in_bmu"]],
                              [cc["in_bmu"], cc["n"] - cc["in_bmu"]]])
    print(f"\ncrude: odds ratio {orr:.2f}   Fisher P = {p:.3e}")

    # Mantel-Haenszel over strata of clinical carrier count
    strata = sorted(set(ac["strata"]) | set(cc["strata"]))
    num = den = 0.0
    print(f"\n{'clinical carriers':>18s} {'AC in/total':>14s} {'CC in/total':>14s} {'OR':>8s}")
    for k in strata:
        a1 = ac.get("focus" + str(k), 0)
        n1 = ac["strata"][k]
        a0 = cc.get("focus" + str(k), 0)
        n0 = cc["strata"][k]
        if n1 == 0 and n0 == 0:
            continue
        N = n1 + n0
        num += a1 * (n0 - a0) / N
        den += a0 * (n1 - a1) / N
        o = ((a1 / (n1 - a1)) / (a0 / (n0 - a0))) if a1*(n1-a1) and a0*(n0-a0) and a0 and (n0-a0) else float("nan")
        print(f"{k:18d} {a1:6d}/{n1:<7d} {a0:6d}/{n0:<7d} {o:8.2f}")
    print(f"\nMantel-Haenszel OR = {num/den:.2f}   (prevalence-stratified)")


if __name__ == "__main__":
    main()
