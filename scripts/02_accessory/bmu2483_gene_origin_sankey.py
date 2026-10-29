#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BMU02483 gene origins: a four-column flow, from the gene pool to the donor's identity.

Motivation
    BMU02483 is the environmental phylogenetic outlier whose accessory genome looks
    clinical, but it has no long-read assembly, so its own accessory chromosomes
    cannot be identified and inferring them from other isolates is unsound. This
    script asks the question the data CAN answer: for every BMU02483 gene, which
    other isolates carry an identical copy (>=99% identity), and what is the donor's
    species, isolation background and chromosomal compartment?

Three dimensions, each read off the same donor set
    1. SPECIES SOURCE -- the species of the donor isolates
         falciformis   only F. falciformis donors
         solani        only F. solani donors
         both          donors of both species -> the species origin is unresolvable
         other         neither species (other named species, "Unknown", or no donor)
    2. ISOLATE SOURCE -- the donor's isolation background
         from nearby env isolates  donors are the F. falciformis ENVIRONMENTAL
                                   isolates within +/-10 positions of BMU02483 in the
                                   core-genome tree order (the same control set as
                                   scripts/ac_bmu2483_gene_content.py -> 17 isolates)
         from clinical isolates    donors are clinical isolates
         both                      donors of BOTH kinds -> the direction is
                                   unresolvable, kept as its own class rather than
                                   forced to one side
         other                     no donor, or donors that are neither nearby-env nor
                                   clinical (i.e. environmental isolates elsewhere)
    3. KNOWN CHROMOSOME SOURCE -- among the 13 chromosome-level assemblies only
         ACs       a named accessory chromosome carries a donor copy (pooled node)
         CCs       a core chromosome carries it
         unknown   every donor is a draft assembly, so its compartment is unknown

Every gene contributes ONE path, so each column sums to the same gene total.

WHY 99% AND NOT 100%
    The figure is drawn from the >=99%-identity clustering (tmp/mmseqs99/out99_cluster),
    re-run locally with the same protocol as the 90/95/100% levels. Only the fully
    self-consistent four-level series is used: the ARCHIVED cluster{90,95,100} files
    CANNOT be reproduced on this machine (fresh 100%: 977,779 clusters vs the archive's
    906,306; only 34.4% identical member sets), so mixing an archived level with a fresh
    one would be silently incomparable. The four-level sensitivity table
    (suppl_table.bmu2483_cutoff_sensitivity.tsv) ships alongside, and the one
    cutoff-robust reading is the "clinical only" fraction (~6-7% at every level).

ID FORM
    tmp/mmseqs99/out99_cluster.tsv uses BARE transcript ids ("BMU02483_001069-T1"),
    while the archived cluster100 file uses "ISOLATE|transcript". Both forms appear in
    this project, so every lookup here is built on the bare form and the isolate is
    recovered by LONGEST-PREFIX match against the 211-name list -- never by
    splitting on "_" (19 isolate names contain an underscore: FSSC_6, IMI_366586,
    NRRL_20438, ...).

Outputs (06.accessory-chrs/BMU02483_gene_content/)
    suppl_table.bmu2483_gene_origin_flows_99.tsv  species x isolate x chromosome, with counts
    suppl_table.bmu2483_ac_donor_detail.tsv       per-AC gene counts
    bmu2483_gene_origin_summary.md                the numbers, incl. the "other" composition

Usage: python3 scripts/bmu2483_gene_origin_sankey.py
"""

import collections
import csv
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TGS = os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies")
ANN = os.path.join(REPO, "02.annotation", "per_isolates")
OUT = os.path.join(REPO, "06.accessory-chrs", "BMU02483_gene_content")
CLUSTERS = os.path.join(REPO, "tmp", "mmseqs99", "out99_cluster.tsv")
FOCUS = "BMU02483"

ISO = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
       "BMU12683", "BMU12684", "BMU12685", "BMU12686",
       "Fu3.1", "Fu6.1", "SB1", "SX13"]


def isolate_list():
    names = [l.strip() for l in
             open(os.path.join(REPO, "05.pan-genome", "sample.solaniSC.211.list"))
             if l.strip()]
    # longest first, so "FSSC_6" wins over any shorter name that is also a prefix
    return sorted(names, key=len, reverse=True)


ISOLATES = isolate_list()


def isolate_of(gene_id):
    """bare transcript id -> its isolate, by longest-prefix match"""
    for n in ISOLATES:
        if gene_id.startswith(n + "_"):
            return n
    return None


def gene_compartment():
    """bare transcript id -> its contig, or 'CC'  (13 chromosome-level assemblies)

    The gff3 mRNA ids are already bare full ids ("BMU01933_000001-T1"), i.e. the same
    form the 99% clustering writes, so no prefix assembly is needed.
    """
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
            tid = None
            for a in f[8].split(";"):
                if a.startswith("ID="):
                    tid = a[3:]
                    break
            if tid:
                tid2c[tid] = f[0] if f[0] in acc else "CC"
    return tid2c


SPEC_ORDER = ["falciformis", "inter-species", "unknown"]
ISO_ORDER = ["from nearby env isolates", "from clinical isolates", "both", "unknown"]
CHROM_ORDER = ["ACs", "CCs", "unknown"]


def main():
    grp = {r["BMU_ID"]: r for r in csv.DictReader(
        open(os.path.join(REPO, "05.pan-genome", "bmu_ref_groups.211.txt")),
        delimiter="\t")}
    order = [l.strip() for l in open(os.path.join(TGS, "tree_order_allSC.txt")) if l.strip()]
    b = order.index(FOCUS)
    near = order[max(0, b - 10):b + 11]
    f_env_all = {k for k, v in grp.items()
                 if v["Species"] == "falciformis" and v["MedEnv"] == "env"}
    NEAR = {s for s in near if s in f_env_all}
    CLIN = {k for k, v in grp.items() if v["MedEnv"] == "med"}
    print(f"nearby environmental F. falciformis control: {len(NEAR)} isolates")
    print(f"clinical isolates: {len(CLIN)}")

    tid2c = gene_compartment()

    clusters = collections.defaultdict(list)
    with open(CLUSTERS) as fh:
        for line in fh:
            rep, _, mem = line.rstrip("\n").partition("\t")
            if mem:
                clusters[rep].append(mem)

    recs = []
    for _rep, members in clusters.items():
        mine = [m for m in members if m.startswith(FOCUS + "_")]
        if not mine:
            continue
        others = [m for m in members if not m.startswith(FOCUS + "_")]
        donors = {isolate_of(m) for m in others}
        donors.discard(None)

        species = {grp[d]["Species"] for d in donors if d in grp}
        # THREE classes: any non-falciformis donor makes the gene inter-species,
        # whether or not a falciformis donor is present too. "unknown" is reserved for
        # genes that cluster with NO other isolate at all.
        spec = ("unknown" if not donors else
                "inter-species" if (species - {"falciformis"}) else
                "falciformis")

        ne = bool(donors & NEAR)
        cli = bool(donors & CLIN)
        src = ("both" if (ne and cli) else
               "from nearby env isolates" if ne else
               "from clinical isolates" if cli else
               "unknown")

        acs = sorted({tid2c[m] for m in others if tid2c.get(m) not in (None, "CC")})
        chrom = ("ACs" if acs else
                 "CCs" if any(tid2c.get(m) == "CC" for m in others) else
                 "unknown")

        recs.extend({"gene": m, "spec": spec, "src": src, "chrom": chrom,
                     "ac": acs[0] if acs else None,
                     "n_donor": len(donors), "donors": donors} for m in mine)

    n = len(recs)
    assert n > 0, "no BMU02483 genes found -- check the cluster file id form"

    spec_c = collections.Counter(r["spec"] for r in recs)
    src_c = collections.Counter(r["src"] for r in recs)
    chrom_c = collections.Counter(r["chrom"] for r in recs)

    # the two residual classes must stay auditable: what is actually inside "unknown"?
    other_src = collections.Counter()
    for r in recs:
        if r["src"] != "unknown":
            continue
        if not r["donors"]:
            other_src["no donor at all"] += 1
        else:
            other_src["environmental donors outside the nearby window"] += 1
    # and how much of the inter-species class still keeps a falciformis donor
    inter = collections.Counter()
    for r in recs:
        if r["spec"] != "inter-species":
            continue
        inter["mixed (also has a falciformis donor)" if any(
            grp.get(d, {}).get("Species") == "falciformis" for d in r["donors"])
            else "non-falciformis only"] += 1

    ac_flow = collections.Counter(r["ac"] for r in recs if r["ac"])

    os.makedirs(OUT, exist_ok=True)
    flows = collections.Counter((r["spec"], r["src"], r["chrom"]) for r in recs)
    with open(os.path.join(OUT, "suppl_table.bmu2483_gene_origin_flows_99.tsv"), "w") as fh:
        fh.write("species_source\tisolate_source\tchromosome_source\tn_genes\n")
        for s in SPEC_ORDER:
            for i in ISO_ORDER:
                for c in CHROM_ORDER:
                    k = flows.get((s, i, c), 0)
                    if k:
                        fh.write(f"{s}\t{i}\t{c}\t{k}\n")

    with open(os.path.join(OUT, "suppl_table.bmu2483_ac_donor_detail.tsv"), "w") as fh:
        fh.write("AC\tgenes_with_a_donor_copy_on_this_AC\n")
        for k, v in ac_flow.most_common():
            fh.write(f"{k}\t{v}\n")

    with open(os.path.join(OUT, "bmu2483_gene_origin_summary.md"), "w") as fh:
        fh.write("# BMU02483 gene origins\n\n")
        fh.write(f"{n} genes; donors from the >=99%-identity clusters over 211 isolates.\n")
        fh.write("Re-clustered locally (90/95/99/100% from one protocol); the archived "
                 "cluster{90,95,100} files are NOT reproducible on this machine, so no "
                 "archived level is mixed in here.\n")
        fh.write(f"Nearby environmental control = the {len(NEAR)} F. falciformis env isolates "
                 "within +/-10 positions of BMU02483 in the core-genome tree order.\n\n")

        fh.write("## Species source\n\n| species | genes | % |\n|---|---|---|\n")
        for s in SPEC_ORDER:
            fh.write(f"| {s} | {spec_c[s]} | {100*spec_c[s]/n:.1f} |\n")
        fh.write("\n'inter-species' = at least one donor belongs to a species other than "
                 "F. falciformis (including donors whose species is Unknown); composition: " +
                 "; ".join(f"{k} {v}" for k, v in inter.most_common()) + ".\n")
        fh.write("'unknown' = the gene clusters with no other isolate at all.\n\n")

        fh.write("## Isolate source\n\n| source | genes | % |\n|---|---|---|\n")
        for s in ISO_ORDER:
            fh.write(f"| {s} | {src_c[s]} | {100*src_c[s]/n:.1f} |\n")
        fh.write("\n'unknown' composition: " +
                 "; ".join(f"{k} {v}" for k, v in other_src.most_common()) + ".\n")
        fh.write("\n'both' = donors of both kinds, so the direction is unresolvable; those "
                 "genes are kept as their own class rather than assigned to one side.\n\n")

        fh.write("## Known chromosome source (13 chromosome-level assemblies only)\n\n")
        fh.write("| compartment | genes | % |\n|---|---|---|\n")
        for c in CHROM_ORDER:
            fh.write(f"| {c} | {chrom_c[c]} | {100*chrom_c[c]/n:.1f} |\n")
        fh.write("\nACs take priority over CCs when a gene has donors on both, so each gene "
                 "contributes one compartment. The %d ACs are pooled into one node; the "
                 "per-AC breakdown is in suppl_table.bmu2483_ac_donor_detail.tsv.\n"
                 % len(ac_flow))

    print(f"\nBMU02483 genes: {n}")
    for title, order_, cnt in (("species", SPEC_ORDER, spec_c),
                                ("isolate source", ISO_ORDER, src_c),
                                ("chromosome", CHROM_ORDER, chrom_c)):
        print(f"{title}:")
        for k in order_:
            print(f"   {k:<26} {cnt[k]:6d}  {100*cnt[k]/n:5.1f}%")
    print(f"\nwrote suppl_table.bmu2483_gene_origin_flows_99.tsv and bmu2483_gene_origin_summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
