#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-gene AC-class and cross-species annotation of transfer-involved genes.

For every gene on a 13-TGS accessory chromosome in a transfer cluster:
  - AC_class: Type1_block / Type2_scattered (from ac_transfer_profile.tsv)
  - cross_species: the gene's cluster contains members from >=2 species
    (solani / falciformis / keratoplastica) among the 13 TGS assemblies
  - partner_species: set of species of partner chromosomes

Output: suppl_table.transfer_gene_ac_class.tsv
"""

import argparse
import os
from collections import defaultdict

import pandas as pd

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]
SPECIES = {"BMU01933": "falciformis", "BMU02028": "solani",
           "BMU02474": "solani", "BMU03267": "falciformis",
           "BMU06673": "keratoplastica", "BMU12683": "solani",
           "BMU12684": "solani", "BMU12685": "solani",
           "BMU12686": "solani", "Fu3.1": "falciformis",
           "Fu6.1": "keratoplastica", "SB1": "solani", "SX13": "solani"}


def parse_args():
    p = argparse.ArgumentParser(description="Per-gene AC class + cross-species")
    p.add_argument("--cluster100", required=True)
    p.add_argument("--affiliation", required=True,
                   help="full transcript->contig mapping; if the small 13-isolate "
                        "affiliation is given, pass --annotation-root to build the "
                        "211-genome mapping instead")
    p.add_argument("--annotation-root", default=None)
    p.add_argument("--accessory-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--ac-profile", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    contig_type = {}
    for iso in ISOLATES:
        acc = set(open(os.path.join(args.accessory_dir,
                                    f"{iso}.accessory_chrs.id")).read().split())
        with open(os.path.join(args.length_dir, iso, "repeatmasker",
                               f"{iso}.assembly.fasta.masked.length")) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    c, L = parts[0], int(parts[1])
                    contig_type[c] = ("accessory" if c in acc
                                      else "core" if L >= 500000 else "scaffold")
    aff = {}
    if args.annotation_root:
        # build from all 211 GFF3s (same as stats_gene_transfer_211.py)
        iso_list = [x.strip() for x in open(
            os.path.join(args.annotation_root, "..", "..", "02.annotation",
                         "sample.solaniSC.211.list"))] if False else None
        import glob
        iso_list = sorted(os.listdir(args.annotation_root))
        for iso in iso_list:
            gff = os.path.join(args.annotation_root, iso, "funannotate",
                               "predict_results", "fusarium.gff3")
            if not os.path.exists(gff):
                continue
            with open(gff) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 9 or parts[2] != "mRNA":
                        continue
                    tid = None
                    for a in parts[8].split(";"):
                        if a.startswith("ID="):
                            tid = a[3:]
                            break
                    if tid:
                        aff[tid] = parts[0]
    else:
        with open(args.affiliation) as f:
            next(f)
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    aff[parts[0]] = parts[2]

    ac_class = {}
    prof = pd.read_csv(args.ac_profile, sep="\t")
    ac_class = dict(zip(prof["AC"], prof["AC_class"]))

    iso_set = set(ISOLATES)
    # The archived clustering writes "ISOLATE|transcript"; a fresh mmseqs run writes the
    # bare transcript id. A parser that assumes one form matches nothing for the other
    # WITHOUT any error, so both are handled (isolate by longest-prefix match, because 19
    # isolate names contain an underscore).
    ISO_NAMES = sorted([l.strip() for l in
                        open("05.pan-genome/sample.solaniSC.211.list") if l.strip()],
                       key=len, reverse=True)

    def split_member(m):
        if "|" in m:
            i, _, t = m.partition("|")
            return i, t
        for n in ISO_NAMES:
            if m.startswith(n + "_"):
                return n, m
        return None, m

    rep_members = defaultdict(list)
    with open(args.cluster100) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            iso, trans = split_member(parts[1])
            c = aff.get(trans)
            if c is None:
                continue
            # keep ALL 211 members so the cluster contig set is complete;
            # only 13-TGS genes are output later
            rep_members[parts[0]].append((iso, trans, c))

    rows = []
    for rep, members in rep_members.items():
        contigs = {c for _, _, c in members}
        if len(contigs) < 2:
            continue
        if not any(contig_type.get(c) == "accessory" for c in contigs):
            continue
        species = {SPECIES[iso] for iso, _, _ in members if iso in iso_set}
        cross = len(species) > 1
        for iso, trans, c in members:
            if iso not in iso_set or contig_type.get(c) != "accessory":
                continue
            partner_species = {SPECIES[piso] for piso, _, pc in members
                               if pc != c and piso in iso_set}
            rows.append({
                "Gene": f"{iso}|{trans}", "Isolate": iso,
                "Transcript": trans, "AC": c,
                "Species": SPECIES[iso],
                "AC_class": ac_class.get(c, "NA"),
                "Cluster_contigs": len(contigs),
                "Cluster_species": sorted(species),
                "Cross_species": cross,
                "Partner_species": ";".join(sorted(partner_species)),
            })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(args.output_dir,
                           "suppl_table.transfer_gene_ac_class.tsv"),
              sep="\t", index=False)
    print(f"transfer genes on ACs: {len(df)}")
    print("AC_class:", df.AC_class.value_counts().to_dict())
    print("cross-species:", df.Cross_species.value_counts().to_dict())


if __name__ == "__main__":
    main()
