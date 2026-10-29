#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AC hub and cross-species analysis over ALL 211 isolates.

The AC nodes are the 13 TGS accessory chromosomes (only these have AC status),
but their gene-sharing REACH is measured over the full 211-genome set using
the 100%-identity clusters and the per-isolate species map.

Hub definitions:
  within-species hub: transfers many genes to many isolates of its own species
                      (own-species genes x own-species isolates reached)
  SC-wide hub:        capable of spreading genes ACROSS SPECIES
                      (genes x isolates reached x species reached)

Q3: cross-species-transferred genes = 13-TGS AC genes whose cluster contains
    members in isolates of >=2 species (211 context).

Outputs (06.accessory-chrs/gene_transfer_211/):
  suppl_table.ac_hub_211.tsv               per-AC hub profile (211 reach)
  suppl_table.transfer_gene_cross211.tsv   per-gene cross-species flag (211)
"""

import argparse
import os
from collections import Counter, defaultdict

import numpy as np
import pandas as pd

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]


def parse_args():
    p = argparse.ArgumentParser(description="AC hub / cross-species on 211")
    p.add_argument("--cluster100", required=True)
    p.add_argument("--affiliation", required=True)
    p.add_argument("--annotation-root", default=None,
                   help="build the full 211-genome transcript->contig map from "
                        "the per-isolate funannotate GFF3s")
    p.add_argument("--species-map", required=True)
    p.add_argument("--accessory-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    species = {}
    with open(args.species_map) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                species[parts[0]] = parts[1]

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
        for iso in sorted(os.listdir(args.annotation_root)):
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

    # clusters over all 211 (contig set complete); track per-AC partner isolates
    rep_members = defaultdict(list)
    with open(args.cluster100) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            iso, _, trans = parts[1].partition("|")
            c = aff.get(trans)
            if c is None:
                continue
            rep_members[parts[0]].append((iso, trans, c))
    print(f"clusters with mapped members: {len(rep_members)}")

    # per-AC: partner isolates (211) + per-isolate gene counts + species
    ac_iso_genes = defaultdict(Counter)   # AC -> iso -> n_genes
    ac_genes_cross = defaultdict(list)    # AC -> list of (gene, partner_species_set)
    ac_species_partners = defaultdict(set)
    iso_set = set(ISOLATES)
    for rep, members in rep_members.items():
        contigs = {c for _, _, c in members}
        if len(contigs) < 2:
            continue
        if not any(contig_type.get(c) == "accessory" for c in contigs):
            continue
        # cluster gene IDs per contig
        genes_by_contig = defaultdict(list)
        for iso, trans, c in members:
            genes_by_contig[c].append(f"{iso}|{trans}")
        for c in genes_by_contig:
            if contig_type.get(c) != "accessory":
                continue
            for iso, trans, _ in members:
                ac_iso_genes[c][iso] += 1
            # species set of partner members (by isolate)
            ps = set()
            for iso, _, pc in members:
                if pc != c:
                    ps.add(species[iso])
            ac_species_partners[c] |= ps
            for gene in genes_by_contig[c]:
                ac_genes_cross[c].append((gene, tuple(sorted(ps))))

    rows = []
    for ac, iso_genes in ac_iso_genes.items():
        n_genes = sum(iso_genes.values())
        partner_iso = set(iso_genes)
        own_sp = species[ac.split("_")[0]]
        n_iso_own = sum(1 for i in partner_iso if species.get(i) == own_sp)
        n_iso_other = len(partner_iso) - n_iso_own
        n_species = len(ac_species_partners.get(ac, set()))
        own_genes = sum(v for i, v in iso_genes.items() if species.get(i) == own_sp)
        other_genes = n_genes - own_genes
        rows.append({
            "AC": ac, "Isolate": ac.split("_")[0], "Species": own_sp,
            "n_shared_genes_211": n_genes,
            "n_partner_isolates_211": len(partner_iso),
            "n_partner_isolates_own_species": n_iso_own,
            "n_partner_isolates_other_species": n_iso_other,
            "n_partner_species_211": n_species,
            "n_genes_own_species": own_genes,
            "n_genes_other_species": other_genes,
            "within_species_hub_score": own_genes * max(n_iso_own, 1),
            "sc_wide_hub_score": n_genes * len(partner_iso) * max(n_species, 1),
            "cross_species": n_species > 1,
        })
    hub = pd.DataFrame(rows)
    hub = hub.sort_values("sc_wide_hub_score", ascending=False)
    hub.to_csv(os.path.join(args.output_dir, "suppl_table.ac_hub_211.tsv"),
               sep="\t", index=False)
    print(f"ACs: {len(hub)}")
    print("\n=== top 12 ACs by SC-wide hub score (211 reach) ===")
    print(hub.head(12)[["AC", "Species", "n_shared_genes_211",
                        "n_partner_isolates_211", "n_partner_species_211",
                        "n_genes_own_species", "n_genes_other_species",
                        "within_species_hub_score", "sc_wide_hub_score"]]
          .round(0).to_string(index=False))
    print("\n=== top 10 by within-species hub score ===")
    print(hub.nlargest(10, "within_species_hub_score")[
        ["AC", "Species", "n_shared_genes_211", "n_partner_isolates_own_species",
         "within_species_hub_score"]].round(0).to_string(index=False))

    # per-gene cross-species flag (211 context)
    rows = []
    for ac, records in ac_genes_cross.items():
        for gene, ps in records:
            rows.append({"Gene": gene, "AC": ac,
                         "Cross_species_211": len(ps) > 1,
                         "Partner_species_211": ";".join(sorted(ps))})
    gdf = pd.DataFrame(rows)
    gdf.to_csv(os.path.join(args.output_dir,
                            "suppl_table.transfer_gene_cross211.tsv"),
               sep="\t", index=False)
    print(f"\nper-gene cross-species (211): {gdf.Cross_species_211.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
