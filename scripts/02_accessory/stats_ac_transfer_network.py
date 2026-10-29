#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AC transfer-network profiling.

For every accessory chromosome (AC) of the 13 TGS assemblies, compute its
transfer profile from the 100%-identity mmseqs clusters (all 211 isolates):
  - n_transferred_genes: genes on this AC whose cluster spans >=2 chromosomes
  - n_partner_chromosomes: distinct 13-TGS chromosomes sharing >=1 cluster
  - genes per partner (concentration of transfers)
  - n_partner_isolates / species spread (cross-species flag)

Partners are restricted to the chromosome-level 13 TGS assemblies (the only
genomes with chromosome structure and AC calls).

Outputs (06.accessory-chrs/gene_transfer_211/):
  suppl_table.ac_transfer_profile.tsv
"""

import argparse
import os
from collections import Counter, defaultdict

import numpy as np
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
    p = argparse.ArgumentParser(description="AC transfer network profiling")
    p.add_argument("--cluster100", required=True)
    p.add_argument("--affiliation", required=True)
    p.add_argument("--accessory-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # contig types + lengths (13 TGS)
    contig_type, contig_len = {}, {}
    for iso in ISOLATES:
        acc = set(open(os.path.join(args.accessory_dir,
                                    f"{iso}.accessory_chrs.id")).read().split())
        with open(os.path.join(args.length_dir, iso, "repeatmasker",
                               f"{iso}.assembly.fasta.masked.length")) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    c, L = parts[0], int(parts[1])
                    contig_len[c] = L
                    contig_type[c] = ("accessory" if c in acc
                                      else "core" if L >= 500000 else "scaffold")

    aff = {}
    with open(args.affiliation) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                aff[parts[0]] = parts[2]

    # build clusters: rep -> members [(iso, contig)]
    print("loading clusters...")
    rep_members = defaultdict(list)
    iso_set = set(ISOLATES)
    with open(args.cluster100) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            iso, _, trans = parts[1].partition("|")
            if iso not in iso_set:
                continue
            contig = aff.get(trans)
            if contig is None:
                continue
            rep_members[parts[0]].append((iso, trans, contig))
    print(f"clusters with 13-TGS members: {len(rep_members)}")

    # per-AC transfer records: for each gene on an AC, its partner chromosomes
    # (13-TGS) in the same cluster
    ac_genes = defaultdict(list)  # AC -> [(gene, partner_contigs, cluster_contigs)]
    for rep, members in rep_members.items():
        contigs = {c for _, _, c in members}
        if len(contigs) < 2:
            continue
        if not any(contig_type.get(c) == "accessory" for c in contigs):
            continue
        for iso, trans, contig in members:
            if contig_type.get(contig) == "accessory":
                partners = sorted({(miso, mcontig) for miso, _, mcontig in members
                                   if mcontig != contig})
                ac_genes[contig].append((f"{iso}|{trans}", partners,
                                         len(contigs)))

    rows = []
    for ac, records in ac_genes.items():
        n_genes = len(records)
        # distinct partner chromosomes (13-TGS)
        partner_chroms = set()
        partner_genes = Counter()
        partner_isolates = set()
        partner_species = set()
        for gene, partners, n_clust_contigs in records:
            for piso, pcontig in partners:
                partner_chroms.add(pcontig)
                partner_genes[pcontig] += 1
                partner_isolates.add(piso)
                partner_species.add(SPECIES[piso])
        n_partners = len(partner_chroms)
        # concentration: fraction of transferred genes going to the top partner
        top_frac = (max(partner_genes.values()) / n_genes
                    if n_genes else 0)
        # block score: genes concentrated on few partners
        conc = n_genes / max(n_partners, 1)
        rows.append({
            "AC": ac,
            "Isolate": ac.split("_")[0],
            "Species": SPECIES[ac.split("_")[0]],
            "Length": contig_len.get(ac, np.nan),
            "n_transferred_genes": n_genes,
            "n_partner_chromosomes": n_partners,
            "n_partner_isolates": len(partner_isolates),
            "n_partner_species": len(partner_species),
            "genes_per_partner": conc,
            "top_partner_frac": top_frac,
            "cross_species": len(partner_species) > 1,
        })
    df = pd.DataFrame(rows).sort_values("n_transferred_genes", ascending=False)
    df.to_csv(os.path.join(args.output_dir,
                           "suppl_table.ac_transfer_profile.tsv"),
              sep="\t", index=False)
    print(f"ACs profiled: {len(df)}")
    print(df[["AC", "Species", "n_transferred_genes", "n_partner_chromosomes",
              "n_partner_species", "genes_per_partner", "top_partner_frac"]]
          .round(2).to_string(index=False))


if __name__ == "__main__":
    main()
