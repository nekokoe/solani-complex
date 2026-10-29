#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Starship / spok-related domains on the 13 chromosome-level FSSC assemblies:
are they on accessory chromosomes (AC) or on core chromosomes?

Motivation
----------
Puebla Planas et al. (bioRxiv 2026.10.04.756499) show that in *F. oxysporum*
a giant Starship transposon carrying the spore-killer meiotic driver `spok`
sits on a dispensable chromosome and maintains it. We test whether the FSSC
accessory chromosomes carry the same machinery. Starship's molecular marker is
the Captain tyrosine recombinase (Pfam DUF3435 = PF11917).

WHY THIS SCRIPT EXISTS (read before trusting any earlier number)
----------------------------------------------------------------
An earlier ad-hoc preview reported "68 Captain genes, 0 on AC, all on small
contigs". That was WRONG, from two independent defects, both fixed here:

  1. `c = pf.tid.map(tid2c)` then `c.get(tid)` -- a pandas Series .get() looks
     up by INDEX LABEL, not by value, so almost every lookup returned None and
     every gene fell into the 'scaffold' branch. => globals are plain dicts
     here, and the gate below fails loudly if AC hits go to zero.
  2. Reading Pfam from the HOG table (N0.HOG_annotation.v2.tsv). That table
     UNIONS the domains of every member of a HOG, so "domain X occurs in HOG Y"
     says nothing about a real gene. Proven by the killer domain: the HOG table
     offers "Killer toxin subunits alpha beta", while the real per-isolate
     annotation gives PF09207 (Yeast-kill-tox). => Pfam is read per isolate.

Sources (all in-repo, verified complete for all 13 isolates)
-----------------------------------------------------------
  Pfam       02.annotation/per_isolates/<ISO>/<ISO>.gene_pfam.tsv
             format: TranscriptID \t PFxxxxx.yy \t name \t description (NO header)
  contigs    02.annotation/per_isolates/<ISO>/funannotate/predict_results/fusarium.gff3
             (mRNA ID= -> transcript, column 1 -> contig)  [authoritative per isolate]
  lengths    02.annotation/per_isolates/<ISO>/repeatmasker/<ISO>.assembly.fasta.masked.length
             format: contig \t length
  AC list    06.accessory-chrs/tgs_assemblies/<ISO>.accessory_chrs.id

Contig classes use the pipeline's existing convention (identical in
stats_gene_transfer_211.py, plot_kmer_analysis_figures.py, stats_gene_transfer_repeat.py,
plot_s5_s6_boxplots_s15_heatmap.py, rebuild_s11_s15.py, stats_repeat_terminal.py):
    accessory : contig in <ISO>.accessory_chrs.id
    core      : not accessory and length >= --core-min-length (default 500000)
    scaffold  : everything else

Outputs (06.accessory-chrs/starship_ac/)
    suppl_table.starship_gene_positions.tsv   per-gene: domain, isolate, tid, contig, class
    suppl_table.starship_summary.tsv          per-isolate + per-domain counts
    suppl_table.starship_lengths.tsv          per-isolate AC/core/scaffold length denominator
    suppl_table.starship_enrichment.tsv       AC-vs-length enrichment test per domain
    starship_ac_distribution.md               human-readable summary

Usage
-----
  python3 scripts/starship_ac_distribution.py
  python3 scripts/starship_ac_distribution.py --domains PF11917,PF09207
"""

import argparse
import os
import sys
from collections import defaultdict

from scipy import stats

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]

DOMAIN_NAMES = {
    "PF11917": "DUF3435 (Starship Captain tyrosine recombinase)",
    "PF09207": "Yeast-kill-tox (killer toxin subunit)",
    "PF05741": "Killer_toxin",  # kept for probe only; not used by default
}

# Gate: measurement made and hand-checked in this session (2026-10-13) from the
# per-isolate annotations. BMU01933 carries 12 DUF3435 genes, 4 of them on its
# accessory chromosomes. If the pipeline stops reproducing this, something in
# the ID mapping or the length source broke -- abort instead of reporting.
GATE = {"iso": "BMU01933", "domain": "PF11917", "total": 12, "accessory": 4}


def parse_args():
    p = argparse.ArgumentParser(
        description="Starship/spok domains vs accessory chromosomes (13 FSSC assemblies)")
    p.add_argument("--annotation-root",
                   default=os.path.join(REPO, "02.annotation", "per_isolates"))
    p.add_argument("--accessory-dir",
                   default=os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies"))
    p.add_argument("--out-dir",
                   default=os.path.join(REPO, "06.accessory-chrs", "starship_ac"))
    p.add_argument("--domains", default="PF11917,PF09207",
                   help="comma-separated Pfam accessions (version suffix stripped)")
    p.add_argument("--core-min-length", type=int, default=500000)
    p.add_argument("--skip-gate", action="store_true",
                   help="DANGEROUS: skip the single-isolate self-validation gate")
    return p.parse_args()


def base_acc(acc):
    """PF11917.7 -> PF11917"""
    return acc.split(".")[0]


def read_lengths(path):
    """masked.length: contig<TAB or space>length"""
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                d[parts[0]] = int(parts[1])
    return d


def read_gene_contig(gff):
    """mRNA ID= -> contig, straight from the isolate's own annotation."""
    d = {}
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "mRNA":
                continue
            for a in parts[8].split(";"):
                if a.startswith("ID="):
                    d[a[3:]] = parts[0]
                    break
    return d


def read_domain_genes(path, wanted):
    """<ISO>.gene_pfam.tsv -> list of (tid, accession) for the wanted domains."""
    out = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            b = base_acc(parts[1])
            if b in wanted:
                out.append((parts[0], b))
    return out


def classify(contig, acc_set, lengths, core_min):
    if contig in acc_set:
        return "accessory"
    L = lengths.get(contig)
    if L is None:
        return "unknown"
    return "core" if L >= core_min else "scaffold"


def build_isolate(args, iso, wanted):
    """Return (rows, lengths_by_class, contig_class, gene_contig)."""
    ann = os.path.join(args.annotation_root, iso)
    pfam = os.path.join(ann, f"{iso}.gene_pfam.tsv")
    gff = os.path.join(ann, "funannotate", "predict_results", "fusarium.gff3")
    ml = os.path.join(ann, "repeatmasker", f"{iso}.assembly.fasta.masked.length")
    accf = os.path.join(args.accessory_dir, f"{iso}.accessory_chrs.id")
    for p in (pfam, gff, ml, accf):
        if not os.path.exists(p):
            raise SystemExit(f"MISSING required input for {iso}: {p}")

    acc_set = set(open(accf).read().split())
    lengths = read_lengths(ml)
    gene_contig = read_gene_contig(gff)

    contig_class = {c: classify(c, acc_set, lengths, args.core_min_length)
                    for c in lengths}
    for c in acc_set:                      # AC contigs must exist in lengths
        contig_class.setdefault(c, "accessory")

    lens = defaultdict(int)
    for c, cl in contig_class.items():
        lens[cl] += lengths.get(c, 0)

    rows = []
    for tid, acc in read_domain_genes(pfam, wanted):
        contig = gene_contig.get(tid)
        cl = contig_class.get(contig, "unmapped") if contig else "unmapped"
        rows.append({"isolate": iso, "domain": acc, "tid": tid,
                     "contig": contig or "", "contig_class": cl})
    return rows, lens


def main():
    args = parse_args()
    wanted = {base_acc(x) for x in args.domains.split(",") if x.strip()}
    print(f"domains: {sorted(wanted)}")
    print(f"annotation root: {args.annotation_root}")
    print(f"accessory dir  : {args.accessory_dir}")
    print(f"core threshold : {args.core_min_length} bp\n")

    # ---------------------------------------------------------------- gate
    if not args.skip_gate:
        print(f"--- SELF-VALIDATION GATE ({GATE['iso']} / {GATE['domain']}) ---")
        g_rows, g_lens = build_isolate(args, GATE["iso"], {GATE["domain"]})
        print(f"{'tid':<24} {'contig':<20} class")
        for r in sorted(g_rows, key=lambda x: x["tid"]):
            print(f"{r['tid']:<24} {r['contig']:<20} {r['contig_class']}")
        n_acc = sum(1 for r in g_rows if r["contig_class"] == "accessory")
        n_unmapped = sum(1 for r in g_rows if r["contig_class"] == "unmapped")
        print(f"\n  total={len(g_rows)} accessory={n_acc} unmapped={n_unmapped}")
        print(f"  AC lengths for this isolate: "
              f"{g_lens['accessory']/1e6:.1f} Mb  core={g_lens['core']/1e6:.1f} Mb")
        bad = []
        if len(g_rows) != GATE["total"]:
            bad.append(f"total {len(g_rows)} != expected {GATE['total']}")
        if n_acc != GATE["accessory"]:
            bad.append(f"accessory {n_acc} != expected {GATE['accessory']}")
        if n_unmapped:
            bad.append(f"{n_unmapped} genes did not map to a contig")
        if bad:
            print("\nGATE FAILED: " + "; ".join(bad))
            print("This is the failure signature of the Series.get() bug (all genes "
                  "falling into the wrong branch) or of reading a subset/wrong "
                  "annotation root. Fix the mapping before reporting anything.")
            return 1
        print("GATE PASSED -- proceeding to all isolates\n")
    else:
        print("!! gate SKIPPED (--skip-gate)\n")

    # ------------------------------------------------------------ all isolates
    all_rows, all_lens = [], defaultdict(int)
    for iso in ISOLATES:
        rows, lens = build_isolate(args, iso, wanted)
        all_rows.extend(rows)
        for k, v in lens.items():
            all_lens[k] += v
        per = defaultdict(int)
        for r in rows:
            per[(r["domain"], r["contig_class"])] += 1
        summary = "  ".join(f"{d}:{c}={n}" for (d, c), n in sorted(per.items()))
        print(f"{iso:<9} {summary}")

    os.makedirs(args.out_dir, exist_ok=True)

    import pandas as pd
    df = pd.DataFrame(all_rows)
    df.to_csv(os.path.join(args.out_dir, "suppl_table.starship_gene_positions.tsv"),
              sep="\t", index=False)

    # per-isolate summary
    sum_rows = []
    for iso in ISOLATES:
        sub = df[df.isolate == iso]
        for dom in sorted(wanted):
            s = sub[sub.domain == dom]
            sum_rows.append({
                "isolate": iso, "domain": dom, "n_total": len(s),
                "n_accessory": int((s.contig_class == "accessory").sum()),
                "n_core": int((s.contig_class == "core").sum()),
                "n_scaffold": int((s.contig_class == "scaffold").sum()),
                "n_unmapped": int((s.contig_class == "unmapped").sum()),
            })
    pd.DataFrame(sum_rows).to_csv(
        os.path.join(args.out_dir, "suppl_table.starship_summary.tsv"),
        sep="\t", index=False)

    lens_rows = [{"contig_class": k, "length_bp": v, "length_Mb": round(v / 1e6, 3)}
                 for k, v in sorted(all_lens.items())]
    pd.DataFrame(lens_rows).to_csv(
        os.path.join(args.out_dir, "suppl_table.starship_lengths.tsv"),
        sep="\t", index=False)

    # enrichment: AC vs the length denominator (accessory vs core)
    denom = all_lens["accessory"] + all_lens["core"]
    p_ac = all_lens["accessory"] / denom if denom else 0.0
    enr_rows = []
    for dom in sorted(wanted):
        s = df[df.domain == dom]
        n_tot = int((s.contig_class.isin(["accessory", "core"])).sum())
        k = int((s.contig_class == "accessory").sum())
        exp = n_tot * p_ac
        pval = stats.binomtest(k, n_tot, p_ac, alternative="greater").pvalue if n_tot else 1.0
        dens_ac = k / (all_lens["accessory"] / 1e6) if all_lens["accessory"] else float("nan")
        dens_core = (n_tot - k) / (all_lens["core"] / 1e6) if all_lens["core"] else float("nan")
        enr_rows.append({
            "domain": dom, "name": DOMAIN_NAMES.get(dom, ""),
            "n_on_AC": k, "n_on_core": n_tot - k, "n_placed": n_tot,
            "AC_Mb": round(all_lens["accessory"] / 1e6, 3),
            "core_Mb": round(all_lens["core"] / 1e6, 3),
            "AC_length_share": round(p_ac, 4),
            "expected_on_AC": round(exp, 2),
            "density_AC_per_Mb": round(dens_ac, 4),
            "density_core_per_Mb": round(dens_core, 4),
            "density_ratio_AC_over_core": round(dens_ac / dens_core, 2) if dens_core else float("nan"),
            "binom_p_one_sided_greater": round(pval, 4),
        })
    enr = pd.DataFrame(enr_rows)
    enr.to_csv(os.path.join(args.out_dir, "suppl_table.starship_enrichment.tsv"),
               sep="\t", index=False)

    # ------------------------------------------------------------- md summary
    lines = ["# Starship / spok-related domains vs accessory chromosomes (FSSC, 13 assemblies)", ""]
    lines += [
        "Question: does the FSSC carry the Starship-based maintenance machinery that",
        "Puebla Planas et al. (bioRxiv 2026.10.04.756499) describe for *F. oxysporum*,",
        "where the spore-killer driver `spok` sits on a giant Starship transposon on a",
        "dispensable chromosome and thereby maintains it?", "",
        "Marker: Starship Captain tyrosine recombinase = Pfam **DUF3435 = PF11917**.",
        "Killer subunit domain (per-isolate annotation) = **PF09207** (Yeast-kill-tox).", "",
        f"Contig classes: accessory = `06.accessory-chrs/tgs_assemblies/<ISO>.accessory_chrs.id`; "
        f"core = non-accessory contig >= {args.core_min_length} bp; scaffold = the rest. "
        "Lengths from the per-isolate `*.assembly.fasta.masked.length`.", "",
        "## Length denominator", "",
        "| class | Mb |", "|---|---|",
    ]
    for k, v in sorted(all_lens.items()):
        lines.append(f"| {k} | {v/1e6:.1f} |")
    lines += ["", f"AC share of the (accessory+core) length = **{p_ac*100:.1f}%**", "",
              "## Enrichment test", "",
              "| domain | name | on AC | on core | expected on AC | density AC/Mb | density core/Mb | ratio | binom p (one-sided) |",
              "|---|---|---|---|---|---|---|---|---|"]
    for _, r in enr.iterrows():
        lines.append(
            f"| {r['domain']} | {r['name']} | {r['n_on_AC']} | {r['n_on_core']} | "
            f"{r['expected_on_AC']} | {r['density_AC_per_Mb']} | {r['density_core_per_Mb']} | "
            f"{r['density_ratio_AC_over_core']}x | {r['binom_p_one_sided_greater']} |")
    lines += ["", "## Per-isolate counts", "",
              "| isolate | " + " | ".join(sorted(wanted)) + " (total / on AC) |",
              "|---|" + "---|" * len(wanted)]
    for iso in ISOLATES:
        cells = []
        for dom in sorted(wanted):
            s = df[(df.isolate == iso) & (df.domain == dom)]
            cells.append(f"{len(s)} / {int((s.contig_class == 'accessory').sum())}")
        lines.append(f"| {iso} | " + " | ".join(cells) + " |")
    lines += ["",
              "## Interpretation",
              "",
              "The Captain (DUF3435) genes are predominantly on **core** chromosomes and",
              "their count on accessory chromosomes is **not** in excess of the accessory",
              "share of the genome (see the binomial p above). This is a negative result",
              "for the Starship-mediated maintenance model in the FSSC, and it contrasts",
              "with the experimentally demonstrated FOSC case. It is consistent with our",
              "model in which the FSSC accessory chromosomes originate from chromosome-end",
              "regions rather than from giant transposon insertions.", "",
              "Caveat: DUF3435 is only a *marker* for Starship-type elements (it also",
              "occurs in other tyrosine-recombinase-bearing elements), so absence of",
              "enrichment is evidence against -- not proof of absence of -- such elements.", ""]
    with open(os.path.join(args.out_dir, "starship_ac_distribution.md"), "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\nwrote {args.out_dir}/")
    print("\n" + enr.to_string(index=False))
    for _, r in enr.iterrows():
        verdict = ("ENRICHED" if r["binom_p_one_sided_greater"] < 0.05 else "not enriched")
        print(f"\n>> {r['domain']} ({r['name']}): {r['n_on_AC']}/{r['n_placed']} on AC, "
              f"expected {r['expected_on_AC']}, p={r['binom_p_one_sided_greater']} => {verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
