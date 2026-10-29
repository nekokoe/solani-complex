#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Position of gene-exchange events along the chromosome, measured from the nearest
chromosome end (0 = the end), for the two accessory-chromosome exchange classes:

    A-A  accessory-to-accessory exchange      -> positions on accessory chromosomes
    C-A  accessory-to-core exchange           -> positions on core chromosomes
                                                  (the receiving, subtelomeric side)

Three variants are drawn so the honest one can be chosen:

    v1_raw         the two density curves alone  (what was asked for)
    v2_with_bg     the two curves plus the gene-density background of the SAME
                   chromosome class (dashed)
    v3_normalized  observed / background, i.e. how much more often an exchange
                   event sits at that distance than a random gene would

WHY v2/v3 EXIST: genes themselves cluster near chromosome ends, and they do so
far more strongly on accessory chromosomes (10.5% of accessory genes lie within
50 kb of an end versus 1.2% of core genes). A raw two-curve plot therefore makes
A-A and C-A look similar, when in fact C-A is ~15x end-enriched over its own
background and A-A is slightly end-DEPLETED over its own. Any statement of the
form "exchange is concentrated at the chromosome ends" must be made against the
gene-density baseline, not against zero.

Method
    distance = min(start, contig_length - end + 1)      [so 0 is the nearest end]
    density  = events per kb in 5 kb bins over 0..--max-kb, normalised so that the
               area under the curve is 1 (i.e. a probability density)
    background = every mRNA of the same contig class, treated the same way

Inputs (in-repo)
    02.annotation/per_isolates/<ISO>/funannotate/predict_results/fusarium.gff3
    02.annotation/per_isolates/<ISO>/repeatmasker/<ISO>.assembly.fasta.masked.length
    06.accessory-chrs/tgs_assemblies/<ISO>.accessory_chrs.id
    06.accessory-chrs/gene_transfer_211/suppl_table.gene_transfer_211.tsv

Outputs (06.accessory-chrs/gene_transfer_211/terminal_density/)
    transfer_terminal_positions.tsv   per-gene: class, contig, class of contig, distance
    transfer_terminal_density.tsv     per-bin: density for each series
    fig_terminal_density_v1_raw.pdf/.png
    fig_terminal_density_v2_with_bg.pdf/.png
    fig_terminal_density_v3_normalized.pdf/.png
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "font.size": 9,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
})

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, "06.accessory-chrs", "gene_transfer_211", "terminal_density")

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]

C_AA, C_CA = "#C44E52", "#3B7DD8"   # A-A red, C-A blue
BG = "#9A9A9A"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--max-kb", type=int, default=200, help="x-axis extent from the end")
    p.add_argument("--bin-kb", type=int, default=5)
    p.add_argument("--core-min-length", type=int, default=500000)
    return p.parse_args()


def build_table(args):
    ann = os.path.join(REPO, "02.annotation", "per_isolates")
    tgs = os.path.join(REPO, "06.accessory-chrs", "tgs_assemblies")

    clen, cls = {}, {}
    for iso in ISOLATES:
        acc = set(open(os.path.join(tgs, f"{iso}.accessory_chrs.id")).read().split())
        for line in open(os.path.join(ann, iso, "repeatmasker",
                                      f"{iso}.assembly.fasta.masked.length")):
            f = line.split()
            if len(f) < 2:
                continue
            clen[f[0]] = int(f[1])
            cls[f[0]] = ("accessory" if f[0] in acc
                         else "core" if int(f[1]) >= args.core_min_length else "scaffold")

    rows = []
    for iso in ISOLATES:
        gff = os.path.join(ann, iso, "funannotate", "predict_results", "fusarium.gff3")
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
            if not tid or f[0] not in clen:
                continue
            start, end = int(f[3]), int(f[4])
            rows.append((iso, tid, f[0], cls[f[0]], clen[f[0]],
                         min(start, clen[f[0]] - end + 1)))
    g = pd.DataFrame(rows, columns=["Isolate", "Transcript", "Contig", "contig_class",
                                    "contig_len", "dist_from_end"])
    t = pd.read_csv(os.path.join(REPO, "06.accessory-chrs", "gene_transfer_211",
                                 "suppl_table.gene_transfer_211.tsv"), sep="\t",
                    usecols=["Isolate", "Transcript", "Transfer_type",
                             "AC_related_transfer"])
    m = g.merge(t, on=["Isolate", "Transcript"], how="left")
    print(f"mRNAs: {len(m)}  | with a transfer class: {m.Transfer_type.notna().sum()}")
    return m


def density(d, n_total, args):
    """probability density per kb in fixed bins over 0..max-kb"""
    edges = np.arange(0, args.max_kb + args.bin_kb, args.bin_kb)
    h, _ = np.histogram(d / 1000.0, bins=edges)
    width = args.bin_kb / args.max_kb            # area under the curve = 1
    return h / (n_total * width)


def main():
    args = parse_args()
    os.makedirs(OUT, exist_ok=True)
    m = build_table(args)

    series = {
        "A-A": m[(m.Transfer_type == "A-A") & (m.contig_class == "accessory")],
        "C-A": m[(m.Transfer_type == "C-A") & (m.contig_class == "core")],
        "bg_accessory": m[m.contig_class == "accessory"],
        "bg_core": m[m.contig_class == "core"],
    }
    print("\n" + "  ".join(f"{k}={len(v)}" for k, v in series.items()))

    bgs = {"A-A": "bg_accessory", "C-A": "bg_core"}
    x = np.arange(args.bin_kb / 2, args.max_kb + args.bin_kb / 2, args.bin_kb)

    dens = {k: density(v.dist_from_end.values, len(v), args) for k, v in series.items()}
    norm = {k: dens[k] / dens[bgs[k]] for k in ("A-A", "C-A")}
    for k, v in series.items():
        d = v.dist_from_end.values
        print(f"  {k:13s} median {np.median(d)/1000:7.1f} kb   <=50kb {100*(d<=5e4).mean():5.1f}%")

    # ---------------------------------------------------------------- products
    m.to_csv(os.path.join(OUT, "transfer_terminal_positions.tsv"), sep="\t", index=False)
    tab = pd.DataFrame({"dist_kb": x, "A-A": dens["A-A"], "C-A": dens["C-A"],
                        "bg_accessory": dens["bg_accessory"], "bg_core": dens["bg_core"],
                        "A-A_norm": norm["A-A"], "C-A_norm": norm["C-A"]})
    tab.to_csv(os.path.join(OUT, "transfer_terminal_density.tsv"), sep="\t", index=False)

    # ----------------------------------------------------------------- figures
    def base(ax, ylab):
        ax.set_xlim(0, args.max_kb)
        ax.set_xlabel("distance from chromosome end (kb)")
        ax.set_ylabel(ylab)
        ax.set_xticks([0, 50, 100, 150, 200])
        ax.axvline(100, color="grey", lw=0.6, ls=":")

    # v1 : the two curves alone
    fig, ax = plt.subplots(figsize=(120 / 25.4, 120 / 25.4))
    for k in ("A-A", "C-A"):
        ax.plot(x, dens[k], color={"A-A": C_AA, "C-A": C_CA}[k], lw=1.6, label=k)
    base(ax, "density (per kb, area = 1)")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("Exchange-event position from the chromosome end", fontsize=9, loc="left")
    fig.subplots_adjust(left=0.185, right=0.955, bottom=0.145, top=0.90)
    fig.savefig(os.path.join(OUT, "fig_terminal_density_v1_raw.pdf"))
    fig.savefig(os.path.join(OUT, "fig_terminal_density_v1_raw.png"), dpi=300)
    plt.close(fig)

    # v2 : the two curves + their own gene-density background
    fig, ax = plt.subplots(figsize=(120 / 25.4, 120 / 25.4))
    for k in ("A-A", "C-A"):
        ax.plot(x, dens[bgs[k]], color=BG, lw=0.9, ls="--",
                label=("accessory genes (background)" if k == "A-A"
                       else "core genes (background)"))
    for k in ("A-A", "C-A"):
        ax.plot(x, dens[k], color={"A-A": C_AA, "C-A": C_CA}[k], lw=1.7, label=k)
    base(ax, "density (per kb, area = 1)")
    ax.legend(loc="upper right", fontsize=6.4, handlelength=1.6, labelspacing=0.25)
    ax.set_title("Exchange events vs the gene-density background", fontsize=9, loc="left")
    fig.subplots_adjust(left=0.185, right=0.955, bottom=0.145, top=0.90)
    fig.savefig(os.path.join(OUT, "fig_terminal_density_v2_with_bg.pdf"))
    fig.savefig(os.path.join(OUT, "fig_terminal_density_v2_with_bg.png"), dpi=300)
    plt.close(fig)

    # v3 : observed / background
    fig, ax = plt.subplots(figsize=(120 / 25.4, 120 / 25.4))
    ax.axhline(1, color="grey", lw=0.8, ls="--")
    for k in ("A-A", "C-A"):
        ax.plot(x, norm[k], color={"A-A": C_AA, "C-A": C_CA}[k], lw=1.7, label=k)
    base(ax, "observed / background gene density (log scale)")
    # the ratio spans ~0.4 to ~50, so a linear axis would flatten A-A onto zero
    ax.set_yscale("log")
    ax.set_ylim(0.1, 100)
    ax.set_yticks([0.1, 0.3, 1, 3, 10, 30, 100])
    ax.set_yticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.text(202, 1, "  no enrichment", fontsize=6, color="grey", va="center", ha="left")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title("End-enrichment of exchange events over gene density",
                 fontsize=9, loc="left")
    fig.subplots_adjust(left=0.205, right=0.955, bottom=0.145, top=0.90)
    fig.savefig(os.path.join(OUT, "fig_terminal_density_v3_normalized.pdf"))
    fig.savefig(os.path.join(OUT, "fig_terminal_density_v3_normalized.png"), dpi=300)
    plt.close(fig)

    print(f"\nwrote 3 figure variants + 2 tables to {os.path.relpath(OUT, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
