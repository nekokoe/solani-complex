#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sliding-window position profile of accessory-chromosome gene-exchange events,
measured from the chromosome end (0 = the nearest end).

The window starts at the chromosome end and slides inward. Inside every window
the two exchange classes are counted AND normalised against the gene density of
the same chromosome class, so the gene-density gradient along the chromosome
cancels out and the curve has an interpretable baseline of 1.

    A-A  accessory-to-accessory exchange   counted on accessory chromosomes
    C-A  accessory-to-core exchange        counted on core chromosomes
                                           (the receiving, subtelomeric side)

    window at position p, width W:
        n_AA(p)  = A-A genes  with distance in [p, p+W)
        n_CA(p)  = C-A genes  with distance in [p, p+W)
        N_acc(p) = ALL accessory genes with distance in [p, p+W)
        N_cor(p) = ALL core genes      with distance in [p, p+W)

        rate_AA(p) = (n_AA / n_AA_total) / (N_acc / N_acc_total)
        rate_CA(p) = (n_CA / n_CA_total) / (N_cor / N_cor_total)

    i.e. how much MORE often an exchange event sits in this window than a random
    gene of the same chromosome class would. 1.0 = no end preference.

WHY THE NORMALISATION IS NOT OPTIONAL: genes cluster near chromosome ends, and
far more strongly on accessory chromosomes (10.5% of accessory genes lie within
50 kb of an end versus 1.2% of core genes). Without it the two classes look
similar; with it, C-A is up to ~48x enriched right at the end while A-A never
exceeds its own baseline by more than ~1.2x.

The window width is a free parameter, so the robustness across widths is written
out and drawn rather than asserted.

Inputs (in-repo)
    02.annotation/per_isolates/<ISO>/funannotate/predict_results/fusarium.gff3
    02.annotation/per_isolates/<ISO>/repeatmasker/<ISO>.assembly.fasta.masked.length
    06.accessory-chrs/tgs_assemblies/<ISO>.accessory_chrs.id
    06.accessory-chrs/gene_transfer_211/suppl_table.gene_transfer_211.tsv

Outputs (06.accessory-chrs/terminal_transfer_density/)
    slidingwindow_positions.tsv          per-gene: class, contig class, distance
    slidingwindow_profile.tsv            per window and width: counts + rates
    slidingwindow_summary.md             the numbers behind the curves
    fig_terminal_sliding_log.pdf/.png    main, log y   (both classes readable)
    fig_terminal_sliding_linear.pdf/.png same, linear y (A-A readable, C-A clipped)
    fig_terminal_sliding_windows.pdf/.png robustness across window widths
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
OUT = os.path.join(REPO, "06.accessory-chrs", "terminal_transfer_density")

ISOLATES = ["BMU01933", "BMU02028", "BMU02474", "BMU03267", "BMU06673",
            "BMU12683", "BMU12684", "BMU12685", "BMU12686",
            "Fu3.1", "Fu6.1", "SB1", "SX13"]
C_AA, C_CA = "#C44E52", "#3B7DD8"
MM = 25.4


def parse_args():
    p = argparse.ArgumentParser(description="sliding-window terminal profile")
    p.add_argument("--window-kb", type=float, default=20.0)
    p.add_argument("--step-kb", type=float, default=5.0)
    p.add_argument("--max-kb", type=float, default=500.0)
    p.add_argument("--widths-kb", default="10,20,40",
                   help="window widths for the robustness panel")
    p.add_argument("--core-min-length", type=int, default=500000)
    return p.parse_args()


def build_positions(args):
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
            L = clen[f[0]]
            rows.append((iso, tid, f[0], cls[f[0]],
                         min(int(f[3]), L - int(f[4]) + 1)))
    g = pd.DataFrame(rows, columns=["Isolate", "Transcript", "Contig",
                                    "contig_class", "dist"])
    t = pd.read_csv(os.path.join(REPO, "06.accessory-chrs", "gene_transfer_211",
                                 "suppl_table.gene_transfer_211.tsv"), sep="\t",
                    usecols=["Isolate", "Transcript", "Transfer_type"])
    return g.merge(t, on=["Isolate", "Transcript"], how="left")


def windows_for(d, window_kb, step_kb, max_kb):
    """sliding windows starting at the chromosome end"""
    starts = np.arange(0, max_kb - window_kb + step_kb, step_kb)
    return starts


def profile(pos, args, window_kb, step_kb):
    """rate for each class in each window (see the module docstring)"""
    acc = pos[pos.contig_class == "accessory"].dist.values
    cor = pos[pos.contig_class == "core"].dist.values
    aa = pos[(pos.Transfer_type == "A-A") & (pos.contig_class == "accessory")].dist.values
    ca = pos[(pos.Transfer_type == "C-A") & (pos.contig_class == "core")].dist.values
    tot = {"AA": len(aa), "CA": len(ca), "acc": len(acc), "cor": len(cor)}
    rec = []
    W = window_kb * 1000
    for p in windows_for(pos, window_kb, step_kb, args.max_kb):
        lo, hi = p * 1000, (p + window_kb) * 1000
        c = {k: int(((v >= lo) & (v < hi)).sum()) for k, v in
             (("AA", aa), ("CA", ca), ("acc", acc), ("cor", cor))}
        e_aa = c["AA"] / tot["AA"] if tot["AA"] else np.nan
        e_acc = c["acc"] / tot["acc"] if tot["acc"] else np.nan
        e_ca = c["CA"] / tot["CA"] if tot["CA"] else np.nan
        e_cor = c["cor"] / tot["cor"] if tot["cor"] else np.nan
        rec.append({"window_kb": window_kb, "start_kb": p,
                    "n_AA": c["AA"], "n_CA": c["CA"],
                    "n_genes_accessory": c["acc"], "n_genes_core": c["cor"],
                    "rate_AA": (e_aa / e_acc) if e_acc else np.nan,
                    "rate_CA": (e_ca / e_cor) if e_cor else np.nan,
                    "frac_AA_of_genes": (c["AA"] / c["acc"]) if c["acc"] else np.nan,
                    "frac_CA_of_genes": (c["CA"] / c["cor"]) if c["cor"] else np.nan})
    return pd.DataFrame(rec)


def main():
    args = parse_args()
    os.makedirs(OUT, exist_ok=True)
    pos = build_positions(args)
    print(f"mRNAs: {len(pos)}")
    pos.to_csv(os.path.join(OUT, "slidingwindow_positions.tsv"), sep="\t", index=False)

    widths = [float(x) for x in args.widths_kb.split(",")]
    all_rec = [profile(pos, args, w, args.step_kb) for w in widths]
    prof = pd.concat(all_rec, ignore_index=True)
    prof.to_csv(os.path.join(OUT, "slidingwindow_profile.tsv"), sep="\t", index=False)

    main_p = profile(pos, args, args.window_kb, args.step_kb)
    pk = main_p.rate_CA.max()
    print(f"\nwindow {args.window_kb:g} kb, step {args.step_kb:g} kb")
    print(f"  peak C-A rate {pk:.1f}x at {main_p.loc[main_p.rate_CA.idxmax(),'start_kb']:g} kb")
    print(f"  peak A-A rate {main_p.rate_AA.max():.2f}x")
    print("  C-A rate by distance:",
          "  ".join(f"{r.start_kb:g}kb={r.rate_CA:.1f}x"
                   for _, r in main_p[main_p.start_kb.isin([0, 5, 10, 20, 50, 100])].iterrows()))

    # ------------------------------------------------------------- summary md
    with open(os.path.join(OUT, "slidingwindow_summary.md"), "w") as f:
        f.write("# Sliding-window profile of exchange events from the chromosome end\n\n")
        f.write("Window starts at the chromosome end (0 = nearest end) and slides inward.\n")
        f.write("Inside each window the two classes are normalised against the gene density of\n")
        f.write("the same chromosome class, so 1.0 means \"no end preference\" and the\n")
        f.write("gene-density gradient along the chromosome cancels out.\n\n")
        f.write(f"Main window {args.window_kb:g} kb, step {args.step_kb:g} kb.\n\n")
        f.write("| window (kb) | A-A | C-A | genes acc | genes core | rate A-A | rate C-A |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for _, r in main_p.iterrows():
            f.write(f"| {r.start_kb:.0f}-{r.start_kb+args.window_kb:.0f} | {r.n_AA:.0f} | "
                    f"{r.n_CA:.0f} | {r.n_genes_accessory:.0f} | {r.n_genes_core:.0f} | "
                    f"{r.rate_AA:.2f} | {r.rate_CA:.2f} |\n")
        f.write("\n## Robustness across window widths\n\n")
        f.write("| width (kb) | peak rate C-A | peak at (kb) | max rate A-A |\n|---|---|---|---|\n")
        for w in widths:
            s = prof[prof.window_kb == w]
            if not len(s):
                continue
            f.write(f"| {w:g} | {s.rate_CA.max():.1f} | "
                    f"{s.loc[s.rate_CA.idxmax(),'start_kb']:.0f} | {s.rate_AA.max():.2f} |\n")
        f.write("\nC-A = accessory-to-core exchange counted on core chromosomes; "
                "A-A = accessory-to-accessory exchange counted on accessory chromosomes.\n")

    # ---------------------------------------------------------------- figures
    def draw(ax, df, logy):
        ax.axhline(1.0, color="grey", lw=0.9, ls="--")
        ax.plot(df.start_kb, df.rate_AA, color=C_AA, lw=1.7, label="A-A")
        ax.plot(df.start_kb, df.rate_CA, color=C_CA, lw=1.7, label="C-A")
        ax.set_xlim(0, args.max_kb)
        ax.set_xticks(list(range(0, int(args.max_kb) + 1, 100)))
        ax.set_xlabel("distance from chromosome end (kb)")
        ax.set_ylabel("observed / expected event rate"
                      + ("  (log)" if logy else ""))
        for s in ("right", "top"):
            ax.spines[s].set_visible(False)
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(0.1, 100)
            ax.set_yticks([0.1, 0.3, 1, 3, 10, 30, 100])
            ax.set_yticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
        else:
            ax.set_ylim(0, min(12, max(3, df.rate_CA.max() * 1.05)))

    for tag, logy in (("log", True), ("linear", False)):
        fig, ax = plt.subplots(figsize=(120 / MM, 120 / MM))
        draw(ax, main_p, logy)
        ax.legend(loc="upper right", fontsize=8)
        ax.set_title("Exchange events near the chromosome end\n"
                     f"({args.window_kb:g} kb sliding window)", fontsize=8.5, loc="left")
        fig.subplots_adjust(left=0.215, right=0.955, bottom=0.15, top=0.845)
        fig.savefig(os.path.join(OUT, f"fig_terminal_sliding_{tag}.pdf"))
        fig.savefig(os.path.join(OUT, f"fig_terminal_sliding_{tag}.png"), dpi=300)
        plt.close(fig)

    # robustness: one line style per window width, C-A solid / A-A dashed
    fig, ax = plt.subplots(figsize=(120 / MM, 120 / MM))
    ax.axhline(1.0, color="grey", lw=0.9, ls="--")
    lw_cycle = [("solid", 1.8), ("dashed", 1.2), ("dashdot", 1.0)]
    for (w, (lsty, lw)) in zip(widths, lw_cycle):
        s = prof[prof.window_kb == w]
        ax.plot(s.start_kb, s.rate_CA, color=C_CA, ls=lsty, lw=lw, label=f"C-A {w:g} kb")
    for (w, (lsty, lw)) in zip(widths, lw_cycle):
        s = prof[prof.window_kb == w]
        ax.plot(s.start_kb, s.rate_AA, color=C_AA, ls=lsty, lw=lw, label=f"A-A {w:g} kb")
    ax.set_xlim(0, args.max_kb)
    ax.set_xticks(list(range(0, int(args.max_kb) + 1, 100)))
    ax.set_xlabel("distance from chromosome end (kb)")
    ax.set_ylabel("observed / expected event rate")
    ax.set_yscale("log")
    ax.set_ylim(0.1, 100)
    ax.set_yticks([0.1, 0.3, 1, 3, 10, 30, 100])
    ax.set_yticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.legend(loc="upper right", fontsize=5.6, ncol=2, handlelength=1.7,
              columnspacing=0.9, labelspacing=0.22)
    ax.set_title("Robustness across window widths", fontsize=8.5, loc="left")
    fig.subplots_adjust(left=0.215, right=0.955, bottom=0.15, top=0.86)
    fig.savefig(os.path.join(OUT, "fig_terminal_sliding_windows.pdf"))
    fig.savefig(os.path.join(OUT, "fig_terminal_sliding_windows.png"), dpi=300)
    plt.close(fig)

    print(f"\nwrote 3 figures + 3 tables to {os.path.relpath(OUT, REPO)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
