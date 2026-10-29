#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Per-repeat-family comparison of core-chromosome terminal regions vs accessory
chromosomes at the optimal terminal window (default: 3% of each chromosome end,
where the overall repeat content of the ends matches the accessory chromosomes).

Outputs:
  suppl_table.repeat_terminal_family_<frac>.tsv
  figure_repeat_terminal.families_<frac>.pdf
"""

import argparse
import importlib.util
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

spec = importlib.util.spec_from_file_location(
    "srt", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "stats_repeat_terminal.py"))
srt = importlib.util.module_from_spec(spec)
spec.loader.exec_module(srt)


def parse_args():
    p = argparse.ArgumentParser(description="Per-family terminal vs accessory at f")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--terminal-fraction", type=float, default=0.03)
    p.add_argument("--terminal-length", type=int, default=None,
                   help="fixed terminal window per end (bp); overrides fraction")
    p.add_argument("--core-min-length", type=int, default=500000)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    f = args.terminal_fraction

    isolates = sorted({x[:-len(".assembly.fasta.out")]
                       for x in os.listdir(args.out_dir)
                       if x.endswith(".assembly.fasta.out")})
    isolates = [i for i in isolates if os.path.exists(
        os.path.join(args.out_dir, f"{i}.accessory_chrs.id"))]

    rows = []
    for iso in isolates:
        lengths = srt.load_lengths(os.path.join(
            args.length_dir, iso, "repeatmasker",
            f"{iso}.assembly.fasta.masked.length"))
        acc = {c for c in open(os.path.join(args.out_dir,
                                            f"{iso}.accessory_chrs.id"))
               .read().split() if c in lengths}
        core = {c for c, l in lengths.items()
                if c not in acc and l >= args.core_min_length}
        hits = srt.parse_out(os.path.join(args.out_dir,
                                          f"{iso}.assembly.fasta.out"),
                             set(lengths))
        for chrom, L in lengths.items():
            if chrom in acc:
                ctype = "accessory"
            elif chrom in core:
                ctype = "core"
            else:
                continue
            fam_bp = defaultdict(int)
            fam_term = defaultdict(int)
            wfix = (min(args.terminal_length, L / 2)
                    if args.terminal_length else f * L)
            for b, e, n, fam in hits.get(chrom, []):
                if n:
                    continue
                fam_g = srt.family_of(fam)
                fam_bp[fam_g] += e - b + 1
                if ctype == "core" and (e < wfix or b > L - wfix):
                    fam_term[fam_g] += e - b + 1
            fams = set(fam_bp) | set(fam_term)
            for g in fams:
                rows.append({
                    "Isolate": iso, "Chromosome": chrom, "Type": ctype,
                    "Family": g,
                    "Repeat_pct": fam_bp.get(g, 0) / L * 100,
                    "Terminal_pct": (fam_term.get(g, 0) / (2 * wfix) * 100
                                     if ctype == "core" else np.nan),
                })
    fam_df = pd.DataFrame(rows)

    # tests per family: terminal (f) vs accessory
    test_rows = []
    for g in sorted(fam_df["Family"].unique()):
        t = fam_df[(fam_df["Family"] == g) & fam_df["Type"].eq("core") &
                   fam_df["Terminal_pct"].notna()]["Terminal_pct"]
        a = fam_df[(fam_df["Family"] == g) & (fam_df["Type"] == "accessory")][
            "Repeat_pct"]
        if len(t) < 3 or len(a) < 3:
            continue
        u = stats.mannwhitneyu(t, a, alternative="two-sided")
        test_rows.append({
            "Family": g, "n_terminal": len(t), "n_accessory": len(a),
            "terminal_mean": t.mean(), "terminal_median": t.median(),
            "accessory_mean": a.mean(), "accessory_median": a.median(),
            "diff_term_acc": t.mean() - a.mean(), "MW_U": u.statistic,
            "MW_p": u.pvalue,
        })
    tests = pd.DataFrame(test_rows).sort_values("terminal_mean", ascending=False)
    suffix = f"{int(f * 1000):03d}" if not args.terminal_length else "100kb"
    win_label = (f"{args.terminal_length // 1000} kb" if args.terminal_length
                 else f"{f * 100:.0f}%")
    tests.to_csv(os.path.join(args.output_dir,
                              f"suppl_table.repeat_terminal_family_{suffix}.tsv"),
                 sep="\t", index=False)
    print(tests[["Family", "terminal_mean", "accessory_mean", "diff_term_acc",
                 "MW_p"]].round(3).to_string())

    # figure: paired dot plot (families with enough data, ordered)
    keep = tests[tests["MW_p"].notna()].head(12)
    pf = keep.iloc[::-1]
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42,
    })
    fig, ax = plt.subplots(figsize=(4.9, 4.6))
    y = np.arange(len(pf))
    ax.plot([pf["terminal_mean"], pf["accessory_mean"]], [y, y],
            color="#B0B0B0", linewidth=0.8, zorder=1)
    ax.scatter(pf["terminal_mean"], y, s=30, color="#DD8452", edgecolor="black",
               linewidths=0.4, zorder=2, label=f"Core terminal {win_label}")
    ax.scatter(pf["accessory_mean"], y, s=30, color="#C44E52", edgecolor="black",
               linewidths=0.4, zorder=2, label="Accessory")
    for i, (_, r) in enumerate(pf.iterrows()):
        if r["MW_p"] < 0.05:
            ax.annotate("*", (max(r["terminal_mean"], r["accessory_mean"]) + 0.3,
                              i), fontsize=8, va="center")
    ax.set_yticks(y)
    ax.set_yticklabels(pf["Family"], fontsize=6.5)
    ax.set_xlabel("Repeat content (%)", fontsize=8)
    ax.set_title(f"Per-family repeat content: core terminal {win_label} vs "
                 "accessory", fontsize=8.5, pad=3)
    ax.tick_params(axis="x", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.output_dir,
                             f"figure_repeat_terminal.families_{suffix}.pdf"),
                dpi=300)
    print(f"tables+figure -> {args.output_dir}")


if __name__ == "__main__":
    main()
