#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fixed-length terminal-window scan: test whether a fixed physical length
(e.g., 250 kb) of each core-chromosome end matches the accessory-chromosome
repeat content better than a percentage-based window (e.g., 3% per end).

For each fixed window length w in {50, 100, 150, 200, 250, 300, 400, 500,
750, 1000} kb (per chromosome end; capped at half the chromosome):
  - terminal repeat content per core chromosome
  - pooled Mann-Whitney vs accessory + per-isolate paired Wilcoxon

Outputs:
  suppl_table.terminal_fixed_scan.tsv
  figure_repeat_terminal.fixed_scan.pdf
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import stats_repeat_terminal as srt

FIXED_LENS_KB = [40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150,
                   160, 180, 200, 250, 300, 400, 500, 750, 1000]


def parse_args():
    p = argparse.ArgumentParser(description="Fixed-length terminal window scan")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--length-dir", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--core-min-length", type=int, default=500000)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    isolates = sorted({f[:-len(".assembly.fasta.out")]
                       for f in os.listdir(args.out_dir)
                       if f.endswith(".assembly.fasta.out")})
    isolates = [i for i in isolates if os.path.exists(
        os.path.join(args.out_dir, f"{i}.accessory_chrs.id"))]

    iso_data = {}
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
        iso_data[iso] = (lengths, acc, core, hits)

    # accessory reference (same for all windows)
    acc_all, iso_acc = [], {}
    for iso, (Ls, acc, core, hits) in iso_data.items():
        ap = [srt.merged_bp([(b, e) for b, e, n, f in hits.get(c, []) if not n])
              / Ls[c] * 100 for c in acc]
        acc_all += ap
        iso_acc[iso] = np.mean(ap) if ap else np.nan
    acc_all = np.array(acc_all)
    print(f"accessory: n={len(acc_all)} mean={acc_all.mean():.2f} "
          f"median={np.median(acc_all):.2f}")

    rows = []
    for w in FIXED_LENS_KB:
        win = w * 1000
        term_all, iso_term = [], {}
        for iso, (Ls, acc, core, hits) in iso_data.items():
            tp = []
            for c in core:
                L = Ls[c]
                ww = min(win, L / 2)
                nn = [(b, e) for b, e, n, f in hits.get(c, []) if not n]
                tiv = [(b, e) for b, e in nn if e < ww or b > L - ww]
                tp.append(srt.merged_bp(tiv) / (2 * ww) * 100)
            term_all += tp
            iso_term[iso] = np.mean(tp)
        term_all = np.array(term_all)
        common = [i for i in iso_term if i in iso_acc and not np.isnan(iso_acc[i])]
        u = stats.mannwhitneyu(term_all, acc_all, alternative="two-sided")
        wtest = stats.wilcoxon([iso_term[i] for i in common],
                               [iso_acc[i] for i in common])
        rows.append({
            "Window_kb_per_end": w,
            "n_terminal": len(term_all),
            "terminal_mean": term_all.mean(),
            "terminal_median": np.median(term_all),
            "accessory_mean": acc_all.mean(),
            "diff_term_acc": term_all.mean() - acc_all.mean(),
            "MW_U": u.statistic, "MW_p": u.pvalue,
            "paired_n": len(common),
            "paired_W": wtest.statistic, "paired_p": wtest.pvalue,
        })
    scan = pd.DataFrame(rows)
    scan.to_csv(os.path.join(args.output_dir,
                             "suppl_table.terminal_fixed_scan.tsv"),
                sep="\t", index=False)
    print(scan[["Window_kb_per_end", "terminal_mean", "accessory_mean",
                "diff_term_acc", "MW_p", "paired_p"]].round(3).to_string())

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.4))
    x = scan["Window_kb_per_end"]

    ax = axes[0]
    ax.plot(x, scan["terminal_mean"], "o-", color="#DD8452", markersize=4,
            linewidth=1.2, label="Core chromosome ends")
    ax.axhline(acc_all.mean(), color="#C44E52", linestyle="--", linewidth=0.9,
               label="Accessory chromosomes")
    ax.set_xlabel("Fixed terminal window (kb, each end)", fontsize=7.5)
    ax.set_ylabel("Mean repeat content (%)", fontsize=7.5)
    ax.set_title("Repeat content vs fixed window size", fontsize=8.5, pad=3)
    ax.legend(fontsize=6.5, frameon=False, loc="center right")
    ax.tick_params(labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    ax.plot(x, -np.log10(scan["paired_p"].clip(lower=1e-300)), "s-",
            color="#4C72B0", markersize=4, linewidth=1.2,
            label="Paired Wilcoxon (per isolate)")
    ax.axhline(-np.log10(0.05), color="black", linestyle=":", linewidth=0.9)
    ax.text(x.iloc[-1], -np.log10(0.05) + 0.12, "p = 0.05", fontsize=6,
            ha="right")
    ax.set_xlabel("Fixed terminal window (kb, each end)", fontsize=7.5)
    ax.set_ylabel("-log10(p)", fontsize=7.5)
    ax.set_title("Terminal vs accessory significance", fontsize=8.5, pad=3)
    ax.legend(fontsize=6.5, frameon=False)
    ax.tick_params(labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    fig.savefig(os.path.join(args.output_dir,
                             "figure_repeat_terminal.fixed_scan.pdf"), dpi=300)
    print(f"tables+figure -> {args.output_dir}")


if __name__ == "__main__":
    main()
