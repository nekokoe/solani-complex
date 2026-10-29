#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures for the accessory-chromosome repeat content analysis.

Inputs (06.accessory-chrs/repeat_terminal/):
  suppl_table.chromosome_repeat_content.tsv
  suppl_table.chromosome_repeat_by_family.tsv
  suppl_table.repeat_terminal_tests.tsv

Outputs:
  figure_repeat_terminal.overall.pdf   - repeat % of core whole chromosomes,
                                         core terminal-10% regions and
                                         accessory chromosomes (box + points)
  figure_repeat_terminal.paired.pdf    - per-isolate paired comparison of
                                         core-terminal-10% vs accessory means
  figure_repeat_terminal.families.pdf  - per-family repeat %: terminal-10% vs
                                         accessory (paired dot plot)
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

CORE_C = "#4C72B0"
TERM_C = "#DD8452"
ACC_C = "#C44E52"
GRID_C = "#E8E8E8"


def parse_args():
    p = argparse.ArgumentParser(description="Accessory repeat content figures")
    p.add_argument("--repeat-table", required=True)
    p.add_argument("--family-table", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def box_points(ax, data, colors, labels, ylabel, title, p_annot=None):
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", linewidth=0.8),
                    whiskerprops=dict(linewidth=0.6),
                    capprops=dict(linewidth=0.6))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
    for i, d in enumerate(data):
        x = np.random.default_rng(i).normal(i + 1, 0.05, size=len(d))
        ax.scatter(x, d, s=6, color="black", alpha=0.35, linewidths=0)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_title(title, fontsize=8.5, pad=3)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=GRID_C, linewidth=0.5)
    ax.set_axisbelow(True)
    if p_annot:
        ax.text(0.02, 0.97, p_annot, transform=ax.transAxes, va="top",
                fontsize=6.5)


def main():
    args = parse_args()
    df = pd.read_csv(args.repeat_table, sep="\t")
    fam = pd.read_csv(args.family_table, sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42,
    })

    # ---------------- 1. overall: core whole / core terminal10 / accessory ----
    core_whole = df.loc[df["Type"] == "core", "Repeat_pct"]
    core_term = df.loc[df["Type"] == "core", "Terminal100kb_pct"].dropna()
    acc = df.loc[df["Type"] == "accessory", "Repeat_pct"]

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    box_points(ax, [core_whole, core_term, acc],
               [CORE_C, TERM_C, ACC_C],
               ["Core chromosomes\n(whole)", "Core chromosomes\nterminal 100 kb",
                "Accessory\nchromosomes"],
               "Repeat content (%)",
               "Repeat content of core vs accessory chromosomes",
               p_annot="n = 166 / 166 / 78 chromosomes\n13 isolates")
    fig.tight_layout()
    fig.savefig(f"{args.output_dir}/figure_repeat_terminal.overall.pdf", dpi=300)
    plt.close(fig)

    # ---------------- 2. paired per isolate: terminal10 vs accessory ----
    term_iso = df[df["Type"] == "core"].groupby("Isolate")["Terminal100kb_pct"].mean()
    acc_iso = df[df["Type"] == "accessory"].groupby("Isolate")["Repeat_pct"].mean()
    common = sorted(set(term_iso.index) & set(acc_iso.index))

    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    x = np.arange(len(common))
    ax.plot([x, x], [term_iso[common].values, acc_iso[common].values],
            color="#B0B0B0", linewidth=0.7, zorder=1)
    ax.scatter(x, term_iso[common].values, s=26, color=TERM_C, edgecolor="black",
               linewidths=0.4, zorder=2, label="Core chromosome terminal 100 kb")
    ax.scatter(x, acc_iso[common].values, s=26, color=ACC_C, edgecolor="black",
               linewidths=0.4, zorder=2, label="Accessory chromosomes")
    ax.set_xticks(x)
    ax.set_xticklabels(common, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("Mean repeat content (%)", fontsize=8)
    ax.set_title("Per-isolate: terminal 100 kb vs accessory repeat content",
                 fontsize=8.5, pad=3)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color=GRID_C, linewidth=0.5)
    ax.set_axisbelow(True)
    w = stats.wilcoxon(term_iso[common], acc_iso[common])
    ax.legend(loc="upper left", fontsize=6.5, frameon=False)
    ax.text(0.98, 0.03, f"paired Wilcoxon p = {w.pvalue:.2e}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5)
    fig.tight_layout()
    fig.savefig(f"{args.output_dir}/figure_repeat_terminal.paired.pdf", dpi=300)
    plt.close(fig)

    print(f"figures -> {args.output_dir}")


if __name__ == "__main__":
    main()
