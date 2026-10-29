#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Figures for the gene-transfer vs repeat-mediation analysis.

Primary control (per user decision): accessory-derived genes found on CORE
chromosomes (core-side copies of C-A clusters) vs other core-resident genes
that are never associated with accessory chromosomes. Accessory-resident
"non-transferred" genes are unreliable as controls because only 13 genomes
are sampled (a gene may be transferred in un-sampled genomes).

Inputs: 06.accessory-chrs/gene_transfer/suppl_table.gene_transfer_repeat.{cluster90,cluster100}.tsv
Outputs (06.accessory-chrs/gene_transfer/):
  figure_transfer_types.pdf          - transfer-type gene counts
  figure_transfer_mediation.pdf      - TE-flanking: accessory-derived-on-core
                                       vs core background (1 kb / 5 kb)
  figure_transfer_nearest_dist.pdf   - ECDF of nearest-repeat distance
"""

import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

C_TR = "#C44E52"
C_CB = "#4C72B0"
FAMS = ["DNA", "LTR", "LINE", "RC", "Unknown", "Simple_repeat",
        "Low_complexity"]


def parse_args():
    p = argparse.ArgumentParser(description="Gene transfer repeat figures")
    p.add_argument("--table90", required=True)
    p.add_argument("--table100", required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    dfs = {"cluster90": pd.read_csv(args.table90, sep="\t"),
           "cluster100": pd.read_csv(args.table100, sep="\t")}

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.linewidth": 0.6, "pdf.fonttype": 42,
    })

    # ---------------- 1. transfer types ----------------
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    x = np.arange(2)
    w = 0.35
    for i, label in enumerate(["cluster90", "cluster100"]):
        df = dfs[label]
        n_aa = (df["Transfer_type"] == "A-A").sum()
        n_ca = (df["Transfer_type"] == "C-A").sum()
        ax.bar(x[i] - w / 2, n_aa, w, color="#E64B35",
               label="A-A (accessory-accessory)" if i == 0 else None)
        ax.bar(x[i] + w / 2, n_ca, w, color="#4C72B0",
               label="C-A (core-accessory)" if i == 0 else None)
    ax.set_xticks(x)
    ax.set_xticklabels(["90% identity", "100% identity"], fontsize=7)
    ax.set_ylabel("Transfer-involved genes", fontsize=8)
    ax.set_title("Transfer-involved genes (accessory-related)", fontsize=8.5,
                 pad=3)
    ax.tick_params(labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{args.output_dir}/figure_transfer_types.pdf", dpi=300)
    plt.close(fig)

    # ---------------- 2. mediation: accessory-derived on core vs core bg ----
    for label, df in dfs.items():
        ca_core = df[(df["Transfer_type"] == "C-A") &
                     (df["Contig_type"] == "core")]
        core_back = df[(df["Transfer_type"].isin(["single", "C-C_ortho"])) &
                       (df["Contig_type"] == "core")]
        fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.4))
        for ax, col, met in [(axes[0], "Flank1kb_families", "1 kb"),
                             (axes[1], "Flank5kb_families", "5 kb")]:
            x = np.arange(len(FAMS))
            w = 0.38
            v_tr = [ca_core[col].fillna("").str.contains(f).mean() for f in FAMS]
            v_cb = [core_back[col].fillna("").str.contains(f).mean() for f in FAMS]
            ax.bar(x - w / 2, v_cb, w, color=C_CB, label="Core background")
            ax.bar(x + w / 2, v_tr, w, color=C_TR,
                   label="Accessory-derived (on core)")
            ax.set_xticks(x)
            ax.set_xticklabels(FAMS, rotation=30, ha="right", fontsize=6)
            ax.set_ylabel("Genes with TE within " + met, fontsize=7.5)
            ax.set_title(f"Within {met} ({label})", fontsize=8, pad=3)
            ax.tick_params(labelsize=6.5)
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="#E8E8E8", linewidth=0.5)
            ax.set_axisbelow(True)
            if met == "1 kb":
                ax.legend(fontsize=6, frameon=False)
        fig.suptitle(f"TE-flanking of accessory-derived genes on core "
                     f"chromosomes ({label})", fontsize=9, y=0.99)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(f"{args.output_dir}/figure_transfer_mediation.{label}.pdf",
                    dpi=300)
        plt.close(fig)

    # ---------------- 3. nearest-repeat distance ECDF ----------------
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    df = dfs["cluster90"]
    ca_core = df[(df["Transfer_type"] == "C-A") &
                 (df["Contig_type"] == "core")]["Nearest_repeat_dist"].dropna()
    core_back = df[(df["Transfer_type"].isin(["single", "C-C_ortho"])) &
                   (df["Contig_type"] == "core")]["Nearest_repeat_dist"].dropna()
    acc = df[(df["Transfer_type"] == "A-A")]["Nearest_repeat_dist"].dropna()
    for v, c, l in [(core_back, C_CB, "Core background"),
                    (ca_core, C_TR, "Accessory-derived (on core)"),
                    (acc, "#E64B35", "A-A (on accessory)")]:
        v = np.sort(v)
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, linewidth=1.1,
                label=l)
    ax.axvline(1000, color="black", linestyle=":", linewidth=0.8)
    ax.text(1050, 0.95, "1 kb", fontsize=6)
    ax.set_xlim(-200, 5000)
    ax.set_xlabel("Distance to nearest repeat (bp)", fontsize=8)
    ax.set_ylabel("Cumulative fraction of genes", fontsize=8)
    ax.set_title("Nearest-repeat distance (90% identity clusters)",
                 fontsize=8.5, pad=3)
    ax.tick_params(labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(color="#E8E8E8", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{args.output_dir}/figure_transfer_nearest_dist.pdf", dpi=300)
    plt.close(fig)
    print(f"figures -> {args.output_dir}")


if __name__ == "__main__":
    main()
