#!/usr/bin/env python3
"""Supplementary Figure S42 - shared k-mers along individual accessory chromosomes.

Six chromosome strips from the accessory_genome working directory
(/share/pasteur/gaosh/Fusarium/phase2/accessory_genome/solaniSC/tgs_assemblies/), each
a sliding-window heatmap of the proportion of k-mers shared between one accessory
chromosome and the chromosome-level assemblies (rows). Block examples on top (a-c),
mosaic examples below (d-f). The original images are embedded as rendered; the script
adds panel letters, the chromosome and source-species label, and the class.

Writes 00.manuscript/supplementary/Supplementary-Fig-S42-block-mosaic-kmer-patterns.{pdf,png}
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = "/share/pasteur/gaosh/Fusarium/phase2/accessory_genome/solaniSC/tgs_assemblies"
BASE = "all_bmu.assembly.combined.accessory_chrs.sliding.all_tgs.common_table.pct_"
OUTDIR = os.path.join(ROOT, "00.manuscript/supplementary")
OUT = os.path.join(OUTDIR, "Supplementary-Fig-S42-block-mosaic-kmer-patterns")

PANELS = [
    ("BMU01933_0021", "F. falciformis", "block"),
    ("BMU02474_0015", "F. solani", "block"),
    ("BMU06673_0013", "F. keratoplastica", "block"),
    ("BMU03267_0019", "F. falciformis", "mosaic"),
    ("BMU01933_0017", "F. falciformis", "mosaic"),
    ("BMU02028_0014", "F. solani", "mosaic"),
]
W_TARGET = 4200  # common display width; strips are never upscaled


def load(path, width):
    im = Image.open(path).convert("RGB")
    width = min(width, im.width)  # never upscale
    scale = width / im.width
    return im.resize((width, int(round(im.height * scale))), Image.LANCZOS)


def main():
    strips = [load(os.path.join(SRC, BASE + name + ".png"), W_TARGET)
              for name, _, _ in PANELS]
    heights = [im.height for im in strips]
    pad = 70
    fig = plt.figure(figsize=(11, sum(heights) / 300 + len(PANELS) * 0.35))
    gs = fig.add_gridspec(len(PANELS), 1, height_ratios=heights,
                          left=0.20, right=0.84, top=0.99, bottom=0.01,
                          hspace=0.22)
    for i, (im, (name, sp, cls)) in enumerate(zip(strips, PANELS)):
        ax = fig.add_subplot(gs[i])
        ax.imshow(im)
        ax.set_xlim(0, im.width)
        ax.set_ylim(im.height, 0)
        ax.axis("off")
        ax.text(-0.010, 0.5, f"{chr(97+i)}  {name}",
                transform=ax.transAxes, ha="right", va="center",
                fontsize=9, fontweight="bold")
        pos = ax.get_position()
        yc = pos.y0 + pos.height / 2
        fig.text(0.855, yc, sp, ha="left", va="center", fontsize=8.5,
                 fontstyle="italic")
        fig.text(0.950, yc, cls, ha="left", va="center", fontsize=8.5)
    fig.savefig(OUT + ".pdf")
    fig.savefig(OUT + ".png", dpi=300)
    plt.close(fig)
    print(f"wrote {os.path.relpath(OUT, ROOT)}.pdf/.png")


if __name__ == "__main__":
    main()
