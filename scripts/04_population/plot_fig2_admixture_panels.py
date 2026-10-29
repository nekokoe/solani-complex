#!/usr/bin/env python3
"""plot_fig2_admixture_panels.py — rebuild the Figure 2 D/E/F ADMIXTURE bar strips
with pop colours keyed to the pop rank used everywhere else in the paper.

Why: scripts/plot_admixture.R colours the stacked bars positionally on the ADMIXTURE
.Q column order, while 04.population/_fig2/build_groups.py numbers populations
pop1..popK by clinical fraction descending. The two orders differ per species, so the
shipped panels label the wrong colour as Pop1/Pop2/Pop3. This script keeps the row
geometry of the shipped panels and only re-keys the colours to the rank order.

Row order = the order of the sample names in the shipped figure, top to bottom. That
order comes from tree_order_<species>.txt where the file exists and is read off the
figure itself otherwise (panel F has no order file, and the per-species .contree /
full-211 .contree / ward.D2 clustering all fail to reproduce it).

Usage:  python3 04.population/plot_fig2_admixture_panels.py
Output: 04.population/fig2_admixture_panels_fixed/fig2_panel<D|E|F>_<species>_admixture_bars.{pdf,png}
"""
import os
import re
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import fitz
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MS = os.path.join(ROOT, "00.manuscript")
FIG = os.path.join(MS, "figure2-v2.pdf")
OUT = os.path.join(HERE, "fig2_admixture_panels_fixed")

# pop rank -> fill colour, sampled from the shipped figure2-v2.pdf panels
POPCOL = {1: "#DB5243", 2: "#4BB2C7", 3: "#2D9D88"}

# panel -> (species, .Q dir, K, bar_left, bar_right, row_pitch) in PDF points
PANELS = {
    "D": ("solani", "solani-short-solani", 3, 140.237, 177.087, 4.33),
    "E": ("falciformis", "solani-short-falciformis", 3, 298.104, 335.444, 4.33),
    "F": ("keratoplastica", "solani-short-keratoplastica", 2, 491.033, 528.156, 4.53),
}
# panel -> bounding box of the panel in figure2-v2.pdf, used to read the sample names
BBOX = {
    "D": (15, 205, 195, 645),
    "E": (190, 205, 370, 560),
    "F": (365, 205, 535, 380),
}
NAME_SIZE_MAX = 5.0     # sample names are ArialMT 3.6 pt; axis ticks/legend are larger


def clinical_source():
    src = {}
    path = os.path.join(ROOT, "05.stats_pan_by_species_shorter/bmu_ref_groups.med_env.txt")
    with open(path) as fh:
        for line in fh:
            a = line.split()
            src[a[0]] = a[1]
    return src


def load_q(species, dname, k):
    pre = os.path.join(HERE, dname, "%s.raw.snps.combined.genotype.maf5pct.recode" % dname)
    fam = [l.split()[0] for l in open(pre + ".fam")]
    q = [[float(x) for x in l.split()] for l in open("%s.%d.Q" % (pre, k))]
    assert len(fam) == len(q), (species, len(fam), len(q))
    return fam, q


def pop_rank(fam, q, k, src):
    """pop1..popK = clinical fraction descending, then n descending (build_groups.py)."""
    hard = [max(range(k), key=lambda i: row[i]) for row in q]
    stat = []
    for kk in range(k):
        mem = [s for s, h in zip(fam, hard) if h == kk]
        n = len(mem)
        stat.append((kk, sum(1 for s in mem if src.get(s) == "med") / n if n else 0.0, n))
    order = [kk for kk, _, _ in sorted(stat, key=lambda t: (-t[1], -t[2]))]
    return hard, {kk: i + 1 for i, kk in enumerate(order)}


def names_from_figure(panel):
    """Sample names down the panel, top to bottom, read out of the shipped figure."""
    x0, y0, x1, y1 = BBOX[panel]
    doc = fitz.open(FIG)
    page = doc[0]
    hits = []
    for block in page.get_text("dict")["blocks"]:
        for line in block.get("lines", []):
            for span in line["spans"]:
                text = span["text"].strip()
                x, y = span["bbox"][0], span["bbox"][1]
                if not (x0 <= x <= x1 and y0 <= y <= y1):
                    continue
                if span["size"] > NAME_SIZE_MAX or re.fullmatch(r"[0-9.]+", text):
                    continue
                hits.append((round(y, 2), text))
    doc.close()
    hits.sort()
    return [t for _, t in hits]


def row_order(panel, species):
    """tree_order_<species>.txt when present, else the order read off the figure."""
    path = os.path.join(HERE, "solani-short-%s" % species, "tree_order_%s.txt" % species)
    if os.path.exists(path):
        return [l.strip() for l in open(path) if l.strip()], path
    return names_from_figure(panel), FIG


def draw(panel, species, dname, k, left, right, pitch, order, fam, q, rank):
    width = round(right - left, 3)
    height = round(len(order) * pitch, 3)
    index = {s: i for i, s in enumerate(fam)}
    missing = [n for n in order if n.replace(" ", "_") not in index]
    assert not missing, (panel, missing[:5])

    fig = plt.figure(figsize=(width / 72.0, height / 72.0))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")
    membership = {}
    for i, name in enumerate(order):
        row = q[index[name.replace(" ", "_")]]
        total = sum(row)
        y = height - (i + 1) * pitch
        x = 0.0
        for kk in range(k):
            w = row[kk] / total * width
            if w > 0:
                ax.add_patch(Rectangle((x, y), w, pitch, facecolor=POPCOL[rank[kk]],
                                       edgecolor="none", linewidth=0, antialiased=False))
            x += w
        r = rank[max(range(k), key=lambda j: row[j])]
        membership[r] = membership.get(r, 0) + 1

    base = os.path.join(OUT, "fig2_panel%s_%s_admixture_bars" % (panel, species))
    fig.savefig(base + ".pdf", format="pdf")
    fig.savefig(base + ".png", dpi=600)
    plt.close(fig)
    return base, width, height, dict(sorted(membership.items()))


def main():
    os.makedirs(OUT, exist_ok=True)
    src = clinical_source()
    summary = {}
    for panel, (species, dname, k, left, right, pitch) in PANELS.items():
        fam, q = load_q(species, dname, k)
        hard, rank = pop_rank(fam, q, k, src)
        order, order_src = row_order(panel, species)
        base, width, height, membership = draw(panel, species, dname, k, left, right, pitch,
                                               order, fam, q, rank)
        summary[panel] = dict(species=species, rows=len(order), k=k,
                              width_pt=width, height_pt=height, pitch_pt=pitch,
                              membership=membership, order_source=order_src,
                              files=[base + ".pdf", base + ".png"])
        print("%s %-15s rows=%2d  %.2f x %.2f pt  pitch=%.2f  pop membership=%s  order=%s"
              % (panel, species, len(order), width, height, pitch, membership,
                 "tree_order file" if order_src.endswith(".txt") else "read off figure"))
    json.dump(summary, open(os.path.join(OUT, "panels_summary.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
