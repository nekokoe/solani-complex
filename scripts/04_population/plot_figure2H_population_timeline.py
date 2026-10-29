#!/usr/bin/env python3
"""Figure 2H - collection years of the ADMIXTURE populations.

The main text states that the clinical isolates "comprise a set of lineages that persist
side by side" and that the 2000-2005 concentration reflects one outbreak episode. This
figure shows that directly: every population of all three species is plotted against a
single year axis, so the spans can be compared across species at a glance and it is
visible that the populations are sampled both inside and outside the outbreak window.

Encoding (all colours derived from Fig. 2 so that this figure can sit next to it):
  * dot fill          clinical vs environmental  -> deepened versions of the pale red / pale
                      mint used for the source annotation track of Fig. 2D-F, drawn translucent
                      so that co-located isolates read without a hard edge
  * y-axis label      ADMIXTURE population       -> the Fig. 2D-F population colours
                      (#DB5243 Pop1 / #4BB2C7 Pop2 / #2D9D88 Pop3)
  * row background    species                    -> lightened versions of the Fig. 2I species
                      panel backgrounds (#CADAE8 / #FAE9C4 / #DACFE7), kept pale so they
                      stay behind the data; the outbreak band is opaque and cuts through them

Population labels follow Fig. 2D-F (argmax Q, components ranked by clinical fraction and
then by size, as in _fig2/build_groups.py), and the species are stacked in the same order
as those panels, F. solani / F. falciformis / F. keratoplastica.

  python3 04.population/plot_figure2H_population_timeline.py
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GROUPS = os.path.join(HERE, "_fig2", "data", "fig2_groups.tsv")
META = os.path.join(ROOT, "00.manuscript", "supplementary", "suppl_table.S1_genome_sources.xlsx")
OUTDIR = os.path.join(HERE, "fig2_panelH_population_timeline")
WIN = (2000, 2005)

# stacked so that the rows run in the same order as Fig. 2D (solani), 2E (falciformis), 2F (keratoplastica)
SPECIES = ["solani", "falciformis", "keratoplastica"]
NICE = {"solani": "F. solani", "falciformis": "F. falciformis", "keratoplastica": "F. keratoplastica"}
# Fig. 2I species panel backgrounds, mixed 45% towards white so the data stays dominant
SPCOL = {"solani": "#E2EBF2", "falciformis": "#FCF3DF", "keratoplastica": "#EBE5F2"}
# Fig. 2D-F source track colours, deepened from the pale track tints (#EEACAA / #9CD4B9)
SRCCOL = {"med": "#C8442E", "env": "#4E8C74", "unknown": "#9E9E9E"}
SRCLAB = {"med": "clinical", "env": "environmental", "unknown": "unknown"}
SRC_ALPHA = 0.80
S_UNIT = 26          # marker area (pt^2) per isolate; area stays proportional to the count
# Fig. 2D-F ADMIXTURE population colours, used here on the y-axis text
POPCOL = {"pop1": "#DB5243", "pop2": "#4BB2C8", "pop3": "#2D9D89"}
# outbreak band: a warm sand deliberately darker than the F. falciformis row tint (#FCF3DF),
# otherwise the opaque band and that species background merge where they meet
WINCOL = "#DCC49A"

FIGSIZE = (7.92, 3.105)   # 7.2 x 3.45 +10% width, -10% height
LAYOUT = dict(left=0.295, right=0.995, top=0.985, bottom=0.265)

plt.rcParams.update({
    "font.family": ["Nimbus Sans", "DejaVu Sans"],
    "font.size": 8.4,
    "axes.linewidth": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    # keep the italic species epithets in the same family as the body text;
    # the default mathtext font is DejaVu, which would not match Nimbus Sans
    "mathtext.fontset": "custom",
    "mathtext.it": "Nimbus Sans:italic",
    "mathtext.rm": "Nimbus Sans",
})


def main():
    g = pd.read_csv(GROUPS, sep="\t")
    # both sheets carry a collection year in Table S1 (the public genomes now have theirs
    # taken from the NCBI BioSample record), so merge both: reading only BMU used to drop
    # the public genomes and made the plotted counts fall short of the Fig. 2 population sizes
    _bmu = pd.read_excel(META, sheet_name="BMU")[["Isolate", "Collection year"]]
    _pub = pd.read_excel(META, sheet_name="public")[["Isolate", "Collection year"]]
    # the Fig. 2 group file keys public genomes by their short name; Table S1 carries the
    # full 'short_GCA_...' identifier, so strip the accession suffix before merging
    _pub["Isolate"] = _pub["Isolate"].str.split("_GCA_").str[0]
    m = pd.concat([_bmu, _pub]).rename(
        columns={"Isolate": "sample", "Collection year": "year"})
    df = g.merge(m[["sample", "year"]], on="sample", how="left").dropna(subset=["year"])
    df["year"] = df.year.astype(int)

    # one row per population, species blocks stacked from the top
    rows, y = [], {}
    for sp in SPECIES:
        pops = sorted(df[df.species == sp]["pop"].unique())
        for p in pops:
            rows.append((sp, p, df[df["pop"] == p]))
            y[p] = -len(rows) + 1

    fig, ax = plt.subplots(figsize=FIGSIZE)

    # species blocks carry the Fig. 2I panel background instead of alternating grey bands
    starts = [i for i, (sp, _, _) in enumerate(rows) if i == 0 or rows[i - 1][0] != sp]
    bounds = starts + [len(rows)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        # rows are laid out downwards, so the lower row has the smaller y: hand axhspan the
        # pair in ascending order or the band collapses to a single row (the last species
        # would get a zero-height band and vanish)
        ax.axhspan(y[rows[b - 1][1]] - 0.5, y[rows[a][1]] + 0.5, color=SPCOL[rows[a][0]], zorder=0)
    # the outbreak band is opaque and sits above the species backgrounds, so it cuts a
    # clean vertical stripe through the three species tints rather than tinting them again
    ax.axvspan(*WIN, color=WINCOL, zorder=0.6)

    # one dot per (population, year, source): area proportional to the number of isolates,
    # clinical above and environmental below the row centre
    SRC_OFF = {"med": 0.16, "env": -0.16}
    max_n = 0
    for sp, pop, s_ in rows:
        yy = y[pop]
        ax.plot([s_.year.min(), s_.year.max()], [yy, yy], color="#b9b9b9", lw=3.0,
                solid_capstyle="round", zorder=1)
        grp = s_.groupby(["year", "source"]).size().reset_index(name="n")
        max_n = max(max_n, int(grp.n.max()))
        for src, off in SRC_OFF.items():
            g2 = grp[grp.source == src]
            if not len(g2):
                continue
            ax.scatter(g2.year, np.full(len(g2), yy + off), s=g2.n * S_UNIT, facecolor=SRCCOL[src],
                       alpha=SRC_ALPHA, edgecolor="none", linewidth=0, zorder=3)
        other = grp[~grp.source.isin(SRC_OFF)]
        if len(other):
            ax.scatter(other.year, np.full(len(other), yy), s=other.n * S_UNIT,
                       facecolor=SRCCOL["unknown"], alpha=SRC_ALPHA, edgecolor="none",
                       linewidth=0, zorder=3)

    # one label per row, prefixed with the species so the blocks read as groups.
    # The species epithet is set in italics via mathtext (a plain string cannot be
    # half-italic), with the space inside the maths escaped.
    def _italic_species(sp):
        return r"$\it{" + NICE[sp].replace(" ", r"\ ") + r"}$"

    labels = [f"{_italic_species(sp)}  {p.split('-')[-1].replace('pop', 'Pop')} n={len(s)}"
              for sp, p, s in rows]
    ax.set_yticks([y[p] for _, p, _ in rows])
    ax.set_yticklabels(labels, fontsize=7.4)
    # the population is carried by the colour of its own axis label, so no Pop legend is needed
    for tick, (_, p, _) in zip(ax.get_yticklabels(), rows):
        tick.set_color(POPCOL.get(p.split("-")[-1], "#333333"))
    # headroom has to fit the largest marker, otherwise the topmost / bottom-most bubble is
    # sliced by the axes edge (F. solani Pop1 has a 13-isolate bubble in 2001). Solve for the
    # margin M in row units so that M * unit_pt >= offset * unit_pt + half the marker height:
    #   M >= (off + k * (n-1)) / (1 - 2k),  k = half_marker_pt / axes_height_pt
    axes_h_pt = (LAYOUT["top"] - LAYOUT["bottom"]) * FIGSIZE[1] * 72.0
    k = (np.sqrt(max_n * S_UNIT) / 2.0) / axes_h_pt
    off = max(abs(o) for o in SRC_OFF.values())
    margin = (off + k * (len(rows) - 1)) / (1 - 2 * k) * 1.30        # 30% slack, so the biggest bubble clears the frame visibly
    ax.set_ylim(-(len(rows) - 1) - margin, margin)
    # the earliest collection year in the plotted set is 1995 (two isolates of
    # F. solani pop3); the axis has to start before it or those two markers are clipped
    # in half by the left spine
    ax.set_xlim(1994, 2026)
    ax.set_xticks([1995, 2000, 2005, 2010, 2015, 2020])
    ax.set_xlabel("Year of isolation", fontsize=8.4)
    ax.tick_params(axis="x", labelsize=7.6)
    ax.tick_params(axis="y", length=2.5)

    # thin separators between species blocks
    for i in starts[1:]:
        ax.axhline(y[rows[i][1]] + 0.5, color="#ffffff", lw=1.0, zorder=2)

    fig.subplots_adjust(**LAYOUT)
    from matplotlib.lines import Line2D
    handles = [
        Line2D([], [], marker="o", ls="", markersize=5.2, markerfacecolor=SRCCOL["med"],
               markeredgecolor="none", alpha=SRC_ALPHA, label=SRCLAB["med"]),
        Line2D([], [], marker="o", ls="", markersize=5.2, markerfacecolor=SRCCOL["env"],
               markeredgecolor="none", alpha=SRC_ALPHA, label=SRCLAB["env"]),
        # the size key is drawn smaller than on the plot and given extra handle room:
        # at full size the 5-isolate circle overran its slot and covered the next label
        Line2D([], [], marker="o", ls="", markersize=np.sqrt(1 * 26) * 0.55, markerfacecolor="white",
               markeredgecolor="#666666", label="1 isolate"),
        Line2D([], [], marker="o", ls="", markersize=np.sqrt(5 * 26) * 0.55, markerfacecolor="white",
               markeredgecolor="#666666", label="5 isolates"),
        plt.Rectangle((0, 0), 1, 1, color=WINCOL, label="2000\u20132005 outbreak"),
    ]
    fig.legend(handles=handles, frameon=False, fontsize=7.2, loc="lower center",
               bbox_to_anchor=(0.63, 0.005), ncol=5, handlelength=1.5, handletextpad=0.8,
               columnspacing=1.4)

    os.makedirs(OUTDIR, exist_ok=True)
    for ext in ("pdf", "png"):
        path = os.path.join(OUTDIR, f"figure2H-population-timeline.{ext}")
        fig.savefig(path, dpi=400 if ext == "png" else None)
        print("wrote", path)

    for sp, p, s in rows:
        print(f"{NICE[sp]:<18} {p.split('-')[-1]}: n={len(s):>2}  {s.year.min()}-{s.year.max()}  "
              f"{int((s.year > WIN[1]).sum())} after 2005")


if __name__ == "__main__":
    main()
