#!/usr/bin/env python3
"""Assemble Supplementary Figure S40 - within-species PCoA of core CNV and non-core PAV.

These panels were in the main figure at an earlier stage and were removed when Figure 3 was
trimmed; they were never carried over to the supplement, which left the within-species paragraph
citing figure panels that no longer existed. The per-species ordinations themselves were always
present on disk (05.pan-genome/accessory_genes/*.pcoa/), so nothing is recomputed here - the six
finished vector panels are placed on one page.

One panel is a re-render: PC1 of the keratoplastica core ordination spans about two hundredths,
so matplotlib's default locator prints six overlapping four-decimal tick labels there. That panel
is regenerated once by `scripts/plot_pcoa.py --max_ticks 5 --width 5.39 --height 4.88` into
05.pan-genome/accessory_genes/S40-panels/, which this script prefers when present. The other five
are used as they are.

Rows are species, columns are the two matrices the Results text contrasts: core orthogroups
ranked by gene copy number, and non-core orthogroups scored by presence/absence.

Writes 00.manuscript/supplementary/Supplementary-Fig-S40-within-species-PCoA.pdf
       00.manuscript/supplementary/Supplementary-Fig-S40-within-species-PCoA.png
"""

import os

import fitz

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "05.pan-genome/accessory_genes")
SRC_OVERRIDE = os.path.join(SRC, "S40-panels")
OUTDIR = os.path.join(ROOT, "00.manuscript/supplementary")
BASE = "Supplementary-Fig-S40-within-species-PCoA"

FONT_R = "/usr/share/fonts/urw-base35/NimbusSans-Regular.otf"
FONT_B = "/usr/share/fonts/urw-base35/NimbusSans-Bold.otf"
FONT_I = "/usr/share/fonts/urw-base35/NimbusSans-Italic.otf"

PW, PH = 388.0, 351.0
LEFT, TOP, GAP_COL, GAP_ROW = 108.0, 78.0, 44.0, 34.0
W = LEFT + PW * 2 + GAP_COL + 8
H = TOP + PH * 3 + GAP_ROW * 2 + 10

ROWS = [("solani", "F. solani", "86"), ("falciformis", "F. falciformis", "66"),
        ("keratoplastica", "F. keratoplastica", "23")]
COLS = [("N0.HOG.hog_count.{sp}_core", "Core orthogroups (gene copy number)"),
        ("N0.HOG.hog_pav.{sp}_uncore", "Non-core orthogroups (presence/absence)")]


def put(page, x, y, text, size=8.5, font="r", color=(0, 0, 0)):
    f = {"r": FONT_R, "b": FONT_B, "i": FONT_I}[font]
    page.insert_text((x, y), text, fontsize=size, fontfile=f, fontname=os.path.basename(f),
                     color=color)


def main():
    doc = fitz.open()
    page = doc.new_page(width=W, height=H)

    # column headers
    for ci, (_, title) in enumerate(COLS):
        x = LEFT + ci * (PW + GAP_COL)
        put(page, x + 6, TOP - 26, title, size=9.5, font="b")

    letters = "abcdef"
    for ri, (sp, label, n) in enumerate(ROWS):
        y = TOP + ri * (PH + GAP_ROW)
        # row label, species epithet in italics
        put(page, 10, y + PH / 2 - 2, label, size=10, font="i")
        put(page, 10, y + PH / 2 + 10, f"(n = {n})", size=8.5, font="r", color=(0.3, 0.3, 0.3))
        for ci, (pattern, _) in enumerate(COLS):
            x = LEFT + ci * (PW + GAP_COL)
            name = pattern.format(sp=sp)
            override = os.path.join(SRC_OVERRIDE, name + ".pdf")
            path = override if os.path.exists(override) else os.path.join(
                SRC, name + ".pcoa_plot.pdf")
            src = fitz.open(path)
            page.show_pdf_page(fitz.Rect(x, y, x + PW, y + PH), src, 0)
            src.close()
            put(page, x + 2, y + 12, letters[ri * 2 + ci], size=12, font="b")

    pdf = os.path.join(OUTDIR, BASE + ".pdf")
    doc.save(pdf, deflate=True, garbage=4)
    doc.close()

    png = os.path.join(OUTDIR, BASE + ".png")
    d = fitz.open(pdf)
    d[0].get_pixmap(dpi=600).save(png)
    d.close()
    print(f"wrote {os.path.relpath(pdf, ROOT)}  ({os.path.getsize(pdf)/1024:.1f} KB)")
    print(f"wrote {os.path.relpath(png, ROOT)}  ({os.path.getsize(png)/1024:.1f} KB)")
    print(f"page {W:.0f} x {H:.0f} pt = {W/72*25.4:.0f} x {H/72*25.4:.0f} mm")


if __name__ == "__main__":
    main()
