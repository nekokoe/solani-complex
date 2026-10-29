#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Render a PNG companion for every delivered supplementary-figure PDF that lacks one.

Why
    The vision-capable main model can only inspect RASTER images -- a PDF handed to `read`
    is dropped by the provider. So every delivered figure needs a PNG twin, both for review
    and for downstream reuse. From now on every figure producer must emit PDF *and* PNG
    (the -dpi default here, 300, is the delivery standard).

Behaviour
    For each 00.manuscript/supplementary/Supplementary-Fig-S*.pdf without a sibling .png,
    render page 1 at --dpi. Existing PNGs are left untouched unless --force.

Usage
    python3 scripts/pdf2png_supplementary.py [--dpi 300] [--force] [--only S13 S16]
"""

import argparse
import glob
import os
import sys

import fitz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="00.manuscript/supplementary")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--only", nargs="*", default=None,
                    help="restrict to these S numbers, e.g. --only S13 S16")
    a = ap.parse_args()

    pdfs = sorted(glob.glob(os.path.join(a.dir, "Supplementary-Fig-S*.pdf")))
    if not pdfs:
        sys.exit(f"no PDFs found under {a.dir}")

    made = skipped = 0
    for pdf in pdfs:
        base = pdf[:-4]
        tag = os.path.basename(base).split("-")[2]          # e.g. S13 / S20g
        if a.only and tag not in a.only:
            continue
        png = base + ".png"
        if os.path.exists(png) and not a.force:
            skipped += 1
            continue
        doc = fitz.open(pdf)
        page = doc[0]
        w_mm = page.rect.width / 72 * 25.4
        h_mm = page.rect.height / 72 * 25.4
        pix = page.get_pixmap(dpi=a.dpi)
        pix.save(png)
        doc.close()
        made += 1
        print(f"  {tag:6s} -> {os.path.basename(png):58s} "
              f"{pix.width}x{pix.height}px  ({w_mm:.0f}x{h_mm:.0f} mm @{a.dpi}dpi)")
    print(f"\nrendered {made}, skipped {skipped} (already had a PNG)")


if __name__ == "__main__":
    main()
