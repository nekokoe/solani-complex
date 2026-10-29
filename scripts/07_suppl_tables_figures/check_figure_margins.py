#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Acceptance check for rendered figures: is there enough ink-free margin on all
four sides?

Project convention (learned the hard way): measure the distance from the CANVAS
EDGE to the MIN/MAX DARK PIXEL -- never test only the border row/column, which
passes even when content is flush against the edge. Ink threshold is
`im < 160` on an 8-bit grayscale render. Default floor 1.2 mm per side.

Usage:
    python3 scripts/check_figure_margins.py FIG.pdf [FIG2.pdf ...]
    python3 scripts/check_figure_margins.py --min-mm 1.2 06.*/figs/*.pdf
"""

import argparse
import os
import sys

import fitz
import numpy as np
from PIL import Image


def render(pdf, dpi):
    doc = fitz.open(pdf)
    page = doc[0]
    pix = page.get_pixmap(dpi=dpi)
    im = Image.frombytes("RGB", (pix.width, pix.height), pix.samples).convert("L")
    arr = np.asarray(im)
    doc.close()
    return arr, pix.width, pix.height


def margins_mm(arr, dpi):
    """Return (left, top, right, bottom) ink-free margins in mm."""
    ink = arr < 160
    if not ink.any():
        return None
    ys, xs = np.where(ink)
    ppm = arr.shape[1] / (arr.shape[1] / dpi * 25.4)   # px per mm
    ppm = dpi / 25.4
    return (xs.min() / ppm,
            ys.min() / ppm,
            (arr.shape[1] - 1 - xs.max()) / ppm,
            (arr.shape[0] - 1 - ys.max()) / ppm)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("pdfs", nargs="+")
    p.add_argument("--min-mm", type=float, default=1.2)
    p.add_argument("--dpi", type=int, default=300)
    a = p.parse_args()

    pad = max(len(os.path.basename(x)) for x in a.pdfs) + 2
    print(f"{'figure':<{pad}}{'size (mm)':<18}{'L':>7}{'T':>7}{'R':>7}{'B':>7}   verdict")
    print("-" * (pad + 18 + 34))
    fails = 0
    for pdf in a.pdfs:
        if not os.path.exists(pdf):
            print(f"{os.path.basename(pdf):<{pad}}MISSING")
            fails += 1
            continue
        arr, w, h = render(pdf, a.dpi)
        m = margins_mm(arr, a.dpi)
        size = f"{w / a.dpi * 25.4:.0f} x {h / a.dpi * 25.4:.0f}"
        if m is None:
            print(f"{os.path.basename(pdf):<{pad}}{size:<18}{'':>28}   BLANK")
            fails += 1
            continue
        # locate the offending band so the failure is diagnosable without guessing
        where = ""
        if min(m) < a.min_mm:
            ppm = a.dpi / 25.4
            ink = arr < 160
            hh = arr.shape[0]
            side = ["left", "top", "right", "bottom"][int(np.argmin(m))]
            if side in ("right", "left"):
                col = ink[:, -1] if side == "right" else ink[:, 0]
            else:
                col = ink[-1, :] if side == "bottom" else ink[0, :]
            band = []
            for i in range(20):
                r0, r1 = int(i * hh / 20), int((i + 1) * hh / 20)
                if ink[r0:r1].any():
                    c = np.where(ink[r0:r1].any(axis=0))[0]
                    gap = ((arr.shape[1] - 1 - c.max()) if side == "right"
                           else c.min()) / ppm
                    if gap < a.min_mm:
                        band.append(f"{r0/ppm:.0f}-{r1/ppm:.0f}mm")
            where = f"   touching: {side}, y-band {'/'.join(band) or '?'}"
        ok = all(v >= a.min_mm for v in m)
        fails += 0 if ok else 1
        print(f"{os.path.basename(pdf):<{pad}}{size:<18}"
              f"{m[0]:>7.2f}{m[1]:>7.2f}{m[2]:>7.2f}{m[3]:>7.2f}   "
              f"{'PASS' if ok else 'FAIL (< %.1f mm)' % a.min_mm}{where}")
    print("-" * (pad + 18 + 34))
    print(f"{len(a.pdfs) - fails}/{len(a.pdfs)} pass the >= {a.min_mm} mm margin floor")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
