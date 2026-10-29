#!/usr/bin/env python3
"""layout_fig2_IJ.py — render the 2I / 2J panels with a guaranteed left whitespace band.

Requirement: for every 2I and 2J panel (full + half matrix, with and without legend) the
left half of the figure is blank, so nothing can be clipped on that side. X-axis labels stay
at 45 degrees and carry the full species names, which makes the label column wide -- so the
canvas is sized FROM the measured content, per panel, not fixed.

Method
  pass A   render each panel on a deliberately over-wide canvas with the final height (the
           panel is height-constrained under coord_equal, so its size is unchanged) and
           generous left/right room, then measure the true ink width
  solve    for each panel: canvas width W = 2 * (ink_width + right_margin) and left margin
           MX_L = W/2 - offset, where `offset` is the measured distance between the ggplot
           plot margin and the leftmost ink (ggplot does not account for rotated markdown
           axis text, so this offset must be measured rather than assumed)
  iterate  until the measured left-blank fraction is 50 +/- 1 % and every gap exceeds 1.5 mm
"""
import os, subprocess, sys, json
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))

MX_R = 11.0          # right margin: the 45-degree top labels overflow into it
MX_T = 11.0
MX_B = 4.0
H_NL = 62.0          # height, no legend
H_LG = 76.0          # height, with legend
XANG = 45
TARGET = 0.50        # left-blank fraction
MIN_GAP = 1.5        # mm, right/top/bottom

TOKENS = ['I', 'J', 'Ihalf', 'Jhalf']
STEMS = {'I': ('fig2I_nolegend', 'fig2I'),
         'J': ('fig2J_nolegend', 'fig2J'),
         'Ihalf': ('fig2I_half_nolegend', 'fig2I_half'),
         'Jhalf': ('fig2J_half_nolegend', 'fig2J_half')}


def render(token, outdir, w, mx_l):
    env = os.environ.copy()
    env.update({'FIG_W_IJ': str(w), 'FIG_H_NL': str(H_NL), 'FIG_H_LG': str(H_LG),
                'FIG_MX_L': str(mx_l), 'FIG_MX_R': str(MX_R),
                'FIG_MX_T': str(MX_T), 'FIG_MX_B': str(MX_B),
                'FIG_XANG': str(XANG), 'FIG_ONLY': token})
    r = subprocess.run(['Rscript', 'plot_fig2_HIJ.R', 'data', outdir],
                       cwd=HERE, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:]); print(r.stderr[-3000:]); sys.exit(1)


def extents(png, w_mm, h_mm):
    im = np.array(Image.open(png).convert('L'))
    dark = im < 160
    H, W = im.shape
    cols = np.nonzero(dark.any(axis=0))[0]
    rows = np.nonzero(dark.any(axis=1))[0]
    px, py = W / w_mm, H / h_mm
    return dict(l=cols.min() / px, r=(W - 1 - cols.max()) / px,
                t=rows.min() / py, b=(H - 1 - rows.max()) / py,
                cw=(cols.max() - cols.min()) / px)


def measure(tok, outdir, w):
    """worst case across the two legend variants"""
    m = {}
    for stem, h in [(STEMS[tok][0], H_NL), (STEMS[tok][1], H_LG)]:
        m[stem] = extents(os.path.join(outdir, stem + '.png'), w, h)
    return m


def main():
    tmp = '/tmp/fig2_layout'
    os.makedirs(tmp, exist_ok=True)
    out = os.path.join(HERE, 'panels')

    state = {}
    # ---- pass A: true ink width on an over-wide canvas --------------------------------
    for tok in TOKENS:
        render(tok, tmp, 320.0, 40.0)
        m = measure(tok, tmp, 320.0)
        cw = max(e['cw'] for e in m.values())
        state[tok] = dict(cw=cw, w=2 * (cw + MX_R), ml=None, off=None)
        print(f'  pass A {tok:6s} ink width {cw:6.1f} mm  ->  canvas {state[tok]["w"]:6.1f} mm')

    # ---- solve the left margin, then iterate ------------------------------------------
    for it in range(4):
        todo = []
        for tok, st in state.items():
            if st['ml'] is None:
                st['ml'] = st['w'] / 2.0          # first guess: ink on the midline
                todo.append(tok)
                continue
            m = measure(tok, out, st['w'])
            l = min(e['l'] for e in m.values())
            r = min(e['r'] for e in m.values())
            cw = max(e['cw'] for e in m.values())
            off = l - st['ml']                     # plot-margin -> ink distance
            w_new = max(2 * (cw + MX_R), st['w'])
            ml_new = TARGET * w_new - off
            # a wider canvas is needed if the right side is tight
            if r < MIN_GAP:
                w_new = 2 * (cw + MX_R + (MIN_GAP - r))
                ml_new = TARGET * w_new - off
            if abs(l / st['w'] - TARGET) > 0.01 or r < MIN_GAP or abs(w_new - st['w']) > 0.5:
                st['w'], st['ml'] = w_new, ml_new
                todo.append(tok)
        if not todo:
            break
        print(f'  solve iteration {it+1}: {", ".join(todo)}')
        for tok in todo:
            render(tok, out, state[tok]['w'], state[tok]['ml'])

    print('\nfinal:')
    bad = 0
    for tok, st in state.items():
        m = measure(tok, out, st['w'])
        for stem, h in [(STEMS[tok][0], H_NL), (STEMS[tok][1], H_LG)]:
            e = m[stem]
            frac = e['l'] / st['w']
            ok = abs(frac - TARGET) <= 0.015 and min(e['r'], e['t'], e['b']) > 1.2
            bad += 0 if ok else 1
            print(f'  {stem:22s} {st["w"]:6.1f} x {h:4.0f} mm  left blank {frac*100:5.1f}%  '
                  f'gaps R={e["r"]:5.1f} T={e["t"]:5.1f} B={e["b"]:5.1f}  [{"OK" if ok else "PROBLEM"}]')
    print('\nproblems:', bad)
    json.dump({k: dict(w=v['w'], ml=v['ml'], cw=v['cw']) for k, v in state.items()},
              open(os.path.join(HERE, 'fig2_ij_layout.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
