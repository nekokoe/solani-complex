#!/usr/bin/env python3
"""panelI_options.py — compare candidate tests for the pi dot matrix (panel I).

(A) plain Mann-Whitney U on all 10-kb windows (5-kb step) + BH-FDR   [literal spec]
(B) Mann-Whitney U on DISJOINT 10-kb windows only (every 2nd window)
(C) circular block permutation of the window series (block = 100 windows = 1 Mb)
    using the difference in mean pi as the statistic, 1000 permutations
(D) reported as effect size: |delta mean pi| with MWU significance as colour

Prints, for each option, how many of the 36 cells are "significant" and the spread of
the significance measure, so the dot-size encoding can be judged.
"""
import os, itertools, math
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

ORDER = ['solani-all', 'solani-env', 'solani-clinical',
         'falciformis-all', 'falciformis-env', 'falciformis-clinical',
         'keratoplastica-all', 'keratoplastica-env', 'keratoplastica-clinical']

BLOCK = 100
N_PERM = 1000


def load():
    d = {}
    for l in open(os.path.join(DATA, 'panelH_pi_windows.tsv')):
        if l.startswith('group_id'):
            continue
        gid, sp, st, binid, pi = l.rstrip('\n').split('\t')
        d.setdefault(gid, []).append(float(pi))
    return {k: np.array(v) for k, v in d.items()}


def bh(p):
    p = np.asarray(p, float); n = p.size
    o = np.argsort(p); r = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        v = min(prev, p[o[i]] * n / (i + 1)); r[o[i]] = v; prev = v
    return np.clip(r, 0, 1)


def main():
    data = load()
    pairs = list(itertools.combinations(ORDER, 2))
    lbl = {k: k.replace('solani', 'So').replace('falciformis', 'Fa').replace('keratoplastica', 'Ke')
              .replace('-clinical', '-cli').replace('-environmental', '-env') for k in ORDER}

    # ---- (A) all windows ----
    qA = []
    for a, b in pairs:
        _, p = stats.mannwhitneyu(data[a], data[b], alternative='two-sided', method='asymptotic')
        qA.append(p)
    qA = bh(qA)

    # ---- (B) disjoint windows only (start offset consistent) ----
    # bins are labelled chrom:start; keep every 2nd window by position within chromosome
    def disjoint(arr, keys):
        return arr[::2]
    idx = {}
    for l in open(os.path.join(DATA, 'panelH_pi_windows.tsv')):
        if l.startswith('group_id'):
            continue
        gid, sp, st, binid, pi = l.rstrip('\n').split('\t')
        idx.setdefault(gid, []).append(binid)
    dA = {k: data[k][::2] for k in ORDER}
    qB = []
    for a, b in pairs:
        _, p = stats.mannwhitneyu(dA[a], dA[b], alternative='two-sided', method='asymptotic')
        qB.append(p)
    qB = bh(qB)

    # ---- (C) block-wise label permutation (blocks of 100 consecutive windows) ----
    rng = np.random.default_rng(11)
    pC, dC = [], []
    for a, b in pairs:
        x, y = data[a], data[b]
        L = min(len(x), len(y))
        x, y = x[:L].copy(), y[:L].copy()
        obs = abs(x.mean() - y.mean())
        dC.append(x.mean() - y.mean())
        nb = int(np.ceil(L / BLOCK))
        cnt = 0
        for _ in range(N_PERM):
            xs, ys = x.copy(), y.copy()
            swap = rng.random(nb) < 0.5
            for bi in np.nonzero(swap)[0]:
                lo, hi = bi * BLOCK, min((bi + 1) * BLOCK, L)
                xs[lo:hi], ys[lo:hi] = y[lo:hi], x[lo:hi]
            if abs(xs.mean() - ys.mean()) >= obs:
                cnt += 1
        pC.append((1 + cnt) / (N_PERM + 1))
    pC = np.array(pC)
    qC = bh(pC)

    def report(name, q):
        sig = (q < 0.05).sum()
        nz = q[q > 0]
        neg = -np.log10(nz) if nz.size else np.array([0.0])
        print('%-52s sig=%2d/36  -log10q: min=%.1f med=%.1f max=%.1f  underflow=%d'
              % (name, sig, neg.min(), np.median(neg), neg.max(), int((q == 0).sum())))

    print('CELL SIGNIFICANCE SUMMARY (36 pairs)')
    report('(A) Mann-Whitney, all 10-kb windows (5-kb step) + FDR', qA)
    report('(B) Mann-Whitney, disjoint 10-kb windows only + FDR', qB)
    report('(C) circular block permutation (1 Mb blocks), BH-FDR', qC)

    print('\nPER-PAIR DETAIL  (values are -log10 of the adjusted P)')
    print('%-24s %-24s %8s %8s %8s %9s' % ('group A', 'group B', 'A', 'B', 'C', 'delta_pi'))
    for (a, b), va, vb, vc, dd in zip(pairs, qA, qB, qC, dC):
        def f(v):
            return 'under' if v == 0 else '%.1f' % -math.log10(v)
        print('%-24s %-24s %8s %8s %8s %9.5f' % (a, b, f(va), f(vb), f(vc), dd))

    np.savetxt(os.path.join(DATA, 'panelI_options.tsv'),
               np.column_stack([[f'{a}|{b}' for a, b in pairs], qA, qB, qC, dC]),
               fmt='%s', delimiter='\t',
               header='pair\tq_mwu_all\tq_mwu_disjoint\tq_blockperm\tdelta_mean_pi', comments='')
    print('\nwrote data/panelI_options.tsv')


if __name__ == '__main__':
    main()
