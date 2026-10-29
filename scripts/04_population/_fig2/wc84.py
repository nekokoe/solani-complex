"""wc84.py — vectorised Weir & Cockerham (1984) Fst on PLINK bed genotypes.

Reproduces the statistic that `vcftools --weir-fst-pop A --weir-fst-pop B
--fst-window-size 10000 --fst-window-step 5000` reports as the mean of the
per-window MEAN_FST column (mean over sites of the per-site Fst).
"""
import numpy as np

MISSING = -1


def wc84_site(g, i1, i2):
    """Per-site Fst for two population index sets. Returns (fst, valid_mask).

    g: (n_snp, n_ind) int16 array with -1 = missing, 0/1/2 = allele dosage.
    """
    a1 = g[:, i1]
    a2 = g[:, i2]
    o1 = a1 != MISSING
    o2 = a2 != MISSING
    n1 = o1.sum(axis=1).astype(np.float64)
    n2 = o2.sum(axis=1).astype(np.float64)
    d1 = np.where(o1, a1, 0).sum(axis=1).astype(np.float64)
    d2 = np.where(o2, a2, 0).sum(axis=1).astype(np.float64)
    h1 = np.where(o1 & (a1 == 1), 1, 0).sum(axis=1).astype(np.float64)
    h2 = np.where(o2 & (a2 == 1), 1, 0).sum(axis=1).astype(np.float64)

    valid = (n1 >= 2) & (n2 >= 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        p1 = d1 / (2 * n1)
        p2 = d2 / (2 * n2)
        h1 = h1 / n1
        h2 = h2 / n2
        nbar = (n1 + n2) / 2.0
        pbar = (n1 * p1 + n2 * p2) / (n1 + n2)
        hbar = (n1 * h1 + n2 * h2) / (n1 + n2)
        s2 = (n1 * (p1 - pbar) ** 2 + n2 * (p2 - pbar) ** 2) / nbar
        nc = 2 * nbar - (n1 ** 2 + n2 ** 2) / (2 * nbar)
        pq = pbar * (1 - pbar)
        A = (nbar / nc) * (s2 - (pq - 0.5 * s2 - hbar / 4.0) / (nbar - 1.0))
        B = (nbar / (nbar - 1.0)) * (pq - 0.5 * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
        C = hbar / 2.0
        fst = A / (A + B + C)
    fst = np.where(valid, fst, np.nan)
    return fst, valid


def window_mean_fst(g, bim, i1, i2, win=10000, step=5000, chrom_col=0, pos_col=3):
    """Mean over windows of the per-window mean per-site Fst (vcftools convention)."""
    fst, valid = wc84_site(g, i1, i2)
    chrom = np.array([r[chrom_col] for r in bim])
    pos = np.array([int(r[pos_col]) for r in bim])
    out = []
    for c in np.unique(chrom):
        sel = chrom == c
        if not sel.any():
            continue
        pc = pos[sel]
        fc = fst[sel]
        vc = valid[sel]
        lo, hi = 1, pc.max()
        for s in range(lo, hi + 1, step):
            m = vc & (pc >= s) & (pc <= s + win - 1)
            if m.sum() == 0:
                continue
            vals = fc[m]
            vals = vals[~np.isnan(vals)]
            if vals.size == 0:
                continue
            out.append((c, s, s + win - 1, vals.size, vals.mean()))
    return out
