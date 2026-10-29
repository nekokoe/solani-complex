"""Fast PLINK .bed genotype access + Fst/pi estimators used for Figure 2 H/I/J."""
import numpy as np

MISSING = -1


def read_bed(prefix):
    fam = [l.split()[0] for l in open(prefix + '.fam')]
    bim = [l.split() for l in open(prefix + '.bim')]
    n, m = len(fam), len(bim)
    raw = np.fromfile(prefix + '.bed', dtype=np.uint8)
    assert raw[0] == 0x6C and raw[1] == 0x1B, 'not a PLINK bed file'
    bps = (n + 3) // 4
    body = raw[3:3 + bps * m].reshape(m, bps)
    lut = np.array([0, MISSING, 1, 2], dtype=np.int16)
    g = np.empty((m, n), dtype=np.int16)
    for j in range(n):
        b = body[:, j // 4].astype(np.int16)
        g[:, j] = lut[(b >> (2 * (j % 4))) & 3]
    return fam, bim, g


def dosage_stats(g):
    """per-SNP: allele dosage sum, number of non-missing, N (2*non-missing)"""
    ok = g != MISSING
    d = np.where(ok, g, 0).sum(axis=1).astype(np.float64)
    n_ind = ok.sum(axis=1).astype(np.float64)
    return d, n_ind


def pi_group(g, idx):
    """Nei's unbiased pi per site (mean over SNPs) for the given sample index."""
    sub = g[:, idx]
    ok = sub != MISSING
    k = ok.sum(axis=1).astype(np.float64)
    p = np.where(ok, sub, 0).sum(axis=1).astype(np.float64) / np.maximum(2 * k, 1)
    with np.errstate(invalid='ignore', divide='ignore'):
        v = 2 * p * (1 - p) * (k / (k - 1))
    v[k < 2] = np.nan
    return v


def hudson_fst(g, i1, i2):
    """Hudson (1992) Fst per SNP; robust for unequal/小 sample sizes.

    returns (num, den) per-SNP components; ratio-of-sums is the genome-wide estimate.
    """
    a = g[:, i1]
    b = g[:, i2]
    oka = a != MISSING
    okb = b != MISSING
    n1 = oka.sum(axis=1).astype(np.float64)
    n2 = okb.sum(axis=1).astype(np.float64)
    p1 = np.where(oka, a, 0).sum(axis=1).astype(np.float64) / np.maximum(2 * n1, 1)
    p2 = np.where(okb, b, 0).sum(axis=1).astype(np.float64) / np.maximum(2 * n2, 1)
    both = (n1 >= 1) & (n2 >= 1)
    p1, p2 = p1[both], p2[both]
    n1, n2 = n1[both], n2[both]
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    ok = den > 0
    return num[ok], den[ok]


def fst_value(g, i1, i2):
    num, den = hudson_fst(g, i1, i2)
    return num.sum() / den.sum()
