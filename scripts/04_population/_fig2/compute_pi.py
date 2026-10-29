#!/usr/bin/env python3
"""compute_pi.py — panel H (pi boxplot) and panel I (pairwise pi significance matrix).

Panel H : per-window nucleotide diversity (10-kb windows, 5-kb step) for 9 groups
          = {F. solani, F. falciformis, F. keratoplastica} x {all, environmental, clinical}
          source = existing `*.windowed.pi` files
          (vcftools --window-pi 10000 --window-pi-step 5000)

Panel I : two significance measures for every pair of groups, BH-FDR corrected:
          * `q_mwu`  Mann-Whitney U on all per-window pi values (the literal spec)
          * `q_block` block-wise label permutation, blocks of 200 consecutive windows
            (~1 Mb), 100,000 permutations. Windows are spatially autocorrelated, so the
            plain MWU is anticonservative; the block test is the defensible one.
          The permutation is exact for the difference in mean pi because only block SUMS
          are needed:  diff(p) = 2 * (S @ (by - bx)) / L.

Outputs (data/)
  panelH_pi_windows.tsv   group_id, species, stratum, bin_id, pi
  panelH_pi_summary.tsv   group_id, ..., n_windows, mean_pi, median_pi, sd_pi
  panelI_pi_pairwise.tsv  group_a, group_b, U, p_mwu, q_mwu, p_block, q_block,
                          mean_a, mean_b, delta_pi, sig_mwu, sig_block
"""
import os, itertools, math
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
POP = os.path.join(ROOT, '04.population')
OUT = os.path.join(HERE, 'data')
os.makedirs(OUT, exist_ok=True)

BLOCK_WIN = 200          # 200 windows x 5 kb step ~ 1 Mb
N_PERM = 100000

SPECIES_DIR = {'solani': 'solani-short-solani',
               'falciformis': 'solani-short-falciformis',
               'keratoplastica': 'solani-short-keratoplastica'}
ORDER = [('solani', 'all', ''), ('solani', 'env', '.env'), ('solani', 'clinical', '.med'),
         ('falciformis', 'all', ''), ('falciformis', 'env', '.env'), ('falciformis', 'clinical', '.med'),
         ('keratoplastica', 'all', ''), ('keratoplastica', 'env', '.env'), ('keratoplastica', 'clinical', '.med')]
LABEL = {'solani': 'F. solani', 'falciformis': 'F. falciformis', 'keratoplastica': 'F. keratoplastica'}
STRAT = {'all': 'all', 'env': 'environmental', 'clinical': 'clinical'}


def pi_path(sp, suffix):
    d = SPECIES_DIR[sp]
    return os.path.join(POP, d, '%s.raw.snps.combined.genotype.maf5pct.recode%s.windowed.pi' % (d, suffix))


def read_pi(path):
    rows = [l.split() for l in open(path)][1:]
    return ([r[0] for r in rows], [int(r[1]) for r in rows],
            np.array([int(r[3]) for r in rows]), np.array([float(r[4]) for r in rows]))


def bh_fdr(p):
    p = np.asarray(p, float)
    n = p.size
    o = np.argsort(p)
    out = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        v = min(prev, p[o[i]] * n / (i + 1))
        out[o[i]] = v
        prev = v
    return np.clip(out, 0, 1)


def _block_means(x, chrom, rng=None):
    """mean pi within each full block of BLOCK_WIN consecutive windows.

    Blocks are built inside contigs (never straddling one) and partial tail blocks are
    dropped so every block carries the same number of windows.
    """
    out = []
    i, n = 0, len(x)
    while i < n:
        c = chrom[i]
        j = i
        while j < n and chrom[j] == c and (j - i) < BLOCK_WIN:
            j += 1
        if (j - i) == BLOCK_WIN:
            out.append(x[i:j].mean())
        i = j
    return np.array(out)


def block_permutation(x, cx, y, cy, rng):
    """two-sided block permutation P for mean(x) - mean(y), UNPAIRED.

    The two series may live in different coordinate systems (different reference
    assemblies), so blocks are pooled and randomly re-assigned to the two groups while
    preserving each group's block count. Vectorised over all permutations.
    """
    bx = _block_means(x, cx)
    by = _block_means(y, cy)
    nx, ny = bx.size, by.size
    if nx < 4 or ny < 4:
        return float(bx.mean() - by.mean()), float('nan')
    pool = np.concatenate([bx, by])
    nb = pool.size
    obs = bx.mean() - by.mean()
    r = rng.random((N_PERM, nb))
    order = np.argsort(r, axis=1)
    A = order[:, :nx]
    sA = pool[A].sum(axis=1)
    tot = pool.sum()
    meanA = sA / nx
    meanB = (tot - sA) / (nb - nx)
    diff = meanA - meanB
    p = (1 + np.count_nonzero(np.abs(diff) >= abs(obs))) / (N_PERM + 1)
    return float(obs), float(p)


def main():
    data, order_keys, summ = {}, {}, []
    fh = open(os.path.join(OUT, 'panelH_pi_windows.tsv'), 'w')
    fh.write('group_id\tspecies\tstratum\tbin_id\tpi\n')
    for sp, stratum, suffix in ORDER:
        gid = '%s-%s' % (sp, stratum)
        chrom, start, nvar, pi = read_pi(pi_path(sp, suffix))
        data[gid] = (pi, chrom)
        order_keys[gid] = pi
        for c, s, v in zip(chrom, start, pi):
            fh.write('%s\t%s\t%s\t%s:%d\t%.8g\n' % (gid, sp, stratum, c, s, v))
        summ.append((gid, sp, stratum, len(pi), pi.mean(), float(np.median(pi)), float(pi.std(ddof=1))))
    fh.close()

    with open(os.path.join(OUT, 'panelH_pi_summary.tsv'), 'w') as fh:
        fh.write('group_id\tspecies\tstratum\tlabel\tn_windows\tmean_pi\tmedian_pi\tsd_pi\n')
        for gid, sp, stratum, n, m, med, sd in summ:
            fh.write('%s\t%s\t%s\t%s (%s)\t%d\t%.6f\t%.6f\t%.6f\n'
                     % (gid, sp, stratum, LABEL[sp], STRAT[stratum], n, m, med, sd))

    print('PANEL H — mean per-window pi (10 kb / 5 kb)')
    print('%-26s %6s %10s %10s' % ('group', 'nWin', 'mean_pi', 'median_pi'))
    for gid, sp, stratum, n, m, med, sd in summ:
        print('%-26s %6d %10.5f %10.5f' % (gid, n, m, med))

    gids = ['%s-%s' % (s, k) for s, k, _ in ORDER]
    pairs = list(itertools.combinations(gids, 2))
    rng = np.random.default_rng(20260916)
    rec = []
    for a, b in pairs:
        xa, ca = data[a]
        xb, cb = data[b]
        u, pm = stats.mannwhitneyu(xa, xb, alternative='two-sided', method='asymptotic')
        obs, pb = block_permutation(xa, ca, xb, cb, rng)
        rec.append([a, b, float(u), float(pm), xa.mean(), xb.mean(), obs, pb,
                    _block_means(xa, ca).size + _block_means(xb, cb).size])
    qm = bh_fdr([r[3] for r in rec])
    qb = bh_fdr([r[7] for r in rec])

    with open(os.path.join(OUT, 'panelI_pi_pairwise.tsv'), 'w') as fh:
        fh.write('group_a\tgroup_b\tU\tp_mwu\tq_mwu\tp_block\tq_block\tmean_a\tmean_b\tdelta_pi\t'
                 'n_blocks\tsig_mwu\tsig_block\n')
        for r, q1, q2 in zip(rec, qm, qb):
            fh.write('%s\t%s\t%.1f\t%.6g\t%.6g\t%.6g\t%.6g\t%.6f\t%.6f\t%.6f\t%d\t%s\t%s\n'
                     % (r[0], r[1], r[2], r[3], q1, r[7], q2, r[4], r[5], r[6], r[8],
                        'yes' if q1 < 0.05 else 'no', 'yes' if q2 < 0.05 else 'no'))

    print('\nPANEL I — %d pairs' % len(rec))
    print('  MWU      : %d/36 significant, median -log10(q)=%.1f' % ((qm < 0.05).sum(), np.median(-np.log10(np.clip(qm, 1e-300, None)))))
    print('  blockperm: %d/36 significant, median -log10(q)=%.1f' % ((qb < 0.05).sum(), np.median(-np.log10(np.clip(qb, 1e-300, None)))))
    print('\n  block-permutation result (q values):')
    print('  %-24s %-24s %10s %10s %10s' % ('group A', 'group B', 'delta_pi', 'q_block', 'sig'))
    for r, q1, q2 in zip(rec, qm, qb):
        print('  %-24s %-24s %10.5f %10.3g %10s' % (r[0], r[1], r[6], q2, 'yes' if q2 < 0.05 else 'NO'))


if __name__ == '__main__':
    main()
