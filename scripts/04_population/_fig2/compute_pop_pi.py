#!/usr/bin/env python3
"""compute_pop_pi.py — nucleotide diversity of each ADMIXTURE population (panel J groups).

Same convention as the panel-H pi files: per-site nucleotide diversity aggregated into
10-kb windows (5-kb step), summarised as the mean over windows. Computed from the PLINK
bed so the populations are defined identically to the Fst matrix.

Output: data/panelJ_pop_pi.tsv
"""
import os, sys, csv
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pg_tools as pt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
POP = os.path.join(ROOT, '04.population')
OUT = os.path.join(HERE, 'data')

WIN_MM, STEP = 10000, 5000
SPECIES_DIR = {'solani': 'solani-short-solani',
               'falciformis': 'solani-short-falciformis',
               'keratoplastica': 'solani-short-keratoplastica'}


def block_ids(bim):
    """assign every site to its 5-kb block, and record the block order per contig.

    vcftools reports 10-kb windows stepping by 5 kb, so each reported window is the
    union of TWO consecutive 5-kb blocks. We therefore aggregate at 5-kb block level and
    combine consecutive blocks below.
    """
    chrom = np.array([r[0] for r in bim])
    pos = np.array([int(r[3]) for r in bim])
    start = 1 + STEP * ((pos - 1) // STEP)
    key, wid = {}, np.empty(len(bim), dtype=np.int64)
    blocks = []          # (chrom, start) in discovery order
    for i, (c, s) in enumerate(zip(chrom, start)):
        k = (c, s)
        if k not in key:
            key[k] = len(key)
            blocks.append(k)
        wid[i] = key[k]
    return wid, blocks


def main():
    groups = {}
    for l in open(os.path.join(OUT, 'fig2_groups.tsv')):
        if l.startswith('sample'):
            continue
        s, sp, src, pop = l.rstrip('\n').split('\t')
        groups.setdefault((sp, pop), []).append(s)

    rows = []
    win_rows = []
    for sp, d in SPECIES_DIR.items():
        pre = os.path.join(POP, d, '%s.raw.snps.combined.genotype.maf5pct.recode' % d)
        fam, bim, g = pt.read_bed(pre)
        wid, blocks = block_ids(bim)
        nblk = len(blocks)
        next_same = np.full(nblk, -1, dtype=np.int64)
        for j in range(nblk - 1):
            if blocks[j + 1][0] == blocks[j][0]:
                next_same[j] = j + 1
        idx_of = {s: i for i, s in enumerate(fam)}
        # whole species + each population, same convention as panel H
        pops = sorted([pp for (ss, pp) in groups if ss == sp],
                      key=lambda x: int(x.rsplit('pop', 1)[1]))
        targets = [('%s-all' % sp, fam)] + [(pp, groups[(sp, pp)]) for pp in pops]
        for label, samples in targets:
            idx = np.array([idx_of[x] for x in samples if x in idx_of])
            if idx.size < 3:
                continue
            v = pt.pi_group(g, idx)
            ok = ~np.isnan(v)
            bs = np.bincount(wid[ok], weights=v[ok], minlength=nblk)
            bc = np.bincount(wid[ok], minlength=nblk)
            # 10-kb window starting at block j  =  block j  +  block j+1 (same contig)
            ws = bs.copy()
            wc = bc.copy()
            has_next = next_same >= 0
            ws[has_next] += bs[next_same[has_next]]
            wc[has_next] += bc[next_same[has_next]]
            m = wc > 0
            pi_win = ws[m] / WIN_MM          # vcftools convention: sum(pi)/window_bp
            rows.append(dict(group=label, species=sp, n_samples=idx.size,
                             n_windows=int(m.sum()),
                             mean_pi=float(pi_win.mean()),
                             median_pi=float(np.median(pi_win))))
            # per-window values, for the boxplot view (panel d of Supplementary Fig. S2)
            valid_j = [j for j in range(nblk) if m[j]]
            for k, j in enumerate(valid_j):
                win_rows.append(dict(group=label, species=sp,
                                     bin_id='%s:%d' % (blocks[j][0], blocks[j][1]),
                                     pi=float(pi_win[k])))
        del g

    with open(os.path.join(OUT, 'panelJ_pop_pi.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['group', 'species', 'n_samples', 'n_windows',
                                           'mean_pi', 'median_pi'], delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(OUT, 'panelJ_pop_pi_windows.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['group', 'species', 'bin_id', 'pi'], delimiter='\t')
        w.writeheader()
        for r in win_rows:
            w.writerow(r)
    print('%-24s %5s %8s %10s %10s' % ('group', 'n', 'nWin', 'mean_pi', 'median_pi'))
    for r in rows:
        print('%-24s %5d %8d %10.5f %10.5f' % (r['group'], r['n_samples'], r['n_windows'],
                                               r['mean_pi'], r['median_pi']))
    print('\nwrote data/panelJ_pop_pi.tsv')


if __name__ == '__main__':
    main()
