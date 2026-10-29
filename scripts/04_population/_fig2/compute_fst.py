#!/usr/bin/env python3
"""compute_fst.py — panel J: 11x11 Fst dot matrix with individual-level permutation tests.

Groups (11, order as requested)
  F. falciformis   : all, pop1, pop2, pop3
  F. solani        : all, pop1, pop2, pop3
  F. keratoplastica: all, pop1, pop2

Estimator
  Weir & Cockerham (1984) per site, summarised as the mean over 10-kb windows (5-kb step)
  of the per-window mean per-site Fst -- i.e. exactly the vcftools
  `--weir-fst-pop A --weir-fst-pop B --fst-window-size 10000 --fst-window-step 5000`
  "MEAN_FST" column averaged over windows. Validated against the existing
  falciformis med_vs_env file (0.278573 vs 0.278579).

SNP set per pair
  same species  -> that species' maf5pct.recode bed
  cross species -> the whole-complex bed (211 isolates)

Permutation test (1000x)
  disjoint pairs : resample the A/B membership within A u B, sizes preserved
  nested pairs   : 'species-all' held fixed; the population is resampled as a random
                   subset of 'species-all' of the same size
  P = (1 + #{stat_perm >= stat_obs}) / (n_perm + 1); stat = the same windowed Fst statistic,
  computed on a SNP subset thinned to <= PERM_SNP markers (observed effect size
  reported on the full set).

Output: data/panelJ_fst_matrix.tsv
"""
import os, sys, itertools, time, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wc84

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
POP = os.path.join(ROOT, '04.population')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
os.makedirs(OUT, exist_ok=True)

WIN, STEP = 10000, 5000
N_PERM = int(os.environ.get('N_PERM', 1000))
PERM_SNP = int(os.environ.get('PERM_SNP', 200000))
ONLY = os.environ.get('ONLY_PAIRS')

SPECIES_DIR = {'solani': 'solani-short-solani',
               'falciformis': 'solani-short-falciformis',
               'keratoplastica': 'solani-short-keratoplastica'}
FULLSC = 'solani-short-fullSC/solani-bmu-ref-short.raw.snps.combined.genotype.maf5pct.recode'

# panel J group order (as requested)
J_ORDER = ['falciformis-all', 'falciformis-pop1', 'falciformis-pop2', 'falciformis-pop3',
           'solani-all', 'solani-pop1', 'solani-pop2', 'solani-pop3',
           'keratoplastica-all', 'keratoplastica-pop1', 'keratoplastica-pop2']


def bed_prefix(species=None):
    if species is None:
        return os.path.join(POP, FULLSC)
    return os.path.join(POP, SPECIES_DIR[species], '%s.raw.snps.combined.genotype.maf5pct.recode' % SPECIES_DIR[species])


def load_groups():
    g = {}
    for l in open(os.path.join(OUT, 'fig2_groups.tsv')):
        if l.startswith('sample'):
            continue
        s, sp, src, pop = l.rstrip('\n').split('\t')
        g[s] = (sp, pop)
    return g


def win_id_of(bim):
    chrom = np.array([r[0] for r in bim])
    pos = np.array([int(r[3]) for r in bim])
    start = 1 + STEP * ((pos - 1) // STEP)
    # map (chrom,start) -> dense id
    key = {}
    wid = np.empty(len(bim), dtype=np.int64)
    for i, (c, s) in enumerate(zip(chrom, start)):
        k = (c, s)
        j = key.get(k)
        if j is None:
            j = len(key)
            key[k] = j
        wid[i] = j
    return wid, len(key)


def windowed_mean_fst_from_site(fst, valid, wid, nwin):
    m = valid & ~np.isnan(fst)
    if not m.any():
        return np.nan
    s = np.bincount(wid[m], weights=fst[m], minlength=nwin)
    c = np.bincount(wid[m], minlength=nwin)
    ok = c > 0
    return float((s[ok] / c[ok]).mean())


def matrices(g):
    """build float32 dosage / non-missing / het matrices (used on the THINNED subset only)"""
    return ((np.where(g != -1, g, 0)).astype(np.float32),
            (g != -1).astype(np.float32),
            (g == 1).astype(np.float32))


def fst_from_vecs(D, O, H, v1, v2):
    """per-site WC84 Fst from indicator vectors (n_ind,)"""
    n1 = O @ v1
    n2 = O @ v2
    d1 = D @ v1
    d2 = D @ v2
    h1 = H @ v1
    h2 = H @ v2
    valid = (n1 >= 2) & (n2 >= 2)
    with np.errstate(invalid='ignore', divide='ignore'):
        p1 = d1 / (2 * n1)
        p2 = d2 / (2 * n2)
        q1 = h1 / np.maximum(n1, 1)
        q2 = h2 / np.maximum(n2, 1)
        nbar = (n1 + n2) / 2.0
        pbar = (n1 * p1 + n2 * p2) / (n1 + n2)
        hbar = (n1 * q1 + n2 * q2) / (n1 + n2)
        s2 = (n1 * (p1 - pbar) ** 2 + n2 * (p2 - pbar) ** 2) / nbar
        nc = 2 * nbar - (n1 ** 2 + n2 ** 2) / (2 * nbar)
        pq = pbar * (1 - pbar)
        A = (nbar / nc) * (s2 - (pq - 0.5 * s2 - hbar / 4.0) / (nbar - 1.0))
        B = (nbar / (nbar - 1.0)) * (pq - 0.5 * s2 - ((2 * nbar - 1) / (4 * nbar)) * hbar)
        C = hbar / 2.0
        fst = A / (A + B + C)
    return fst, valid


def main():
    groups = load_groups()
    gid2samples = {}
    for s, (sp, pop) in groups.items():
        gid2samples.setdefault('%s-all' % sp, []).append(s)
        gid2samples.setdefault(pop, []).append(s)
    for k in J_ORDER:
        print('%-26s n=%d' % (k, len(gid2samples.get(k, []))))

    beds = {}
    for tag in ['solani', 'falciformis', 'keratoplastica', None]:
        pre = bed_prefix(tag)
        fam = [l.split()[0] for l in open(pre + '.fam')]
        bim = [l.split() for l in open(pre + '.bim')]
        beds[tag] = {'pre': pre, 'fam': fam, 'bim': bim,
                     'idx': {s: i for i, s in enumerate(fam)}}
        print('bed %-16s n=%d snps=%d' % (tag or 'FULLSC', len(fam), len(bim)))

    pairs = [p for p in itertools.combinations(J_ORDER, 2)]
    if ONLY:
        want = {tuple(x.split('|')) for x in ONLY.split(',')}
        pairs = [p for p in pairs if p in want]
    print('\n%d pairs to evaluate\n' % len(pairs))

    results = []
    t0 = time.time()
    for n, (A, B) in enumerate(pairs, 1):
        spA = groups[gid2samples[A][0]][0]
        spB = groups[gid2samples[B][0]][0]
        same = spA == spB
        tag = spA if same else None
        bed = beds[tag]
        idxA = np.array([bed['idx'][s] for s in gid2samples[A] if s in bed['idx']])
        idxB = np.array([bed['idx'][s] for s in gid2samples[B] if s in bed['idx']])
        nested = same and (A.endswith('-all') or B.endswith('-all'))

        fam, bim, g = read_bed_subset(bed['pre'], bed['fam'], np.union1d(idxA, idxB) if not nested
                                     else np.arange(len(bed['fam'])))
        remap = {old: i for i, old in enumerate(np.union1d(idxA, idxB) if not nested else np.arange(len(bed['fam'])))}
        a = np.array([remap[i] for i in idxA])
        b = np.array([remap[i] for i in idxB])

        wid, nwin = win_id_of(bim)

        # ---- observed, full SNP set (computed from int16 genotypes: no big float matrices) ----
        fst_site, valid = wc84.wc84_site(g, a, b)
        fst_full = windowed_mean_fst_from_site(fst_site, valid, wid, nwin)

        # ---- permutation on a thinned SNP subset ----
        m = g.shape[0]
        step = max(1, int(np.ceil(m / PERM_SNP)))
        sel = np.arange(0, m, step)
        gsub = g[sel]
        del g
        Dt, Ot, Ht = matrices(gsub)
        del gsub
        widt = wid[sel]
        v1 = np.zeros(Dt.shape[1], dtype=np.float32); v1[a] = 1
        v2 = np.zeros(Dt.shape[1], dtype=np.float32); v2[b] = 1
        f_obs_t = windowed_mean_fst_from_site(*fst_from_vecs(Dt, Ot, Ht, v1, v2), widt, nwin)

        rng = np.random.default_rng(20260916)
        cnt = 0
        null_stats = []
        na, nb = len(a), len(b)
        for _ in range(N_PERM):
            if nested:
                # 'all' fixed; population resampled from 'all'
                if A.endswith('-all'):
                    fixed, target = a, b
                else:
                    fixed, target = b, a
                pool = fixed
                vA = np.zeros(Dt.shape[1], dtype=np.float32)
                vA[fixed] = 1
                pick = rng.choice(pool, size=len(target), replace=False)
                vB = np.zeros(Dt.shape[1], dtype=np.float32)
                vB[pick] = 1
                stat = windowed_mean_fst_from_site(*fst_from_vecs(Dt, Ot, Ht, vA, vB), widt, nwin)
            else:
                pool = np.concatenate([a, b])
                perm = rng.permutation(pool)
                vA = np.zeros(Dt.shape[1], dtype=np.float32)
                vB = np.zeros(Dt.shape[1], dtype=np.float32)
                vA[perm[:na]] = 1
                vB[perm[na:]] = 1
                stat = windowed_mean_fst_from_site(*fst_from_vecs(Dt, Ot, Ht, vA, vB), widt, nwin)
            if np.isfinite(stat):
                null_stats.append(stat)
                if stat >= f_obs_t:
                    cnt += 1
        p = (1 + cnt) / (N_PERM + 1)
        ns = np.asarray(null_stats, dtype=float)
        null_mean = float(ns.mean()) if ns.size else float('nan')
        null_sd = float(ns.std(ddof=1)) if ns.size > 1 else float('nan')
        z = (f_obs_t - null_mean) / null_sd if null_sd and np.isfinite(null_sd) and null_sd > 0 else float('nan')

        results.append(dict(group_a=A, group_b=B, fst=fst_full, fst_test_set=f_obs_t,
                            p_perm=p, n_snp_full=int(m), n_snp_test=int(sel.size),
                            nested=bool(nested), bed=tag or 'FULLSC',
                            n_a=int(na), n_b=int(nb),
                            null_mean=null_mean, null_sd=null_sd, z=z))
        print('[%2d/%d] %-24s %-24s Fst=%.4f  P=%.4f  %s  (%.0fs)'
              % (n, len(pairs), A, B, fst_full, p, 'nested' if nested else 'disjoint', time.time() - t0))
        flush(results)

    write_out(results)
    print('\ndone in %.1f min' % ((time.time() - t0) / 60))


def read_bed_subset(prefix, fam, wanted):
    """read the bed and keep only the columns in `wanted` (indices into fam)"""
    n, m = len(fam), None
    with open(prefix + '.bim') as fh:
        m = sum(1 for _ in fh)
    raw = np.fromfile(prefix + '.bed', dtype=np.uint8)
    bps = (n + 3) // 4
    body = raw[3:3 + bps * m].reshape(m, bps)
    lut = np.array([0, -1, 1, 2], dtype=np.int16)
    g = np.empty((m, len(wanted)), dtype=np.int16)
    for out_j, j in enumerate(wanted):
        bb = body[:, j // 4].astype(np.int16)
        g[:, out_j] = lut[(bb >> (2 * (j % 4))) & 3]
    bim = [l.split() for l in open(prefix + '.bim')]
    return [fam[i] for i in wanted], bim, g


OUTNAME = os.environ.get('OUTNAME', 'panelJ_fst_matrix')


def flush(results):
    with open(os.path.join(OUT, OUTNAME + '.partial.tsv'), 'w') as fh:
        fh.write('group_a\tgroup_b\tfst\tfst_test_set\tp_perm\tn_snp_full\tn_snp_test\tnested\tbed\tn_a\tn_b'
                 '\tnull_mean\tnull_sd\tz\n')
        for r in results:
            fh.write('%(group_a)s\t%(group_b)s\t%(fst).6f\t%(fst_test_set).6f\t%(p_perm).5f\t'
                     '%(n_snp_full)d\t%(n_snp_test)d\t%(nested)s\t%(bed)s\t%(n_a)d\t%(n_b)d\t'
                     '%(null_mean).6f\t%(null_sd).6f\t%(z).3f\n' % r)


def write_out(results):
    flush(results)
    os.replace(os.path.join(OUT, OUTNAME + '.partial.tsv'),
               os.path.join(OUT, OUTNAME + '.tsv'))
    print('wrote data/%s.tsv (%d pairs)' % (OUTNAME, len(results)))


if __name__ == '__main__':
    main()
