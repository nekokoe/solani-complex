#!/usr/bin/env python3
"""build_groups.py — assemble the Figure 2 H/I/J group definitions.

Outputs (into data/):
  fig2_groups.tsv   sample, species, source, pop_label, in_species_vcf
  pi_groups.tsv     group_id, species, source_stratum, n_samples   (9 groups for panel H/I)
  fst_groups.tsv    group_id, species, kind, n_samples            (11 groups for panel J)

Conventions
-----------
* source            : med = clinical, env = environmental, unknown dropped from strata
* populations       : ADMIXTURE (maf5pct core SNP set), argmax Q assignment
                      solani K=3, falciformis K=3, keratoplastica K=2
* pop label order   : clinical fraction DESC, tie-break n DESC  -> pop1..popK
"""
import os, collections, csv

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
POP = os.path.join(ROOT, '04.population')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
os.makedirs(OUT, exist_ok=True)

SPECIES = {
    'solani': ('solani-short-solani', 3),
    'falciformis': ('solani-short-falciformis', 3),
    'keratoplastica': ('solani-short-keratoplastica', 2),
}
PREFIX = {sp: '%s/%s/%s.raw.snps.combined.genotype.maf5pct.recode' % (POP, d, d)
          for sp, (d, _) in SPECIES.items()}


def load_meta():
    src, spc = {}, {}
    with open(os.path.join(ROOT, '05.stats_pan_by_species_shorter/bmu_ref_groups.med_env.txt')) as fh:
        for l in fh:
            a = l.split()
            src[a[0]] = a[1]
    with open(os.path.join(ROOT, '05.stats_pan_by_species_shorter/bmu_ref_groups.species.txt')) as fh:
        for l in fh:
            a = l.split()
            spc[a[0]] = a[1]
    return src, spc


def read_fam(prefix):
    return [l.split()[0] for l in open(prefix + '.fam')]


def read_q(prefix, k):
    return [[float(x) for x in l.split()] for l in open('%s.%d.Q' % (prefix, k))]


def main():
    src, spc = load_meta()
    rows, pi_groups, fst_groups = [], [], []
    for sp, (d, K) in SPECIES.items():
        pre = PREFIX[sp]
        fam = read_fam(pre)
        q = read_q(pre, K)
        assert len(fam) == len(q), (sp, len(fam), len(q))
        # argmax assignment
        hard = [max(range(K), key=lambda i: r[i]) for r in q]
        # order components by clinical fraction desc, then n desc
        stat = []
        for k in range(K):
            mem = [s for s, h in zip(fam, hard) if h == k]
            n = len(mem)
            med = sum(1 for s in mem if src.get(s) == 'med')
            stat.append((k, med / n if n else 0.0, n))
        order = [k for k, _, _ in sorted(stat, key=lambda t: (-t[1], -t[2]))]
        rank = {k: i + 1 for i, k in enumerate(order)}
        for s, h in zip(fam, hard):
            rows.append(dict(sample=s, species=sp, source=src.get(s, 'unknown'),
                             pop='%s-pop%d' % (sp, rank[h]),
                             qmax='%.4f' % max(q[fam.index(s)])))
        # panel H/I strata
        counts = collections.Counter(src.get(s, 'unknown') for s in fam)
        pi_groups.append(dict(group_id='%s-all' % sp, species=sp, source_stratum='all', n_samples=len(fam)))
        for stratum in ('env', 'med'):
            pi_groups.append(dict(group_id='%s-%s' % (sp, 'clinical' if stratum == 'med' else 'env'),
                                  species=sp, source_stratum=stratum, n_samples=counts.get(stratum, 0)))
        # panel J
        fst_groups.append(dict(group_id='%s-all' % sp, species=sp, kind='all', n_samples=len(fam)))
        for k in order:
            n = sum(1 for h in hard if h == k)
            fst_groups.append(dict(group_id='%s-pop%d' % (sp, rank[k]), species=sp, kind='pop', n_samples=n))

    with open(os.path.join(OUT, 'fig2_groups.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['sample', 'species', 'source', 'pop'], delimiter='\t')
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in w.fieldnames})
    with open(os.path.join(OUT, 'pi_groups.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['group_id', 'species', 'source_stratum', 'n_samples'], delimiter='\t')
        w.writeheader()
        for r in pi_groups:
            w.writerow(r)
    with open(os.path.join(OUT, 'fst_groups.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['group_id', 'species', 'kind', 'n_samples'], delimiter='\t')
        w.writeheader()
        for r in fst_groups:
            w.writerow(r)

    print('fig2_groups.tsv : %d samples' % len(rows))
    print('pi_groups.tsv   : %d groups' % len(pi_groups))
    for r in pi_groups:
        print('   %-24s n=%d' % (r['group_id'], r['n_samples']))
    print('fst_groups.tsv  : %d groups' % len(fst_groups))
    for r in fst_groups:
        print('   %-24s n=%d' % (r['group_id'], r['n_samples']))
    # pop composition sanity
    print('\npop composition (clinical fraction):')
    comp = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        c = comp[(r['species'], r['pop'])]
        if r['source'] == 'med':
            c[0] += 1
        c[1] += 1
    for k in sorted(comp):
        m, n = comp[k]
        print('   %-24s n=%-3d clinical=%-3d frac=%.2f' % (k[1], n, m, m / n))


if __name__ == '__main__':
    main()
