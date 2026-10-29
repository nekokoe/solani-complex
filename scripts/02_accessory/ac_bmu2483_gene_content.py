#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AC-T1: BMU02483 基因成分富集分析。

科学问题: BMU02483 (临床株, 无 chromosome-level 组装, 核心系统发育居环境枝)
与相邻 non-clinical falciformis 株相比, 其特异携带的 HOG 是否富集在
clinical AC (BMU03267_0019 等) 的基因集上?
若是, 支持 "BMU02483 的独特基因恰好落在 clinical 来源 AC 上" —
即其可能通过 AC 获得临床相关基因成分 (non-clinical -> clinical 转化候选)。

分析:
1. BMU02483 特异 HOG = BMU02483 携带 (count>0) 且 tree-order 相邻 falciformis env 株携带率低
2. clinical AC 基因集 = BMU03267 + BMU01933 的临床 AC 上基因对应的 HOG
3. Fisher 检验: 特异 HOG 在 clinical AC 基因集 vs 全基因组背景的富集
4. 补充: BMU02483 全部携带 HOG 中 clinical AC 基因占比 vs 背景
"""
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, mannwhitneyu
import os

os.makedirs("06.accessory-chrs/BMU02483_gene_content", exist_ok=True)
OUT = "06.accessory-chrs/BMU02483_gene_content"

# ---------- 1. 数据加载 ----------
hog = pd.read_csv('05.pan-genome/orthofinder-shorter/N0.HOG.hog_count.tsv', sep='\t')
grp = pd.read_csv('05.pan-genome/bmu_ref_groups.211.txt', sep='\t')
gene2hog = pd.read_csv('05.pan-genome/orthofinder-shorter/N0.gene_vs_HOG.tsv',
                       sep='\t', usecols=['GeneID', 'SampleID', 'HOG_ID'])

# ---------- 2. BMU02483 特异 HOG ----------
b2483 = hog[hog['SampleID'] == 'BMU02483']
b2483_hogs = set(b2483.loc[b2483['GeneCount'] > 0, 'HOG_ID'])
# env 对照 = tree order 上与 BMU02483 相邻 (±10 位置) 的 falciformis env 株
# (排除远处/疑误标 env 株, 如 CR12 在 +21)
tree_order = [l.strip() for l in
              open('06.accessory-chrs/tgs_assemblies/tree_order_allSC.txt')]
b_idx = tree_order.index('BMU02483')
near = tree_order[max(0, b_idx - 10): b_idx + 11]
f_env_all = set(grp[(grp['Species'] == 'falciformis') &
                    (grp['MedEnv'] == 'env')]['BMU_ID'])
f_env = [s for s in near if s in f_env_all]
f_med = grp[(grp['Species'] == 'falciformis') & (grp['MedEnv'] == 'med')]['BMU_ID'].tolist()
print(f"BMU02483 carried HOGs: {len(b2483_hogs)}")
print(f"near env control: {len(f_env)} | falc med: {len(f_med)}")
print(f"near env: {f_env}")

# env carriage per HOG
env_hog = hog[(hog['SampleID'].isin(f_env)) & (hog['GeneCount'] > 0)]
env_counts = env_hog.groupby('HOG_ID').size()
# all carried HOGs across env
all_carried = set(hog[hog['SampleID'].isin(f_env)]['HOG_ID'])  # includes count=0 rows? filter
all_carried = set(hog[(hog['SampleID'].isin(f_env)) & (hog['GeneCount'] > 0)]['HOG_ID'])
print(f"falc env carried HOG union: {len(all_carried)}")

# BMU02483-specific: carried by 2483 but in 0 env isolates
specific = b2483_hogs - all_carried
print(f"BMU02483-specific HOGs (absent in near env): {len(specific)}")

# ---------- 3. Clinical AC 基因集 ----------
def build_ac_hog_sets(iso, ac_ids, gene2hog):
    """返回 AC 染色体 -> HOG set 映射"""
    gff = f'02.annotation/per_isolates/{iso}/funannotate/predict_results/fusarium.gff3'
    chr_genes = {}  # chr -> set of gene ids
    with open(gff) as f:
        for line in f:
            if '\tgene\t' not in line:
                continue
            parts = line.split('\t')
            chrom = parts[0]
            if chrom in ac_ids:
                import re
                m = re.search(r'ID=([^;]+)', parts[8])
                if m:
                    chr_genes.setdefault(chrom, set()).add(m.group(1))
    # map to HOGs
    iso_genes = gene2hog[gene2hog['SampleID'] == iso]
    g2h = dict(zip(iso_genes['GeneID'].str.replace('-T1', '', regex=False),
                   iso_genes['HOG_ID']))
    ac_hog = {}
    for chrom, genes in chr_genes.items():
        hogs = set()
        for g in genes:
            h = g2h.get(g)
            if h:
                hogs.add(h)
        ac_hog[chrom] = hogs
    return ac_hog

# BMU03267 clinical ACs (from kmer_clinical_ac.assessment/burden_controlled + Q4)
# 已知临床 AC: BMU03267_0013/0016/0018/0019/0022/0024 (falc), BMU01933_0015/0017/0021/0022
bmu03267_ac_ids = set(open('06.accessory-chrs/tgs_assemblies/BMU03267.accessory_chrs.id').read().split())
bmu01933_ac_ids = set(open('06.accessory-chrs/tgs_assemblies/BMU01933.accessory_chrs.id').read().split())

ac_hog_03267 = build_ac_hog_sets('BMU03267', bmu03267_ac_ids, gene2hog)
ac_hog_01933 = build_ac_hog_sets('BMU01933', bmu01933_ac_ids, gene2hog)

# known clinical ACs (Q4/kmer_clinical_ac.assessment)
clinical_ac_03267 = ['BMU03267_0013', 'BMU03267_0016', 'BMU03267_0018',
                     'BMU03267_0019', 'BMU03267_0022', 'BMU03267_0024']
clinical_ac_01933 = ['BMU01933_0015', 'BMU01933_0017', 'BMU01933_0021', 'BMU01933_0022']
clinical_ac_hogs = set()
for c in clinical_ac_03267:
    clinical_ac_hogs |= ac_hog_03267.get(c, set())
for c in clinical_ac_01933:
    clinical_ac_hogs |= ac_hog_01933.get(c, set())
print(f"clinical AC HOG set: {len(clinical_ac_hogs)} HOGs "
      f"(from {clinical_ac_03267} + {clinical_ac_01933})")

# ---------- 4. Fisher 富集检验 ----------
# 背景 = 全部 HOG (29080)
all_hogs = set(hog['HOG_ID'].unique())
# 特异 HOG 中 clinical AC 占比 vs 背景中 clinical AC 占比
a = len(specific & clinical_ac_hogs)          # specific & clinical-AC
b = len(specific - clinical_ac_hogs)          # specific & not
c = len((all_hogs - specific) & clinical_ac_hogs)  # non-specific & clinical-AC
d = len(all_hogs - specific - clinical_ac_hogs)    # non-specific & not
table = [[a, b], [c, d]]
odds, p = fisher_exact(table, alternative='greater')
print(f"\n=== Fisher: BMU02483-specific HOGs enriched in clinical AC gene set ===")
print(f"specific & clinicalAC: {a} | specific & non: {b}")
print(f"non-spec & clinicalAC: {c} | non-spec & non: {d}")
print(f"odds ratio = {odds:.2f}, P (one-sided) = {p:.3e}")

# 全携带 HOG (非特异) 的富集
carried = b2483_hogs
a2 = len(carried & clinical_ac_hogs)
b2 = len(carried - clinical_ac_hogs)
c2 = len((all_hogs - carried) & clinical_ac_hogs)
d2 = len(all_hogs - carried - clinical_ac_hogs)
odds2, p2 = fisher_exact([[a2, b2], [c2, d2]], alternative='greater')
print(f"\n=== Fisher: ALL carried HOGs enriched in clinical AC gene set ===")
print(f"carried & clinicalAC: {a2} | carried & non: {b2}")
print(f"odds ratio = {odds2:.2f}, P = {p2:.3e}")

# 对照: falc env 株携带 HOG 的临床 AC 富集 (预期不富集)
env_carried = all_carried
a3 = len(env_carried & clinical_ac_hogs)
b3 = len(env_carried - clinical_ac_hogs)
c3 = len((all_hogs - env_carried) & clinical_ac_hogs)
d3 = len(all_hogs - env_carried - clinical_ac_hogs)
odds3, p3 = fisher_exact([[a3, b3], [c3, d3]], alternative='greater')
print(f"\n=== Fisher: env-carried HOGs (对照) ===")
print(f"env & clinicalAC: {a3} | env & non: {b3}")
print(f"odds ratio = {odds3:.2f}, P = {p3:.3e}")

# ---------- 5. 结果保存 ----------
res = {
    'bm2483_carried_hogs': len(b2483_hogs),
    'falc_env_count': len(f_env),
    'bm2483_specific_hogs': len(specific),
    'clinical_ac_hog_set': len(clinical_ac_hogs),
    'specific_in_clinicalAC': a, 'specific_not_in_clinicalAC': b,
    'fisher_specific_odds': odds, 'fisher_specific_p': p,
    'carried_in_clinicalAC': a2, 'fisher_carried_odds': odds2, 'fisher_carried_p': p2,
    'env_in_clinicalAC': a3, 'fisher_env_odds': odds3, 'fisher_env_p': p3,
}
with open(f'{OUT}/bm2483_gene_content_stats.tsv', 'w') as f:
    for k, v in res.items():
        f.write(f"{k}\t{v}\n")

# specific HOGs with clinical AC annotation
spec_df = pd.DataFrame({'HOG_ID': sorted(specific)})
spec_df['in_clinical_AC'] = spec_df['HOG_ID'].isin(clinical_ac_hogs)
spec_df.to_csv(f'{OUT}/bm2483_specific_hogs.tsv', sep='\t', index=False)
print(f"\nsaved -> {OUT}/bm2483_gene_content_stats.tsv + bm2483_specific_hogs.tsv")
