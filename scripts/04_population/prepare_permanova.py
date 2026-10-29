#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""T2 审稿回应 #505/#61 — 预处理 PAV 矩阵供 PERMANOVA 使用。
pivot hog_pav.tsv (长表 SampleID/HOG_ID/Present) -> 样本 x HOG 二进制矩阵。
与正文 PCoA 口径一致: 保留可变 HOG (presence rate < 100% 且有变异)，
core 恒定 HOG (全部 211 株携带或全部缺失) 无判别信息, 但为与正文完整 PCoA
一致, 此处保留全部 HOG; 输出 tmp/permanova_pav_matrix.tsv。
"""
import pandas as pd

pav = pd.read_csv('05.pan-genome/orthofinder-shorter/N0.HOG.hog_pav.tsv', sep='\t')
mat = pav.pivot(index='SampleID', columns='HOG_ID', values='Present')
mat = mat.astype(int)
mat.to_csv('tmp/permanova_pav_matrix.tsv', sep='\t')
print("matrix:", mat.shape)
