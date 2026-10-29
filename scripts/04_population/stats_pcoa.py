#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCoA Matrix Calculation Script

This script reads a TSV file containing SampleID, Feature_ID, and Count columns,
performs PCoA analysis, and saves the results to the output directory.

Usage:
    python stats_pcoa.py -i input.tsv -g groups.tsv -o output_dir [-f feature_list.txt]

Features:
    - PCoA calculation with distance metrics
    - Optional feature subset analysis using a feature list
    - Saves both full and subset data to output directory
"""

import argparse
import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import pdist, squareform
import logging
import warnings
warnings.filterwarnings('ignore')

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='PCoA Matrix Calculation Script')
    parser.add_argument('-i', '--input', required=True, 
                        help='Input long-format TSV file (SampleID, Feature_ID, Count)')
    parser.add_argument('-g', '--groups', required=True, 
                        help='Sample group information TSV file (SampleID, Group)')
    parser.add_argument('-s', '--species', required=False, default=None, 
                        help='Sample species information TSV file (SampleID, Species)')
    parser.add_argument('-o', '--output_dir', required=True, 
                        help='Output directory for PCoA results')
    parser.add_argument('-f', '--feature_list', required=False, default=None, 
                        help='Optional feature list file (one feature ID per line, no header)')
    parser.add_argument('-l', '--sample_list', required=False, default=None, 
                        help='Optional sample list file (one sample ID per line, no header). If not provided, all samples will be used.')
    parser.add_argument('--distance', 
                        choices=['braycurtis','canberra','chebyshev','cityblock','correlation','cosine',
                                'dice','euclidean','hamming','jaccard','jensenshannon','kulczynski1',
                                'mahalanobis','matching','minkowski','rogerstanimoto','russellrao',
                                'seuclidean','sokalmichener','sokalsneath','sqeuclidean','yule'], 
                        default='braycurtis', 
                        help='Distance metric to use (default: braycurtis)')
    args = parser.parse_args()

    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    
    # 读取数据
    try:
        data_df = pd.read_csv(args.input, sep='\t')
        groups_df = pd.read_csv(args.groups, sep='\t', header=None, names=['SampleID', 'Group'])
        
        # 检查是否提供了物种信息文件
        if args.species:
            species_df = pd.read_csv(args.species, sep='\t', header=None, names=['SampleID', 'Species'])
            # 合并分组和物种信息
            groups_df = pd.merge(groups_df, species_df, on='SampleID', how='outer')
        else:
            logging.warning("No species information file provided. All samples will be labeled as 'Unknown'.")
            groups_df['Species'] = 'Unknown'
    except Exception as e:
        logging.error(f"Error reading file: {e}")
        return

    # 转换为宽格式矩阵
    if len(data_df.columns) < 3:
        logging.error(f"Input data must have at least 3 columns. Found {len(data_df.columns)} columns.")
        return

    # 获取列名
    sample_col = data_df.columns[0]
    feature_col = data_df.columns[1]
    value_col = data_df.columns[2]

    logging.info(f"Detected columns: {sample_col} (SampleID), {feature_col} (FeatureID), {value_col} (Value)")

    # 转换为宽格式矩阵
    try:
        matrix_df = data_df.pivot(index=sample_col, columns=feature_col, values=value_col).fillna(0)
    except Exception as e:
        logging.error(f"Error transforming data: {e}")
        return
    
    # 保存原始数据信息
    logging.info(f"Original matrix: {matrix_df.shape[0]} samples x {matrix_df.shape[1]} features")
    
    # 处理feature list
    if args.feature_list:
        try:
            feature_list_df = pd.read_csv(args.feature_list, sep='\t', header=None)
            selected_features = feature_list_df[0].tolist()
            logging.info(f"Loaded {len(selected_features)} features from feature list")
            
            # 筛选数据
            common_features = list(set(matrix_df.columns) & set(selected_features))
            logging.info(f"Found {len(common_features)} features in both matrix and feature list")
            
            if len(common_features) == 0:
                logging.error("No common features found between matrix and feature list!")
                return
            
            # 保存feature list到输出目录
            feature_list_output = os.path.join(args.output_dir, 'used_feature_list.txt')
            with open(feature_list_output, 'w') as f:
                for feat in common_features:
                    f.write(f"{feat}\n")
            logging.info(f"Saved used feature list to {feature_list_output}")
            
            # 提取子集数据
            subset_matrix_df = matrix_df[common_features]
            
            # 保存子集矩阵
            subset_matrix_output = os.path.join(args.output_dir, 'subset_feature_matrix.tsv')
            subset_matrix_df.to_csv(subset_matrix_output, sep='\t')
            logging.info(f"Saved subset matrix to {subset_matrix_output}")
            
            # 使用子集数据进行后续分析
            matrix_df = subset_matrix_df
            
        except Exception as e:
            logging.error(f"Error processing feature list: {e}")
            return

    # 处理sample list - 筛选特定样本
    if args.sample_list:
        try:
            sample_list_df = pd.read_csv(args.sample_list, sep='\t', header=None)
            selected_samples = sample_list_df[0].tolist()
            logging.info(f"Loaded {len(selected_samples)} samples from sample list")
            
            # 筛选在矩阵和样本列表中都存在的样本
            common_samples = list(set(matrix_df.index) & set(selected_samples))
            logging.info(f"Found {len(common_samples)} samples in both matrix and sample list")
            
            if len(common_samples) == 0:
                logging.error("No common samples found between matrix and sample list!")
                return
            
            # 保存使用的sample list到输出目录
            sample_list_output = os.path.join(args.output_dir, 'used_sample_list.txt')
            with open(sample_list_output, 'w') as f:
                for samp in common_samples:
                    f.write(f"{samp}\n")
            logging.info(f"Saved used sample list to {sample_list_output}")
            
            # 筛选样本
            matrix_df = matrix_df.loc[common_samples]
            logging.info(f"Filtered matrix to {len(common_samples)} samples")
            
        except Exception as e:
            logging.error(f"Error processing sample list: {e}")
            return

    # 确保所有样本都有分组和物种信息
    missing_samples = set(matrix_df.index) - set(groups_df['SampleID'])
    if missing_samples:
        logging.warning(f"The following samples were not found in the group file: {', '.join(missing_samples)}")
        matrix_df = matrix_df.loc[groups_df['SampleID']]

    # 填充缺失的物种信息
    groups_df['Species'] = groups_df['Species'].fillna('Unknown')

    # 数据标准化
    if args.distance in ['braycurtis', 'jaccard']:
        logging.info(f"Using {args.distance} distance, skipping data normalization")
        scaled_data = matrix_df.values
    else:
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(matrix_df)

    # 计算距离矩阵
    logging.info(f"Calculating distance matrix using {args.distance} metric")
    distance_matrix = squareform(pdist(scaled_data, metric=args.distance))

    # 执行PCoA (使用MDS)
    logging.info(f"Performing PCoA using MDS with {args.distance} distance")
    mds = MDS(n_components=2, dissimilarity='precomputed', random_state=42, n_init=4, max_iter=300)
    principal_components = mds.fit_transform(distance_matrix)
    logging.info(f"PCoA completed, generated {principal_components.shape[0]} samples x {principal_components.shape[1]} dimensions")

    # 计算PCA以估计方差解释率
    logging.info("Calculating variance explained by PCA")
    pca = PCA(n_components=2)
    pca.fit(scaled_data)
    explained_variance_ratio = pca.explained_variance_ratio_ * 100
    logging.info(f"PC1 explained variance: {explained_variance_ratio[0]:.2f}%")
    logging.info(f"PC2 explained variance: {explained_variance_ratio[1]:.2f}%")

    # 创建结果DataFrame
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'])
    pca_df['SampleID'] = matrix_df.index
    pca_df = pd.merge(pca_df, groups_df, on='SampleID')
    
    # 添加方差解释率信息
    pca_df['PC1_Variance'] = explained_variance_ratio[0]
    pca_df['PC2_Variance'] = explained_variance_ratio[1]

    # 处理标签：隐藏GCA（含）之后的内容
    pca_df['Label'] = pca_df['SampleID'].apply(lambda x: (
        x[3:] if x.startswith('BMU') else 
        x.split('GCA')[0].rstrip('_') if 'GCA' in x else 
        x
    ))
    
    # 对human和med分组的样本进行聚类
    human_med_data = pca_df[pca_df['Group'].str.lower().isin(['human', 'med'])].copy()
    
    if len(human_med_data) > 0:
        if len(human_med_data) >= 3:
            dbscan = DBSCAN(eps=0.02, min_samples=5, metric='euclidean', n_jobs=-1)
            cluster_labels = dbscan.fit_predict(human_med_data[['PC1', 'PC2']])
            
            noise_count = (cluster_labels == -1).sum()
            valid_labels = [label for label in cluster_labels if label != -1]
            num_clusters = len(set(valid_labels)) if valid_labels else 0
            
            logging.info(f"DBSCAN clustering: {num_clusters} clusters found, {noise_count} noise points")
            
            cluster_map = {}
            for sample_id, group, label in zip(human_med_data['SampleID'], human_med_data['Group'], cluster_labels):
                if label != -1:
                    cluster_map[sample_id] = f"{group}_{label+1}"
                else:
                    cluster_map[sample_id] = f"{group}_Noise"
            
            pca_df['Clustered_Group'] = pca_df['SampleID'].map(cluster_map).fillna(pca_df['Group'])
            logging.info("Added Clustered_Group column for human/med samples")
        else:
            logging.warning("Too few samples for DBSCAN clustering")
            pca_df['Clustered_Group'] = pca_df['Group']
    else:
        pca_df['Clustered_Group'] = pca_df['Group']

    # 保存PCoA结果
    pcoa_output = os.path.join(args.output_dir, 'pcoa_results.tsv')
    pca_df.to_csv(pcoa_output, sep='\t', index=False)
    logging.info(f"Saved PCoA results to {pcoa_output}")
    
    # 保存距离矩阵
    distance_output = os.path.join(args.output_dir, 'distance_matrix.tsv')
    distance_df = pd.DataFrame(distance_matrix, index=matrix_df.index, columns=matrix_df.index)
    distance_df.to_csv(distance_output, sep='\t')
    logging.info(f"Saved distance matrix to {distance_output}")
    
    # 保存元数据信息
    metadata = {
        'distance_metric': args.distance,
        'n_samples': int(principal_components.shape[0]),
        'n_features': int(matrix_df.shape[1]),
        'pc1_variance': float(explained_variance_ratio[0]),
        'pc2_variance': float(explained_variance_ratio[1])
    }
    
    if args.feature_list:
        metadata['feature_list_used'] = 'yes'
        metadata['n_features_in_list'] = len(common_features)
    else:
        metadata['feature_list_used'] = 'no'
    
    if args.sample_list:
        metadata['sample_list_used'] = 'yes'
        metadata['n_samples_in_list'] = len(common_samples)
    else:
        metadata['sample_list_used'] = 'no'
    
    metadata_df = pd.DataFrame([metadata])
    metadata_output = os.path.join(args.output_dir, 'metadata.tsv')
    metadata_df.to_csv(metadata_output, sep='\t', index=False)
    logging.info(f"Saved metadata to {metadata_output}")
    
    logging.info("="*60)
    logging.info("PCoA calculation completed successfully!")
    logging.info(f"Results saved to: {args.output_dir}")
    if args.sample_list:
        logging.info(f"Note: Sample list was used, saved to: {os.path.join(args.output_dir, 'used_sample_list.txt')}")
    logging.info(f"To generate plot, run: python plot_pcoa.py -d {args.output_dir}")
    logging.info("="*60)


if __name__ == '__main__':
    main()
