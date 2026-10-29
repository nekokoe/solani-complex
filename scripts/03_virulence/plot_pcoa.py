#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCoA Plotting Script

This script reads pre-calculated PCoA results and generates plots.

Usage:
    python plot_pcoa.py -d data_dir -o output_file
    
Requirements:
    - pcoa_results.tsv in data directory
    - Optional: metadata.tsv, distance_matrix.tsv
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os
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
    parser = argparse.ArgumentParser(description='PCoA Plotting Script')
    parser.add_argument('-d', '--data_dir', required=True, 
                        help='Directory containing PCoA results (pcoa_results.tsv)')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output image file path (without extension)')
    parser.add_argument('--title', default='PCoA Plot', 
                        help='Plot title')
    parser.add_argument('--width', type=float, default=16, 
                        help='Image width (inches, default: 16)')
    parser.add_argument('--height', type=float, default=12, 
                        help='Image height (inches, default: 12)')
    parser.add_argument('--dpi', type=int, default=300, 
                        help='Image resolution (DPI, default: 300)')
    parser.add_argument('--fontsize', type=int, default=6, 
                        help='Font size for sample labels (default: 6)')
    parser.add_argument('--show_labels', action='store_true', 
                        help='Show sample labels on plot')
    parser.add_argument('--show_ellipse', action='store_true', 
                        help='Show confidence ellipses for each group')
    parser.add_argument('--ellipse_std', type=float, default=2.0, 
                        help='Standard deviation for ellipse (default: 2.0)')
    parser.add_argument('--max_ticks', type=int, default=0,
                        help='Cap the number of major ticks per axis; 0 keeps the matplotlib '
                             'default. Needed where an axis spans a few hundredths, because the '
                             'default locator then prints overlapping four-decimal labels')
    args = parser.parse_args()

    # 检查数据目录
    pcoa_file = os.path.join(args.data_dir, 'pcoa_results.tsv')
    if not os.path.exists(pcoa_file):
        logging.error(f"PCoA results file not found: {pcoa_file}")
        logging.info("Please run stats_pcoa.py first to generate PCoA results.")
        return

    # 读取PCoA结果
    try:
        pca_df = pd.read_csv(pcoa_file, sep='\t')
        logging.info(f"Loaded {len(pca_df)} samples from PCoA results")
    except Exception as e:
        logging.error(f"Error reading PCoA results: {e}")
        return

    # 检查元数据
    metadata_file = os.path.join(args.data_dir, 'metadata.tsv')
    if os.path.exists(metadata_file):
        try:
            metadata = pd.read_csv(metadata_file, sep='\t')
            if 'pc1_variance' in metadata.columns and 'pc2_variance' in metadata.columns:
                pc1_var = metadata['pc1_variance'].iloc[0]
                pc2_var = metadata['pc2_variance'].iloc[0]
                logging.info(f"PC1 variance: {pc1_var:.2f}%, PC2 variance: {pc2_var:.2f}%")
            else:
                pc1_var = None
                pc2_var = None
        except:
            pc1_var = None
            pc2_var = None
    else:
        pc1_var = None
        pc2_var = None

    # 设置matplotlib参数
    plt.rcParams['axes.unicode_minus'] = False

    # 绘制PCoA图
    fig, ax = plt.subplots(figsize=(args.width, args.height))
    
    # 获取唯一的分组和物种
    groups = pca_df['Group'].unique()
    species = pca_df['Species'].unique() if 'Species' in pca_df.columns else ['Unknown']
    
    # 只保留样本最多的10个物种，其他归为Others
    if len(species) > 10:
        species_counts = pca_df['Species'].value_counts()
        top_10_species = species_counts.head(10).index.tolist()
        other_species = species_counts.index[10:].tolist()
        
        logging.info(f"Total {len(species)} species, keeping top 10: {top_10_species}")
        logging.info(f"Other species merged into 'Others': {other_species}")
        
        # 将其他物种的样本合并为Others
        pca_df.loc[pca_df['Species'].isin(other_species), 'Species'] = 'Others'
        
        # 更新species列表，Others放最后
        species = [sp for sp in top_10_species if sp in pca_df['Species'].unique()]
        if 'Others' in pca_df['Species'].values:
            species.append('Others')
    
    # 定义颜色映射
    color_map = {}
    np.random.seed(42)
    used_colors = []
    
    # 颜色定义
    red_colors = ['#E41A1C', '#FF6B6B', '#FF8C8C', '#FF5252', '#F44336', '#D32F2F']
    orange_colors = ['#FF8C00', '#FFA726', '#FFB74D', '#FF9800', '#F57C00', '#EF6C00']
    soil_colors = ['#DAA520', '#E6BE8A', '#CCA01D', '#CD853F', '#DEB887', '#BC8F8F']
    plant_colors = ['#228B22', '#3CB371', '#2E8B57', '#4CAF50', '#66BB6A', '#81C784']
    env_colors = ['#1E90FF', '#6495ED', '#87CEEB', '#4682B4', '#5F9EA0', '#708090']
    unknown_color = '#555555'
    
    color_indices = {'red': 0, 'orange': 0, 'soil': 0, 'plant': 0, 'env': 0}
    
    other_colors = ['#377EB8', '#984EA3', '#FFFF33', '#A65628', '#F781BF',
                   '#999999', '#7570B3', '#E7298A', '#E6AB02']
    
    for group in groups:
        group_lower = group.lower()
        
        if group_lower in ['human', 'med']:
            idx = color_indices['red'] % len(red_colors)
            color = red_colors[idx]
            color_map[group] = color
            used_colors.append(color)
            color_indices['red'] += 1
        elif group_lower == 'mammal':
            idx = color_indices['orange'] % len(orange_colors)
            color = orange_colors[idx]
            color_map[group] = color
            used_colors.append(color)
            color_indices['orange'] += 1
        elif group_lower == 'unknown':
            color_map[group] = unknown_color
            used_colors.append(unknown_color)
        elif group_lower == 'soil':
            idx = color_indices['soil'] % len(soil_colors)
            color = soil_colors[idx]
            color_map[group] = color
            used_colors.append(color)
            color_indices['soil'] += 1
        elif group_lower == 'plant':
            idx = color_indices['plant'] % len(plant_colors)
            color = plant_colors[idx]
            color_map[group] = color
            used_colors.append(color)
            color_indices['plant'] += 1
        elif group_lower == 'env':
            idx = color_indices['env'] % len(env_colors)
            color = env_colors[idx]
            color_map[group] = color
            used_colors.append(color)
            color_indices['env'] += 1
        else:
            for color in other_colors:
                if color not in used_colors:
                    color_map[group] = color
                    used_colors.append(color)
                    break
    
    # 定义形状映射
    shape_map = ['o', 's', 'D', '^', 'v', '<', '>', 'p', '*', 'h', 
                'H', '+', 'x', 'd', '|', '_', 'P', 'X', '8', '1',
                '2', '3', '4']
    
    # 特殊物种的固定形状
    fixed_shape_species = {
        'solani': 'o',
        'falciformis': 'D',
        'keratoplastica': 's'
    }
    
    # 为每个物种分配形状
    species_to_shape = {}
    available_shapes = shape_map.copy()
    for sp in species:
        if sp in fixed_shape_species:
            species_to_shape[sp] = fixed_shape_species[sp]
            if fixed_shape_species[sp] in available_shapes:
                available_shapes.remove(fixed_shape_species[sp])
        else:
            if available_shapes:
                species_to_shape[sp] = available_shapes.pop(0)
            else:
                species_to_shape[sp] = shape_map[species.index(sp) % len(shape_map)]
    
    # 获取聚类分组
    if 'Clustered_Group' in pca_df.columns:
        clustered_groups = pca_df['Clustered_Group'].unique()
    else:
        clustered_groups = pca_df['Group'].unique()
    
    # 为聚类分组创建颜色映射
    clustered_color_map = {}
    for group in clustered_groups:
        if 'Clustered_Group' in pca_df.columns and group in pca_df['Clustered_Group'].unique():
            original_group = group.split('_')[0] if '_' in group else group
            if original_group in color_map:
                clustered_color_map[group] = color_map[original_group]
            else:
                clustered_color_map[group] = color_map.get(group, '#555555')
        else:
            clustered_color_map[group] = color_map.get(group, '#555555')
    
    # 分离group和species的图例
    legend_handles_group = {}
    legend_handles_species = {}
    
    # 绘制散点
    for group in groups:
        group_data = pca_df[pca_df['Group'] == group]
        
        if len(species) > 1:
            for sp in species:
                sp_data = group_data[group_data['Species'] == sp]
                if not sp_data.empty:
                    handle = ax.scatter(
                        sp_data['PC1'], 
                        sp_data['PC2'], 
                        alpha=0.8, 
                        s=50, 
                        color=color_map[group],
                        marker=species_to_shape.get(sp, 'o'),
                        edgecolors='none'
                    )
                    # 只添加不重复的图例项
                    if group not in legend_handles_group:
                        legend_handles_group[group] = (handle, color_map[group])
                    if sp not in legend_handles_species:
                        legend_handles_species[sp] = species_to_shape.get(sp, 'o')
        else:
            handle = ax.scatter(
                group_data['PC1'], 
                group_data['PC2'], 
                alpha=0.8, 
                s=50, 
                color=color_map[group],
                label=group
            )
            if group not in legend_handles_group:
                legend_handles_group[group] = (handle, color_map[group])
    
    # 绘制置信椭圆
    if args.show_ellipse:
        for group in clustered_groups:
            group_data = pca_df[pca_df['Clustered_Group'] == group]
            if len(group_data) >= 2:
                plot_confidence_ellipse(
                    group_data['PC1'].values,
                    group_data['PC2'].values,
                    ax,
                    n_std=args.ellipse_std,
                    edgecolor=clustered_color_map.get(group, '#555555'),
                    alpha=0.3
                )
    
    # 设置标题和标签
    if pc1_var and pc2_var:
        ax.set_xlabel(f'PC1 ({pc1_var:.1f}%)', fontsize=9)
        ax.set_ylabel(f'PC2 ({pc2_var:.1f}%)', fontsize=9)
        ax.set_title(f'{args.title}', fontsize=11)
    else:
        ax.set_xlabel('PC1', fontsize=9)
        ax.set_ylabel('PC2', fontsize=9)
        ax.set_title(f'{args.title}', fontsize=11)
    
    # 显示样本标签
    if args.show_labels and 'Label' in pca_df.columns:
        texts = []
        for idx, row in pca_df.iterrows():
            group = row['Group']
            label_color = color_map.get(group, '#555555')
            texts.append(ax.text(row['PC1'], row['PC2'], row['Label'], fontsize=args.fontsize, color=label_color))
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))
    
    # 创建独立的Group图例（颜色）
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    
    # Group图例
    group_handles = []
    for group, (handle, color) in legend_handles_group.items():
        group_handles.append(Patch(facecolor=color, label=group))
    
    # Species图例（形状）
    species_handles = []
    for sp, marker in legend_handles_species.items():
        species_handles.append(Line2D([0], [0], marker=marker, color='none', markerfacecolor='gray',
                                     markersize=7, label=sp, markeredgecolor='none'))
    
    # 创建两个独立的图例
    legend1 = ax.legend(handles=group_handles, title='Source', loc='upper left', 
                        fontsize=6, title_fontsize=7, framealpha=0.9)
    ax.add_artist(legend1)
    
    if len(species_handles) > 0:
        legend2 = ax.legend(handles=species_handles, title='Species', loc='upper right', 
                            fontsize=6, title_fontsize=7, framealpha=0.9)
        ax.add_artist(legend1)
    
    ax.grid(True, alpha=0.3)

    if args.max_ticks:
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(args.max_ticks))
        ax.yaxis.set_major_locator(MaxNLocator(args.max_ticks))

    plt.tight_layout()
    
    # 保存PNG
    plt.savefig(f'{args.output}.png', dpi=args.dpi, bbox_inches='tight')
    logging.info(f"Saved PNG to {args.output}.png")
    
    # 保存PDF
    plt.savefig(f'{args.output}.pdf', bbox_inches='tight')
    logging.info(f"Saved PDF to {args.output}.pdf")
    
    plt.close()
    
    logging.info("Plotting completed successfully!")


def plot_confidence_ellipse(x, y, ax, n_std=2.0, facecolor='none', edgecolor='red', alpha=0.3, **kwargs):
    """绘制数据点的置信椭圆"""
    if len(x) < 2:
        return
    
    cov = np.cov(x, y)
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    order = eigenvalues.argsort()[::-1]
    eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
    
    vx, vy = eigenvectors[:, 0]
    theta = np.degrees(np.arctan2(vy, vx))
    width, height = 2 * n_std * np.sqrt(eigenvalues)
    
    from matplotlib.patches import Ellipse
    ellipse = Ellipse(
        xy=(np.mean(x), np.mean(y)),
        width=width,
        height=height,
        angle=theta,
        facecolor=facecolor,
        edgecolor=edgecolor,
        alpha=alpha,
        **kwargs
    )
    ax.add_patch(ellipse)
    
    return ellipse


if __name__ == '__main__':
    main()
