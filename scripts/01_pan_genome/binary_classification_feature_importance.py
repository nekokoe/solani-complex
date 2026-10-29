#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
二分类机器学习训练与特征重要性分析

该脚本实现对指定样本的2分类机器学习训练，并找出与分类最相关的特征。
支持长格式特征矩阵输入，使用易于解释的机器学习模型，并输出特征重要性排序。
"""

import os
import argparse
import pandas as pd
import numpy as np
import joblib
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, accuracy_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from scipy import stats


def load_and_transform_feature_matrix(feature_matrix_file, all_samples_file):
    """
    加载特征矩阵并转换为宽格式
    
    参数:
        feature_matrix_file: 特征矩阵文件路径 (长格式tsv)
        all_samples_file: 总体样本ID文件路径
        
    返回:
        features_df: 转换后的宽格式特征矩阵 (样本×特征)
    """
    print(f"读取特征矩阵文件: {feature_matrix_file}")
    # 读取特征矩阵（长格式）
    feature_df = pd.read_csv(feature_matrix_file, sep='\t')
    
    # 检查列数是否符合要求
    if feature_df.shape[1] < 3:
        raise ValueError("特征矩阵文件必须至少包含3列: 样本ID, 特征ID, 特征值")
    
    # 读取总体样本ID
    print(f"读取总体样本ID文件: {all_samples_file}")
    all_samples_df = pd.read_csv(all_samples_file, sep='\t', header=None, names=['SampleID'])
    all_sample_ids = set(all_samples_df['SampleID'].tolist())
    
    # 过滤特征矩阵，只保留总体样本ID中的样本
    filtered_feature_df = feature_df[feature_df.iloc[:, 0].isin(all_sample_ids)]
    
    # 转换为宽格式 (样本×特征)
    # 使用第一列为样本ID，第二列为特征ID，第三列为特征值
    features_df = filtered_feature_df.pivot(index=filtered_feature_df.columns[0], 
                                           columns=filtered_feature_df.columns[1], 
                                           values=filtered_feature_df.columns[2])
    
    print(f"特征矩阵转换完成。原始样本数: {len(feature_df.iloc[:, 0].unique())}, "
          f"过滤后样本数: {len(features_df)}, "
          f"特征数量: {features_df.shape[1]}")
    
    return features_df


def create_labels(features_df, positive_samples_file):
    """
    创建分类标签
    
    参数:
        features_df: 特征矩阵
        positive_samples_file: 阳性样本ID文件路径
        
    返回:
        y: 标签数组
        label_counts: 标签统计信息
    """
    print(f"读取阳性样本ID文件: {positive_samples_file}")
    positive_samples_df = pd.read_csv(positive_samples_file, sep='\t', header=None, names=['SampleID'])
    positive_sample_ids = set(positive_samples_df['SampleID'].tolist())
    
    # 创建标签：阳性样本为1，其他为0
    y = np.zeros(len(features_df), dtype=int)
    for i, sample_id in enumerate(features_df.index):
        if sample_id in positive_sample_ids:
            y[i] = 1
    
    # 统计标签分布
    label_counts = pd.Series(y).value_counts()
    print(f"标签创建完成。类别分布: {dict(label_counts)}")
    
    # 检查类别是否平衡
    if min(label_counts) < 5:
        print("警告: 样本数量较少，可能会影响模型性能。建议阳性和阴性样本数量均不少于5个。")
    
    return y, label_counts


def train_model(X, y, random_state=42):
    """
    训练随机森林分类模型
    
    参数:
        X: 特征矩阵
        y: 标签
        random_state: 随机种子
        
    返回:
        model: 训练好的模型
        X_train, X_test, y_train, y_test: 训练集和测试集
    """
    print("分割训练集和测试集...")
    # 分割训练集和测试集（70%训练，30%测试）
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=random_state, stratify=y
    )
    
    print(f"训练集大小: {len(X_train)}, 测试集大小: {len(X_test)}")
    
    # 初始化随机森林分类器
    # 设置n_jobs=-1使用所有可用的CPU核心加速训练
    model = RandomForestClassifier(
        n_estimators=100,          # 树的数量
        max_depth=10,              # 树的最大深度，防止过拟合
        min_samples_split=2,       # 内部节点再划分所需最小样本数
        min_samples_leaf=1,        # 叶子节点最少样本数
        random_state=random_state, # 随机种子
        class_weight='balanced',   # 自动平衡类别权重
        n_jobs=-1                  # 使用所有可用CPU核心
    )
    
    print("开始训练随机森林模型...")
    model.fit(X_train, y_train)
    
    # 交叉验证评估模型性能
    print("进行5折交叉验证...")
    cv_scores = cross_val_score(model, X, y, cv=5)
    print(f"交叉验证准确率: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")
    
    return model, X_train, X_test, y_train, y_test


def evaluate_model(model, X_test, y_test, output_dir):
    """
    评估模型性能
    
    参数:
        model: 训练好的模型
        X_test: 测试集特征
        y_test: 测试集标签
        output_dir: 输出目录
        
    返回:
        metrics: 模型评估指标
    """
    print("评估模型性能...")
    
    # 在测试集上进行预测
    y_pred = model.predict(X_test)
    y_pred_proba = model.predict_proba(X_test)[:, 1]  # 预测为阳性类别的概率
    
    # 计算评估指标
    accuracy = accuracy_score(y_test, y_pred)
    auc_score = roc_auc_score(y_test, y_pred_proba)
    
    # 生成分类报告
    report = classification_report(y_test, y_pred, output_dict=True)
    
    # 生成混淆矩阵
    conf_matrix = confusion_matrix(y_test, y_pred)
    
    # 打印评估结果
    print(f"测试集准确率: {accuracy:.4f}")
    print(f"ROC-AUC 分数: {auc_score:.4f}")
    print("分类报告:")
    print(classification_report(y_test, y_pred))
    
    # 可视化混淆矩阵
    plt.figure(figsize=(8, 6))
    sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Negative', 'Positive'], yticklabels=['Negative', 'Positive'])
    plt.xlabel('Predicted Labels')
    plt.ylabel('True Labels')
    plt.title('Confusion Matrix')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'confusion_matrix.png'), dpi=300)
    plt.close()
    
    # 保存评估指标
    metrics = {
        'accuracy': accuracy,
        'roc_auc': auc_score,
        'classification_report': report,
        'confusion_matrix': conf_matrix.tolist()
    }
    
    # 绘制ROC-AUC曲线
    fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {auc_score:.3f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'roc_curve.png'), dpi=300)
    plt.close()
    
    return metrics


def get_feature_importance(model, X, y, feature_names, output_dir, top_n=300):
    """
    获取并可视化特征重要性，并计算每个特征的显著性p值和分布统计信息
    
    参数:
        model: 训练好的模型
        X: 特征矩阵
        y: 标签
        feature_names: 特征名称列表
        output_dir: 输出目录
        top_n: 保留的最重要特征数量
        
    返回:
        top_features_df: top_n个最重要的特征及其重要性、p值和分布统计信息
    """
    # 获取特征重要性（基于树模型的内置重要性）
    feature_importance = model.feature_importances_
    
    print("计算每个特征的显著性p值和分布统计信息...")
    # 为每个特征计算Mann-Whitney U检验的p值（适合两组样本数量不一致的情况）
    p_values = []
    # 新增统计信息列表
    positive_samples_total = []  # 阳性组样本总数
    positive_means = []
    positive_vars = []
    positive_counts = []  # count > 0 的样本数
    negative_samples_total = []  # 阴性组样本总数
    negative_means = []
    negative_vars = []
    negative_counts = []  # count > 0 的样本数
    
    for feature_idx, feature_name in enumerate(feature_names):
        # 从矩阵中提取每个特征包含的样本的原始数据
        feature_values = X.iloc[:, feature_idx].values
        
        # 按照positive和negative的分组
        positive_group = feature_values[y == 1]
        negative_group = feature_values[y == 0]
        
        # 计算分布统计信息
        # 阳性组统计
        pos_total = len(positive_group)  # 阳性组样本总数
        pos_mean = np.mean(positive_group)
        pos_var = np.var(positive_group)
        pos_count = np.sum(positive_group > 0)  # count > 0 的样本数
        
        # 阴性组统计
        neg_total = len(negative_group)  # 阴性组样本总数
        neg_mean = np.mean(negative_group)
        neg_var = np.var(negative_group)
        neg_count = np.sum(negative_group > 0)  # count > 0 的样本数
        
        # 添加到列表
        positive_samples_total.append(pos_total)
        positive_means.append(pos_mean)
        positive_vars.append(pos_var)
        positive_counts.append(pos_count)
        negative_samples_total.append(neg_total)
        negative_means.append(neg_mean)
        negative_vars.append(neg_var)
        negative_counts.append(neg_count)
        
        # 使用Mann-Whitney U检验计算两组之间的显著性
        try:
            # 检查数组是否为空
            if len(positive_group) == 0 or len(negative_group) == 0:
                print(f"警告: 特征 {feature_name} 的阳性或阴性样本组为空，无法计算p值。")
                p_values.append(np.nan)
            else:
                # 直接使用Mann-Whitney U检验
                _, p_value = stats.mannwhitneyu(positive_group, negative_group)
                p_values.append(p_value)
        except Exception as e:
            print(f"计算特征 {feature_name} 的p值时出错: {str(e)}")
            p_values.append(np.nan)
    
    # 创建特征重要性DataFrame
    importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': feature_importance,
        # 阳性组统计信息
        'Positive_Exist': positive_counts,
        'Positive_Total': positive_samples_total,
        'Positive_Mean': positive_means,
        'Positive_Variance': positive_vars,
        # 阴性组统计信息
        'Negative_Exist': negative_counts,
        'Negative_Total': negative_samples_total,
        'Negative_Mean': negative_means,
        'Negative_Variance': negative_vars,
        # 显著性p值
        'P_value': p_values
    })
    
    # 保留浮点数据的8位有效数字
    float_columns = ['Importance', 'Positive_Exist', 'Positive_Total', 'Positive_Mean', 'Positive_Variance', 
                     'Negative_Exist', 'Negative_Total', 'Negative_Mean', 'Negative_Variance', 'P_value']
    for col in float_columns:
        if col in importance_df.columns:
            importance_df[col] = importance_df[col].apply(lambda x: '{0:.8g}'.format(x) if pd.notnull(x) else x)
    
    # 将格式化后的字符串转换回浮点数类型
    for col in float_columns:
        if col in importance_df.columns:
            importance_df[col] = pd.to_numeric(importance_df[col], errors='coerce')
    
    # 计算Bonferroni校正的p值
    num_tests = len(feature_names)
    importance_df['Bonferroni_p_value'] = importance_df['P_value'] * num_tests
    # 确保校正后的p值不超过1
    importance_df['Bonferroni_p_value'] = importance_df['Bonferroni_p_value'].clip(upper=1.0)
    
    # 为Bonferroni_p_value列也保留8位有效数字
    if 'Bonferroni_p_value' in importance_df.columns:
        importance_df['Bonferroni_p_value'] = importance_df['Bonferroni_p_value'].apply(lambda x: '{0:.8g}'.format(x) if pd.notnull(x) else x)
        importance_df['Bonferroni_p_value'] = pd.to_numeric(importance_df['Bonferroni_p_value'], errors='coerce')
    
    # 计算Positive_Rate和Negative_Rate
    importance_df['Positive_Rate'] = importance_df['Positive_Exist'] / importance_df['Positive_Total']
    importance_df['Negative_Rate'] = importance_df['Negative_Exist'] / importance_df['Negative_Total']
    
    # 删除原始的Exist和Total列
    importance_df = importance_df.drop(['Positive_Exist', 'Positive_Total', 'Negative_Exist', 'Negative_Total'], axis=1)
    
    # 按重要性排序
    importance_df = importance_df.sort_values('Importance', ascending=False)
    
    # 输出所有重要性不为0的features
    all_non_zero_features_df = importance_df[importance_df['Importance'] > 0].copy()
    print(f"获取 {len(all_non_zero_features_df)} 个重要性不为0的特征...")
    
    # 保存所有重要性不为0的特征到all_features.txt
    all_features_file = os.path.join(output_dir, 'all_features.txt')
    all_non_zero_features_df.to_csv(all_features_file, sep='\t', index=False)
    print(f"所有重要性不为0的特征已保存至 {all_features_file}")
    
    # 筛选p_value小于0.05且重要性不为0的特征，按重要性排序
    significant_features_df = importance_df[(importance_df['P_value'] < 0.05) & (importance_df['Importance'] > 0)].copy()
    print(f"筛选出 {len(significant_features_df)} 个p_value < 0.05且重要性不为0的显著特征...")
    
    # 保存p值显著的特征到top_significant_features.txt
    top_significant_file = os.path.join(output_dir, 'top_significant_features.txt')
    significant_features_df.to_csv(top_significant_file, sep='\t', index=False)
    print(f"p值显著的特征已保存至 {top_significant_file}")
    
    # 可视化特征重要性（p_value < 0.05的前100个）
    viz_features_count = min(100, len(significant_features_df))
    significant_features_viz = significant_features_df.head(viz_features_count)
    
    plt.figure(figsize=(12, 20))  # 减小高度使特征更紧凑
    significant_viz_df = significant_features_viz.sort_values('Importance', ascending=True)
    plt.barh(significant_viz_df['Feature'], significant_viz_df['Importance'], height=0.6)  # 设置条形高度，控制间距
    plt.xlabel('Importance')
    plt.ylabel('Features')
    plt.title(f'Top {viz_features_count} Significant Features (P_value < 0.05)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'top_100_significant_feature_importance.png'), dpi=300)
    plt.close()
    
    return all_non_zero_features_df


def main():
    """
    主函数
    """
    # 解析命令行参数
    parser = argparse.ArgumentParser(description='二分类机器学习训练与特征重要性分析')
    parser.add_argument('--feature_matrix', required=True, 
                        help='特征矩阵文件路径（长格式tsv，有标题行）')
    parser.add_argument('--all_samples', required=True, 
                        help='总体样本ID文件路径（tsv，无标题行）')
    parser.add_argument('--positive_samples', required=True, 
                        help='阳性样本ID文件路径（tsv，无标题行）')
    parser.add_argument('--output_dir', default='./output', 
                        help='输出目录（默认：./output）')
    parser.add_argument('--random_state', type=int, default=42, 
                        help='随机种子（默认：42）')
    parser.add_argument('--top_n', type=int, default=300, 
                        help='保留的最重要特征数量（默认：300）')
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"输出将保存到目录: {args.output_dir}")
    
    try:
        # 1. 加载并转换特征矩阵
        features_df = load_and_transform_feature_matrix(args.feature_matrix, args.all_samples)
        
        # 2. 创建分类标签
        y, label_counts = create_labels(features_df, args.positive_samples)
        
        # 3. 处理缺失值
        print("处理特征矩阵中的缺失值...")
        # 填充缺失值（使用中位数）
        X = features_df.fillna(features_df.median())
        
        # 4. 训练模型
        model, X_train, X_test, y_train, y_test = train_model(X, y, random_state=args.random_state)
        
        # 5. 评估模型性能
        metrics = evaluate_model(model, X_test, y_test, args.output_dir)
        
        # 6. 获取特征重要性
        feature_names = X.columns.tolist()
        all_features_df = get_feature_importance(model, X, y, feature_names, args.output_dir, top_n=args.top_n)
        
        # 7. 保存模型
        model_file = os.path.join(args.output_dir, 'trained_model.joblib')
        joblib.dump(model, model_file)
        print(f"训练好的模型已保存至 {model_file}")
        
        # 8. 保存数据处理后的特征矩阵和标签，供后续分析使用
        data_file = os.path.join(args.output_dir, 'processed_data.joblib')
        joblib.dump({
            'features': X,
            'labels': y,
            'feature_names': feature_names,
            'X_train': X_train,
            'X_test': X_test,
            'y_train': y_train,
            'y_test': y_test,
            'metrics': metrics
        }, data_file)
        print(f"处理后的数据已保存至 {data_file}")
        
        print("\n分析完成！")
        print(f"1. 模型保存在: {model_file}")
        print(f"2. 所有重要性不为0的特征（包含p值和分布统计）保存在: {os.path.join(args.output_dir, 'all_features.txt')}")
        print(f"3. p值显著的特征（p_value < 0.05，按重要性排序）保存在: {os.path.join(args.output_dir, 'top_significant_features.txt')}")
        print(f"4. 评估指标和可视化结果保存在: {args.output_dir}")
        print(f"   - 混淆矩阵: confusion_matrix.png")
        print(f"   - ROC曲线: roc_curve.png")
        print(f"   - 特征重要性可视化: top_100_significant_feature_importance.png")
    except Exception as e:
        print(f"分析过程中出错: {str(e)}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == "__main__":
    main()