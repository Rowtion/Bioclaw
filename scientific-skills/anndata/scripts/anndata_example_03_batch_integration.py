#!/usr/bin/env python3
"""
anndata_example_03_batch_integration.py
批量数据整合示例

本脚本演示:
1. 加载多个批次的数据
2. 批次效应校正
3. 跨批次整合分析
"""

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from pathlib import Path

sc.settings.verbosity = 2


def create_batch_data(n_cells_per_batch=200, n_genes=1500, n_batches=3):
    """
    创建模拟的多批次单细胞数据
    每个批次有轻微的表达差异 (模拟批次效应)
    """
    adatas = []
    
    cell_types = ['Type_A', 'Type_B', 'Type_C']
    
    for batch_idx in range(n_batches):
        np.random.seed(batch_idx)
        
        n_cells = n_cells_per_batch
        
        # 基础表达矩阵
        X = np.random.poisson(lam=2, size=(n_cells, n_genes)).astype(np.float32)
        
        # 添加批次效应 (系统偏移)
        batch_effect = np.random.normal(0, 0.5, n_genes)
        X = X + batch_effect * 0.3  # 批次效应强度
        X = np.clip(X, 0, None)  # 确保非负
        
        # 创建细胞类型 (各批次略有不同比例)
        type_probs = [0.4 - batch_idx*0.05, 0.35, 0.25 + batch_idx*0.05]
        cell_labels = np.random.choice(cell_types, n_cells, p=type_probs)
        
        obs = pd.DataFrame({
            'cell_type': cell_labels,
            'batch': f'Batch_{batch_idx+1}',
            'sample_id': np.random.choice([f'S{i}' for i in range(1, 4)], n_cells),
        }, index=[f'B{batch_idx+1}_C{i:04d}' for i in range(n_cells)])
        
        var = pd.DataFrame({
            'gene_name': [f'Gene_{i}' for i in range(n_genes)],
        }, index=[f'ENSG{i:06d}' for i in range(n_genes)])
        
        adata = ad.AnnData(X=X, obs=obs, var=var)
        
        # 添加观测矩阵
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
        adata.obs['total_counts'] = adata.X.sum(axis=1)
        
        adatas.append(adata)
        print(f"创建 {adata.obs['batch'][0]}: {adata.n_obs} 细胞")
    
    return adatas


def merge_batches(adatas):
    """合并多个批次的数据"""
    print("\n" + "="*60)
    print("合并批次数据")
    print("="*60)
    
    # 使用 concat 合并
    adata_merged = ad.concat(
        adatas,
        axis=0,
        join='inner',
        label='batch',
        index_unique='-'
    )
    
    print(f"\n合并后: {adata_merged.n_obs} 细胞 × {adata_merged.n_vars} 基因")
    print("批次分布:")
    print(adata_merged.obs['batch'].value_counts())
    
    return adata_merged


def preprocess_for_integration(adata):
    """为整合预处理"""
    print("\n" + "="*60)
    print("预处理")
    print("="*60)
    
    # 基础过滤
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # 标准化和对数转换
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 保存原始数据
    adata.raw = adata.copy()
    
    # 高变基因
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=1000,
        batch_key='batch'  # 考虑批次信息
    )
    
    print(f"高变基因: {adata.var['highly_variable'].sum()}")
    
    return adata


def analyze_before_correction(adata):
    """分析校正前的批次效应"""
    print("\n" + "="*60)
    print("校正前分析 (显示批次效应)")
    print("="*60)
    
    adata_temp = adata.copy()
    
    sc.pp.scale(adata_temp)
    sc.pp.pca(adata_temp, n_comps=30)
    sc.pp.neighbors(adata_temp)
    sc.tl.umap(adata_temp)
    
    # 保存结果用于对比
    adata.obsm['X_umap_before'] = adata_temp.obsm['X_umap'].copy()
    adata.obsm['X_pca_before'] = adata_temp.obsm['X_pca'].copy()
    
    print("✓ 已计算校正前的 UMAP")
    
    return adata


def batch_correction_combat(adata):
    """使用 Combat 进行批次校正"""
    print("\n" + "="*60)
    print("批次校正 (Combat)")
    print("="*60)
    
    # 使用 scanpy 内置的 Combat
    sc.pp.combat(adata, key='batch')
    
    print("✓ Combat 校正完成")
    
    return adata


def analyze_after_correction(adata):
    """分析校正后的数据"""
    print("\n" + "="*60)
    print("校正后分析")
    print("="*60)
    
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=30)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    
    print(f"✓ 发现 {adata.obs['leiden'].nunique()} 个聚类")
    
    return adata


def compare_batch_effects(adata, output_dir='./output'):
    """对比批次校正效果"""
    print("\n" + "="*60)
    print("对比校正效果")
    print("="*60)
    
    import matplotlib.pyplot as plt
    
    Path(output_dir).mkdir(exist_ok=True)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 校正前 - 按批次着色
    sc.pl.embedding(
        adata, basis='umap_before', color='batch', ax=axes[0, 0],
        show=False, title='Before Correction (by Batch)'
    )
    
    # 校正前 - 按细胞类型着色
    sc.pl.embedding(
        adata, basis='umap_before', color='cell_type', ax=axes[0, 1],
        show=False, title='Before Correction (by Cell Type)'
    )
    
    # 校正后 - 按批次着色
    sc.pl.umap(adata, color='batch', ax=axes[1, 0],
               show=False, title='After Correction (by Batch)')
    
    # 校正后 - 按细胞类型着色
    sc.pl.umap(adata, color='cell_type', ax=axes[1, 1],
               show=False, title='After Correction (by Cell Type)')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'batch_correction_comparison.png',
                dpi=150, bbox_inches='tight')
    print(f"✓ 保存对比图: {output_dir}/batch_correction_comparison.png")
    plt.close()


def analyze_batch_mixing(adata):
    """分析批次混合程度"""
    print("\n" + "="*60)
    print("批次混合分析")
    print("="*60)
    
    # 计算每个聚类中的批次分布
    batch_mixing = pd.crosstab(
        adata.obs['leiden'],
        adata.obs['batch'],
        normalize='index'
    )
    
    print("\n各聚类的批次分布 (比例):")
    print(batch_mixing.round(3))
    
    # 计算混合分数 (熵)
    from scipy.stats import entropy
    
    mixing_scores = []
    for cluster in batch_mixing.index:
        probs = batch_mixing.loc[cluster].values
        mix_score = entropy(probs) / np.log(len(probs))  # 归一化熵
        mixing_scores.append(mix_score)
    
    batch_mixing['mixing_score'] = mixing_scores
    
    print("\n混合分数 (1.0 = 完全混合):")
    for cluster, score in zip(batch_mixing.index, mixing_scores):
        status = "✓" if score > 0.8 else "⚠"
        print(f"  聚类 {cluster}: {score:.3f} {status}")
    
    return batch_mixing


def cross_batch_consistency(adata):
    """评估跨批次一致性"""
    print("\n" + "="*60)
    print("跨批次一致性分析")
    print("="*60)
    
    # 计算细胞类型在不同批次中的比例一致性
    consistency = pd.crosstab(
        adata.obs['cell_type'],
        adata.obs['batch'],
        normalize='columns'
    )
    
    print("\n各批次中的细胞类型比例:")
    print(consistency.round(3))
    
    # 计算变异系数
    cv = consistency.std(axis=1) / consistency.mean(axis=1)
    print(f"\n细胞类型比例变异系数 (越低越一致):")
    for cell_type, coeff in cv.items():
        status = "✓" if coeff < 0.3 else "⚠"
        print(f"  {cell_type}: {coeff:.3f} {status}")


def save_integration_results(adata, output_dir='./output'):
    """保存整合结果"""
    print("\n" + "="*60)
    print("保存结果")
    print("="*60)
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # 保存整合后的数据
    output_path = Path(output_dir) / 'integrated_data.h5ad'
    adata.write_h5ad(output_path)
    print(f"✓ 保存: {output_path}")
    
    # 导出整合后的元数据
    meta_path = Path(output_dir) / 'integration_metadata.csv'
    adata.obs.to_csv(meta_path)
    print(f"✓ 导出元数据: {meta_path}")
    
    # 导出批次信息
    batch_stats = adata.obs.groupby(['batch', 'cell_type']).size().unstack(fill_value=0)
    batch_path = Path(output_dir) / 'batch_celltype_counts.csv'
    batch_stats.to_csv(batch_path)
    print(f"✓ 导出批次统计: {batch_path}")


def main():
    """主流程"""
    print("AnnData 批次整合示例")
    print("="*60)
    
    # 1. 创建模拟的多批次数据
    adatas = create_batch_data(
        n_cells_per_batch=200,
        n_genes=1500,
        n_batches=3
    )
    
    # 2. 合并数据
    adata = merge_batches(adatas)
    
    # 3. 预处理
    adata = preprocess_for_integration(adata)
    
    # 4. 校正前分析
    adata = analyze_before_correction(adata)
    
    # 5. 批次校正
    adata = batch_correction_combat(adata)
    
    # 6. 校正后分析
    adata = analyze_after_correction(adata)
    
    # 7. 对比可视化
    compare_batch_effects(adata)
    
    # 8. 混合分析
    analyze_batch_mixing(adata)
    
    # 9. 一致性分析
    cross_batch_consistency(adata)
    
    # 10. 保存结果
    save_integration_results(adata)
    
    print("\n" + "="*60)
    print("批次整合完成!")
    print("="*60)


if __name__ == '__main__':
    main()
