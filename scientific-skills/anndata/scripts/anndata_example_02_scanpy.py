#!/usr/bin/env python3
"""
anndata_example_02_scanpy.py
使用 AnnData 与 Scanpy 进行单细胞分析

本脚本演示:
1. 使用 scanpy 进行预处理
2. 降维和聚类
3. 可视化
4. 差异表达分析
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from pathlib import Path

# 设置 scanpy 参数
sc.settings.verbosity = 3  # 输出信息级别
sc.settings.set_figure_params(dpi=80, facecolor='white')


def create_realistic_data(n_obs=500, n_vars=3000):
    """创建模拟的真实单细胞数据"""
    np.random.seed(42)
    
    # 模拟稀疏表达矩阵 (UMI counts)
    X = np.random.poisson(lam=1.5, size=(n_obs, n_vars))
    X = X.astype(np.float32)
    
    # 创建细胞类型标签 (有生物学意义的分组)
    cell_types = ['CD4_T', 'CD8_T', 'B_cell', 'Monocyte', 'NK']
    cell_type_probs = [0.3, 0.2, 0.25, 0.15, 0.1]
    cell_labels = np.random.choice(cell_types, n_obs, p=cell_type_probs)
    
    obs = pd.DataFrame({
        'cell_type': cell_labels,
        'sample': np.random.choice(['S1', 'S2', 'S3'], n_obs),
    }, index=[f'cell_{i:04d}' for i in range(n_obs)])
    
    # 基因信息
    var = pd.DataFrame({
        'gene_name': [f'Gene_{i}' for i in range(n_vars)],
    }, index=[f'g{i}' for i in range(n_vars)])
    
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # 计算质量控制指标
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.obs['total_counts'] = adata.X.sum(axis=1)
    
    return adata


def quality_control(adata):
    """质量控制"""
    print("\n" + "="*60)
    print("质量控制")
    print("="*60)
    
    print("\n原始数据:", adata.shape)
    
    # 过滤低质量细胞
    print("\n过滤低质量细胞...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=500)
    sc.pp.filter_genes(adata, min_cells=3)
    
    print("过滤后:", adata.shape)
    
    # 计算线粒体基因比例
    # 这里模拟线粒体基因 (通常以 MT- 开头)
    mito_genes = adata.var_names.str.startswith('g0')  # 模拟
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1
    ) / np.sum(adata.X, axis=1) * 100
    
    # 过滤高线粒体比例的细胞
    adata = adata[adata.obs['percent_mito'] < 20, :]
    
    print(f"过滤线粒体后: {adata.shape}")
    
    return adata


def preprocess(adata):
    """数据预处理"""
    print("\n" + "="*60)
    print("数据预处理")
    print("="*60)
    
    # 标准化
    print("\n1. 标准化 (每细胞总计数 1e4)")
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # 对数转换
    print("2. 对数转换 (log1p)")
    sc.pp.log1p(adata)
    
    # 识别高变基因
    print("3. 识别高变基因")
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    
    print(f"高变基因数量: {adata.n_vars}")
    
    return adata


def dimensionality_reduction(adata):
    """降维"""
    print("\n" + "="*60)
    print("降维分析")
    print("="*60)
    
    # PCA
    print("\n1. PCA 降维")
    sc.pp.scale(adata)  # 标准化到零均值单位方差
    sc.pp.pca(adata, n_comps=50)
    
    print(f"PCA 形状: {adata.obsm['X_pca'].shape}")
    
    # 计算邻居图
    print("2. 计算邻居图")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    
    # UMAP
    print("3. UMAP 降维")
    sc.tl.umap(adata)
    
    # t-SNE (可选)
    print("4. t-SNE 降维")
    sc.tl.tsne(adata)
    
    return adata


def clustering(adata):
    """聚类分析"""
    print("\n" + "="*60)
    print("聚类分析")
    print("="*60)
    
    # Leiden 聚类
    print("\nLeiden 聚类...")
    sc.tl.leiden(adata, resolution=0.8)
    
    print(f"发现 {adata.obs['leiden'].nunique()} 个聚类")
    print("聚类分布:")
    print(adata.obs['leiden'].value_counts().sort_index())
    
    return adata


def differential_expression(adata):
    """差异表达分析"""
    print("\n" + "="*60)
    print("差异表达分析")
    print("="*60)
    
    # 使用 Wilcoxon 秩和检验
    print("\n计算差异表达基因...")
    sc.tl.rank_genes_groups(
        adata,
        groupby='leiden',
        method='wilcoxon',
        use_raw=False
    )
    
    # 显示每个聚类的 top 基因
    print("\n各聚类 top 基因:")
    for cluster in adata.obs['leiden'].cat.categories[:3]:  # 只显示前3个
        genes = adata.uns['rank_genes_groups']['names'][cluster][:5]
        scores = adata.uns['rank_genes_groups']['scores'][cluster][:5]
        print(f"\n聚类 {cluster}:")
        for gene, score in zip(genes, scores):
            print(f"  {gene}: {score:.2f}")
    
    return adata


def visualize(adata, output_dir='./output'):
    """可视化"""
    print("\n" + "="*60)
    print("生成可视化")
    print("="*60)
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # 1. UMAP 图 - 按聚类着色
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    sc.pl.umap(adata, color='leiden', ax=axes[0], show=False, 
               title='Leiden Clusters', legend_loc='on data')
    sc.pl.umap(adata, color='cell_type', ax=axes[1], show=False,
               title='True Cell Types')
    
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'umap_clusters.png', dpi=150, bbox_inches='tight')
    print(f"✓ 保存: {output_dir}/umap_clusters.png")
    plt.close()
    
    # 2. t-SNE 图
    fig, ax = plt.subplots(figsize=(6, 5))
    sc.pl.tsne(adata, color='leiden', ax=ax, show=False,
               title='t-SNE Clusters')
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'tsne_clusters.png', dpi=150, bbox_inches='tight')
    print(f"✓ 保存: {output_dir}/tsne_clusters.png")
    plt.close()
    
    # 3. 差异表达基因热图
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.rank_genes_groups_heatmap(
        adata,
        n_genes=5,
        groupby='leiden',
        show_gene_labels=True,
        show=False,
        ax=ax
    )
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'de_heatmap.png', dpi=150, bbox_inches='tight')
    print(f"✓ 保存: {output_dir}/de_heatmap.png")
    plt.close()
    
    # 4. 小提琴图 - QC 指标
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    sc.pl.violin(adata, 'n_genes', groupby='leiden', ax=axes[0], show=False)
    sc.pl.violin(adata, 'total_counts', groupby='leiden', ax=axes[1], show=False)
    plt.tight_layout()
    plt.savefig(Path(output_dir) / 'qc_violin.png', dpi=150, bbox_inches='tight')
    print(f"✓ 保存: {output_dir}/qc_violin.png")
    plt.close()


def save_results(adata, output_dir='./output'):
    """保存分析结果"""
    print("\n" + "="*60)
    print("保存结果")
    print("="*60)
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # 保存 AnnData 对象
    output_path = Path(output_dir) / 'scanpy_analysis_result.h5ad'
    adata.write_h5ad(output_path)
    print(f"✓ 保存 AnnData: {output_path}")
    
    # 导出细胞元数据
    obs_path = Path(output_dir) / 'cell_metadata.csv'
    adata.obs.to_csv(obs_path)
    print(f"✓ 导出细胞元数据: {obs_path}")
    
    # 导出聚类结果
    cluster_path = Path(output_dir) / 'clustering_results.csv'
    adata.obs[['cell_type', 'leiden', 'sample']].to_csv(cluster_path)
    print(f"✓ 导出聚类结果: {cluster_path}")
    
    # 导出差异表达结果
    de_results = []
    for cluster in adata.obs['leiden'].cat.categories:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df['cluster'] = cluster
        de_results.append(df.head(10))  # 每个聚类 top 10
    
    de_df = pd.concat(de_results, ignore_index=True)
    de_path = Path(output_dir) / 'differential_expression.csv'
    de_df.to_csv(de_path, index=False)
    print(f"✓ 导出差异表达结果: {de_path}")


def main():
    """主分析流程"""
    print("Scanpy 单细胞分析流程")
    print("="*60)
    
    # 1. 创建数据
    adata = create_realistic_data(n_obs=500, n_vars=3000)
    
    # 2. 质量控制
    adata = quality_control(adata)
    
    # 3. 预处理
    adata = preprocess(adata)
    
    # 4. 降维
    adata = dimensionality_reduction(adata)
    
    # 5. 聚类
    adata = clustering(adata)
    
    # 6. 差异表达
    adata = differential_expression(adata)
    
    # 7. 可视化
    visualize(adata)
    
    # 8. 保存结果
    save_results(adata)
    
    print("\n" + "="*60)
    print("分析完成!")
    print("="*60)
    print("\n输出文件:")
    print("  - output/scanpy_analysis_result.h5ad")
    print("  - output/umap_clusters.png")
    print("  - output/tsne_clusters.png")
    print("  - output/de_heatmap.png")
    print("  - output/qc_violin.png")
    print("  - output/*.csv")


if __name__ == '__main__':
    main()
