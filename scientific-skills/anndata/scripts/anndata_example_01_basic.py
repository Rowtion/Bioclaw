#!/usr/bin/env python3
"""
anndata_example_01_basic.py
单细胞数据分析基础示例 - 创建和操作 AnnData 对象

本脚本演示:
1. 创建 AnnData 对象
2. 添加观察值和变量元数据
3. 保存和加载 h5ad 文件
4. 基本数据操作
"""

import numpy as np
import pandas as pd
import anndata as ad
from pathlib import Path


def create_sample_anndata(n_obs=100, n_vars=2000, seed=42):
    """
    创建示例 AnnData 对象
    
    Args:
        n_obs: 细胞数量 (observations)
        n_vars: 基因数量 (variables)
        seed: 随机种子
    
    Returns:
        AnnData 对象
    """
    np.random.seed(seed)
    
    # 创建表达矩阵 (稀疏数据模拟)
    X = np.random.poisson(lam=2, size=(n_obs, n_vars)).astype(np.float32)
    
    # 创建细胞元数据 (obs)
    cell_types = ['T_cell', 'B_cell', 'Monocyte', 'NK_cell']
    obs = pd.DataFrame({
        'cell_type': np.random.choice(cell_types, n_obs),
        'sample': np.random.choice(['sample_A', 'sample_B', 'sample_C'], n_obs),
        'n_genes': (X > 0).sum(axis=1),
        'total_counts': X.sum(axis=1),
        'percent_mito': np.random.uniform(0, 20, n_obs),
    }, index=[f'cell_{i:04d}' for i in range(n_obs)])
    
    # 创建基因元数据 (var)
    var = pd.DataFrame({
        'gene_name': [f'Gene_{i}' for i in range(n_vars)],
        'highly_variable': np.random.choice([True, False], n_vars, p=[0.1, 0.9]),
        'gene_length': np.random.randint(500, 10000, n_vars),
    }, index=[f'ENSG{i:08d}' for i in range(n_vars)])
    
    # 创建 AnnData 对象
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    # 添加非结构化数据
    adata.uns['study'] = 'Example Study'
    adata.uns['date'] = '2024-01-15'
    
    # 添加观测矩阵 (obsm) - 降维结果示例
    adata.obsm['X_pca'] = np.random.randn(n_obs, 50)
    adata.obsm['X_umap'] = np.random.randn(n_obs, 2)
    
    print(f"✓ 创建 AnnData 对象: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    print(f"  - 内存使用: {adata.n_obs * adata.n_vars * 4 / 1024 / 1024:.2f} MB")
    
    return adata


def demonstrate_basic_operations(adata):
    """演示 AnnData 基本操作"""
    print("\n" + "="*60)
    print("基本操作演示")
    print("="*60)
    
    # 1. 查看基本信息
    print("\n1. 基本信息:")
    print(f"   Shape: {adata.shape}")
    print(f"   观测值数量: {adata.n_obs}")
    print(f"   变量数量: {adata.n_vars}")
    
    # 2. 访问 obs 和 var
    print("\n2. 细胞类型分布:")
    print(adata.obs['cell_type'].value_counts())
    
    # 3. 子集选择
    print("\n3. 子集选择 (仅 T 细胞):")
    t_cells = adata[adata.obs['cell_type'] == 'T_cell']
    print(f"   选择 {t_cells.n_obs} 个 T 细胞")
    
    # 4. 选择高变基因
    print("\n4. 高变基因:")
    hv_genes = adata[:, adata.var['highly_variable']]
    print(f"   选择 {hv_genes.n_vars} 个高变基因")
    
    # 5. 添加新列
    print("\n5. 添加新元数据列:")
    adata.obs['quality'] = np.random.choice(['good', 'medium', 'poor'], adata.n_obs)
    print(f"   已添加 'quality' 列")
    print(adata.obs['quality'].value_counts())


def demonstrate_concatenation():
    """演示数据集合并"""
    print("\n" + "="*60)
    print("数据集合并演示")
    print("="*60)
    
    # 创建多个批次的数据
    adata1 = create_sample_anndata(n_obs=50, seed=1)
    adata1.obs['batch'] = 'batch1'
    
    adata2 = create_sample_anndata(n_obs=50, seed=2)
    adata2.obs['batch'] = 'batch2'
    
    adata3 = create_sample_anndata(n_obs=50, seed=3)
    adata3.obs['batch'] = 'batch3'
    
    # 合并数据
    adata_merged = ad.concat(
        [adata1, adata2, adata3],
        axis=0,  # 沿观测值方向合并
        join='inner',  # 只保留共同基因
        label='batch',
        keys=['batch1', 'batch2', 'batch3']
    )
    
    print(f"\n合并后的数据: {adata_merged.n_obs} 细胞 × {adata_merged.n_vars} 基因")
    print("批次分布:")
    print(adata_merged.obs['batch'].value_counts())
    
    return adata_merged


def save_and_load(adata, output_dir='./output'):
    """保存和加载 AnnData 对象"""
    print("\n" + "="*60)
    print("文件 I/O 演示")
    print("="*60)
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # 保存为 h5ad
    h5ad_path = Path(output_dir) / 'example_data.h5ad'
    adata.write_h5ad(h5ad_path)
    print(f"✓ 已保存: {h5ad_path}")
    print(f"  文件大小: {h5ad_path.stat().st_size / 1024:.2f} KB")
    
    # 使用压缩保存
    h5ad_compressed = Path(output_dir) / 'example_data_compressed.h5ad'
    adata.write_h5ad(h5ad_compressed, compression='gzip')
    print(f"✓ 已保存 (gzip 压缩): {h5ad_compressed}")
    print(f"  文件大小: {h5ad_compressed.stat().st_size / 1024:.2f} KB")
    
    # 加载数据
    print("\n加载数据...")
    adata_loaded = ad.read_h5ad(h5ad_path)
    print(f"✓ 成功加载: {adata_loaded.n_obs} 细胞 × {adata_loaded.n_vars} 基因")
    
    # 验证数据一致性
    assert np.allclose(adata.X.sum(), adata_loaded.X.sum())
    print("✓ 数据验证通过")
    
    return adata_loaded


def demonstrate_backed_mode(h5ad_path):
    """演示 backed 模式处理大文件"""
    print("\n" + "="*60)
    print("Backed 模式演示 (内存映射)")
    print("="*60)
    
    # 以 backed 模式打开 (只读)
    adata_backed = ad.read_h5ad(h5ad_path, backed='r')
    print(f"✓ 以 backed 模式打开: {adata_backed.n_obs} 细胞 × {adata_backed.n_vars} 基因")
    print(f"  注意: 数据仍保存在磁盘，不占用内存")
    
    # 基于元数据筛选 (不加载数据)
    print("\n筛选高表达细胞...")
    mask = adata_backed.obs['n_genes'] > 1000
    selected_cells = adata_backed[mask]
    print(f"  筛选出 {selected_cells.n_obs} 个细胞")
    
    # 将筛选后的数据加载到内存
    adata_in_memory = selected_cells.to_memory()
    print(f"✓ 已加载到内存: {adata_in_memory.n_obs} 细胞")
    
    adata_backed.file.close()  # 关闭 backed 文件
    return adata_in_memory


def main():
    """主函数"""
    print("AnnData 基础示例")
    print("="*60)
    
    # 1. 创建示例数据
    adata = create_sample_anndata(n_obs=100, n_vars=2000)
    
    # 2. 基本操作
    demonstrate_basic_operations(adata)
    
    # 3. 合并数据
    adata_merged = demonstrate_concatenation()
    
    # 4. 保存和加载
    adata_loaded = save_and_load(adata)
    
    # 5. Backed 模式
    output_path = Path('./output/example_data.h5ad')
    if output_path.exists():
        demonstrate_backed_mode(output_path)
    
    print("\n" + "="*60)
    print("示例完成!")
    print("="*60)


if __name__ == '__main__':
    main()
