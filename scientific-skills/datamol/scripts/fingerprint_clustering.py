#!/usr/bin/env python3
"""
Datamol 分子指纹相似性搜索与聚类脚本
功能：
1. 生成分子指纹（ECFP, MACCS等）
2. 计算分子间相似性
3. 基于相似性进行Butina聚类
4. 选择多样化分子子集

依赖：datamol, numpy, pandas
"""

import datamol as dm
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional
from scipy.spatial.distance import squareform


def generate_fingerprints(
    mols: List[dm.Mol],
    fp_type: str = "ecfp",
    radius: int = 2,
    n_bits: int = 2048
) -> np.ndarray:
    """
    批量生成分子指纹
    
    Args:
        mols: 分子对象列表
        fp_type: 指纹类型 ('ecfp', 'maccs', 'topological', 'atompair')
        radius: ECFP指纹半径
        n_bits: 指纹位数
    
    Returns:
        指纹数组 (n_mols, n_bits)
    """
    fingerprints = []
    for mol in mols:
        if mol is not None:
            fp = dm.to_fp(
                mol,
                fp_type=fp_type,
                radius=radius,
                n_bits=n_bits
            )
            fingerprints.append(fp)
        else:
            fingerprints.append(np.zeros(n_bits))
    
    return np.array(fingerprints)


def calculate_similarity_matrix(
    mols: List[dm.Mol],
    n_jobs: int = -1
) -> np.ndarray:
    """
    计算分子间Tanimoto距离矩阵
    
    Args:
        mols: 分子对象列表
        n_jobs: 并行处理核心数
    
    Returns:
        距离矩阵 (n_mols, n_mols)
    """
    # 使用datamol的并行距离计算
    condensed_dist = dm.pdist(mols, n_jobs=n_jobs)
    dist_matrix = squareform(condensed_dist)
    return dist_matrix


def cluster_molecules(
    mols: List[dm.Mol],
    cutoff: float = 0.35,
    n_jobs: int = -1
) -> Tuple[List[List[int]], int]:
    """
    Butina聚类分析
    
    Args:
        mols: 分子对象列表
        cutoff: Tanimoto距离阈值（越小越严格）
        n_jobs: 并行处理核心数
    
    Returns:
        聚类列表和聚类数量
    """
    print(f"执行Butina聚类 (cutoff={cutoff})...")
    clusters = dm.cluster_mols(mols, cutoff=cutoff, n_jobs=n_jobs)
    
    # 统计聚类大小
    sizes = [len(c) for c in clusters]
    print(f"  聚类数量: {len(clusters)}")
    print(f"  平均聚类大小: {np.mean(sizes):.1f}")
    print(f"  最大聚类: {max(sizes)}")
    print(f"  单例聚类: {sum(1 for s in sizes if s == 1)}")
    
    return clusters, len(clusters)


def select_diverse_subset(
    mols: List[dm.Mol],
    n_pick: int,
    method: str = "diverse"
) -> List[dm.Mol]:
    """
    选择多样化分子子集
    
    Args:
        mols: 分子对象列表
        n_pick: 选择的分子数量
        method: 选择方法 ('diverse'或'centroids')
    
    Returns:
        选择的分子列表
    """
    if method == "diverse":
        return dm.pick_diverse(mols, npick=n_pick)
    elif method == "centroids":
        return dm.pick_centroids(mols, npick=n_pick)
    else:
        raise ValueError(f"未知方法: {method}")


def find_similar_molecules(
    query_mol: dm.Mol,
    library_mols: List[dm.Mol],
    top_n: int = 10,
    n_jobs: int = -1
) -> List[Tuple[int, float]]:
    """
    在分子库中搜索相似分子
    
    Args:
        query_mol: 查询分子
        library_mols: 分子库
        top_n: 返回前N个相似分子
        n_jobs: 并行处理核心数
    
    Returns:
        (索引, 相似度)列表
    """
    distances = dm.cdist([query_mol], library_mols, n_jobs=n_jobs)
    similarities = 1 - distances[0]  # 转换距离为相似度
    
    # 获取前N个
    top_indices = np.argsort(similarities)[::-1][:top_n]
    results = [(int(idx), float(similarities[idx])) for idx in top_indices]
    
    return results


def main(
    input_file: str,
    output_dir: str = "clustering_results",
    fp_type: str = "ecfp",
    cluster_cutoff: float = 0.35,
    n_diverse: int = 100,
    n_jobs: int = -1
):
    """
    主函数：指纹相似性分析与聚类
    
    Args:
        input_file: 输入分子文件
        output_dir: 输出目录
        fp_type: 指纹类型
        cluster_cutoff: 聚类距离阈值
        n_diverse: 选择的多样化分子数
        n_jobs: 并行处理核心数
    """
    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"=" * 60)
    print(f"分子指纹相似性分析与聚类")
    print(f"输入文件: {input_file}")
    print(f"=" * 60)
    
    # 1. 加载分子
    print("\n[1/5] 加载分子...")
    mols = []
    path = Path(input_file)
    
    if path.suffix == ".sdf":
        df = dm.read_sdf(str(path), mol_column="mol")
    elif path.suffix == ".csv":
        df = dm.read_csv(str(path), smiles_column="SMILES", mol_column="mol")
    else:
        df = dm.open_df(str(path), mol_column="mol")
    
    mols = df["mol"].dropna().tolist()
    print(f"  加载了 {len(mols)} 个分子")
    
    # 2. 生成分子指纹
    print(f"\n[2/5] 生成{fp_type.upper()}指纹...")
    fps = generate_fingerprints(mols, fp_type=fp_type)
    print(f"  指纹维度: {fps.shape}")
    
    # 3. 分子聚类
    print(f"\n[3/5] 执行分子聚类...")
    clusters, n_clusters = cluster_molecules(
        mols,
        cutoff=cluster_cutoff,
        n_jobs=n_jobs
    )
    
    # 保存聚类结果
    cluster_data = []
    for cluster_id, cluster_indices in enumerate(clusters):
        for mol_idx in cluster_indices:
            smiles = dm.to_smiles(mols[mol_idx])
            cluster_data.append({
                "mol_index": mol_idx,
                "cluster_id": cluster_id,
                "cluster_size": len(cluster_indices),
                "is_centroid": mol_idx == cluster_indices[0],
                "SMILES": smiles
            })
    
    cluster_df = pd.DataFrame(cluster_data)
    cluster_file = Path(output_dir) / "clusters.csv"
    cluster_df.to_csv(cluster_file, index=False)
    print(f"  聚类结果已保存: {cluster_file}")
    
    # 4. 选择多样化子集
    print(f"\n[4/5] 选择{n_diverse}个多样化分子...")
    if n_diverse < len(mols):
        diverse_mols = select_diverse_subset(mols, n_diverse, method="diverse")
        diverse_smiles = [dm.to_smiles(m) for m in diverse_mols]
        
        diverse_df = pd.DataFrame({
            "SMILES": diverse_smiles,
            "selection_method": "maxmin_diversity"
        })
        diverse_file = Path(output_dir) / "diverse_subset.csv"
        diverse_df.to_csv(diverse_file, index=False)
        print(f"  多样化子集已保存: {diverse_file}")
    else:
        print(f"  警告: 请求的多样化分子数({n_diverse})超过总分子数({len(mols)})")
    
    # 5. 相似性分析（仅分析前10个分子作为示例）
    print(f"\n[5/5] 相似性搜索示例 (前10个分子)...")
    query_results = []
    n_queries = min(10, len(mols))
    
    for i in range(n_queries):
        similar = find_similar_molecules(mols[i], mols, top_n=5, n_jobs=n_jobs)
        for rank, (idx, sim) in enumerate(similar, 1):
            query_results.append({
                "query_index": i,
                "query_smiles": dm.to_smiles(mols[i]),
                "rank": rank,
                "match_index": idx,
                "match_smiles": dm.to_smiles(mols[idx]),
                "similarity": sim
            })
    
    sim_df = pd.DataFrame(query_results)
    sim_file = Path(output_dir) / "similarity_search.csv"
    sim_df.to_csv(sim_file, index=False)
    print(f"  相似性结果已保存: {sim_file}")
    
    # 可视化（如果matplotlib可用）
    try:
        import matplotlib.pyplot as plt
        
        # 聚类大小分布
        cluster_sizes = [len(c) for c in clusters]
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.hist(cluster_sizes, bins=30, edgecolor='black')
        plt.xlabel("Cluster Size")
        plt.ylabel("Count")
        plt.title(f"Cluster Size Distribution (n={n_clusters})")
        
        # 相似性分布
        plt.subplot(1, 2, 2)
        sample_distances = dm.pdist(mols[:100], n_jobs=n_jobs)  # 采样前100个
        plt.hist(1 - sample_distances, bins=30, edgecolor='black')
        plt.xlabel("Tanimoto Similarity")
        plt.ylabel("Count")
        plt.title("Pairwise Similarity Distribution")
        
        plt.tight_layout()
        viz_file = Path(output_dir) / "clustering_analysis.png"
        plt.savefig(viz_file, dpi=300, bbox_inches='tight')
        print(f"  可视化图表已保存: {viz_file}")
    except ImportError:
        print("  跳过可视化 (matplotlib未安装)")
    
    print("\n" + "=" * 60)
    print("分析完成!")
    print("=" * 60)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="分子指纹相似性分析与聚类"
    )
    parser.add_argument(
        "input_file",
        help="输入分子文件 (.sdf, .csv)"
    )
    parser.add_argument(
        "-o", "--output",
        default="clustering_results",
        help="输出目录"
    )
    parser.add_argument(
        "--fp-type",
        choices=["ecfp", "maccs", "topological", "atompair"],
        default="ecfp",
        help="指纹类型"
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.35,
        help="聚类距离阈值 (默认0.35)"
    )
    parser.add_argument(
        "--n-diverse",
        type=int,
        default=100,
        help="选择的多样化分子数"
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="并行处理核心数 (-1=所有核心)"
    )
    
    args = parser.parse_args()
    
    main(
        input_file=args.input_file,
        output_dir=args.output,
        fp_type=args.fp_type,
        cluster_cutoff=args.cutoff,
        n_diverse=args.n_diverse,
        n_jobs=args.n_jobs
    )
