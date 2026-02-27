#!/usr/bin/env python3
"""
Datamol 分子支架分析与SAR研究脚本
功能：
1. 提取Bemis-Murcko分子支架
2. 按支架对分子进行分组
3. 分析支架-活性关系(SAR)
4. 支架多样性评估

依赖：datamol, pandas, matplotlib
"""

import datamol as dm
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple
from collections import Counter


def extract_scaffolds(mols: List[dm.Mol]) -> Tuple[List[dm.Mol], List[str]]:
    """
    提取分子的Bemis-Murcko支架
    
    Args:
        mols: 分子对象列表
    
    Returns:
        支架对象列表和支架SMILES列表
    """
    scaffolds = []
    scaffold_smiles = []
    
    for mol in mols:
        if mol is not None:
            try:
                scaffold = dm.to_scaffold_murcko(mol)
                if scaffold is not None:
                    scaffolds.append(scaffold)
                    smi = dm.to_smiles(scaffold)
                    scaffold_smiles.append(smi)
                else:
                    scaffolds.append(None)
                    scaffold_smiles.append(None)
            except Exception:
                scaffolds.append(None)
                scaffold_smiles.append(None)
        else:
            scaffolds.append(None)
            scaffold_smiles.append(None)
    
    return scaffolds, scaffold_smiles


def analyze_scaffold_frequency(scaffold_smiles: List[str]) -> pd.DataFrame:
    """
    分析支架频率分布
    
    Args:
        scaffold_smiles: 支架SMILES列表
    
    Returns:
        支架频率统计DataFrame
    """
    # 过滤None值
    valid_smiles = [s for s in scaffold_smiles if s is not None]
    
    # 统计频率
    counter = Counter(valid_smiles)
    
    df = pd.DataFrame([
        {"scaffold_smiles": smi, "count": count}
        for smi, count in counter.most_common()
    ])
    
    df["percentage"] = 100 * df["count"] / len(valid_smiles)
    df["cumulative_percentage"] = df["percentage"].cumsum()
    
    return df


def group_by_scaffold(
    mols: List[dm.Mol],
    activities: List[float],
    scaffold_smiles: List[str]
) -> Dict[str, Dict]:
    """
    按支架对分子进行分组，分析活性
    
    Args:
        mols: 分子对象列表
        activities: 活性值列表
        scaffold_smiles: 支架SMILES列表
    
    Returns:
        支架分组信息字典
    """
    scaffold_groups = {}
    
    for i, (mol, smi, act) in enumerate(zip(mols, scaffold_smiles, activities)):
        if smi is None:
            continue
        
        if smi not in scaffold_groups:
            scaffold_groups[smi] = {
                "molecules": [],
                "activities": [],
                "indices": []
            }
        
        scaffold_groups[smi]["molecules"].append(mol)
        scaffold_groups[smi]["activities"].append(act)
        scaffold_groups[smi]["indices"].append(i)
    
    # 计算统计信息
    for smi, data in scaffold_groups.items():
        acts = np.array(data["activities"])
        data["mean_activity"] = np.mean(acts)
        data["std_activity"] = np.std(acts)
        data["min_activity"] = np.min(acts)
        data["max_activity"] = np.max(acts)
        data["n_compounds"] = len(acts)
    
    return scaffold_groups


def scaffold_train_test_split(
    mols: List[dm.Mol],
    scaffold_smiles: List[str],
    test_ratio: float = 0.2,
    random_seed: int = 42
) -> Tuple[List[int], List[int]]:
    """
    基于支架的训练/测试集划分
    确保训练集和测试集包含不同的支架
    
    Args:
        mols: 分子对象列表
        scaffold_smiles: 支架SMILES列表
        test_ratio: 测试集比例
        random_seed: 随机种子
    
    Returns:
        训练集和测试集的索引列表
    """
    import random
    random.seed(random_seed)
    
    # 按支架分组
    scaffold_to_indices = {}
    for i, smi in enumerate(scaffold_smiles):
        if smi not in scaffold_to_indices:
            scaffold_to_indices[smi] = []
        scaffold_to_indices[smi].append(i)
    
    # 随机划分支架
    scaffolds = list(scaffold_to_indices.keys())
    random.shuffle(scaffolds)
    
    split_idx = int(len(scaffolds) * (1 - test_ratio))
    train_scaffolds = scaffolds[:split_idx]
    test_scaffolds = scaffolds[split_idx:]
    
    # 获取分子索引
    train_indices = []
    for scaf in train_scaffolds:
        train_indices.extend(scaffold_to_indices[scaf])
    
    test_indices = []
    for scaf in test_scaffolds:
        test_indices.extend(scaffold_to_indices[scaf])
    
    return train_indices, test_indices


def visualize_scaffolds(
    scaffold_mols: List[dm.Mol],
    scaffold_smiles: List[str],
    counts: List[int],
    output_file: str,
    top_n: int = 20
):
    """
    可视化最常见的支架
    
    Args:
        scaffold_mols: 支架分子对象列表
        scaffold_smiles: 支架SMILES列表
        counts: 出现次数列表
        output_file: 输出图像文件
        top_n: 显示前N个支架
    """
    # 去重并按频率排序
    unique_data = []
    seen = set()
    for mol, smi, count in zip(scaffold_mols, scaffold_smiles, counts):
        if smi not in seen and mol is not None:
            unique_data.append((mol, smi, count))
            seen.add(smi)
    
    unique_data.sort(key=lambda x: x[2], reverse=True)
    top_data = unique_data[:top_n]
    
    if not top_data:
        print("没有有效的支架用于可视化")
        return
    
    top_mols = [d[0] for d in top_data]
    top_counts = [d[2] for d in top_data]
    legends = [f"Count: {c}" for c in top_counts]
    
    # 生成图像
    img = dm.viz.to_image(
        top_mols,
        legends=legends,
        n_cols=5,
        mol_size=(250, 250)
    )
    
    # 保存
    img.save(output_file)
    print(f"  支架可视化已保存: {output_file}")


def main(
    input_file: str,
    output_dir: str = "scaffold_analysis",
    activity_column: str = "activity",
    smiles_column: str = "SMILES",
    top_n_visualize: int = 20
):
    """
    主函数：支架分析流程
    
    Args:
        input_file: 输入CSV文件（包含SMILES和活性列）
        output_dir: 输出目录
        activity_column: 活性值列名
        smiles_column: SMILES列名
        top_n_visualize: 可视化的支架数量
    """
    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"=" * 60)
    print(f"分子支架分析")
    print(f"输入文件: {input_file}")
    print(f"=" * 60)
    
    # 1. 加载数据
    print("\n[1/6] 加载分子数据...")
    df = pd.read_csv(input_file)
    print(f"  加载了 {len(df)} 行数据")
    
    # 解析SMILES
    smiles_list = df[smiles_column].tolist()
    mols = [dm.to_mol(smi) for smi in smiles_list]
    valid_mols = [m for m in mols if m is not None]
    print(f"  成功解析 {len(valid_mols)}/{len(mols)} 个分子")
    
    # 获取活性值（如果有）
    has_activity = activity_column in df.columns
    if has_activity:
        activities = df[activity_column].tolist()
        print(f"  检测到活性列: {activity_column}")
    
    # 2. 提取支架
    print("\n[2/6] 提取Bemis-Murcko支架...")
    scaffolds, scaffold_smiles = extract_scaffolds(mols)
    valid_scaffolds = [s for s in scaffolds if s is not None]
    print(f"  提取了 {len(valid_scaffolds)} 个支架")
    print(f"  唯一支架数: {len(set(s for s in scaffold_smiles if s is not None))}")
    
    # 3. 支架频率分析
    print("\n[3/6] 分析支架频率...")
    scaffold_freq_df = analyze_scaffold_frequency(scaffold_smiles)
    freq_file = Path(output_dir) / "scaffold_frequency.csv"
    scaffold_freq_df.to_csv(freq_file, index=False)
    print(f"  支架频率统计已保存: {freq_file}")
    
    print(f"\n  前10个最常见支架:")
    for i, row in scaffold_freq_df.head(10).iterrows():
        print(f"    {i+1}. Count={row['count']}, "
              f"Percentage={row['percentage']:.1f}%, "
              f"SMILES={row['scaffold_smiles'][:50]}...")
    
    # 4. SAR分析（如果有活性数据）
    if has_activity:
        print("\n[4/6] 支架-活性关系(SAR)分析...")
        scaffold_groups = group_by_scaffold(mols, activities, scaffold_smiles)
        
        sar_data = []
        for smi, data in scaffold_groups.items():
            if data["n_compounds"] >= 3:  # 至少3个化合物
                sar_data.append({
                    "scaffold_smiles": smi,
                    "n_compounds": data["n_compounds"],
                    "mean_activity": data["mean_activity"],
                    "std_activity": data["std_activity"],
                    "activity_range": data["max_activity"] - data["min_activity"]
                })
        
        sar_df = pd.DataFrame(sar_data)
        sar_df = sar_df.sort_values("activity_range", ascending=False)
        sar_file = Path(output_dir) / "scaffold_sar_analysis.csv"
        sar_df.to_csv(sar_file, index=False)
        print(f"  SAR分析结果已保存: {sar_file}")
        print(f"  分析了 {len(sar_df)} 个系列（每个系列>=3个化合物）")
        
        # 5. 支架划分
        print("\n[5/6] 支架划分（训练/测试集）...")
        train_idx, test_idx = scaffold_train_test_split(
            mols, scaffold_smiles, test_ratio=0.2
        )
        
        split_df = pd.DataFrame({
            "index": range(len(mols)),
            "SMILES": smiles_list,
            "scaffold": scaffold_smiles,
            "split": ["train" if i in train_idx else "test" for i in range(len(mols))]
        })
        split_file = Path(output_dir) / "scaffold_split.csv"
        split_df.to_csv(split_file, index=False)
        print(f"  划分结果已保存: {split_file}")
        print(f"  训练集: {len(train_idx)} 个分子")
        print(f"  测试集: {len(test_idx)} 个分子")
        print(f"  训练集支架数: {len(set(scaffold_smiles[i] for i in train_idx))}")
        print(f"  测试集支架数: {len(set(scaffold_smiles[i] for i in test_idx))}")
    else:
        print("\n[4-5/6] 跳过SAR分析 (未检测到活性列)")
    
    # 6. 可视化
    print(f"\n[6/6] 生成支架可视化...")
    
    # 准备数据
    freq_list = []
    for smi in scaffold_smiles:
        if smi is not None:
            count = scaffold_freq_df[scaffold_freq_df["scaffold_smiles"] == smi]["count"].values
            freq_list.append(count[0] if len(count) > 0 else 1)
        else:
            freq_list.append(0)
    
    viz_file = Path(output_dir) / "top_scaffolds.png"
    visualize_scaffolds(scaffolds, scaffold_smiles, freq_list, str(viz_file), top_n_visualize)
    
    print("\n" + "=" * 60)
    print("支架分析完成!")
    print("=" * 60)
    print(f"\n输出文件:")
    print(f"  - {freq_file.name}")
    if has_activity:
        print(f"  - {sar_file.name}")
        print(f"  - {split_file.name}")
    print(f"  - {viz_file.name}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="分子支架分析与SAR研究"
    )
    parser.add_argument(
        "input_file",
        help="输入CSV文件 (包含SMILES列)"
    )
    parser.add_argument(
        "-o", "--output",
        default="scaffold_analysis",
        help="输出目录"
    )
    parser.add_argument(
        "--activity-column",
        default="activity",
        help="活性值列名"
    )
    parser.add_argument(
        "--smiles-column",
        default="SMILES",
        help="SMILES列名"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=20,
        help="可视化的支架数量"
    )
    
    args = parser.parse_args()
    
    main(
        input_file=args.input_file,
        output_dir=args.output,
        activity_column=args.activity_column,
        smiles_column=args.smiles_column,
        top_n_visualize=args.top_n
    )
