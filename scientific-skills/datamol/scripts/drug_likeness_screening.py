#!/usr/bin/env python3
"""
Datamol 分子属性计算与药物筛选脚本
功能：
1. 从SMILES字符串或SDF文件加载分子
2. 计算分子描述符（分子量、LogP、HBD/HBA等）
3. 应用Lipinski规则筛选药物样分子
4. 输出筛选结果到CSV

依赖：datamol, pandas
"""

import datamol as dm
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional


def load_molecules(input_path: str, smiles_column: str = "SMILES") -> List[dm.Mol]:
    """
    从文件加载分子
    
    Args:
        input_path: 输入文件路径（.sdf, .csv, .smi等）
        smiles_column: CSV文件中的SMILES列名
    
    Returns:
        分子对象列表
    """
    path = Path(input_path)
    mols = []
    
    if path.suffix == ".sdf":
        df = dm.read_sdf(str(path), mol_column="mol")
        mols = df["mol"].dropna().tolist()
    elif path.suffix == ".csv":
        df = dm.read_csv(str(path), smiles_column=smiles_column, mol_column="mol")
        mols = df["mol"].dropna().tolist()
    elif path.suffix == ".smi":
        df = dm.read_smi(str(path), smiles_column="smiles", mol_column="mol")
        mols = df["mol"].dropna().tolist()
    else:
        # 尝试通用读取
        df = dm.open_df(str(path), mol_column="mol")
        mols = df["mol"].dropna().tolist()
    
    return mols


def standardize_molecules(mols: List[dm.Mol]) -> List[dm.Mol]:
    """标准化分子结构"""
    standardized = []
    for mol in mols:
        if mol is not None:
            try:
                std_mol = dm.standardize_mol(
                    mol,
                    disconnect_metals=True,
                    normalize=True,
                    reionize=True
                )
                if std_mol is not None:
                    standardized.append(std_mol)
            except Exception as e:
                print(f"标准化失败: {e}")
    return standardized


def compute_descriptors(mols: List[dm.Mol], n_jobs: int = -1) -> pd.DataFrame:
    """
    批量计算分子描述符
    
    Args:
        mols: 分子对象列表
        n_jobs: 并行处理核心数，-1表示使用所有核心
    
    Returns:
        描述符DataFrame
    """
    desc_df = dm.descriptors.batch_compute_many_descriptors(
        mols,
        n_jobs=n_jobs,
        progress=True
    )
    return desc_df


def lipinski_filter(descriptors: pd.DataFrame) -> pd.Series:
    """
    Lipinski's Rule of Five 药物筛选
    
    规则：
    - 分子量 <= 500
    - LogP <= 5
    - 氢键供体数 <= 5
    - 氢键受体数 <= 10
    
    违反不超过1条的分子被认为是药物样分子
    """
    rules = {
        'mw': descriptors['mw'] <= 500,
        'logp': descriptors['logp'] <= 5,
        'hbd': descriptors['hbd'] <= 5,
        'hba': descriptors['hba'] <= 10
    }
    
    violations = sum(~rule for rule in rules.values())
    return violations <= 1


def main(
    input_file: str,
    output_file: str = "druglike_compounds.csv",
    smiles_column: str = "SMILES",
    standardize: bool = True
):
    """
    主函数：分子筛选流程
    
    Args:
        input_file: 输入分子文件
        output_file: 输出CSV文件
        smiles_column: CSV文件中的SMILES列名
        standardize: 是否标准化分子
    """
    print(f"=" * 60)
    print(f"药物样分子筛选流程")
    print(f"输入文件: {input_file}")
    print(f"=" * 60)
    
    # 1. 加载分子
    print("\n[1/5] 加载分子...")
    mols = load_molecules(input_file, smiles_column)
    print(f"  加载了 {len(mols)} 个分子")
    
    # 2. 标准化
    if standardize:
        print("\n[2/5] 标准化分子...")
        mols = standardize_molecules(mols)
        print(f"  标准化后: {len(mols)} 个分子")
    
    # 3. 计算描述符
    print("\n[3/5] 计算分子描述符...")
    desc_df = compute_descriptors(mols)
    print(f"  计算了 {len(desc_df.columns)} 个描述符")
    
    # 4. Lipinski筛选
    print("\n[4/5] 应用Lipinski规则筛选...")
    druglike_mask = lipinski_filter(desc_df)
    druglike_count = druglike_mask.sum()
    print(f"  药物样分子: {druglike_count}/{len(mols)} ({100*druglike_count/len(mols):.1f}%)")
    
    # 5. 生成输出
    print("\n[5/5] 生成输出文件...")
    
    # 添加SMILES
    smiles_list = [dm.to_smiles(mol) for mol in mols]
    desc_df["SMILES"] = smiles_list
    desc_df["is_druglike"] = druglike_mask
    
    # 保存全部结果
    desc_df.to_csv(output_file, index=False)
    print(f"  已保存全部结果: {output_file}")
    
    # 保存药物样分子
    druglike_file = output_file.replace(".csv", "_druglike.csv")
    druglike_df = desc_df[druglike_mask]
    druglike_df.to_csv(druglike_file, index=False)
    print(f"  已保存药物样分子: {druglike_file}")
    
    # 统计摘要
    print("\n" + "=" * 60)
    print("统计摘要")
    print("=" * 60)
    print(f"总分子数: {len(mols)}")
    print(f"药物样分子: {druglike_count}")
    print(f"通过率: {100*druglike_count/len(mols):.2f}%")
    print("\n描述符统计 (药物样分子):")
    for col in ['mw', 'logp', 'hbd', 'hba', 'tpsa']:
        if col in druglike_df.columns:
            print(f"  {col}: mean={druglike_df[col].mean():.2f}, "
                  f"std={druglike_df[col].std():.2f}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="分子药物筛选 - Lipinski规则筛选药物样分子"
    )
    parser.add_argument(
        "input_file",
        help="输入分子文件 (.sdf, .csv, .smi)"
    )
    parser.add_argument(
        "-o", "--output",
        default="druglike_compounds.csv",
        help="输出CSV文件路径"
    )
    parser.add_argument(
        "--smiles-column",
        default="SMILES",
        help="CSV文件中的SMILES列名"
    )
    parser.add_argument(
        "--no-standardize",
        action="store_true",
        help="跳过分子标准化"
    )
    
    args = parser.parse_args()
    
    main(
        input_file=args.input_file,
        output_file=args.output,
        smiles_column=args.smiles_column,
        standardize=not args.no_standardize
    )
