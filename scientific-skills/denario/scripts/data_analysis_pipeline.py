#!/usr/bin/env python3
"""
Denario 数据分析与可视化生成脚本
功能：
1. 执行数据集探索性分析
2. 生成统计图表和可视化
3. 输出分析结果供Denario使用
4. 创建可重用的分析模板

依赖：denario, pandas, matplotlib, seaborn, numpy
"""

import os
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import matplotlib.pyplot as plt
import seaborn as sns


def exploratory_data_analysis(
    data_path: str,
    output_dir: str,
    target_column: Optional[str] = None
) -> Dict:
    """
    执行探索性数据分析(EDA)
    
    Args:
        data_path: 数据文件路径(CSV, Excel等)
        output_dir: 输出目录
        target_column: 目标变量列名
    
    Returns:
        分析结果字典
    """
    print(f"加载数据: {data_path}")
    
    # 读取数据
    path = Path(data_path)
    if path.suffix == ".csv":
        df = pd.read_csv(data_path)
    elif path.suffix in [".xlsx", ".xls"]:
        df = pd.read_excel(data_path)
    elif path.suffix == ".json":
        df = pd.read_json(data_path)
    else:
        raise ValueError(f"不支持的文件格式: {path.suffix}")
    
    print(f"数据形状: {df.shape}")
    
    results = {
        "dataset_name": path.stem,
        "shape": df.shape,
        "columns": list(df.columns),
        "numeric_columns": list(df.select_dtypes(include=[np.number]).columns),
        "categorical_columns": list(df.select_dtypes(include=['object', 'category']).columns),
        "missing_values": df.isnull().sum().to_dict(),
        "missing_percentage": (df.isnull().sum() / len(df) * 100).to_dict()
    }
    
    # 数值列统计
    if results["numeric_columns"]:
        results["numeric_stats"] = df[results["numeric_columns"]].describe().to_dict()
    
    # 目标变量分析
    if target_column and target_column in df.columns:
        results["target_column"] = target_column
        if df[target_column].dtype in [np.number]:
            results["target_stats"] = {
                "mean": float(df[target_column].mean()),
                "std": float(df[target_column].std()),
                "min": float(df[target_column].min()),
                "max": float(df[target_column].max())
            }
        else:
            results["target_distribution"] = df[target_column].value_counts().to_dict()
    
    # 生成可视化
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # 1. 缺失值热图
    if df.isnull().sum().sum() > 0:
        plt.figure(figsize=(12, 6))
        sns.heatmap(df.isnull(), cbar=True, yticklabels=False, cmap='viridis')
        plt.title("Missing Values Pattern")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / "missing_values.png", dpi=300)
        plt.close()
    
    # 2. 数值列分布
    numeric_cols = results["numeric_columns"][:6]  # 最多6个
    if numeric_cols:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        
        for i, col in enumerate(numeric_cols):
            if i < len(axes):
                df[col].hist(ax=axes[i], bins=30, edgecolor='black')
                axes[i].set_title(col)
                axes[i].set_xlabel("Value")
                axes[i].set_ylabel("Frequency")
        
        # 隐藏多余的子图
        for i in range(len(numeric_cols), len(axes)):
            axes[i].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(Path(output_dir) / "numeric_distributions.png", dpi=300)
        plt.close()
    
    # 3. 相关性热图
    if len(results["numeric_columns"]) > 1:
        plt.figure(figsize=(12, 10))
        corr_matrix = df[results["numeric_columns"]].corr()
        sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0,
                    fmt='.2f', square=True)
        plt.title("Feature Correlation Matrix")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / "correlation_matrix.png", dpi=300)
        plt.close()
    
    # 4. 目标变量分布
    if target_column and target_column in df.columns:
        plt.figure(figsize=(10, 6))
        if df[target_column].dtype in [np.number]:
            df[target_column].hist(bins=30, edgecolor='black')
            plt.title(f"Target Variable Distribution: {target_column}")
            plt.xlabel(target_column)
        else:
            df[target_column].value_counts().plot(kind='bar')
            plt.title(f"Target Variable Distribution: {target_column}")
            plt.xticks(rotation=45)
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(Path(output_dir) / "target_distribution.png", dpi=300)
        plt.close()
    
    # 保存结果
    results_file = Path(output_dir) / "eda_results.json"
    with open(results_file, 'w') as f:
        # 转换numpy类型为Python原生类型
        json_results = json.dumps(results, indent=2, default=str)
        f.write(json_results)
    
    print(f"分析结果已保存: {results_file}")
    
    return results


def generate_analysis_summary(
    eda_results: Dict,
    output_file: str
):
    """
    生成分析摘要文档(Markdown格式)
    
    Args:
        eda_results: EDA结果字典
        output_file: 输出Markdown文件路径
    """
    lines = [
        "# 数据分析摘要\n",
        f"## 数据集: {eda_results['dataset_name']}\n",
        f"**数据形状:** {eda_results['shape'][0]} 行 × {eda_results['shape'][1]} 列\n",
        "\n## 列信息\n",
        f"- **数值列:** {len(eda_results['numeric_columns'])} 个",
        f"- **类别列:** {len(eda_results['categorical_columns'])} 个",
        "\n### 数值列",
    ]
    
    for col in eda_results['numeric_columns'][:10]:
        lines.append(f"- {col}")
    if len(eda_results['numeric_columns']) > 10:
        lines.append(f"- ... 等共 {len(eda_results['numeric_columns'])} 个")
    
    lines.extend(["\n### 类别列"])
    for col in eda_results['categorical_columns'][:10]:
        lines.append(f"- {col}")
    if len(eda_results['categorical_columns']) > 10:
        lines.append(f"- ... 等共 {len(eda_results['categorical_columns'])} 个")
    
    # 缺失值信息
    lines.extend(["\n## 数据质量"])
    missing_cols = [(k, v) for k, v in eda_results['missing_percentage'].items() if v > 0]
    if missing_cols:
        lines.append(f"\n发现 {len(missing_cols)} 列存在缺失值:")
        for col, pct in sorted(missing_cols, key=lambda x: x[1], reverse=True)[:10]:
            lines.append(f"- {col}: {pct:.2f}%")
    else:
        lines.append("\n数据完整，无缺失值。")
    
    # 目标变量
    if 'target_column' in eda_results:
        lines.extend([
            f"\n## 目标变量: {eda_results['target_column']}",
        ])
        if 'target_stats' in eda_results:
            stats = eda_results['target_stats']
            lines.extend([
                f"\n- 均值: {stats['mean']:.4f}",
                f"- 标准差: {stats['std']:.4f}",
                f"- 最小值: {stats['min']:.4f}",
                f"- 最大值: {stats['max']:.4f}",
            ])
    
    # 建议
    lines.extend([
        "\n## 分析建议",
        "\n基于初步探索，建议进行以下分析：",
        "1. **特征工程**: 考虑对数值特征进行标准化或归一化",
        "2. **缺失值处理**: 根据缺失比例选择填充或删除策略",
        "3. **特征选择**: 利用相关性矩阵识别重要特征",
        "4. **模型选择**: 根据数据特性选择合适的机器学习模型",
    ])
    
    with open(output_file, 'w') as f:
        f.write('\n'.join(lines))
    
    print(f"分析摘要已保存: {output_file}")


def prepare_denario_input(
    data_path: str,
    output_dir: str,
    research_goal: str,
    target_column: Optional[str] = None
) -> str:
    """
    准备Denario输入文件
    
    Args:
        data_path: 数据文件路径
        output_dir: 输出目录
        research_goal: 研究目标描述
        target_column: 目标变量
    
    Returns:
        数据描述文本
    """
    # 执行EDA
    eda_results = exploratory_data_analysis(data_path, output_dir, target_column)
    
    # 生成摘要
    summary_file = Path(output_dir) / "analysis_summary.md"
    generate_analysis_summary(eda_results, str(summary_file))
    
    # 构建数据描述
    data_description = f"""
# 研究数据集

## 研究目标
{research_goal}

## 数据集信息
- **数据集名称:** {eda_results['dataset_name']}
- **样本数量:** {eda_results['shape'][0]} 
- **特征数量:** {eda_results['shape'][1]}

## 可用工具
- pandas: 数据处理和分析
- numpy: 数值计算
- matplotlib/seaborn: 数据可视化
- scikit-learn: 机器学习建模

## 数据特征
"""
    
    if eda_results['numeric_columns']:
        data_description += f"\n### 数值特征 ({len(eda_results['numeric_columns'])}个)\n"
        for col in eda_results['numeric_columns'][:10]:
            data_description += f"- {col}\n"
        if len(eda_results['numeric_columns']) > 10:
            data_description += f"- ... 等共{len(eda_results['numeric_columns'])}个特征\n"
    
    if target_column:
        data_description += f"\n### 目标变量\n- {target_column}\n"
    
    # 保存数据描述
    desc_file = Path(output_dir) / "data_description.txt"
    with open(desc_file, 'w') as f:
        f.write(data_description)
    
    print(f"\n数据描述已保存: {desc_file}")
    print(f"可将此描述用于Denario研究流程")
    
    return data_description


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="数据集分析与Denario输入准备"
    )
    parser.add_argument(
        "data_file",
        help="数据文件路径 (CSV, Excel, JSON)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="输出目录"
    )
    parser.add_argument(
        "--target",
        help="目标变量列名"
    )
    parser.add_argument(
        "--goal",
        default="分析数据集中的模式和关系，构建预测模型",
        help="研究目标描述"
    )
    
    args = parser.parse_args()
    
    # 准备Denario输入
    data_description = prepare_denario_input(
        data_path=args.data_file,
        output_dir=args.output,
        research_goal=args.goal,
        target_column=args.target
    )
    
    print("\n" + "=" * 60)
    print("数据准备完成!")
    print("=" * 60)
    print(f"\n输出目录: {args.output}")
    print("\n下一步:")
    print(f"  python automated_research.py \\")
    print(f"    --project my_research \\")
    print(f"    --data-desc {args.output}/data_description.txt")


if __name__ == "__main__":
    main()
