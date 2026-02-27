#!/usr/bin/env python3
"""
ENA数据库元数据提取与报表生成脚本
功能：
1. 提取研究的完整元数据
2. 生成样本信息报表
3. 导出为多种格式（CSV, Excel, JSON）
4. 统计分析

依赖：requests, pandas
"""

import requests
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional
from collections import Counter


ENA_PORTAL_API = "https://www.ebi.ac.uk/ena/portal/api"


def fetch_study_metadata(study_accession: str) -> Dict:
    """
    获取研究的完整元数据
    
    Args:
        study_accession: 研究accession (如 PRJEB1234)
    
    Returns:
        元数据字典
    """
    url = f"{ENA_PORTAL_API}/search"
    params = {
        "result": "study",
        "query": f"accession={study_accession}",
        "format": "json"
    }
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    results = response.json()
    if results:
        return results[0]
    return {}


def fetch_sample_metadata(
    study_accession: str,
    limit: int = 10000
) -> pd.DataFrame:
    """
    获取研究的所有样本元数据
    
    Args:
        study_accession: 研究accession
        limit: 最大样本数
    
    Returns:
        样本元数据DataFrame
    """
    url = f"{ENA_PORTAL_API}/search"
    params = {
        "result": "sample",
        "query": f"study_accession={study_accession}",
        "format": "json",
        "limit": limit
    }
    
    print(f"获取样本元数据: {study_accession}")
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    samples = response.json()
    
    if not samples:
        return pd.DataFrame()
    
    df = pd.DataFrame(samples)
    print(f"  获取了 {len(df)} 个样本")
    
    return df


def fetch_run_metadata(
    study_accession: str,
    limit: int = 10000
) -> pd.DataFrame:
    """
    获取研究的所有run元数据
    
    Args:
        study_accession: 研究accession
        limit: 最大run数
    
    Returns:
        Run元数据DataFrame
    """
    url = f"{ENA_PORTAL_API}/search"
    params = {
        "result": "run",
        "query": f"study_accession={study_accession}",
        "format": "json",
        "limit": limit
    }
    
    print(f"获取Run元数据: {study_accession}")
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    runs = response.json()
    
    if not runs:
        return pd.DataFrame()
    
    df = pd.DataFrame(runs)
    print(f"  获取了 {len(df)} 个runs")
    
    return df


def analyze_sample_characteristics(
    sample_df: pd.DataFrame
) -> Dict:
    """
    分析样本特征
    
    Args:
        sample_df: 样本元数据DataFrame
    
    Returns:
        统计分析结果
    """
    analysis = {
        'total_samples': len(sample_df),
        'columns': list(sample_df.columns)
    }
    
    # 分类变量统计
    categorical_cols = sample_df.select_dtypes(include=['object']).columns
    
    for col in categorical_cols:
        if col in sample_df.columns:
            value_counts = sample_df[col].value_counts()
            analysis[f'{col}_distribution'] = value_counts.to_dict()
    
    # 数值变量统计
    numeric_cols = sample_df.select_dtypes(include=['number']).columns
    
    for col in numeric_cols:
        if col in sample_df.columns:
            analysis[f'{col}_stats'] = {
                'mean': sample_df[col].mean(),
                'std': sample_df[col].std(),
                'min': sample_df[col].min(),
                'max': sample_df[col].max()
            }
    
    return analysis


def generate_study_report(
    study_accession: str,
    output_dir: str = "./ena_reports"
) -> Dict:
    """
    生成完整的研究报告
    
    Args:
        study_accession: 研究accession
        output_dir: 输出目录
    
    Returns:
        报告信息字典
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"=" * 60)
    print(f"生成研究报告: {study_accession}")
    print(f"=" * 60)
    
    # 1. 获取研究元数据
    study_meta = fetch_study_metadata(study_accession)
    
    # 2. 获取样本元数据
    sample_df = fetch_sample_metadata(study_accession)
    
    # 3. 获取run元数据
    run_df = fetch_run_metadata(study_accession)
    
    # 4. 保存原始数据
    study_file = Path(output_dir) / f"{study_accession}_study.json"
    import json
    with open(study_file, 'w') as f:
        json.dump(study_meta, f, indent=2)
    
    if not sample_df.empty:
        sample_file = Path(output_dir) / f"{study_accession}_samples.csv"
        sample_df.to_csv(sample_file, index=False)
    
    if not run_df.empty:
        run_file = Path(output_dir) / f"{study_accession}_runs.csv"
        run_df.to_csv(run_file, index=False)
    
    # 5. 生成分析
    report = {
        'study_accession': study_accession,
        'study_title': study_meta.get('study_title', 'N/A'),
        'study_description': study_meta.get('description', 'N/A'),
        'num_samples': len(sample_df),
        'num_runs': len(run_df)
    }
    
    if not sample_df.empty:
        report['sample_analysis'] = analyze_sample_characteristics(sample_df)
    
    # 6. 生成Markdown报告
    report_md = f"""# ENA研究报告: {study_accession}

## 研究信息

- **Accession**: {report['study_accession']}
- **标题**: {report['study_title']}
- **样本数**: {report['num_samples']}
- **Run数**: {report['num_runs']}

## 描述

{report['study_description']}

## 样本分析

"""
    
    if 'sample_analysis' in report:
        sa = report['sample_analysis']
        report_md += f"\n**总样本数**: {sa['total_samples']}\n\n"
        
        # 添加关键分布
        for key, value in sa.items():
            if key.endswith('_distribution') and isinstance(value, dict):
                col_name = key.replace('_distribution', '')
                report_md += f"\n### {col_name}分布\n\n"
                for val, count in list(value.items())[:10]:
                    report_md += f"- {val}: {count}\n"
    
    report_md += f"""
## 数据文件

- 研究元数据: `{study_file.name}`
- 样本元数据: `{study_accession}_samples.csv`
- Run元数据: `{study_accession}_runs.csv`

---
*报告生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}*
"""
    
    report_file = Path(output_dir) / f"{study_accession}_report.md"
    with open(report_file, 'w') as f:
        f.write(report_md)
    
    print(f"\n报告已生成: {report_file}")
    
    return report


def compare_studies(
    study_accessions: List[str],
    output_file: str = "study_comparison.csv"
) -> pd.DataFrame:
    """
    比较多个研究
    
    Args:
        study_accessions: 研究accession列表
        output_file: 输出文件
    
    Returns:
        比较结果DataFrame
    """
    comparisons = []
    
    for accession in study_accessions:
        try:
            study_meta = fetch_study_metadata(accession)
            sample_df = fetch_sample_metadata(accession)
            run_df = fetch_run_metadata(accession)
            
            comparisons.append({
                'accession': accession,
                'title': study_meta.get('study_title', 'N/A')[:100],
                'num_samples': len(sample_df),
                'num_runs': len(run_df)
            })
        except Exception as e:
            print(f"获取 {accession} 失败: {e}")
            comparisons.append({
                'accession': accession,
                'title': 'ERROR',
                'num_samples': 0,
                'num_runs': 0
            })
    
    df = pd.DataFrame(comparisons)
    df.to_csv(output_file, index=False)
    
    print(f"\n比较结果已保存: {output_file}")
    print(f"\n{'='*60}")
    print(df.to_string(index=False))
    
    return df


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="ENA数据库元数据提取与报表生成"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="可用命令")
    
    # 生成报告
    report_parser = subparsers.add_parser("report", help="生成研究报告")
    report_parser.add_argument("accession", help="研究accession")
    report_parser.add_argument("--output", "-o", default="./ena_reports")
    
    # 比较研究
    compare_parser = subparsers.add_parser("compare", help="比较多个研究")
    compare_parser.add_argument("accessions", nargs="+", help="研究accessions")
    compare_parser.add_argument("--output", "-o", default="study_comparison.csv")
    
    # 获取样本
    sample_parser = subparsers.add_parser("samples", help="获取样本元数据")
    sample_parser.add_argument("accession", help="研究accession")
    sample_parser.add_argument("--output", "-o", help="输出CSV文件")
    
    args = parser.parse_args()
    
    if args.command == "report":
        generate_study_report(args.accession, args.output)
    elif args.command == "compare":
        compare_studies(args.accessions, args.output)
    elif args.command == "samples":
        df = fetch_sample_metadata(args.accession)
        if args.output:
            df.to_csv(args.output, index=False)
            print(f"已保存: {args.output}")
        else:
            print(df.head(20).to_string())
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
