#!/usr/bin/env python3
"""
GWAS Catalog 示例脚本 2: 批量数据提取和分析
演示如何批量下载和分析 GWAS 数据
"""

import requests
import pandas as pd
import json
from time import sleep
from pathlib import Path


def extract_trait_associations(efo_id, p_threshold=5e-8, output_file=None):
    """
    提取特定性状的所有显著关联
    
    Args:
        efo_id: EFO 性状 ID (例如: EFO_0001360 for type 2 diabetes)
        p_threshold: p 值阈值
        output_file: 输出文件路径
    
    Returns:
        DataFrame 包含所有关联
    """
    print(f"\n提取 EFO ID {efo_id} 的关联 (p < {p_threshold})...")
    
    base_url = "https://www.ebi.ac.uk/gwas/rest/api"
    endpoint = f"efoTraits/{efo_id}/associations"
    
    headers = {"Content-Type": "application/json"}
    results = []
    page = 0
    max_pages = 50  # 安全限制
    
    while page < max_pages:
        params = {"page": page, "size": 100}
        
        try:
            response = requests.get(
                f"{base_url}/{endpoint}",
                params=params,
                headers=headers,
                timeout=30
            )
            response.raise_for_status()
            data = response.json()
            
            associations = data.get('_embedded', {}).get('associations', [])
            
            if not associations:
                break
            
            for assoc in associations:
                pvalue = assoc.get('pvalue')
                if pvalue and float(pvalue) <= p_threshold:
                    results.append({
                        'variant': assoc.get('rsId'),
                        'chromosome': assoc.get('chromosomeName'),
                        'position': assoc.get('chromosomePosition'),
                        'pvalue': float(pvalue),
                        'risk_allele': assoc.get('strongestAllele'),
                        'effect_size': assoc.get('orPerCopyNum') or assoc.get('betaNum'),
                        'trait': assoc.get('efoTrait'),
                        'pmid': assoc.get('pubmedId'),
                        'study': assoc.get('studyAccessions', [None])[0]
                    })
            
            print(f"  已处理页面 {page + 1}, 累计结果: {len(results)}")
            
            # 检查是否有更多页面
            if '_links' not in data or 'next' not in data['_links']:
                break
            
            page += 1
            sleep(0.2)  # 速率限制
            
        except Exception as e:
            print(f"  错误: {e}")
            break
    
    df = pd.DataFrame(results)
    
    if output_file and not df.empty:
        df.to_csv(output_file, index=False)
        print(f"✓ 结果已保存到 {output_file}")
    
    return df


def analyze_association_summary(df):
    """分析关联数据的摘要统计"""
    if df.empty:
        print("没有数据可供分析")
        return
    
    print("\n" + "=" * 60)
    print("关联数据摘要")
    print("=" * 60)
    
    print(f"\n总关联数: {len(df)}")
    print(f"唯一变异数: {df['variant'].nunique()}")
    
    if 'chromosome' in df.columns:
        print(f"\n按染色体分布:")
        chrom_counts = df['chromosome'].value_counts().head(10)
        for chrom, count in chrom_counts.items():
            print(f"  {chrom}: {count}")
    
    if 'pvalue' in df.columns:
        print(f"\nP 值分布:")
        print(f"  最小: {df['pvalue'].min():.2e}")
        print(f"  最大: {df['pvalue'].max():.2e}")
        print(f"  中位数: {df['pvalue'].median():.2e}")
        
        # 极显著关联
        very_significant = df[df['pvalue'] < 1e-20]
        print(f"\n  极显著关联 (p < 1e-20): {len(very_significant)}")
    
    if 'effect_size' in df.columns:
        # 过滤掉缺失值
        valid_effects = df['effect_size'].dropna()
        if len(valid_effects) > 0:
            print(f"\n效应值分布:")
            print(f"  范围: {valid_effects.min():.4f} - {valid_effects.max():.4f}")


def create_summary_statistics_report(study_accessions, output_dir="gwas_reports"):
    """
    为多个研究创建汇总统计报告
    
    Args:
        study_accessions: 研究 GCST ID 列表
        output_dir: 输出目录
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    report_data = []
    
    print(f"\n正在处理 {len(study_accessions)} 个研究...")
    
    for accession in study_accessions:
        print(f"\n  查询研究 {accession}...")
        
        study_info = query_study_info(accession)
        
        if study_info:
            report_data.append({
                'accession': accession,
                'pubmed_id': study_info.get('pubmedId'),
                'author': study_info.get('author'),
                'publication_date': study_info.get('publicationDate'),
                'trait': study_info.get('trait'),
                'initial_sample_size': study_info.get('initialSampleSize'),
                'replication_sample_size': study_info.get('replicationSampleSize'),
            })
        
        sleep(0.5)
    
    # 保存报告
    report_df = pd.DataFrame(report_data)
    report_file = Path(output_dir) / "study_summary_report.csv"
    report_df.to_csv(report_file, index=False)
    
    print(f"\n✓ 研究报告已保存到 {report_file}")
    print(f"\n报告预览:")
    print(report_df.head().to_string(index=False))
    
    return report_df


def query_study_info(accession):
    """查询研究详细信息"""
    url = f"https://www.ebi.ac.uk/gwas/rest/api/studies/{accession}"
    
    try:
        response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
        if response.status_code == 200:
            return response.json()
    except Exception as e:
        print(f"  查询失败: {e}")
    
    return None


def main():
    print("=" * 70)
    print("GWAS Catalog 批量数据提取")
    print("=" * 70)
    
    # 示例 1: 提取特定性状数据
    print("\n示例 1: 提取性状关联数据")
    print("-" * 70)
    
    # EFO IDs for common diseases
    efo_examples = {
        "EFO_0001360": "type 2 diabetes",
        # "EFO_0002508": "Parkinson disease",
        # "EFO_0000612": "myocardial infarction",
    }
    
    for efo_id, trait_name in efo_examples.items():
        print(f"\n查询: {trait_name} ({efo_id})")
        
        df = extract_trait_associations(
            efo_id,
            p_threshold=5e-8,
            output_file=f"{trait_name.replace(' ', '_')}_associations.csv"
        )
        
        if not df.empty:
            analyze_association_summary(df)
        else:
            print(f"  未找到数据或 API 限制")
        
        sleep(1)
    
    # 示例 2: 创建研究报告
    print("\n" + "=" * 70)
    print("示例 2: 研究信息汇总")
    print("-" * 70)
    
    # 示例研究 ID
    sample_studies = ["GCST001795", "GCST001234"]
    
    # 注意: 这些 ID 可能不存在，仅作演示
    print("\n演示: 查询研究元数据")
    for study_id in sample_studies:
        study_info = query_study_info(study_id)
        if study_info:
            print(f"\n  研究 {study_id}:")
            print(f"    作者: {study_info.get('author', 'N/A')}")
            print(f"    PubMed: {study_info.get('pubmedId', 'N/A')}")
        else:
            print(f"\n  研究 {study_id}: 未找到 (演示目的)")
    
    # 示例 3: 数据质量检查
    print("\n" + "=" * 70)
    print("示例 3: 数据质量检查清单")
    print("-" * 70)
    
    quality_checks = """
GWAS 数据质量检查清单:

1. 显著性阈值:
   ✓ 全基因组显著性: p < 5×10⁻⁸
   ✓ 建议阈值: p < 1×10⁻⁵ (建议性关联)

2. 研究质量指标:
   ✓ 样本量 > 1000 (发现队列)
   ✓ 有独立重复验证
   ✓ 多祖先群体分析

3. 数据完整性检查:
   ✓ 效应方向一致性
   ✓ 等位基因频率合理性
   ✓ 无显著异质性 (I² < 75%)

4. 下游分析准备:
   ✓ 按 p 值排序
   ✓ 检查连锁不平衡
   ✓ 准备多基因风险评分输入
"""
    print(quality_checks)
    
    print("\n" + "=" * 70)
    print("批量提取完成!")
    print("=" * 70)
    
    # 清理临时文件
    for f in Path(".").glob("*_associations.csv"):
        print(f"\n注意: 临时文件已保存: {f}")


if __name__ == "__main__":
    main()
