#!/usr/bin/env python3
"""
GWAS Catalog 示例脚本 1: 查询 SNP-性状关联
演示如何使用 REST API 查询遗传变异与疾病/性状的关联
"""

import requests
import pandas as pd
from time import sleep


def query_gwas_api(endpoint, params=None, base_url="https://www.ebi.ac.uk/gwas/rest/api"):
    """
    查询 GWAS Catalog REST API
    
    Args:
        endpoint: API 端点路径
        params: 查询参数
        base_url: API 基础 URL
    
    Returns:
        JSON 响应数据
    """
    url = f"{base_url}/{endpoint}"
    headers = {"Content-Type": "application/json"}
    
    try:
        response = requests.get(url, params=params, headers=headers, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"请求错误: {e}")
        return None


def search_variant_associations(variant_id):
    """
    查询特定 SNP 的所有性状关联
    
    Args:
        variant_id: rs ID (例如: rs7903146)
    """
    print(f"\n查询变异 {variant_id} 的关联...")
    
    # 获取变异详细信息
    variant_info = query_gwas_api(f"singleNucleotidePolymorphisms/{variant_id}")
    
    if variant_info:
        print(f"✓ 变异位置: {variant_info.get('chromosomeName')}:{variant_info.get('chromosomePosition')}")
        print(f"  基因: {variant_info.get('geneNames', 'N/A')}")
    
    # 获取关联信息
    associations = query_gwas_api(
        f"singleNucleotidePolymorphisms/{variant_id}/associations",
        params={"projection": "associationBySnp"}
    )
    
    results = []
    if associations and '_embedded' in associations:
        for assoc in associations['_embedded'].get('associations', []):
            results.append({
                'variant': variant_id,
                'trait': assoc.get('efoTrait', 'N/A'),
                'pvalue': assoc.get('pvalue'),
                'risk_allele': assoc.get('strongestAllele'),
                'or_beta': assoc.get('orPerCopyNum') or assoc.get('betaNum'),
                'pmid': assoc.get('pubmedId')
            })
    
    return results


def search_trait_associations(trait_name, p_threshold=5e-8, max_results=50):
    """
    按性状搜索关联
    
    Args:
        trait_name: 性状名称 (例如: type 2 diabetes)
        p_threshold: p 值阈值 (默认: 5e-8, 全基因组显著性)
        max_results: 最大结果数
    """
    print(f"\n搜索性状 '{trait_name}' 的关联...")
    
    # 注意: 这里使用简化的搜索方法
    # 实际应用中可能需要先查询 EFO ID
    
    # 使用摘要统计 API
    base_url = "https://www.ebi.ac.uk/gwas/summary-statistics/api"
    
    # 搜索关联 (需要实际的 EFO ID)
    # 这里展示基本模式
    
    results = []
    
    # 示例: 使用 REST API 搜索
    associations = query_gwas_api(
        "associations/search/findByEfoTrait",
        params={"trait": trait_name, "page": 0, "size": max_results}
    )
    
    return results


def main():
    print("=" * 70)
    print("GWAS Catalog API 查询示例")
    print("=" * 70)
    
    # 示例 1: 查询著名变异
    print("\n示例 1: 查询 TCF7L2 变异 rs7903146 (与 2 型糖尿病相关)")
    print("-" * 70)
    
    variant_results = search_variant_associations("rs7903146")
    
    if variant_results:
        df = pd.DataFrame(variant_results)
        print(f"\n找到 {len(df)} 个关联:")
        print(df.head(10).to_string(index=False))
    else:
        print("未找到关联信息或 API 请求失败")
    
    # 示例 2: 批量查询多个变异
    print("\n" + "=" * 70)
    print("示例 2: 批量查询多个疾病相关变异")
    print("-" * 70)
    
    disease_variants = [
        "rs429358",   # APOE - Alzheimer's disease
        "rs7412",     # APOE - Alzheimer's disease
        "rs1801133",  # MTHFR - various traits
    ]
    
    all_results = []
    for variant in disease_variants:
        results = search_variant_associations(variant)
        all_results.extend(results)
        sleep(0.5)  # 遵守速率限制
    
    if all_results:
        df = pd.DataFrame(all_results)
        print(f"\n总计找到 {len(df)} 个关联")
        
        # 按变异统计
        print("\n按变异统计关联数量:")
        print(df.groupby('variant').size())
    
    # 示例 3: 查询特定区域的变异
    print("\n" + "=" * 70)
    print("示例 3: 查询特定基因组区域的变异")
    print("-" * 70)
    
    # TCF7L2 基因区域 (chr10:114000000-115000000)
    region_variants = query_gwas_api(
        "singleNucleotidePolymorphisms/search/findByChromBpLocationRange",
        params={"chrom": "10", "bpStart": 114000000, "bpEnd": 115000000}
    )
    
    if region_variants:
        print("✓ 成功查询区域变异")
        # 处理结果...
    else:
        print("区域查询演示 (需要有效参数)")
    
    # 示例 4: 获取研究信息
    print("\n" + "=" * 70)
    print("示例 4: 获取研究详细信息")
    print("-" * 70)
    
    study_id = "GCST001234"  # 示例研究 ID
    study_info = query_gwas_api(f"studies/{study_id}")
    
    if study_info:
        print(f"研究 ID: {study_info.get('accessionId')}")
        print(f"PubMed ID: {study_info.get('pubmedId')}")
        print(f"作者: {study_info.get('author')}")
        print(f"发表日期: {study_info.get('publicationDate')}")
        print(f"样本量: {study_info.get('initialSampleSize', 'N/A')}")
    else:
        print(f"研究 {study_id} 查询演示 (可能不存在)")
    
    print("\n" + "=" * 70)
    print("查询完成!")
    print("=" * 70)
    print("\n提示:")
    print("- 访问 https://www.ebi.ac.uk/gwas/ 进行交互式搜索")
    print("- API 文档: https://www.ebi.ac.uk/gwas/rest/docs/api")
    print("- 引用 GWAS Catalog: Sollis et al. (2023) Nucleic Acids Research")


if __name__ == "__main__":
    main()
