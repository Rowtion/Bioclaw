#!/usr/bin/env python3
"""
HMDB 示例脚本 2: 代谢物识别和通路分析
演示如何使用 HMDB 数据进行代谢物识别和通路分析
"""

import csv
import json
from pathlib import Path
from collections import defaultdict


def identify_metabolite_from_ms(mz, ppm_tolerance=10, hmdb_csv=None):
    """
    根据质荷比(m/z)从 HMDB 中识别可能的代谢物
    
    Args:
        mz: 实验测得的质荷比
        ppm_tolerance: ppm 容差
        hmdb_csv: HMDB CSV 文件路径
    
    Returns:
        候选代谢物列表
    """
    print(f"\n根据 m/z = {mz} 识别代谢物...")
    print(f"容差: ±{ppm_tolerance} ppm")
    
    # 计算 m/z 范围
    mz_range = mz * ppm_tolerance / 1e6
    mz_min = mz - mz_range
    mz_max = mz + mz_range
    
    print(f"搜索范围: {mz_min:.4f} - {mz_max:.4f}")
    
    # 如果没有 CSV 文件，返回模拟结果
    if not hmdb_csv or not Path(hmdb_csv).exists():
        print("  (未提供 HMDB 数据文件，返回示例结果)")
        return [
            {'name': 'Example Metabolite 1', 'mw': mz, 'formula': 'C6H12O6'},
            {'name': 'Example Metabolite 2', 'mw': mz * 0.99, 'formula': 'C5H10O5'},
        ]
    
    # 从 CSV 文件搜索
    candidates = []
    
    with open(hmdb_csv, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                mw = float(row.get('monoisotopic_molecular_weight', 0))
                if mz_min <= mw <= mz_max:
                    candidates.append({
                        'hmdb_id': row.get('accession'),
                        'name': row.get('name'),
                        'formula': row.get('chemical_formula'),
                        'mw': mw,
                        'ppm_diff': abs(mz - mw) / mz * 1e6
                    })
            except (ValueError, TypeError):
                continue
    
    # 按 ppm 差异排序
    candidates.sort(key=lambda x: x['ppm_diff'])
    
    return candidates


def pathway_analysis(metabolite_list, hmdb_data=None):
    """
    对代谢物列表进行通路富集分析
    
    Args:
        metabolite_list: 代谢物 HMDB ID 列表
        hmdb_data: HMDB 数据字典
    
    Returns:
        通路分析结果
    """
    print(f"\n通路分析: {len(metabolite_list)} 个代谢物")
    
    # 模拟通路数据
    pathway_db = {
        'SMP00044': {'name': 'Histidine Metabolism', 'metabolites': ['HMDB0000001', 'HMDB0000002']},
        'SMP00001': {'name': 'Amino Acid Metabolism', 'metabolites': ['HMDB0000001', 'HMDB0000003']},
        'SMP00020': {'name': 'Glycolysis', 'metabolites': ['HMDB0000122', 'HMDB0000123']},
        'SMP00021': {'name': 'Citric Acid Cycle', 'metabolites': ['HMDB0000094', 'HMDB0000158']},
        'SMP00022': {'name': 'Fatty Acid Metabolism', 'metabolites': ['HMDB0000220', 'HMDB0000221']},
    }
    
    # 统计每个通路的命中数
    pathway_counts = defaultdict(int)
    pathway_hits = defaultdict(list)
    
    for met_id in metabolite_list:
        for pathway_id, pathway_info in pathway_db.items():
            if met_id in pathway_info['metabolites']:
                pathway_counts[pathway_id] += 1
                pathway_hits[pathway_id].append(met_id)
    
    # 计算富集度 (简化版)
    results = []
    for pathway_id, count in pathway_counts.items():
        pathway_info = pathway_db[pathway_id]
        pathway_size = len(pathway_info['metabolites'])
        
        # 简化富集计算
        enrichment = count / pathway_size if pathway_size > 0 else 0
        
        results.append({
            'pathway_id': pathway_id,
            'pathway_name': pathway_info['name'],
            'hits': count,
            'pathway_size': pathway_size,
            'enrichment': enrichment,
            'metabolites': pathway_hits[pathway_id]
        })
    
    # 按富集度排序
    results.sort(key=lambda x: x['enrichment'], reverse=True)
    
    return results


def biomarker_search(disease_name, hmdb_data=None):
    """
    搜索与疾病相关的生物标志物
    
    Args:
        disease_name: 疾病名称
        hmdb_data: HMDB 数据
    
    Returns:
        相关代谢物列表
    """
    print(f"\n搜索疾病 '{disease_name}' 的生物标志物...")
    
    # 模拟生物标志物数据
    biomarker_db = {
        'prostate cancer': [
            {'hmdb_id': 'HMDB0000001', 'name': '1-Methylhistidine', 'change': 'up', 'fold_change': 2.5},
            {'hmdb_id': 'HMDB0000002', 'name': 'Creatine', 'change': 'down', 'fold_change': 0.6},
        ],
        'diabetes': [
            {'hmdb_id': 'HMDB0000122', 'name': 'Glucose', 'change': 'up', 'fold_change': 3.0},
            {'hmdb_id': 'HMDB0000043', 'name': 'Insulin', 'change': 'variable', 'fold_change': None},
        ],
        'alzheimer': [
            {'hmdb_id': 'HMDB0000159', 'name': 'Homovanillic acid', 'change': 'down', 'fold_change': 0.5},
        ],
    }
    
    disease_lower = disease_name.lower()
    
    if disease_lower in biomarker_db:
        return biomarker_db[disease_lower]
    else:
        print(f"  未找到 '{disease_name}' 的数据")
        return []


def cross_reference_databases(hmdb_id):
    """
    获取代谢物在其他数据库中的链接
    
    Args:
        hmdb_id: HMDB ID
    
    Returns:
        外部数据库链接
    """
    external_db = {
        'HMDB0000122': {  # Glucose
            'kegg': 'C00031',
            'pubchem': '5793',
            'chebi': '4167',
            'metacyc': 'GLC',
            'lipidmaps': None,
        },
        'HMDB0000001': {  # 1-Methylhistidine
            'kegg': 'C01152',
            'pubchem': '92865',
            'chebi': '27596',
            'metacyc': '1-METHYLHISTIDINE',
            'lipidmaps': None,
        },
    }
    
    urls = {
        'kegg': lambda x: f"https://www.genome.jp/dbget-bin/www_bget?cpd:{x}",
        'pubchem': lambda x: f"https://pubchem.ncbi.nlm.nih.gov/compound/{x}",
        'chebi': lambda x: f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:{x}",
        'metacyc': lambda x: f"https://biocyc.org/compound?orgid=META&id={x}",
    }
    
    if hmdb_id in external_db:
        refs = external_db[hmdb_id]
        links = {}
        
        for db, db_id in refs.items():
            if db_id and db in urls:
                links[db] = urls[db](db_id)
        
        return links
    
    return {}


def main():
    print("=" * 70)
    print("HMDB 代谢物识别和通路分析")
    print("=" * 70)
    
    # 示例 1: 代谢物识别
    print("\n示例 1: 基于质谱数据的代谢物识别")
    print("-" * 70)
    
    # 模拟实验测得的 m/z 值
    experimental_mz = [
        180.0634,  # Glucose
        169.0748,  # 1-Methylhistidine
    ]
    
    for mz in experimental_mz:
        candidates = identify_metabolite_from_ms(mz, ppm_tolerance=10)
        
        print(f"\nm/z = {mz} 的候选代谢物:")
        for i, candidate in enumerate(candidates[:5], 1):
            print(f"  {i}. {candidate['name']}")
            if 'formula' in candidate:
                print(f"     分子式: {candidate['formula']}")
            if 'ppm_diff' in candidate:
                print(f"     ppm 差异: {candidate['ppm_diff']:.2f}")
    
    # 示例 2: 通路分析
    print("\n" + "=" * 70)
    print("示例 2: 代谢通路富集分析")
    print("=" * 70)
    
    # 差异代谢物列表
    diff_metabolites = ['HMDB0000001', 'HMDB0000122', 'HMDB0000159']
    
    pathway_results = pathway_analysis(diff_metabolites)
    
    print(f"\n差异代谢物: {len(diff_metabolites)} 个")
    print("\n富集通路:")
    for result in pathway_results:
        print(f"  - {result['pathway_name']}")
        print(f"    命中: {result['hits']}/{result['pathway_size']}")
        print(f"    富集度: {result['enrichment']:.2%}")
    
    # 示例 3: 生物标志物搜索
    print("\n" + "=" * 70)
    print("示例 3: 疾病生物标志物搜索")
    print("=" * 70)
    
    diseases = ['prostate cancer', 'diabetes', 'alzheimer']
    
    for disease in diseases:
        biomarkers = biomarker_search(disease)
        
        if biomarkers:
            print(f"\n{disease.title()} 的生物标志物:")
            for marker in biomarkers:
                print(f"  - {marker['name']} ({marker['hmdb_id']})")
                print(f"    变化: {marker['change'].upper()}")
                if marker['fold_change']:
                    print(f"    倍数变化: {marker['fold_change']:.2f}x")
    
    # 示例 4: 数据库交叉引用
    print("\n" + "=" * 70)
    print("示例 4: 外部数据库交叉引用")
    print("=" * 70)
    
    test_ids = ['HMDB0000122', 'HMDB0000001']
    
    for hmdb_id in test_ids:
        links = cross_reference_databases(hmdb_id)
        
        print(f"\n{hmdb_id}:")
        for db, url in links.items():
            print(f"  {db.upper()}: {url}")
    
    # 示例 5: 数据导出
    print("\n" + "=" * 70)
    print("示例 5: 分析结果导出")
    print("=" * 70)
    
    # 保存通路分析结果
    output = {
        'analysis_type': 'pathway_enrichment',
        'input_metabolites': diff_metabolites,
        'results': pathway_results,
        'timestamp': '2024-01-01T00:00:00'
    }
    
    output_file = 'pathway_analysis_results.json'
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n✓ 结果已保存到 {output_file}")
    
    print("\n" + "=" * 70)
    print("提示:")
    print("-" * 70)
    print("1. 下载 HMDB CSV 文件以获得完整功能")
    print("2. 使用 MSP 或 MGF 格式进行光谱匹配")
    print("3. 结合 KEGG 进行更全面的通路分析")
    print("4. 使用 MetaboAnalyst 进行在线代谢组学分析")


if __name__ == "__main__":
    main()
