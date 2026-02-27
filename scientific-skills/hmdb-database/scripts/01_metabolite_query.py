#!/usr/bin/env python3
"""
HMDB 示例脚本 1: 代谢物查询和数据提取
演示如何查询和解析 HMDB 数据
"""

import requests
import xml.etree.ElementTree as ET
from pathlib import Path
import json


def search_metabolite_by_name(name):
    """
    通过名称搜索代谢物 (使用 HMDB 网站搜索)
    
    Args:
        name: 代谢物名称
    
    Returns:
        搜索结果
    """
    # HMDB 网站搜索 URL
    search_url = f"https://www.hmdb.ca/unearth/q?utf8=%E2%9C%93&query={name}&searcher=metabolites&button="
    
    print(f"搜索代谢物: {name}")
    print(f"搜索 URL: {search_url}")
    
    # 注意: 实际搜索需要网页解析或 API 访问
    # 这里展示搜索概念
    
    return search_url


def parse_hmdb_xml(xml_file):
    """
    解析 HMDB XML 数据文件
    
    Args:
        xml_file: HMDB XML 文件路径
    
    Returns:
        解析后的代谢物数据列表
    """
    print(f"\n解析 XML 文件: {xml_file}")
    
    if not Path(xml_file).exists():
        print(f"  文件不存在，跳过解析")
        return []
    
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        
        # 定义命名空间
        ns = {'hmdb': 'http://www.hmdb.ca'}
        
        metabolites = []
        
        # 解析每个代谢物
        for metabolite in root.findall('.//hmdb:metabolite', ns):
            data = {
                'accession': get_xml_text(metabolite, 'hmdb:accession', ns),
                'name': get_xml_text(metabolite, 'hmdb:name', ns),
                'formula': get_xml_text(metabolite, 'hmdb:chemical_formula', ns),
                'mw': get_xml_text(metabolite, 'hmdb:average_molecular_weight', ns),
                'smiles': get_xml_text(metabolite, 'hmdb:smiles', ns),
                'inchi': get_xml_text(metabolite, 'hmdb:inchi', ns),
                'inchikey': get_xml_text(metabolite, 'hmdb:inchikey', ns),
            }
            metabolites.append(data)
        
        return metabolites
        
    except Exception as e:
        print(f"  解析错误: {e}")
        return []


def get_xml_text(element, tag, ns):
    """安全地获取 XML 元素文本"""
    elem = element.find(tag, ns)
    return elem.text if elem is not None else None


def download_hmdb_dataset(output_dir="./hmdb_data"):
    """
    下载 HMDB 数据集
    
    Args:
        output_dir: 输出目录
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    # HMDB 下载链接
    downloads = {
        'metabolites_xml': 'https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
        'metabolites_sdf': 'https://hmdb.ca/system/downloads/current/structures.zip',
        'proteins_xml': 'https://hmdb.ca/system/downloads/current/hmdb_proteins.zip',
    }
    
    print("\nHMDB 数据集下载链接:")
    print("-" * 60)
    for name, url in downloads.items():
        print(f"{name}:")
        print(f"  {url}")
        print(f"  保存到: {output_dir}/{name}.zip")
    
    print("\n下载命令示例:")
    print("-" * 60)
    print("""
# 使用 wget 下载
wget https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip -P ./hmdb_data/

# 或使用 Python
import urllib.request
urllib.request.urlretrieve(
    'https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip',
    './hmdb_data/hmdb_metabolites.zip'
)

# 解压
unzip ./hmdb_data/hmdb_metabolites.zip -d ./hmdb_data/
""")


def query_metabolite_web(metabolite_id):
    """
    通过网页查询单个代谢物信息
    
    Args:
        metabolite_id: HMDB ID (例如: HMDB0000001)
    
    Returns:
        代谢物信息 URL
    """
    url = f"https://hmdb.ca/metabolites/{metabolite_id}"
    
    print(f"\n代谢物页面: {url}")
    
    # 尝试获取页面
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            print(f"  ✓ 页面可访问")
            return url
        else:
            print(f"  × 页面访问失败 (状态码: {response.status_code})")
    except Exception as e:
        print(f"  × 请求错误: {e}")
    
    return None


def demonstrate_metabolite_fields():
    """展示 HMDB 代谢物数据字段"""
    print("\n" + "=" * 70)
    print("HMDB 代谢物数据字段示例")
    print("=" * 70)
    
    # 示例代谢物数据结构
    example_metabolite = {
        'accession': 'HMDB0000001',
        'name': '1-Methylhistidine',
        'description': '1-Methylhistidine is a...',
        'chemistry': {
            'formula': 'C7H9N3O2',
            'average_mass': '167.164',
            'monisotopic_mass': '167.0695',
            'iupac_name': '2-amino-3-(1-methylimidazol-4-yl)propanoic acid',
            'smiles': 'CN1C=NC(C=C1)CC(C(=O)O)N',
            'inchi': 'InChI=1S/C7H9N3O2/c1-10-4-9-3-5(10)2-6(8)7(11)12/h3-4,6H,2,8H2,1H3,(H,11,12)',
            'inchikey': 'BRMWTNUJHUMWMS-UHFFFAOYSA-N',
        },
        'biological': {
            'biospecimen_locations': ['Blood', 'Urine', 'Saliva'],
            'tissue_locations': ['Liver', 'Kidney', 'Muscle'],
            'pathways': [
                {'name': 'Histidine Metabolism', 'smpdb_id': 'SMP00044'},
                {'name': 'Amino Acid Metabolism', 'smpdb_id': 'SMP00001'},
            ],
        },
        'clinical': {
            'normal_concentrations': [
                {'biospecimen': 'Blood', 'value': '0.5-2.0', 'units': 'uM'},
                {'biospecimen': 'Urine', 'value': '50-200', 'units': 'uM'},
            ],
            'disease_associations': [
                {'disease': 'Prostate Cancer', 'pmid': '12345678'},
                {'disease': 'Muscle Atrophy', 'pmid': '23456789'},
            ],
        },
        'spectra': {
            'nmr': [{'type': '1H NMR', 'solvent': 'Water'}],
            'ms': [{'type': 'LC-MS', 'ion_mode': 'Positive'}],
        },
        'external_ids': {
            'kegg': 'C01152',
            'pubchem': '92865',
            'chebi': '27596',
            'metacyc': 'CPD-12345',
        }
    }
    
    print(f"\n代谢物: {example_metabolite['name']}")
    print(f"HMDB ID: {example_metabolite['accession']}")
    
    print(f"\n化学信息:")
    for key, value in example_metabolite['chemistry'].items():
        print(f"  {key}: {value}")
    
    print(f"\n生物学位置:")
    print(f"  生物样本: {', '.join(example_metabolite['biological']['biospecimen_locations'])}")
    print(f"  组织: {', '.join(example_metabolite['biological']['tissue_locations'])}")
    
    print(f"\n临床关联:")
    for disease in example_metabolite['clinical']['disease_associations']:
        print(f"  - {disease['disease']} (PMID: {disease['pmid']})")
    
    print(f"\n外部数据库链接:")
    for db, id in example_metabolite['external_ids'].items():
        print(f"  {db.upper()}: {id}")


def main():
    print("=" * 70)
    print("HMDB (Human Metabolome Database) 查询示例")
    print("=" * 70)
    
    # 示例 1: 搜索代谢物
    print("\n示例 1: 代谢物搜索")
    print("-" * 70)
    search_metabolite_by_name("glucose")
    search_metabolite_by_name("caffeine")
    
    # 示例 2: 查询特定代谢物
    print("\n示例 2: 查询特定代谢物")
    print("-" * 70)
    query_metabolite_web("HMDB0000122")  # Glucose
    query_metabolite_web("HMDB0001847")  # Caffeine
    
    # 示例 3: 数据字段示例
    demonstrate_metabolite_fields()
    
    # 示例 4: 数据下载
    print("\n" + "=" * 70)
    print("示例 4: 批量数据下载")
    print("=" * 70)
    download_hmdb_dataset()
    
    # 示例 5: XML 解析示例
    print("\n" + "=" * 70)
    print("示例 5: XML 数据解析")
    print("=" * 70)
    
    xml_example = """
# 下载并解析 HMDB XML 文件的示例

from lxml import etree

def parse_hmdb_metabolites(xml_file):
    context = etree.iterparse(xml_file, tag='{http://www.hmdb.ca}metabolite')
    
    for event, elem in context:
        accession = elem.findtext('{http://www.hmdb.ca}accession')
        name = elem.findtext('{http://www.hmdb.ca}name')
        formula = elem.findtext('{http://www.hmdb.ca}chemical_formula')
        
        print(f"{accession}: {name} ({formula})")
        
        # 清理内存
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

# 使用示例
# parse_hmdb_metabolites('hmdb_metabolites.xml')
"""
    print(xml_example)
    
    print("\n" + "=" * 70)
    print("提示:")
    print("-" * 70)
    print("1. 访问 https://www.hmdb.ca/ 进行交互式搜索")
    print("2. 下载完整数据集: https://www.hmdb.ca/downloads")
    print("3. 引用 HMDB: Wishart et al. (2022) Nucleic Acids Research")
    print("4. 商业使用需要联系: samackay@ualberta.ca")


if __name__ == "__main__":
    main()
