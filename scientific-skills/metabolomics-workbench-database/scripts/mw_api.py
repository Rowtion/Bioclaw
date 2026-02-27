#!/usr/bin/env python3
"""
Metabolomics Workbench 数据库查询工具
NIH 代谢组学数据库 REST API 封装
"""

import requests
import json
from typing import Optional, Dict, List, Union
from dataclasses import dataclass
from urllib.parse import quote


BASE_URL = "https://www.metabolomicsworkbench.org/rest"


@dataclass
class Compound:
    """代谢物化合物数据类"""
    regno: str
    name: str
    formula: Optional[str] = None
    exact_mass: Optional[float] = None
    pubchem_cid: Optional[str] = None
    inchi_key: Optional[str] = None
    kegg_id: Optional[str] = None
    hmdb_id: Optional[str] = None
    
    def to_dict(self) -> Dict:
        return {
            'regno': self.regno,
            'name': self.name,
            'formula': self.formula,
            'exact_mass': self.exact_mass,
            'pubchem_cid': self.pubchem_cid,
            'inchi_key': self.inchi_key,
            'kegg_id': self.kegg_id,
            'hmdb_id': self.hmdb_id
        }


class MetabolomicsWorkbenchAPI:
    """Metabolomics Workbench REST API 客户端"""
    
    def __init__(self):
        self.base_url = BASE_URL
        self.session = requests.Session()
    
    def _get(self, endpoint: str, format: str = "json") -> Dict:
        """
        发送 GET 请求
        
        Args:
            endpoint: API 端点
            format: 返回格式 (json 或 txt)
        
        Returns:
            API 响应数据
        """
        url = f"{self.base_url}/{endpoint}/{format}"
        response = self.session.get(url)
        response.raise_for_status()
        
        if format == "json":
            return response.json()
        return response.text
    
    def get_compound_by_pubchem(self, pubchem_cid: str) -> Dict:
        """
        通过 PubChem CID 获取化合物信息
        
        Args:
            pubchem_cid: PubChem CID
        
        Returns:
            化合物详细信息
        """
        return self._get(f"compound/pubchem_cid/{pubchem_cid}/all")
    
    def get_compound_by_regno(self, regno: str) -> Dict:
        """
        通过注册号获取化合物信息
        
        Args:
            regno: Metabolomics Workbench 注册号
        
        Returns:
            化合物信息
        """
        return self._get(f"compound/regno/{regno}/all")
    
    def get_compound_structure_png(self, regno: str) -> bytes:
        """
        获取化合物结构 PNG 图像
        
        Args:
            regno: 注册号
        
        Returns:
            PNG 图像字节数据
        """
        url = f"{self.base_url}/compound/regno/{regno}/png"
        response = self.session.get(url)
        response.raise_for_status()
        return response.content
    
    def save_compound_structure(self, regno: str, output_path: str):
        """
        保存化合物结构图像
        
        Args:
            regno: 注册号
            output_path: 输出文件路径
        """
        png_data = self.get_compound_structure_png(regno)
        with open(output_path, 'wb') as f:
            f.write(png_data)
        print(f"结构已保存到: {output_path}")
    
    def search_studies_by_metabolite(self, refmet_name: str) -> List[Dict]:
        """
        通过代谢物名称搜索相关研究
        
        Args:
            refmet_name: RefMet 标准化名称
        
        Returns:
            研究列表
        """
        encoded_name = quote(refmet_name)
        result = self._get(f"study/refmet_name/{encoded_name}/summary")
        
        # 处理不同格式的响应
        if isinstance(result, dict):
            studies = result.get('study', [])
            if not isinstance(studies, list):
                studies = [studies]
            return studies
        return []
    
    def get_study_summary(self, study_id: str) -> Dict:
        """
        获取研究摘要信息
        
        Args:
            study_id: 研究 ID (如 ST000001)
        
        Returns:
            研究摘要
        """
        return self._get(f"study/study_id/{study_id}/summary")
    
    def get_study_data(self, study_id: str) -> Dict:
        """
        获取研究实验数据
        
        Args:
            study_id: 研究 ID
        
        Returns:
            实验数据
        """
        return self._get(f"study/study_id/{study_id}/data")
    
    def list_available_studies(self) -> List[str]:
        """
        列出所有可用公共研究
        
        Returns:
            研究 ID 列表
        """
        result = self._get("study/study_id/ST/available")
        studies = result.get('study_id', [])
        if isinstance(studies, str):
            return [studies]
        return studies if studies else []
    
    def refmet_match(self, metabolite_name: str) -> Dict:
        """
        匹配代谢物名称到 RefMet 标准名称
        
        Args:
            metabolite_name: 代谢物名称
        
        Returns:
            RefMet 匹配结果
        """
        encoded_name = quote(metabolite_name)
        return self._get(f"refmet/match/{encoded_name}/name")
    
    def search_by_formula(self, formula: str) -> List[Dict]:
        """
        通过分子式搜索代谢物
        
        Args:
            formula: 分子式 (如 C6H12O6)
        
        Returns:
            匹配代谢物列表
        """
        result = self._get(f"refmet/formula/{formula}/all")
        metabolites = result.get('RefMet', [])
        if not isinstance(metabolites, list):
            metabolites = [metabolites]
        return metabolites
    
    def search_by_exact_mass(self, mass: float, tolerance: float = 0.5) -> List[Dict]:
        """
        通过精确质量搜索代谢物
        
        Args:
            mass: 目标质量
            tolerance: 质量容差 (Da)
        
        Returns:
            匹配代谢物列表
        """
        result = self._get(f"refmet/exact_mass/{mass}/{tolerance}")
        metabolites = result.get('RefMet', [])
        if not isinstance(metabolites, list):
            metabolites = [metabolites]
        return metabolites
    
    def search_mz(self, mz: float, adduct: str = "M+H", tolerance: float = 0.5, 
                  database: str = "MB") -> List[Dict]:
        """
        通过 m/z 值搜索化合物
        
        Args:
            mz: m/z 值
            adduct: 离子加合物类型 (M+H, M-H, M+Na, M+NH4, M+2H 等)
            tolerance: 质量容差
            database: 搜索数据库 (MB, LIPIDS, REFMET)
        
        Returns:
            匹配化合物列表
        """
        result = self._get(f"moverz/{database}/{mz}/{adduct}/{tolerance}")
        compounds = result.get('compound', [])
        if not isinstance(compounds, list):
            compounds = [compounds]
        return compounds
    
    def calculate_exact_mass(self, metabolite_name: str, adduct: str = "M+H") -> Optional[float]:
        """
        计算代谢物的精确质量
        
        Args:
            metabolite_name: 代谢物名称
            adduct: 离子加合物
        
        Returns:
            精确质量或 None
        """
        encoded_name = quote(metabolite_name)
        result = self._get(f"moverz/exactmass/{encoded_name}/{adduct}")
        
        # 解析质量值
        if isinstance(result, dict) and 'mass' in result:
            try:
                return float(result['mass'])
            except (ValueError, TypeError):
                pass
        return None
    
    def metstat_search(self, analytical_method: str = "", polarity: str = "",
                       chromatography: str = "", species: str = "",
                       sample_source: str = "", disease: str = "") -> List[Dict]:
        """
        MetStat 高级搜索
        
        Args:
            analytical_method: 分析方法 (LCMS, GCMS, NMR)
            polarity: 离子化极性 (POSITIVE, NEGATIVE)
            chromatography: 色谱类型 (HILIC, RP, GC)
            species: 物种
            sample_source: 样本来源
            disease: 疾病
        
        Returns:
            匹配研究列表
        """
        # 构建查询字符串
        params = [analytical_method, polarity, chromatography, species, sample_source, "", disease]
        query = ";".join(params)
        
        result = self._get(f"metstat/{query}")
        studies = result.get('study', [])
        if not isinstance(studies, list):
            studies = [studies]
        return studies
    
    def get_gene_info(self, gene_symbol: str) -> Dict:
        """
        获取基因信息
        
        Args:
            gene_symbol: 基因符号
        
        Returns:
            基因信息
        """
        return self._get(f"gene/gene_symbol/{gene_symbol}/all")
    
    def get_protein_info(self, uniprot_id: str) -> Dict:
        """
        获取蛋白质信息
        
        Args:
            uniprot_id: UniProt ID
        
        Returns:
            蛋白质信息
        """
        return self._get(f"protein/uniprot_id/{uniprot_id}/all")


class MetabolomicsAnalyzer:
    """代谢组学分析工具"""
    
    def __init__(self):
        self.api = MetabolomicsWorkbenchAPI()
    
    def identify_compound_from_ms(self, mz: float, adduct: str = "M+H", 
                                   tolerance: float = 0.5) -> Dict:
        """
        从质谱数据鉴定化合物
        
        Args:
            mz: 观测到的 m/z 值
            adduct: 离子加合物
            tolerance: 质量容差
        
        Returns:
            包含候选化合物和分析结果的字典
        """
        print(f"搜索 m/z = {mz} ± {tolerance} Da ({adduct})")
        
        # 在多个数据库中搜索
        candidates = []
        for db in ["MB", "REFMET", "LIPIDS"]:
            try:
                results = self.api.search_mz(mz, adduct, tolerance, db)
                for r in results:
                    r['database'] = db
                candidates.extend(results)
            except Exception as e:
                print(f"  {db} 数据库搜索失败: {e}")
        
        # 去重
        seen = set()
        unique_candidates = []
        for c in candidates:
            key = c.get('regno', '') + c.get('name', '')
            if key and key not in seen:
                seen.add(key)
                unique_candidates.append(c)
        
        return {
            'query_mz': mz,
            'adduct': adduct,
            'tolerance': tolerance,
            'candidate_count': len(unique_candidates),
            'candidates': unique_candidates
        }
    
    def find_studies_for_metabolite(self, metabolite_name: str) -> Dict:
        """
        查找包含特定代谢物的研究
        
        Args:
            metabolite_name: 代谢物名称
        
        Returns:
            包含研究列表的字典
        """
        print(f"查找代谢物: {metabolite_name}")
        
        # 首先标准化名称
        refmet_match = self.api.refmet_match(metabolite_name)
        standardized_name = refmet_match.get('name', metabolite_name)
        
        if standardized_name != metabolite_name:
            print(f"  标准化名称: {standardized_name}")
        
        # 搜索研究
        studies = self.api.search_studies_by_metabolite(standardized_name)
        
        return {
            'query_name': metabolite_name,
            'standardized_name': standardized_name,
            'study_count': len(studies),
            'studies': studies
        }
    
    def batch_search_formulas(self, formulas: List[str]) -> Dict[str, List[Dict]]:
        """
        批量搜索分子式
        
        Args:
            formulas: 分子式列表
        
        Returns:
            每个分子式的搜索结果
        """
        results = {}
        for formula in formulas:
            print(f"搜索分子式: {formula}")
            results[formula] = self.api.search_by_formula(formula)
        return results


def demo():
    """演示功能"""
    print("=" * 60)
    print("Metabolomics Workbench API 演示")
    print("=" * 60)
    
    api = MetabolomicsWorkbenchAPI()
    analyzer = MetabolomicsAnalyzer()
    
    # 1. 搜索化合物
    print("\n1. 搜索葡萄糖相关研究:")
    studies = api.search_studies_by_metabolite("Glucose")
    print(f"   找到 {len(studies)} 个研究")
    if studies:
        print(f"   示例研究: {studies[0].get('study_id', 'N/A')}")
    
    # 2. RefMet 名称匹配
    print("\n2. RefMet 名称匹配:")
    match = api.refmet_match("citrate")
    print(f"   输入: citrate")
    print(f"   匹配: {match.get('name', 'N/A')}")
    
    # 3. 分子式搜索
    print("\n3. 搜索分子式 C6H12O6:")
    metabolites = api.search_by_formula("C6H12O6")
    print(f"   找到 {len(metabolites)} 个代谢物")
    for m in metabolites[:3]:
        print(f"   - {m.get('name', 'N/A')} ({m.get('formula', 'N/A')})")
    
    # 4. m/z 搜索
    print("\n4. m/z 搜索 (180.06, M+H):")
    candidates = api.search_mz(180.06, "M+H", 0.5, "MB")
    print(f"   找到 {len(candidates)} 个候选化合物")
    for c in candidates[:3]:
        print(f"   - {c.get('name', 'N/A')}")
    
    # 5. MetStat 搜索
    print("\n5. MetStat 搜索 (人类血液研究):")
    studies = api.metstat_search(species="Human", sample_source="Blood")
    print(f"   找到 {len(studies)} 个研究")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
