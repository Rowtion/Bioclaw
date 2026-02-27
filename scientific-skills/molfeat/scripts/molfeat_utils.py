#!/usr/bin/env python3
"""
Molfeat 分子特征化工具
用于计算分子指纹和描述符
"""

import numpy as np
from typing import List, Union, Optional, Dict, Any, Tuple
from dataclasses import dataclass
from pathlib import Path
import json

# 尝试导入 molfeat
try:
    from molfeat.calc import FPCalculator, RDKitDescriptors2D
    from molfeat.trans import MoleculeTransformer, FeatConcat
    from molfeat.trans.pretrained import PretrainedMolTransformer
    from molfeat.store.modelstore import ModelStore
    from molfeat.calc import MordredDescriptors
    MOLFEAT_AVAILABLE = True
except ImportError:
    MOLFEAT_AVAILABLE = False
    print("警告: molfeat 未安装。运行: uv pip install molfeat")

try:
    import datamol as dm
    DATAMOL_AVAILABLE = True
except ImportError:
    DATAMOL_AVAILABLE = False


@dataclass
class FeaturizationResult:
    """特征化结果"""
    smiles: str
    features: np.ndarray
    featurizer_name: str
    success: bool
    error_message: Optional[str] = None
    
    def to_dict(self) -> Dict:
        return {
            'smiles': self.smiles,
            'features': self.features.tolist() if self.success else None,
            'featurizer_name': self.featurizer_name,
            'success': self.success,
            'error_message': self.error_message
        }


class MolecularFeaturizer:
    """分子特征化器"""
    
    # 可用的指纹类型
    FP_TYPES = {
        'ecfp': 'Extended Connectivity Fingerprint (Morgan)',
        'fcfp': 'Functional Class Fingerprint',
        'maccs': 'MACCS Keys (167 bits)',
        'map4': 'MinHash Atom Pair Fingerprint',
        'erg': 'Extended Reduced Graph',
        'gobbi2D': 'Gobbi Pharmacophore',
        'gobbi3D': 'Gobbi 3D Pharmacophore',
        'atompair': 'Atom Pair Fingerprint',
        'torsion': 'Topological Torsion',
        'rdkit': 'RDKit Topological Fingerprint',
        'pattern': 'Pattern Fingerprint',
        'layered': 'Layered Fingerprint',
        'secfp': 'SMILES Extended Connectivity',
    }
    
    # 预训练模型
    PRETRAINED_MODELS = {
        'chemberta': 'ChemBERTa-77M-MLM',
        'chemberta_mtr': 'ChemBERTa-77M-MTR',
        'chemgpt': 'ChemGPT-1.2B',
        'molt5': 'MolT5',
        'gin_masking': 'gin-supervised-masking',
        'gin_infomax': 'gin-supervised-infomax',
        'gin_edgepred': 'gin-supervised-edgepred',
        'gin_contextpred': 'gin-supervised-contextpred',
        'graphormer': 'Graphormer-pcqm4mv2',
    }
    
    def __init__(self, featurizer_type: str = 'ecfp', **kwargs):
        """
        初始化分子特征化器
        
        Args:
            featurizer_type: 特征化器类型
            **kwargs: 特征化器参数
        """
        if not MOLFEAT_AVAILABLE:
            raise ImportError("molfeat 未安装。运行: uv pip install molfeat")
        
        self.featurizer_type = featurizer_type
        self.kwargs = kwargs
        self.transformer = None
        self.calculator = None
        self._build_featurizer()
    
    def _build_featurizer(self):
        """构建特征化器"""
        if self.featurizer_type in ['ecfp', 'fcfp', 'maccs', 'map4', 'erg', 
                                     'gobbi2D', 'gobbi3D', 'atompair', 'torsion',
                                     'rdkit', 'pattern', 'layered', 'secfp']:
            # 指纹计算器
            radius = self.kwargs.get('radius', 3)
            fp_size = self.kwargs.get('fp_size', 2048)
            
            self.calculator = FPCalculator(
                self.featurizer_type,
                radius=radius,
                fpSize=fp_size
            )
            self.transformer = MoleculeTransformer(self.calculator, n_jobs=-1)
            
        elif self.featurizer_type == 'rdkit2d':
            # RDKit 2D 描述符
            self.calculator = RDKitDescriptors2D()
            self.transformer = MoleculeTransformer(self.calculator, n_jobs=-1)
            
        elif self.featurizer_type == 'mordred':
            # Mordred 描述符
            self.calculator = MordredDescriptors()
            self.transformer = MoleculeTransformer(self.calculator, n_jobs=-1)
            
        elif self.featurizer_type in self.PRETRAINED_MODELS.values():
            # 预训练模型
            self.transformer = PretrainedMolTransformer(
                self.featurizer_type,
                n_jobs=-1
            )
        else:
            raise ValueError(f"未知的特征化器类型: {self.featurizer_type}")
    
    def featurize(self, molecules: Union[str, List[str]], 
                  ignore_errors: bool = True) -> Union[FeaturizationResult, List[FeaturizationResult]]:
        """
        计算分子特征
        
        Args:
            molecules: SMILES 字符串或列表
            ignore_errors: 是否忽略错误
        
        Returns:
            特征化结果
        """
        single_input = isinstance(molecules, str)
        if single_input:
            molecules = [molecules]
        
        results = []
        
        # 批量特征化
        try:
            features_list = self.transformer(molecules, ignore_errors=ignore_errors)
            
            for smiles, features in zip(molecules, features_list):
                if features is not None:
                    result = FeaturizationResult(
                        smiles=smiles,
                        features=np.array(features),
                        featurizer_name=self.featurizer_type,
                        success=True
                    )
                else:
                    result = FeaturizationResult(
                        smiles=smiles,
                        features=None,
                        featurizer_name=self.featurizer_type,
                        success=False,
                        error_message="特征化失败"
                    )
                results.append(result)
                
        except Exception as e:
            # 单个处理
            for smiles in molecules:
                try:
                    features = self.transformer([smiles])[0]
                    result = FeaturizationResult(
                        smiles=smiles,
                        features=np.array(features) if features is not None else None,
                        featurizer_name=self.featurizer_type,
                        success=features is not None
                    )
                except Exception as ex:
                    result = FeaturizationResult(
                        smiles=smiles,
                        features=None,
                        featurizer_name=self.featurizer_type,
                        success=False,
                        error_message=str(ex)
                    )
                results.append(result)
        
        return results[0] if single_input else results
    
    def featurize_numpy(self, molecules: List[str], 
                        ignore_errors: bool = True) -> np.ndarray:
        """
        计算分子特征并返回 numpy 数组
        
        Args:
            molecules: SMILES 列表
            ignore_errors: 是否忽略错误
        
        Returns:
            特征矩阵
        """
        results = self.featurize(molecules, ignore_errors=ignore_errors)
        
        if isinstance(results, FeaturizationResult):
            results = [results]
        
        # 过滤成功的结果
        valid_features = [r.features for r in results if r.success]
        
        if not valid_features:
            return np.array([])
        
        return np.vstack(valid_features)
    
    def get_feature_names(self) -> Optional[List[str]]:
        """
        获取特征名称列表
        
        Returns:
            特征名称列表
        """
        if hasattr(self.calculator, 'columns'):
            return self.calculator.columns
        return None
    
    def save_config(self, filepath: str):
        """
        保存配置
        
        Args:
            filepath: 配置文件路径
        """
        config = {
            'featurizer_type': self.featurizer_type,
            'kwargs': self.kwargs
        }
        with open(filepath, 'w') as f:
            json.dump(config, f, indent=2)
    
    @classmethod
    def load_config(cls, filepath: str) -> 'MolecularFeaturizer':
        """
        从配置文件加载
        
        Args:
            filepath: 配置文件路径
        
        Returns:
            MolecularFeaturizer 实例
        """
        with open(filepath, 'r') as f:
            config = json.load(f)
        return cls(config['featurizer_type'], **config['kwargs'])


class MultiFeaturizer:
    """多特征化器组合"""
    
    def __init__(self, featurizers: Dict[str, MolecularFeaturizer]):
        """
        初始化多特征化器
        
        Args:
            featurizers: 特征化器字典 {名称: 特征化器}
        """
        self.featurizers = featurizers
    
    def featurize(self, molecules: List[str]) -> Dict[str, np.ndarray]:
        """
        使用多个特征化器计算特征
        
        Args:
            molecules: SMILES 列表
        
        Returns:
            特征字典 {名称: 特征矩阵}
        """
        results = {}
        for name, featurizer in self.featurizers.items():
            results[name] = featurizer.featurize_numpy(molecules)
        return results
    
    def featurize_concatenated(self, molecules: List[str]) -> np.ndarray:
        """
        连接所有特征化器的输出
        
        Args:
            molecules: SMILES 列表
        
        Returns:
            连接后的特征矩阵
        """
        features_list = []
        for featurizer in self.featurizers.values():
            features = featurizer.featurize_numpy(molecules)
            if features.size > 0:
                features_list.append(features)
        
        if not features_list:
            return np.array([])
        
        return np.hstack(features_list)


def compare_featurizers(smiles_list: List[str], 
                       featurizer_types: List[str] = None) -> Dict[str, Dict]:
    """
    比较不同特征化器的性能
    
    Args:
        smiles_list: SMILES 列表
        featurizer_types: 要比较的特征化器类型
    
    Returns:
        比较结果
    """
    if featurizer_types is None:
        featurizer_types = ['ecfp', 'maccs', 'rdkit2d']
    
    results = {}
    
    for ft in featurizer_types:
        try:
            featurizer = MolecularFeaturizer(ft)
            import time
            start = time.time()
            features = featurizer.featurize_numpy(smiles_list)
            elapsed = time.time() - start
            
            results[ft] = {
                'shape': features.shape,
                'time_seconds': elapsed,
                'success_rate': features.shape[0] / len(smiles_list),
                'error': None
            }
        except Exception as e:
            results[ft] = {
                'shape': None,
                'time_seconds': None,
                'success_rate': 0,
                'error': str(e)
            }
    
    return results


def list_available_featurizers() -> Dict[str, str]:
    """
    列出所有可用的特征化器
    
    Returns:
        特征化器名称和描述
    """
    return {
        **MolecularFeaturizer.FP_TYPES,
        'rdkit2d': 'RDKit 2D Descriptors (200+ properties)',
        'mordred': 'Mordred Descriptors (1800+ properties)',
        **{k: f"Pretrained: {v}" for k, v in MolecularFeaturizer.PRETRAINED_MODELS.items()}
    }


def demo():
    """演示功能"""
    print("=" * 60)
    print("Molfeat 分子特征化演示")
    print("=" * 60)
    
    if not MOLFEAT_AVAILABLE:
        print("错误: molfeat 未安装")
        print("请运行: uv pip install molfeat")
        return
    
    # 示例分子
    smiles_list = [
        "CCO",  # 乙醇
        "CC(=O)O",  # 乙酸
        "c1ccccc1",  # 苯
        "CC(C)O",  # 异丙醇
        "C1CCCCC1",  # 环己烷
    ]
    
    print("\n示例分子:")
    for i, smi in enumerate(smiles_list, 1):
        print(f"  {i}. {smi}")
    
    # 1. ECFP 指纹
    print("\n1. ECFP 指纹 (radius=3, 2048 bits):")
    featurizer = MolecularFeaturizer('ecfp', radius=3, fp_size=2048)
    results = featurizer.featurize(smiles_list)
    print(f"   成功: {sum(r.success for r in results)}/{len(results)}")
    print(f"   特征维度: {results[0].features.shape if results[0].success else 'N/A'}")
    
    # 2. MACCS 指纹
    print("\n2. MACCS 指纹 (167 bits):")
    featurizer_maccs = MolecularFeaturizer('maccs')
    results_maccs = featurizer_maccs.featurize(smiles_list)
    print(f"   特征维度: {results_maccs[0].features.shape if results_maccs[0].success else 'N/A'}")
    
    # 3. 比较不同特征化器
    print("\n3. 特征化器比较:")
    comparison = compare_featurizers(smiles_list[:3], ['ecfp', 'maccs'])
    for name, result in comparison.items():
        print(f"   {name}: shape={result['shape']}, time={result['time_seconds']:.3f}s")
    
    # 4. 列出可用特征化器
    print("\n4. 可用特征化器类型:")
    available = list_available_featurizers()
    for name, desc in list(available.items())[:5]:
        print(f"   - {name}: {desc}")
    print(f"   ... (共 {len(available)} 种)")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
