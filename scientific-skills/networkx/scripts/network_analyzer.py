#!/usr/bin/env python3
"""
NetworkX 网络分析工具
用于复杂网络分析和可视化
"""

import networkx as nx
import numpy as np
from typing import List, Dict, Tuple, Optional, Union, Any
from dataclasses import dataclass
from pathlib import Path
import json

# 尝试导入可视化库
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False


@dataclass
class NetworkMetrics:
    """网络度量指标"""
    num_nodes: int
    num_edges: int
    density: float
    avg_clustering: float
    avg_shortest_path: Optional[float]
    diameter: Optional[int]
    is_connected: bool
    
    def to_dict(self) -> Dict:
        return {
            'num_nodes': self.num_nodes,
            'num_edges': self.num_edges,
            'density': self.density,
            'avg_clustering': self.avg_clustering,
            'avg_shortest_path': self.avg_shortest_path,
            'diameter': self.diameter,
            'is_connected': self.is_connected
        }


class NetworkAnalyzer:
    """网络分析器"""
    
    def __init__(self, graph: nx.Graph):
        """
        初始化网络分析器
        
        Args:
            graph: NetworkX 图对象
        """
        self.graph = graph
    
    def get_basic_metrics(self) -> NetworkMetrics:
        """
        获取基本网络度量
        
        Returns:
            NetworkMetrics 对象
        """
        is_connected = nx.is_connected(self.graph) if not self.graph.is_directed() else False
        
        avg_path = None
        diameter = None
        if is_connected:
            try:
                avg_path = nx.average_shortest_path_length(self.graph)
                diameter = nx.diameter(self.graph)
            except:
                pass
        
        return NetworkMetrics(
            num_nodes=self.graph.number_of_nodes(),
            num_edges=self.graph.number_of_edges(),
            density=nx.density(self.graph),
            avg_clustering=nx.average_clustering(self.graph),
            avg_shortest_path=avg_path,
            diameter=diameter,
            is_connected=is_connected
        )
    
    def get_centrality_measures(self) -> Dict[str, Dict]:
        """
        计算中心性指标
        
        Returns:
            中心性指标字典
        """
        measures = {}
        
        # 度中心性
        measures['degree'] = nx.degree_centrality(self.graph)
        
        # 介数中心性
        try:
            measures['betweenness'] = nx.betweenness_centrality(self.graph)
        except:
            pass
        
        # 接近中心性
        try:
            measures['closeness'] = nx.closeness_centrality(self.graph)
        except:
            pass
        
        # 特征向量中心性
        try:
            measures['eigenvector'] = nx.eigenvector_centrality(self.graph, max_iter=1000)
        except:
            pass
        
        # PageRank (适用于有向图)
        if self.graph.is_directed():
            try:
                measures['pagerank'] = nx.pagerank(self.graph)
            except:
                pass
        
        return measures
    
    def get_communities(self, method: str = 'greedy') -> List[set]:
        """
        检测社区
        
        Args:
            method: 社区检测方法 ('greedy', 'louvain', 'label_propagation')
        
        Returns:
            社区集合列表
        """
        if method == 'greedy':
            from networkx.algorithms import community
            return list(community.greedy_modularity_communities(self.graph))
        elif method == 'label_propagation':
            from networkx.algorithms import community
            return list(community.label_propagation_communities(self.graph))
        else:
            # 默认使用 greedy
            from networkx.algorithms import community
            return list(community.greedy_modularity_communities(self.graph))
    
    def get_shortest_paths(self, source: Any, target: Any = None) -> Dict:
        """
        计算最短路径
        
        Args:
            source: 源节点
            target: 目标节点 (None 则计算到所有节点的路径)
        
        Returns:
            路径信息
        """
        if target is not None:
            try:
                path = nx.shortest_path(self.graph, source, target)
                length = nx.shortest_path_length(self.graph, source, target)
                return {'path': path, 'length': length}
            except nx.NetworkXNoPath:
                return {'path': None, 'length': float('inf')}
        else:
            try:
                paths = dict(nx.shortest_path_length(self.graph, source))
                return {'paths': paths}
            except:
                return {'paths': {}}
    
    def analyze_node(self, node: Any) -> Dict:
        """
        分析单个节点
        
        Args:
            node: 节点标识符
        
        Returns:
            节点分析结果
        """
        if node not in self.graph:
            return {'error': f'Node {node} not in graph'}
        
        return {
            'degree': self.graph.degree(node),
            'neighbors': list(self.graph.neighbors(node)),
            'clustering': nx.clustering(self.graph, node),
            'eccentricity': nx.eccentricity(self.graph, node) if self.get_basic_metrics().is_connected else None
        }
    
    def export_metrics(self, filepath: str):
        """
        导出所有度量到文件
        
        Args:
            filepath: 输出文件路径 (.json)
        """
        metrics = {
            'basic': self.get_basic_metrics().to_dict(),
            'centrality': self.get_centrality_measures()
        }
        
        with open(filepath, 'w') as f:
            json.dump(metrics, f, indent=2)
        
        print(f"度量已导出到: {filepath}")


class NetworkGenerator:
    """网络生成器"""
    
    @staticmethod
    def create_random_graph(n: int, p: float, seed: int = 42) -> nx.Graph:
        """
        创建 Erdős-Rényi 随机图
        
        Args:
            n: 节点数
            p: 边概率
            seed: 随机种子
        
        Returns:
            随机图
        """
        return nx.erdos_renyi_graph(n, p, seed=seed)
    
    @staticmethod
    def create_scale_free_graph(n: int, m: int, seed: int = 42) -> nx.Graph:
        """
        创建 Barabási-Albert 无标度网络
        
        Args:
            n: 节点数
            m: 每个新节点的边数
            seed: 随机种子
        
        Returns:
            无标度网络
        """
        return nx.barabasi_albert_graph(n, m, seed=seed)
    
    @staticmethod
    def create_small_world_graph(n: int, k: int, p: float, seed: int = 42) -> nx.Graph:
        """
        创建 Watts-Strogatz 小世界网络
        
        Args:
            n: 节点数
            k: 每个节点的邻居数
            p: 重连概率
            seed: 随机种子
        
        Returns:
            小世界网络
        """
        return nx.watts_strogatz_graph(n, k, p, seed=seed)
    
    @staticmethod
    def create_grid_graph(m: int, n: int) -> nx.Graph:
        """
        创建二维网格图
        
        Args:
            m: 行数
            n: 列数
        
        Returns:
            网格图
        """
        return nx.grid_2d_graph(m, n)
    
    @staticmethod
    def create_social_network(model: str = 'karate', **kwargs) -> nx.Graph:
        """
        创建社交网络模型
        
        Args:
            model: 模型名称 ('karate', 'davis_southern_women', 'florentine_families')
        
        Returns:
            社交网络图
        """
        if model == 'karate':
            return nx.karate_club_graph()
        elif model == 'davis_southern_women':
            return nx.davis_southern_women_graph()
        elif model == 'florentine_families':
            return nx.florentine_families_graph()
        else:
            return nx.karate_club_graph()


class NetworkVisualizer:
    """网络可视化器"""
    
    def __init__(self, graph: nx.Graph):
        self.graph = graph
    
    def plot_basic(self, 
                   layout: str = 'spring',
                   figsize: Tuple[int, int] = (10, 8),
                   node_color: str = 'lightblue',
                   with_labels: bool = True,
                   save_path: Optional[str] = None):
        """
        基本网络可视化
        
        Args:
            layout: 布局算法 ('spring', 'circular', 'kamada_kawai', 'spectral')
            figsize: 图形大小
            node_color: 节点颜色
            with_labels: 是否显示标签
            save_path: 保存路径
        """
        if not MATPLOTLIB_AVAILABLE:
            print("错误: matplotlib 未安装")
            return
        
        plt.figure(figsize=figsize)
        
        # 选择布局
        if layout == 'spring':
            pos = nx.spring_layout(self.graph, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(self.graph)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(self.graph)
        elif layout == 'spectral':
            pos = nx.spectral_layout(self.graph)
        else:
            pos = nx.spring_layout(self.graph, seed=42)
        
        # 绘制
        nx.draw(self.graph, pos, 
                node_color=node_color,
                with_labels=with_labels,
                node_size=500,
                font_size=10,
                font_weight='bold',
                edge_color='gray')
        
        plt.title(f'Network: {self.graph.number_of_nodes()} nodes, {self.graph.number_of_edges()} edges')
        plt.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"图形已保存到: {save_path}")
        else:
            plt.show()
    
    def plot_centrality(self,
                       centrality_type: str = 'degree',
                       figsize: Tuple[int, int] = (12, 10),
                       save_path: Optional[str] = None):
        """
        基于中心性的可视化
        
        Args:
            centrality_type: 中心性类型
            figsize: 图形大小
            save_path: 保存路径
        """
        if not MATPLOTLIB_AVAILABLE:
            print("错误: matplotlib 未安装")
            return
        
        # 计算中心性
        analyzer = NetworkAnalyzer(self.graph)
        centralities = analyzer.get_centrality_measures()
        
        if centrality_type not in centralities:
            print(f"中心性类型 {centrality_type} 不可用")
            return
        
        centrality = centralities[centrality_type]
        
        # 根据中心性设置节点大小
        node_sizes = [3000 * centrality.get(node, 0) + 100 for node in self.graph.nodes()]
        
        plt.figure(figsize=figsize)
        pos = nx.spring_layout(self.graph, seed=42)
        
        # 绘制
        nodes = nx.draw_networkx_nodes(self.graph, pos, 
                                       node_size=node_sizes,
                                       node_color=list(centrality.values()),
                                       cmap=plt.cm.viridis)
        nx.draw_networkx_edges(self.graph, pos, alpha=0.3)
        nx.draw_networkx_labels(self.graph, pos, font_size=8)
        
        plt.colorbar(nodes, label=f'{centrality_type.capitalize()} Centrality')
        plt.title(f'Network Colored by {centrality_type.capitalize()} Centrality')
        plt.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
        else:
            plt.show()


def convert_to_dataframe(graph: nx.Graph) -> 'pd.DataFrame':
    """
    将图转换为 DataFrame
    
    Args:
        graph: NetworkX 图
    
    Returns:
        边列表 DataFrame
    """
    if not PANDAS_AVAILABLE:
        raise ImportError("pandas 未安装")
    
    edges_data = []
    for u, v, data in graph.edges(data=True):
        edge_info = {'source': u, 'target': v}
        edge_info.update(data)
        edges_data.append(edge_info)
    
    return pd.DataFrame(edges_data)


def find_bridges(graph: nx.Graph) -> List[Tuple]:
    """
    查找图中的桥接边
    
    Args:
        graph: 无向图
    
    Returns:
        桥接边列表
    """
    if graph.is_directed():
        graph = graph.to_undirected()
    
    return list(nx.bridges(graph))


def find_articulation_points(graph: nx.Graph) -> set:
    """
    查找关节点
    
    Args:
        graph: 无向图
    
        Returns:
        关节点集合
    """
    if graph.is_directed():
        graph = graph.to_undirected()
    
    return set(nx.articulation_points(graph))


def demo():
    """演示功能"""
    print("=" * 60)
    print("NetworkX 网络分析演示")
    print("=" * 60)
    
    # 1. 创建网络
    print("\n1. 创建 Karate Club 社交网络:")
    G = NetworkGenerator.create_social_network('karate')
    print(f"   节点数: {G.number_of_nodes()}")
    print(f"   边数: {G.number_of_edges()}")
    
    # 2. 基本分析
    print("\n2. 基本网络度量:")
    analyzer = NetworkAnalyzer(G)
    metrics = analyzer.get_basic_metrics()
    print(f"   密度: {metrics.density:.4f}")
    print(f"   平均聚类系数: {metrics.avg_clustering:.4f}")
    print(f"   是否连通: {metrics.is_connected}")
    if metrics.avg_shortest_path:
        print(f"   平均最短路径: {metrics.avg_shortest_path:.4f}")
    
    # 3. 中心性分析
    print("\n3. 中心性分析 (Top 5):")
    centralities = analyzer.get_centrality_measures()
    if 'degree' in centralities:
        degree_cent = centralities['degree']
        top_nodes = sorted(degree_cent.items(), key=lambda x: x[1], reverse=True)[:5]
        for node, cent in top_nodes:
            print(f"   节点 {node}: {cent:.4f}")
    
    # 4. 社区检测
    print("\n4. 社区检测:")
    communities = analyzer.get_communities('greedy')
    print(f"   发现 {len(communities)} 个社区")
    for i, comm in enumerate(communities[:3]):
        print(f"   社区 {i+1}: {len(comm)} 个节点")
    
    # 5. 创建随机网络
    print("\n5. 创建随机网络模型:")
    er_graph = NetworkGenerator.create_random_graph(100, 0.05)
    sf_graph = NetworkGenerator.create_scale_free_graph(100, 3)
    sw_graph = NetworkGenerator.create_small_world_graph(100, 6, 0.1)
    
    print(f"   随机网络 (Erdős-Rényi): {er_graph.number_of_edges()} 条边")
    print(f"   无标度网络 (Barabási-Albert): {sf_graph.number_of_edges()} 条边")
    print(f"   小世界网络 (Watts-Strogatz): {sw_graph.number_of_edges()} 条边")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
