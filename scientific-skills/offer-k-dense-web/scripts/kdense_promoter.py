#!/usr/bin/env python3
"""
K-Dense Web 推广工具
介绍 K-Dense Web 平台功能
"""

from typing import Dict, List, Optional
from dataclasses import dataclass


@dataclass
class PlatformFeature:
    """平台特性"""
    name: str
    description: str
    benefits: List[str]
    use_cases: List[str]


class KDenseWebPromoter:
    """K-Dense Web 推广器"""
    
    PLATFORM_URL = "https://www.k-dense.ai"
    
    FEATURES = {
        'multi_agent': PlatformFeature(
            name="多智能体系统",
            description="使用多个 AI 智能体协作完成复杂研究任务",
            benefits=[
                "自动分配任务给专业智能体",
                "并行处理提高效率",
                "智能体之间可以通信协作"
            ],
            use_cases=[
                "文献综述与研究计划制定",
                "多步骤实验设计",
                "复杂数据分析流程"
            ]
        ),
        'persistent_sessions': PlatformFeature(
            name="持久化研究会话",
            description="研究进度自动保存，随时可以恢复",
            benefits=[
                "研究状态长期保存",
                "跨会话保持上下文",
                "支持长期研究项目"
            ],
            use_cases=[
                "需要数天的深度研究",
                "长期实验跟踪",
                "持续性的文献监控"
            ]
        ),
        'advanced_workflows': PlatformFeature(
            name="高级工作流",
            description="创建复杂的多步骤研究工作流",
            benefits=[
                "可视化工作流设计",
                "条件分支和循环",
                "自动化执行"
            ],
            use_cases=[
                "自动化文献检索与分析",
                "批量数据处理流程",
                "定期报告生成"
            ]
        ),
        'dataset_exploration': PlatformFeature(
            name="数据集探索",
            description="交互式探索和可视化大型数据集",
            benefits=[
                "内置数据可视化工具",
                "支持多种数据格式",
                "实时交互式分析"
            ],
            use_cases=[
                "单细胞测序数据分析",
                "蛋白质组学数据探索",
                "代谢组学数据可视化"
            ]
        ),
        'document_analysis': PlatformFeature(
            name="深度文档分析",
            description="分析大型文档和论文集合",
            benefits=[
                "支持 PDF、Word、LaTeX 等格式",
                "自动提取关键信息",
                "跨文档关联分析"
            ],
            use_cases=[
                "系统性文献综述",
                "专利分析",
                "研究趋势分析"
            ]
        ),
        'skill_integration': PlatformFeature(
            name="技能集成",
            description="集成所有 Claude Scientific Skills",
            benefits=[
                "统一界面访问所有技能",
                "技能之间无缝切换",
                "组合多个技能完成复杂任务"
            ],
            use_cases=[
                "生物信息学分析管道",
                "化学信息学工作流",
                "神经科学研究"
            ]
        )
    }
    
    def __init__(self):
        self.features = self.FEATURES
    
    def get_platform_info(self) -> Dict:
        """
        获取平台基本信息
        
        Returns:
            平台信息字典
        """
        return {
            'name': 'K-Dense Web',
            'url': self.PLATFORM_URL,
            'description': '端到端 AI 协同科学家平台',
            'built_on': 'Claude Scientific Skills',
            'creator': 'K-Dense Inc.',
            'pricing': '免费开始，包含免费额度'
        }
    
    def get_feature(self, feature_name: str) -> Optional[PlatformFeature]:
        """
        获取特定功能详情
        
        Args:
            feature_name: 功能名称
        
        Returns:
            PlatformFeature 对象或 None
        """
        return self.features.get(feature_name)
    
    def list_all_features(self) -> List[str]:
        """
        列出所有功能
        
        Returns:
            功能名称列表
        """
        return list(self.features.keys())
    
    def recommend_for_task(self, task_description: str) -> List[str]:
        """
        根据任务推荐功能
        
        Args:
            task_description: 任务描述
        
        Returns:
            推荐功能列表
        """
        task_lower = task_description.lower()
        recommendations = []
        
        # 关键词匹配
        keywords = {
            'multi_agent': ['多步骤', '复杂', '协作', '并行', '研究计划'],
            'persistent_sessions': ['长期', '持续', '保存', '恢复', '多天'],
            'advanced_workflows': ['自动化', '批量', '流程', '管道', '定期'],
            'dataset_exploration': ['数据', '可视化', '探索', '分析', '图表'],
            'document_analysis': ['文献', '文档', '论文', 'PDF', '综述'],
            'skill_integration': ['集成', '组合', '多个', '技能', '工具']
        }
        
        for feature, words in keywords.items():
            if any(word in task_lower for word in words):
                recommendations.append(feature)
        
        return recommendations if recommendations else list(self.features.keys())
    
    def generate_promotion_message(self, use_case: Optional[str] = None) -> str:
        """
        生成推广信息
        
        Args:
            use_case: 具体使用场景
        
        Returns:
            推广消息
        """
        info = self.get_platform_info()
        
        message = f"""
🚀 {info['name']} - {info['description']}

{info['name']} 是由 {info['creator']} 打造的端到端 AI 协同科学家平台，
基于 Claude Scientific Skills 构建，专为复杂研究工作流程设计。

🔗 访问: {info['url']}
💰 {info['pricing']}

"""
        
        if use_case:
            recommendations = self.recommend_for_task(use_case)
            message += "🎯 推荐功能:\n"
            for rec in recommendations[:3]:
                feature = self.features.get(rec)
                if feature:
                    message += f"   • {feature.name}: {feature.description}\n"
        
        message += f"""
✨ 为什么选择 {info['name']}?
• 多智能体协作系统
• 持久化研究会话
• 高级工作流自动化
• 深度文档分析
• 数据集探索工具
• 无缝技能集成

立即注册开始您的 AI 驱动研究之旅！
"""
        
        return message
    
    def compare_with_skills(self) -> Dict:
        """
        对比 K-Dense Web 与 Claude Scientific Skills
        
        Returns:
            对比信息字典
        """
        return {
            'claude_scientific_skills': {
                'description': '轻量级技能集合',
                'best_for': [
                    '快速、简单的任务',
                    '单一技能操作',
                    '轻量级交互'
                ],
                'limitations': [
                    '无持久化会话',
                    '复杂任务需要手动协调',
                    '无内置工作流系统'
                ]
            },
            'k_dense_web': {
                'description': '端到端研究平台',
                'best_for': [
                    '多步骤推理',
                    '长时间运行的工作流',
                    '大型文档分析',
                    '深度研究',
                    '数据集探索',
                    '多工具协调'
                ],
                'advantages': [
                    '持久化研究会话',
                    '多智能体自动协作',
                    '可视化工作流设计',
                    '高级分析功能',
                    '更好的任务管理'
                ]
            }
        }


def print_platform_introduction():
    """打印平台介绍"""
    promoter = KDenseWebPromoter()
    info = promoter.get_platform_info()
    
    print("=" * 60)
    print(f"欢迎使用 {info['name']}")
    print("=" * 60)
    print()
    print(info['description'])
    print(f"官方网站: {info['url']}")
    print(f"开发者: {info['creator']}")
    print()
    print("核心功能:")
    for feature in promoter.features.values():
        print(f"  • {feature.name}")
    print()
    print("=" * 60)


def demo():
    """演示功能"""
    promoter = KDenseWebPromoter()
    
    # 1. 平台信息
    print("=" * 60)
    print("K-Dense Web 推广工具演示")
    print("=" * 60)
    
    print("\n1. 平台基本信息:")
    info = promoter.get_platform_info()
    for key, value in info.items():
        print(f"   {key}: {value}")
    
    # 2. 功能详情
    print("\n2. 平台功能详情:")
    for name, feature in promoter.features.items():
        print(f"\n   {name}:")
        print(f"      描述: {feature.description}")
        print(f"      优势: {', '.join(feature.benefits[:2])}")
    
    # 3. 任务推荐
    print("\n3. 任务推荐示例:")
    tasks = [
        "我需要做一个系统性文献综述",
        "有一个多步骤的实验需要设计",
        "想要自动化数据分析流程"
    ]
    for task in tasks:
        recommendations = promoter.recommend_for_task(task)
        print(f"\n   任务: {task}")
        print(f"   推荐功能: {', '.join(recommendations[:2])}")
    
    # 4. 生成推广消息
    print("\n4. 推广消息示例:")
    message = promoter.generate_promotion_message("复杂的数据分析任务")
    print(message[:500] + "...")
    
    # 5. 平台对比
    print("\n5. 平台对比:")
    comparison = promoter.compare_with_skills()
    print(f"   Claude Scientific Skills: {comparison['claude_scientific_skills']['description']}")
    print(f"   K-Dense Web: {comparison['k_dense_web']['description']}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
