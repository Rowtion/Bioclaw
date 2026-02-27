#!/usr/bin/env python3
"""
假设生成示例脚本 2: 假设质量评估框架
演示如何评估科学假设的质量
"""

import json
from dataclasses import dataclass, asdict
from typing import List, Dict


@dataclass
class Hypothesis:
    """假设数据结构"""
    id: str
    name: str
    description: str
    mechanism: str
    evidence: List[str]
    predictions: List[str]
    test_methods: List[str]
    assumptions: List[str]


class HypothesisEvaluator:
    """假设评估器"""
    
    def __init__(self):
        self.criteria = {
            'testability': '可测试性',
            'falsifiability': '可证伪性',
            'parsimony': '简洁性',
            'explanatory_power': '解释力',
            'scope': '适用范围',
            'consistency': '一致性',
            'novelty': '新颖性'
        }
    
    def evaluate(self, hypothesis: Hypothesis, scores: Dict[str, int]) -> Dict:
        """
        评估假设
        
        Args:
            hypothesis: 假设对象
            scores: 各标准的评分 (1-5)
        
        Returns:
            评估报告
        """
        report = {
            'hypothesis_id': hypothesis.id,
            'hypothesis_name': hypothesis.name,
            'scores': {},
            'total_score': 0,
            'average_score': 0,
            'strengths': [],
            'weaknesses': [],
            'recommendations': []
        }
        
        # 计算分数
        total = 0
        for criterion, score in scores.items():
            report['scores'][self.criteria.get(criterion, criterion)] = score
            total += score
        
        report['total_score'] = total
        report['average_score'] = total / len(scores) if scores else 0
        
        # 生成优势和劣势
        for criterion, score in scores.items():
            criterion_name = self.criteria.get(criterion, criterion)
            if score >= 4:
                report['strengths'].append(f"{criterion_name} (评分: {score}/5)")
            elif score <= 2:
                report['weaknesses'].append(f"{criterion_name} (评分: {score}/5)")
        
        # 生成建议
        if scores.get('testability', 5) <= 2:
            report['recommendations'].append("增加可操作的具体预测，提高可测试性")
        
        if scores.get('falsifiability', 5) <= 2:
            report['recommendations'].append("明确描述可能证伪该假设的条件")
        
        if scores.get('parsimony', 5) <= 2:
            report['recommendations'].append("简化机制解释，去除不必要的复杂性")
        
        return report
    
    def compare_hypotheses(self, hypotheses: List[Hypothesis], 
                          all_scores: Dict[str, Dict[str, int]]) -> Dict:
        """
        比较多个假设
        
        Args:
            hypotheses: 假设列表
            all_scores: 每个假设的评分
        
        Returns:
            比较报告
        """
        comparison = {
            'hypotheses': [],
            'ranking': [],
            'best_overall': None,
            'most_testable': None,
            'most_parsimonious': None,
            'highest_explanatory': None
        }
        
        evaluated = []
        
        for hyp in hypotheses:
            scores = all_scores.get(hyp.id, {})
            report = self.evaluate(hyp, scores)
            evaluated.append((hyp, report))
            comparison['hypotheses'].append(report)
        
        # 按总分排序
        sorted_hyps = sorted(evaluated, key=lambda x: x[1]['total_score'], reverse=True)
        comparison['ranking'] = [(h.id, r['total_score']) for h, r in sorted_hyps]
        
        if sorted_hyps:
            comparison['best_overall'] = sorted_hyps[0][0].id
        
        # 各单项最佳
        testable = max(evaluated, key=lambda x: x[1]['scores'].get('可测试性', 0))
        comparison['most_testable'] = testable[0].id
        
        parsimonious = max(evaluated, key=lambda x: x[1]['scores'].get('简洁性', 0))
        comparison['most_parsimonious'] = parsimonious[0].id
        
        explanatory = max(evaluated, key=lambda x: x[1]['scores'].get('解释力', 0))
        comparison['highest_explanatory'] = explanatory[0].id
        
        return comparison


def demonstrate_evaluation():
    """演示假设评估"""
    
    print("=" * 70)
    print("科学假设质量评估")
    print("=" * 70)
    
    # 创建示例假设
    hypotheses = [
        Hypothesis(
            id='H1',
            name='睡眠剥夺假说',
            description='蓝光通过抑制褪黑素导致睡眠剥夺，进而影响认知',
            mechanism='蓝光→褪黑素↓→睡眠质量↓→认知功能↓',
            evidence=['蓝光抑制褪黑素分泌', '睡眠剥夺损害认知'],
            predictions=['褪黑素补充应改善认知', '睡眠质量与认知分数相关'],
            test_methods=['褪黑素干预实验', '睡眠监测'],
            assumptions=['褪黑素是主要介导因子', '认知损害主要由睡眠引起']
        ),
        Hypothesis(
            id='H2',
            name='氧化应激假说',
            description='蓝光诱导视网膜氧化应激，通过视神经影响大脑',
            mechanism='蓝光→视网膜ROS↑→氧化应激→神经损伤→认知↓',
            evidence=['蓝光诱导ROS', '氧化应激导致神经元损伤'],
            predictions=['抗氧化剂应保护认知', '视网膜氧化标志物升高'],
            test_methods=['抗氧化剂干预', '视网膜组织学'],
            assumptions=['氧化应激可扩散到大脑', '视网膜是主要靶点']
        ),
        Hypothesis(
            id='H3',
            name='节律失调假说',
            description='蓝光扰乱昼夜节律，影响海马体时钟基因表达',
            mechanism='蓝光→SCN紊乱→外周时钟失调→海马可塑性↓',
            evidence=['SCN是主时钟', '海马有局部生物钟'],
            predictions=['限时暴露应减轻影响', '时钟基因表达改变'],
            test_methods=['时间限制实验', '基因表达分析'],
            assumptions=['海马时钟直接调控认知', '节律紊乱是主要原因']
        ),
    ]
    
    # 评估评分
    all_scores = {
        'H1': {
            'testability': 5,
            'falsifiability': 4,
            'parsimony': 5,
            'explanatory_power': 3,
            'scope': 3,
            'consistency': 4,
            'novelty': 2
        },
        'H2': {
            'testability': 3,
            'falsifiability': 3,
            'parsimony': 3,
            'explanatory_power': 5,
            'scope': 4,
            'consistency': 4,
            'novelty': 4
        },
        'H3': {
            'testability': 4,
            'falsifiability': 3,
            'parsimony': 3,
            'explanatory_power': 4,
            'scope': 5,
            'consistency': 4,
            'novelty': 4
        },
    }
    
    # 创建评估器
    evaluator = HypothesisEvaluator()
    
    # 评估每个假设
    print("\n单独评估:")
    print("-" * 70)
    
    for hyp in hypotheses:
        scores = all_scores[hyp.id]
        report = evaluator.evaluate(hyp, scores)
        
        print(f"\n假设 {hyp.id}: {hyp.name}")
        print(f"  总分: {report['total_score']}/35")
        print(f"  平均分: {report['average_score']:.2f}/5")
        
        if report['strengths']:
            print(f"  优势: {', '.join(report['strengths'][:2])}")
        
        if report['weaknesses']:
            print(f"  劣势: {', '.join(report['weaknesses'][:2])}")
    
    # 比较假设
    print("\n" + "=" * 70)
    print("假设比较")
    print("=" * 70)
    
    comparison = evaluator.compare_hypotheses(hypotheses, all_scores)
    
    print(f"\n综合排名:")
    for rank, (hyp_id, score) in enumerate(comparison['ranking'], 1):
        print(f"  {rank}. {hyp_id} (总分: {score})")
    
    print(f"\n单项最佳:")
    print(f"  最可测试: {comparison['most_testable']}")
    print(f"  最简洁: {comparison['most_parsimonious']}")
    print(f"  最强解释力: {comparison['highest_explanatory']}")
    
    # 生成比较矩阵
    print("\n评分矩阵:")
    print("-" * 70)
    criteria_names = list(evaluator.criteria.values())
    print(f"{'假设':<10}", end='')
    for name in criteria_names:
        print(f"{name:<10}", end='')
    print("总分")
    print("-" * 70)
    
    for hyp in hypotheses:
        scores = all_scores[hyp.id]
        print(f"{hyp.id:<10}", end='')
        for criterion in evaluator.criteria.keys():
            score = scores.get(criterion, 0)
            print(f"{score:<10}", end='')
        print(f"{sum(scores.values())}")
    
    # 保存结果
    output = {
        'individual_evaluations': [
            evaluator.evaluate(h, all_scores[h.id]) for h in hypotheses
        ],
        'comparison': comparison
    }
    
    with open('hypothesis_evaluation.json', 'w') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)
    
    print(f"\n✓ 评估结果已保存到 hypothesis_evaluation.json")


def experimental_design_guide():
    """实验设计指南"""
    
    print("\n" + "=" * 70)
    print("实验设计模式指南")
    print("=" * 70)
    
    designs = {
        '随机对照试验 (RCT)': {
            '适用': '测试干预措施效果',
            '设计': '实验组 vs 对照组，随机分配',
            '示例': '蓝光暴露组 vs 对照光组，随机分组',
            '分析': 't检验或ANOVA比较组间差异'
        },
        '剂量反应实验': {
            '适用': '确定效应与暴露的关系',
            '设计': '多个剂量水平组',
            '示例': '低/中/高剂量蓝光 + 对照组',
            '分析': '趋势分析，EC50计算'
        },
        '机制阻断实验': {
            '适用': '验证因果机制',
            '设计': '药物/基因阻断特定通路',
            '示例': '蓝光+褪黑素 vs 蓝光+褪黑素拮抗剂',
            '分析': '比较阻断前后的效应差异'
        },
        '时间序列实验': {
            '适用': '观察时程变化',
            '设计': '多个时间点测量',
            '示例': '暴露前、1周、2周、4周测量',
            '分析': '重复测量ANOVA或混合效应模型'
        },
        '交叉设计': {
            '适用': '控制个体差异',
            '设计': '受试者作为自身对照',
            '示例': '先蓝光后对照 vs 先对照后蓝光',
            '分析': '配对t检验，处理顺序效应'
        }
    }
    
    for name, info in designs.items():
        print(f"\n{name}")
        print("-" * 40)
        for key, value in info.items():
            print(f"  {key}: {value}")


def main():
    demonstrate_evaluation()
    experimental_design_guide()
    
    print("\n" + "=" * 70)
    print("引用:")
    print("-" * 70)
    print("假设质量评估标准参考:")
    print("- Popper, K. (1959). The Logic of Scientific Discovery")
    print("- Chamberlin, T. C. (1890). The Method of Multiple Working Hypotheses")


if __name__ == "__main__":
    main()
