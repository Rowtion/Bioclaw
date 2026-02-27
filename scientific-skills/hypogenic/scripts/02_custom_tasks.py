#!/usr/bin/env python3
"""
Hypogenic 示例脚本 2: 自定义任务实现
演示如何为特定领域创建自定义假设生成任务
"""

import json
import re
from pathlib import Path


def create_deception_detection_task():
    """
    创建欺骗检测任务的示例
    
    应用场景: 酒店评论的真假检测
    """
    
    print("=" * 70)
    print("自定义任务: 欺骗检测 (Deception Detection)")
    print("=" * 70)
    
    # 数据集格式
    print("\n1. 数据集格式")
    print("-" * 70)
    
    data_format = {
        'text_features_1': [
            'This hotel was amazing! The staff was so helpful and the room was perfect.',
            'Terrible experience. The room was dirty and the AC did not work.',
        ],
        'label': [
            'fake',      # 虚假评论
            'genuine',   # 真实评论
        ]
    }
    
    print("数据集格式示例:")
    print(json.dumps(data_format, indent=2))
    
    # 配置文件
    print("\n2. 配置文件 (config.yaml)")
    print("-" * 70)
    
    config = '''
task_name: deception_detection

train_data_path: ./data/deception_train.json
val_data_path: ./data/deception_val.json
test_data_path: ./data/deception_test.json

prompt_templates:
  observations: |
    Review: ${text_features_1}
    Label: ${label}
  
  batched_generation:
    system: |
      You are an expert in linguistic analysis and deception detection.
      Generate hypotheses about linguistic cues that distinguish fake from genuine reviews.
    user: |
      Based on the following review examples, generate ${num_hypotheses} hypotheses
      about linguistic patterns that indicate deception in hotel reviews.
      
      Examples:
      ${observations}
      
      Generate specific, testable hypotheses about:
      - Word choice patterns
      - Emotional expression
      - Specificity of details
      - Writing style indicators
  
  inference:
    system: |
      Determine if a hotel review is fake or genuine based on linguistic cues.
    user: |
      Hypothesis: ${hypothesis}
      
      Review: ${text_features_1}
      
      Is this review FAKE or GENUINE?
      Provide your answer in the format: "Final answer: [fake/genuine]"
'''
    
    print(config)
    
    # 标签提取函数
    print("\n3. 自定义标签提取函数")
    print("-" * 70)
    
    code = '''
def extract_deception_label(llm_output: str) -> str:
    """
    从 LLM 输出中提取欺骗检测标签
    
    Args:
        llm_output: LLM 的原始输出文本
    
    Returns:
        提取的标签 ('fake' 或 'genuine')
    """
    # 匹配 "Final answer: fake/genuine" 格式
    match = re.search(
        r'final answer:\s*(fake|genuine)',
        llm_output,
        re.IGNORECASE
    )
    
    if match:
        return match.group(1).lower()
    
    # 回退: 检查关键词
    text_lower = llm_output.lower()
    if 'fake' in text_lower or 'false' in text_lower:
        return 'fake'
    elif 'genuine' in text_lower or 'real' in text_lower or 'authentic' in text_lower:
        return 'genuine'
    
    # 默认返回
    return 'unknown'

# 使用自定义提取函数
from hypogenic import BaseTask

task = BaseTask(
    config_path='./config/deception_config.yaml',
    extract_label=extract_deception_label
)
'''
    
    print(code)
    
    # 示例假设
    print("\n4. 预期生成的假设示例")
    print("-" * 70)
    
    hypotheses = [
        {
            'hypothesis_id': 'hyp_001',
            'text': 'Fake reviews use more extreme positive adjectives (amazing, perfect, incredible) than genuine reviews',
            'performance': 0.72,
            'rationale': 'Deceptive writers tend to overcompensate with exaggerated positivity'
        },
        {
            'hypothesis_id': 'hyp_002',
            'text': 'Genuine reviews contain more specific details (room number, staff names, specific incidents) than fake reviews',
            'performance': 0.78,
            'rationale': 'Real experiences include concrete, memorable details'
        },
        {
            'hypothesis_id': 'hyp_003',
            'text': 'Fake reviews use fewer first-person pronouns (I, my, we) compared to genuine reviews',
            'performance': 0.65,
            'rationale': 'Psychological distance in deception leads to reduced self-reference'
        },
        {
            'hypothesis_id': 'hyp_004',
            'text': 'Reviews mentioning specific timeframes (check-in date, length of stay) are more likely to be genuine',
            'performance': 0.70,
            'rationale': 'Specific temporal details indicate authentic experience'
        },
    ]
    
    for hyp in hypotheses:
        print(f"\n  {hyp['hypothesis_id']}:")
        print(f"    假设: {hyp['text']}")
        print(f"    性能: {hyp['performance']:.2%}")
        print(f"    原理: {hyp['rationale']}")


def create_health_indicators_task():
    """创建健康指标检测任务示例"""
    
    print("\n" + "=" * 70)
    print("自定义任务: 心理健康指标检测")
    print("=" * 70)
    
    print("""
任务描述:
---------
从社交媒体帖子中检测心理健康指标

数据格式:
---------
{
  "post_content": [
    "今天感觉特别好，完成了所有工作任务！",
    "又一个失眠的夜晚，不知道还能坚持多久...",
  ],
  "label": [
    "positive",
    "distress",
  ]
}

生成假设的方向:
--------------
1. 情感词汇的使用模式
2. 时间表达 (过去/现在/未来)
3. 社交语言特征
4. 自我关注程度
5. 活动描述的具体性

预期假设示例:
------------
- 使用未来时态较少的帖子与负面情绪相关
- 第一人称单数代词频率升高可能表示心理困扰
- 包含活动描述 (运动、社交) 的帖子通常更积极
""")


def demonstrate_evaluation():
    """演示假设评估方法"""
    
    print("\n" + "=" * 70)
    print("假设评估和选择")
    print("=" * 70)
    
    print("""
评估指标:
---------
1. 准确率 (Accuracy): 假设正确预测的比例
2. 精确率 (Precision): 正例预测中真正例的比例
3. 召回率 (Recall): 真正例中被正确预测的比例
4. F1 分数: 精确率和召回率的调和平均
5. 多样性 (Diversity): 假设之间的非重复性

选择策略:
---------
- 迭代优化: 基于错误样例生成新假设
- 多数投票: 多个假设共同决策
- 加权组合: 根据性能给假设赋权重
- 自适应选择: 根据输入特征动态选择假设

代码示例:
---------
# 多假设推理
from examples.multi_hyp_inference import run_multi_hypothesis_inference

results = run_multi_hypothesis_inference(
    config_path='./config/config.yaml',
    hypothesis_bank='./output/hypotheses.json',
    test_data='./data/test.json',
    aggregation_method='majority_vote'  # 或 'weighted'
)
""")


def demonstrate_optimization():
    """演示性能优化"""
    
    print("\n" + "=" * 70)
    print("性能优化技巧")
    print("=" * 70)
    
    print("""
1. 使用 Redis 缓存:
   # 安装 Redis 并运行在端口 6832
   # 自动缓存 LLM 响应，减少 API 调用成本

2. 并行处理:
   - 批量生成假设
   - 并行验证多个假设
   - 使用多进程加速推理

3. 自适应细化:
   - 专注于困难样例
   - 根据错误模式生成针对性假设
   - 动态调整假设库

4. 提示工程优化:
   - 使用 few-shot 示例
   - 明确输出格式要求
   - 添加角色设定

配置示例:
---------
prompt_templates:
  batched_generation:
    system: "You are an expert in [domain]. Generate testable hypotheses."
    user: |
      Based on these ${num_observations} examples:
      ${observations}
      
      Generate ${num_hypotheses} specific hypotheses about [phenomenon].
      Each hypothesis should be:
      - Specific and testable
      - Based on observable patterns
      - Distinct from others
      
      Format: One hypothesis per line, numbered.
""")


def main():
    print("=" * 70)
    print("Hypogenic 自定义任务实现指南")
    print("=" * 70)
    
    # 检查依赖
    try:
        import hypogenic
        print("\n✓ hypogenic 已安装")
    except ImportError:
        print("\n× hypogenic 未安装")
        print("  安装命令: pip install hypogenic")
    
    try:
        import yaml
        print("✓ pyyaml 已安装")
    except ImportError:
        print("  安装命令: pip install pyyaml")
    
    # 运行演示
    create_deception_detection_task()
    create_health_indicators_task()
    demonstrate_evaluation()
    demonstrate_optimization()
    
    print("\n" + "=" * 70)
    print("最佳实践:")
    print("-" * 70)
    print("1. 从小数据集开始测试配置")
    print("2. 仔细设计标签提取函数")
    print("3. 迭代优化提示模板")
    print("4. 监控 API 使用成本")
    print("5. 保存中间结果以便复现")
    print("\n论文引用:")
    print("-" * 70)
    print("Zhou et al. (2024). Hypothesis Generation with Large Language Models.")
    print("Liu et al. (2024). Literature Meets Data: A Synergistic Approach to")
    print("  Hypothesis Generation.")


if __name__ == "__main__":
    main()
