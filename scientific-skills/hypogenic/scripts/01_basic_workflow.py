#!/usr/bin/env python3
"""
Hypogenic 示例脚本 1: 基本假设生成流程
演示如何使用 HypoGeniC 方法从数据生成假设
"""

import json
import yaml
from pathlib import Path


def create_sample_config():
    """创建示例配置文件"""
    
    config = {
        'task_name': 'sentiment_analysis',
        
        'train_data_path': './data/train.json',
        'val_data_path': './data/val.json',
        'test_data_path': './data/test.json',
        
        'prompt_templates': {
            'observations': '''
文本: ${text}
标签: ${label}
''',
            'batched_generation': {
                'system': '你是一个假设生成专家。基于观察数据生成解释性假设。',
                'user': '''
基于以下观察，生成 ${num_hypotheses} 个解释为什么文本会被标记为特定标签的假设。

观察:
${observations}

请生成具体的、可测试的假设。每个假设应该是简洁的陈述。
'''
            },
            'inference': {
                'system': '基于以下假设判断文本标签。',
                'user': '''
假设: ${hypothesis}

文本: ${text}

这个文本的标签是什么？请以 "Final answer: [标签]" 的格式回答。
'''
            }
        }
    }
    
    return config


def create_sample_data():
    """创建示例数据集"""
    
    train_data = {
        'text': [
            'This movie is absolutely fantastic!',
            'Terrible waste of time, do not watch.',
            'Best film I have seen this year.',
            'Boring and predictable storyline.',
            'Amazing performances by all actors.',
            'Worst acting I have ever seen.',
        ],
        'label': [
            'positive',
            'negative',
            'positive',
            'negative',
            'positive',
            'negative',
        ]
    }
    
    val_data = {
        'text': [
            'Great cinematography and direction.',
            'Poor script and weak dialogues.',
        ],
        'label': [
            'positive',
            'negative',
        ]
    }
    
    return train_data, val_data


def demonstrate_hypogenic_workflow():
    """演示 HypoGeniC 工作流程"""
    
    print("=" * 70)
    print("HypoGeniC 假设生成工作流程")
    print("=" * 70)
    
    print("""
HypoGeniC 是一种数据驱动的假设生成方法，使用大语言模型(LLM)
从观察数据中自动发现和生成可测试的科学假设。

工作流程:
---------
1. 准备数据集 (train/val/test)
2. 配置提示模板 (prompt templates)
3. 运行假设生成
4. 验证和评分假设
5. 使用生成的假设进行推理
""")
    
    # 步骤 1: 创建配置
    print("\n步骤 1: 创建配置文件")
    print("-" * 70)
    
    config = create_sample_config()
    
    # 保存配置
    Path('./config').mkdir(exist_ok=True)
    with open('./config/config.yaml', 'w') as f:
        yaml.dump(config, f, default_flow_style=False, allow_unicode=True)
    
    print("✓ 配置文件已创建: ./config/config.yaml")
    print("\n配置内容:")
    print(yaml.dump(config, default_flow_style=False, allow_unicode=True))
    
    # 步骤 2: 创建数据
    print("\n步骤 2: 准备数据集")
    print("-" * 70)
    
    train_data, val_data = create_sample_data()
    
    Path('./data').mkdir(exist_ok=True)
    
    with open('./data/train.json', 'w') as f:
        json.dump(train_data, f, indent=2)
    
    with open('./data/val.json', 'w') as f:
        json.dump(val_data, f, indent=2)
    
    with open('./data/test.json', 'w') as f:
        json.dump(val_data, f, indent=2)
    
    print("✓ 数据集已创建")
    print(f"  - 训练集: {len(train_data['text'])} 条")
    print(f"  - 验证集: {len(val_data['text'])} 条")
    
    # 步骤 3: 展示生成的假设示例
    print("\n步骤 3: 假设生成示例")
    print("-" * 70)
    
    example_hypotheses = [
        {
            'hypothesis': '文本包含积极形容词(如 fantastic, amazing, great)表示正面情感',
            'score': 0.85,
            'supporting_examples': 5
        },
        {
            'hypothesis': '文本包含消极词汇(如 terrible, boring, poor, worst)表示负面情感',
            'score': 0.82,
            'supporting_examples': 4
        },
        {
            'hypothesis': '比较级词汇(best, worst)的出现强烈指示情感极性',
            'score': 0.78,
            'supporting_examples': 3
        },
        {
            'hypothesis': '文本长度与情感极性无显著关联',
            'score': 0.45,
            'supporting_examples': 2
        },
    ]
    
    print("生成的假设示例:")
    for i, hyp in enumerate(example_hypotheses, 1):
        print(f"\n  假设 {i}:")
        print(f"    内容: {hyp['hypothesis']}")
        print(f"    置信度: {hyp['score']:.2%}")
        print(f"    支持样本数: {hyp['supporting_examples']}")
    
    # 步骤 4: 推理示例
    print("\n步骤 4: 使用假设进行推理")
    print("-" * 70)
    
    test_text = "The acting was superb and the plot was engaging."
    
    print(f"测试文本: {test_text}")
    print("\n各假设预测:")
    
    predictions = [
        ('假设 1', 'positive', 0.85),
        ('假设 2', 'negative', 0.20),
        ('假设 3', 'positive', 0.75),
    ]
    
    for hyp_name, pred, conf in predictions:
        print(f"  {hyp_name}: {pred} (置信度: {conf:.2%})")
    
    # 综合预测
    print("\n综合预测: positive")
    print("推理: 多个假设支持正面情感，特别是关于积极形容词的假设")


def demonstrate_command_line_usage():
    """演示命令行用法"""
    
    print("\n" + "=" * 70)
    print("命令行使用示例")
    print("=" * 70)
    
    print("""
# 1. 安装 hypogenic
pip install hypogenic

# 2. 准备数据集和配置
# 参见上面的示例代码

# 3. 运行假设生成
hypogenic_generation \\
    --config ./config/config.yaml \\
    --method hypogenic \\
    --num_hypotheses 20 \\
    --output ./output/hypotheses.json

# 参数说明:
#   --config: 配置文件路径
#   --method: 生成方法 (hypogenic, hyporefine, union)
#   --num_hypotheses: 生成假设数量
#   --output: 输出文件路径

# 4. 运行推理
hypogenic_inference \\
    --config ./config/config.yaml \\
    --hypotheses ./output/hypotheses.json \\
    --test_data ./data/test.json \\
    --output ./output/results.json

# 5. 使用 Python API
python -c "
from hypogenic import BaseTask

task = BaseTask(config_path='./config/config.yaml')
task.generate_hypotheses(method='hypogenic', num_hypotheses=20)
results = task.inference(hypothesis_bank='./output/hypotheses.json')
"
""")


def demonstrate_hyporefine():
    """演示 HypoRefine (文献+数据整合)"""
    
    print("\n" + "=" * 70)
    print("HypoRefine: 文献与数据整合")
    print("=" * 70)
    
    print("""
HypoRefine 结合了文献洞察和经验数据，生成更全面的假设。

工作流程:
---------
1. 准备文献 PDF 文件
2. 使用 GROBID 处理文献
3. 从文献生成理论假设
4. 从数据生成经验假设
5. 整合并优化假设库

目录结构:
---------
literature/
└── your_task/
    └── raw/
        ├── paper1.pdf
        ├── paper2.pdf
        └── paper3.pdf

运行命令:
---------
# 设置 GROBID
bash ./modules/setup_grobid.sh

# 处理文献
bash ./modules/run_grobid.sh

# 生成假设 (整合文献和数据)
hypogenic_generation \\
    --config ./config/config.yaml \\
    --method hyporefine \\
    --num_hypotheses 15 \\
    --literature ./literature/your_task/

输出:
-----
- hypotheses_hyporefine.json (整合方法)
- hypotheses_literature.json (仅文献)
- hypotheses_literature_union_hyporefine.json (并集方法)
""")


def main():
    print("=" * 70)
    print("Hypogenic 假设生成框架示例")
    print("=" * 70)
    
    # 检查安装
    try:
        import hypogenic
        print("\n✓ hypogenic 已安装")
    except ImportError:
        print("\n× hypogenic 未安装")
        print("  安装命令: pip install hypogenic")
    
    # 运行演示
    demonstrate_hypogenic_workflow()
    demonstrate_command_line_usage()
    demonstrate_hyporefine()
    
    print("\n" + "=" * 70)
    print("资源链接:")
    print("-" * 70)
    print("- GitHub: https://github.com/ChicagoHAI/hypothesis-generation")
    print("- PyPI: https://pypi.org/project/hypogenic/")
    print("- 论文: https://arxiv.org/abs/2410.17309")
    
    # 清理
    import shutil
    for d in ['./config', './data']:
        if Path(d).exists():
            shutil.rmtree(d)


if __name__ == "__main__":
    main()
