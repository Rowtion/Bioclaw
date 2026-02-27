#!/usr/bin/env python3
"""
Histolab 示例脚本 2: 批量处理多个切片
演示如何批量处理多个 WSI 文件
"""

import os
from pathlib import Path
import json


def batch_process_slides(slide_dir, output_dir, config=None):
    """
    批量处理多个全切片图像
    
    Args:
        slide_dir: 切片文件目录
        output_dir: 输出目录
        config: 处理配置
    """
    try:
        from histolab.slide import Slide
        from histolab.tiler import RandomTiler
        from histolab.masks import BiggestTissueBoxMask
    except ImportError:
        print("需要安装 histolab: pip install histolab")
        return
    
    # 默认配置
    default_config = {
        'tile_size': (512, 512),
        'n_tiles': 50,
        'level': 0,
        'tissue_percent': 80.0,
        'seed': 42
    }
    
    if config:
        default_config.update(config)
    
    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 配置 tiler
    tiler = RandomTiler(
        tile_size=default_config['tile_size'],
        n_tiles=default_config['n_tiles'],
        level=default_config['level'],
        seed=default_config['seed'],
        check_tissue=True,
        tissue_percent=default_config['tissue_percent']
    )
    
    # 支持的文件格式
    slide_extensions = ['.svs', '.tif', '.tiff', '.ndpi', '.vms', '.vmu', '.scn']
    
    # 查找所有切片文件
    slide_files = []
    for ext in slide_extensions:
        slide_files.extend(Path(slide_dir).glob(f"*{ext}"))
    
    print(f"找到 {len(slide_files)} 个切片文件")
    
    # 处理统计
    stats = {
        'processed': 0,
        'failed': 0,
        'total_tiles': 0,
        'slides': []
    }
    
    # 批量处理
    for slide_file in slide_files:
        print(f"\n处理: {slide_file.name}")
        
        try:
            # 创建切片特定的输出目录
            slide_output = output_path / slide_file.stem
            slide_output.mkdir(exist_ok=True)
            
            # 加载幻灯片
            slide = Slide(str(slide_file), processed_path=str(slide_output))
            
            # 保存缩略图
            slide.save_thumbnail()
            
            # 提取瓦片
            tiler.extract(slide)
            
            # 统计
            n_tiles = len(list(slide_output.glob("*.png")))
            
            stats['processed'] += 1
            stats['total_tiles'] += n_tiles
            stats['slides'].append({
                'name': slide_file.name,
                'tiles': n_tiles,
                'dimensions': slide.dimensions,
                'levels': slide.levels
            })
            
            print(f"  ✓ 提取了 {n_tiles} 个瓦片")
            
        except Exception as e:
            stats['failed'] += 1
            print(f"  × 失败: {e}")
            continue
    
    # 保存统计报告
    report_file = output_path / "processing_report.json"
    with open(report_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\n{'='*60}")
    print("批量处理完成!")
    print(f"  成功: {stats['processed']}")
    print(f"  失败: {stats['failed']}")
    print(f"  总瓦片数: {stats['total_tiles']}")
    print(f"  报告: {report_file}")
    
    return stats


def create_processing_pipeline():
    """创建一个完整的处理流程配置"""
    pipeline = {
        'name': '标准组织学处理流程',
        'description': '用于深度学习训练数据的预处理',
        'steps': [
            {
                'name': '质量控制',
                'description': '检查切片质量，排除模糊或损坏的切片',
                'thresholds': {
                    'min_tissue_percent': 50,
                    'max_blur_score': 0.3
                }
            },
            {
                'name': '组织检测',
                'description': '使用 BiggestTissueBoxMask 提取主要组织区域',
                'params': {
                    'mask_type': 'BiggestTissueBoxMask',
                    'padding': 50
                }
            },
            {
                'name': '瓦片提取',
                'description': '随机提取高质量瓦片',
                'tiler': {
                    'type': 'RandomTiler',
                    'tile_size': [512, 512],
                    'n_tiles': 100,
                    'level': 0,
                    'tissue_percent': 80
                }
            },
            {
                'name': '数据增强',
                'description': '可选的数据增强步骤',
                'augmentation': {
                    'rotation': True,
                    'flip': True,
                    'color_jitter': False
                }
            }
        ]
    }
    
    return pipeline


def demo_pipeline():
    """演示完整流程"""
    print("=" * 70)
    print("Histolab 批量处理示例")
    print("=" * 70)
    
    # 显示流程配置
    pipeline = create_processing_pipeline()
    
    print("\n处理流程配置:")
    print("-" * 70)
    print(f"名称: {pipeline['name']}")
    print(f"描述: {pipeline['description']}")
    
    for i, step in enumerate(pipeline['steps'], 1):
        print(f"\n步骤 {i}: {step['name']}")
        print(f"  {step['description']}")
        if 'params' in step:
            print(f"  参数: {step['params']}")
    
    # 演示批量处理
    print("\n" + "=" * 70)
    print("批量处理演示")
    print("=" * 70)
    
    # 检查示例数据目录
    sample_dirs = ["./slides", "./data/slides", "/path/to/slides"]
    
    for slide_dir in sample_dirs:
        if Path(slide_dir).exists():
            print(f"\n找到切片目录: {slide_dir}")
            
            # 配置
            config = {
                'tile_size': (512, 512),
                'n_tiles': 20,
                'level': 0,
                'tissue_percent': 80.0
            }
            
            # 运行批量处理
            stats = batch_process_slides(slide_dir, "./output/batch", config)
            break
    else:
        print("\n未找到切片目录")
        print("\n示例用法:")
        print("-" * 70)
        print("""
# 批量处理脚本
from pathlib import Path

# 配置
config = {
    'tile_size': (512, 512),
    'n_tiles': 50,
    'level': 0,
    'tissue_percent': 80.0,
    'seed': 42
}

# 运行批量处理
stats = batch_process_slides(
    slide_dir="/path/to/your/slides",
    output_dir="./output",
    config=config
)

# 查看统计
print(f"处理完成: {stats['processed']} 个切片")
print(f"总瓦片数: {stats['total_tiles']}")
""")


def main():
    demo_pipeline()
    
    print("\n" + "=" * 70)
    print("提示:")
    print("-" * 70)
    print("1. 确保有足够的磁盘空间 (WSI 文件通常很大)")
    print("2. 根据 GPU 显存调整 batch_size")
    print("3. 使用 ScoreTiler 提取高质量瓦片")
    print("4. 定期检查输出目录空间")
    print("5. 考虑使用多进程加速处理")


if __name__ == "__main__":
    main()
