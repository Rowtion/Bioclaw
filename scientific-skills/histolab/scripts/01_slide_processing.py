#!/usr/bin/env python3
"""
Histolab 示例脚本 1: 全切片图像处理和瓦片提取
演示如何加载 WSI、检测组织区域并提取图像瓦片
"""

import os
from pathlib import Path


def check_installation():
    """检查必要的包是否安装"""
    try:
        from histolab.slide import Slide
        from histolab.tiler import RandomTiler
        print("✓ histolab 已安装")
        return True
    except ImportError:
        print("× histolab 未安装")
        print("  安装命令: pip install histolab")
        return False


def demo_slide_management():
    """演示幻灯片加载和管理"""
    print("\n" + "=" * 60)
    print("1. 幻灯片加载和管理")
    print("=" * 60)
    
    try:
        from histolab.slide import Slide
        from histolab.data import prostate_tissue
        
        # 使用内置示例数据
        print("\n加载示例前列腺组织...")
        prostate_svs, prostate_path = prostate_tissue()
        
        # 初始化幻灯片
        slide = Slide(prostate_path, processed_path="./output")
        
        print(f"✓ 幻灯片名称: {slide.name}")
        print(f"  尺寸: {slide.dimensions}")
        print(f"  金字塔层级: {slide.levels}")
        
        # 访问属性
        if 'openslide.objective-power' in slide.properties:
            mag = slide.properties['openslide.objective-power']
            print(f"  物镜倍率: {mag}x")
        
        # 保存缩略图
        slide.save_thumbnail()
        print(f"✓ 缩略图已保存")
        
        return slide
        
    except Exception as e:
        print(f"演示错误: {e}")
        return None


def demo_tissue_detection():
    """演示组织检测"""
    print("\n" + "=" * 60)
    print("2. 组织检测和掩码")
    print("=" * 60)
    
    try:
        from histolab.slide import Slide
        from histolab.masks import TissueMask, BiggestTissueBoxMask
        from histolab.data import prostate_tissue
        
        # 加载幻灯片
        prostate_svs, prostate_path = prostate_tissue()
        slide = Slide(prostate_path, processed_path="./output")
        
        # 创建组织掩码 (所有组织区域)
        print("\n创建 TissueMask...")
        tissue_mask = TissueMask()
        
        # 可视化掩码
        try:
            slide.locate_mask(tissue_mask)
            print("✓ 组织掩码已可视化")
        except:
            print("  (可视化需要显示环境)")
        
        # 获取掩码数组
        mask_array = tissue_mask(slide)
        print(f"✓ 掩码数组形状: {mask_array.shape}")
        print(f"  组织像素比例: {mask_array.sum() / mask_array.size:.2%}")
        
        # 使用最大组织区域掩码
        print("\n使用 BiggestTissueBoxMask...")
        biggest_mask = BiggestTissueBoxMask()
        biggest_mask_array = biggest_mask(slide)
        print(f"✓ 最大组织区域掩码形状: {biggest_mask_array.shape}")
        
        return tissue_mask, biggest_mask
        
    except Exception as e:
        print(f"演示错误: {e}")
        return None, None


def demo_tile_extraction():
    """演示瓦片提取"""
    print("\n" + "=" * 60)
    print("3. 图像瓦片提取")
    print("=" * 60)
    
    try:
        from histolab.slide import Slide
        from histolab.tiler import RandomTiler, GridTiler, ScoreTiler
        from histolab.scorer import NucleiScorer
        from histolab.data import prostate_tissue
        
        # 加载幻灯片
        prostate_svs, prostate_path = prostate_tissue()
        slide = Slide(prostate_path, processed_path="./output")
        
        # 方法 1: 随机瓦片提取
        print("\n方法一: RandomTiler (随机采样)")
        random_tiler = RandomTiler(
            tile_size=(512, 512),
            n_tiles=10,
            level=0,
            seed=42,
            check_tissue=True,
            tissue_percent=80.0
        )
        
        print("  配置:")
        print(f"    - 瓦片大小: 512x512")
        print(f"    - 提取数量: 10")
        print(f"    - 组织阈值: 80%")
        
        # 预览瓦片位置
        try:
            random_tiler.locate_tiles(slide, n_tiles=5)
            print("  ✓ 瓦片位置预览已生成")
        except:
            print("  (预览需要显示环境)")
        
        # 提取瓦片
        print("  正在提取随机瓦片...")
        random_tiler.extract(slide)
        print("  ✓ 随机瓦片提取完成")
        
        # 方法 2: 网格瓦片提取
        print("\n方法二: GridTiler (系统网格)")
        grid_tiler = GridTiler(
            tile_size=(512, 512),
            level=1,  # 使用较低分辨率
            pixel_overlap=0,
            check_tissue=True,
            tissue_percent=70.0
        )
        
        print("  提取网格瓦片...")
        # grid_tiler.extract(slide)  # 提取所有网格瓦片
        print("  ✓ 网格配置完成 (实际提取需要更多时间)")
        
        # 方法 3: 基于评分的瓦片提取
        print("\n方法三: ScoreTiler (质量驱动)")
        score_tiler = ScoreTiler(
            tile_size=(512, 512),
            n_tiles=5,
            scorer=NucleiScorer(),
            level=0,
            check_tissue=True
        )
        
        print("  使用 NucleiScorer 评分...")
        # score_tiler.extract(slide, report_path="tiles_report.csv")
        print("  ✓ 评分配置完成 (实际提取需要更多时间)")
        
    except Exception as e:
        print(f"演示错误: {e}")
        import traceback
        traceback.print_exc()


def demo_filters():
    """演示图像滤波器"""
    print("\n" + "=" * 60)
    print("4. 图像滤波和预处理")
    print("=" * 60)
    
    try:
        from histolab.filters.compositions import Compose
        from histolab.filters.image_filters import RgbToGrayscale, OtsuThreshold
        from histolab.filters.morphological_filters import (
            BinaryDilation, RemoveSmallHoles, RemoveSmallObjects
        )
        
        # 标准组织检测流程
        print("\n创建组织检测滤波器链...")
        tissue_detection = Compose([
            RgbToGrayscale(),
            OtsuThreshold(),
            BinaryDilation(disk_size=5),
            RemoveSmallHoles(area_threshold=1000),
            RemoveSmallObjects(area_threshold=500)
        ])
        
        print("✓ 滤波器链创建完成:")
        print("  1. RGB 转灰度")
        print("  2. Otsu 阈值分割")
        print("  3. 二值膨胀")
        print("  4. 移除小孔洞")
        print("  5. 移除小物体")
        
        # 可以在自定义掩码中使用
        from histolab.masks import TissueMask
        custom_mask = TissueMask(filters=tissue_detection)
        print("\n✓ 自定义掩码已创建")
        
    except Exception as e:
        print(f"演示错误: {e}")


def main():
    print("=" * 70)
    print("Histolab 全切片图像处理示例")
    print("=" * 70)
    
    # 检查安装
    if not check_installation():
        print("\n请先安装 histolab:")
        print("  pip install histolab")
        print("  pip install matplotlib pillow")  # 用于可视化
        return
    
    # 创建输出目录
    Path("./output").mkdir(exist_ok=True)
    
    try:
        # 运行各个演示
        slide = demo_slide_management()
        tissue_mask, biggest_mask = demo_tissue_detection()
        demo_tile_extraction()
        demo_filters()
        
        print("\n" + "=" * 70)
        print("演示完成!")
        print("=" * 70)
        print("\n输出文件保存在 ./output 目录")
        print("\n建议下一步:")
        print("1. 查看提取的瓦片")
        print("2. 调整 tissue_percent 参数优化结果")
        print("3. 尝试不同的 tiler 策略")
        print("4. 使用 ScoreTiler 进行质量评估")
        
    except Exception as e:
        print(f"\n运行错误: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
