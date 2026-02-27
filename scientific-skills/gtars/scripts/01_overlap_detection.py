#!/usr/bin/env python3
"""
gtars 示例脚本 1: 基因组区间重叠检测
演示如何使用 IGD 索引和重叠检测功能
"""

import gtars
import numpy as np
from pathlib import Path


def create_sample_bed_file(filepath, regions):
    """创建示例 BED 文件"""
    with open(filepath, 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")


def main():
    # 示例 1: 创建示例数据
    print("=" * 60)
    print("示例 1: 创建基因组区间数据")
    print("=" * 60)
    
    # 定义一些基因组区间 (ChIP-seq peaks)
    peaks_data = [
        ("chr1", 1000, 2000),
        ("chr1", 5000, 6000),
        ("chr2", 3000, 4000),
        ("chr1", 1500, 2500),  # 与第一个区间重叠
        ("chr3", 10000, 11000),
    ]
    
    # 定义启动子区域
    promoter_data = [
        ("chr1", 1200, 2200),
        ("chr1", 5500, 6500),
        ("chr2", 3500, 4500),
        ("chr3", 10500, 11500),
    ]
    
    # 创建 BED 文件
    create_sample_bed_file("peaks.bed", peaks_data)
    create_sample_bed_file("promoters.bed", promoter_data)
    
    print("✓ 已创建 peaks.bed 和 promoters.bed")
    
    # 示例 2: 加载基因组区间
    print("\n" + "=" * 60)
    print("示例 2: 加载和检查 RegionSet")
    print("=" * 60)
    
    try:
        peaks = gtars.RegionSet.from_bed("peaks.bed")
        promoters = gtars.RegionSet.from_bed("promoters.bed")
        
        print(f"✓ Peaks 数量: {len(peaks)}")
        print(f"✓ Promoters 数量: {len(promoters)}")
        
        # 遍历区间
        print("\n前3个 peaks:")
        for i, region in enumerate(peaks[:3]):
            print(f"  {region.chromosome}:{region.start}-{region.end}")
            
    except Exception as e:
        print(f"注意: 需要安装 gtars 包才能运行: pip install gtars")
        print(f"错误信息: {e}")
        return
    
    # 示例 3: 重叠检测
    print("\n" + "=" * 60)
    print("示例 3: 检测区间重叠")
    print("=" * 60)
    
    try:
        # 查找与启动子重叠的 peaks
        overlapping_peaks = peaks.filter_overlapping(promoters)
        print(f"✓ 与启动子重叠的 Peaks: {len(overlapping_peaks)}")
        
        for region in overlapping_peaks:
            print(f"  - {region.chromosome}:{region.start}-{region.end}")
        
        # 统计每个染色体的区间数量
        chrom_counts = {}
        for region in peaks:
            chrom_counts[region.chromosome] = chrom_counts.get(region.chromosome, 0) + 1
        
        print(f"\n按染色体统计 Peaks:")
        for chrom, count in sorted(chrom_counts.items()):
            print(f"  {chrom}: {count}")
            
    except Exception as e:
        print(f"重叠检测示例跳过: {e}")
    
    # 示例 4: 导出结果
    print("\n" + "=" * 60)
    print("示例 4: 导出重叠区间")
    print("=" * 60)
    
    try:
        if len(overlapping_peaks) > 0:
            overlapping_peaks.to_bed("overlapping_peaks.bed")
            print("✓ 已导出到 overlapping_peaks.bed")
    except Exception as e:
        print(f"导出跳过: {e}")
    
    # 清理临时文件
    for f in ["peaks.bed", "promoters.bed", "overlapping_peaks.bed"]:
        if Path(f).exists():
            Path(f).unlink()
    
    print("\n✓ 示例完成!")


if __name__ == "__main__":
    main()
