#!/usr/bin/env python3
"""
gtars 示例脚本 2: 覆盖度分析和 WIG 文件生成
演示如何使用 uniwig 模块生成覆盖度轨道
"""

import subprocess
import sys
from pathlib import Path


def create_fragment_file(filepath):
    """创建示例片段文件 (ATAC-seq 片段)"""
    fragments = [
        # chr, start, end, barcode, count
        ("chr1", 1000, 1100, "cell_A", 1),
        ("chr1", 1020, 1120, "cell_B", 1),
        ("chr1", 1050, 1150, "cell_A", 1),
        ("chr1", 1500, 1600, "cell_C", 1),
        ("chr2", 2000, 2100, "cell_A", 1),
        ("chr2", 2050, 2150, "cell_D", 1),
    ]
    
    with open(filepath, 'w') as f:
        for chrom, start, end, barcode, count in fragments:
            f.write(f"{chrom}\t{start}\t{end}\t{barcode}\t{count}\n")


def main():
    print("=" * 60)
    print("gtars 覆盖度分析示例")
    print("=" * 60)
    
    # 创建示例片段文件
    print("\n1. 创建示例 ATAC-seq 片段文件...")
    create_fragment_file("fragments.bed")
    print("✓ 已创建 fragments.bed")
    
    # 显示文件内容
    print("\n片段文件内容:")
    with open("fragments.bed", 'r') as f:
        for line in f:
            print(f"  {line.strip()}")
    
    # 使用 CLI 生成覆盖度轨道
    print("\n" + "=" * 60)
    print("2. 使用 gtars CLI 生成覆盖度")
    print("=" * 60)
    
    # 注意: 这需要在系统上安装 gtars-cli
    # cargo install gtars-cli
    
    print("\n命令示例 (需要安装 gtars-cli):")
    print("-" * 40)
    print("# 生成 WIG 格式覆盖度")
    print("gtars uniwig generate --input fragments.bed --output coverage.wig --resolution 10")
    print()
    print("# 生成 BigWig 格式覆盖度")
    print("gtars uniwig generate --input fragments.bed --output coverage.bw --format bigwig")
    print("-" * 40)
    
    # 尝试运行 (如果安装了 gtars-cli)
    try:
        result = subprocess.run(
            ["gtars", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            print(f"\n✓ gtars CLI 已安装: {result.stdout.strip()}")
            
            # 生成 WIG 文件
            print("\n3. 生成覆盖度轨道...")
            cmd = [
                "gtars", "uniwig", "generate",
                "--input", "fragments.bed",
                "--output", "coverage.wig",
                "--resolution", "10"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print("✓ 已生成 coverage.wig")
                
                # 显示输出文件
                if Path("coverage.wig").exists():
                    print("\nWIG 文件预览:")
                    with open("coverage.wig", 'r') as f:
                        lines = f.readlines()[:20]
                        for line in lines:
                            print(f"  {line.strip()}")
            else:
                print(f"生成失败: {result.stderr}")
        else:
            print("\n注意: gtars CLI 未安装")
            print("安装命令: cargo install gtars-cli --features 'uniwig'")
            
    except FileNotFoundError:
        print("\n注意: gtars CLI 未安装或不在 PATH 中")
        print("安装命令: cargo install gtars-cli --features 'uniwig'")
    except Exception as e:
        print(f"\nCLI 测试跳过: {e}")
    
    # 示例: Python API 方式 (如果可用)
    print("\n" + "=" * 60)
    print("3. Python API 示例 (需要 gtars Python 包)")
    print("=" * 60)
    
    try:
        import gtars
        
        # 加载片段
        fragments = gtars.RegionSet.from_bed("fragments.bed")
        print(f"✓ 加载了 {len(fragments)} 个片段")
        
        # 计算每个位置的覆盖度 (简化示例)
        coverage = {}
        for region in fragments:
            if region.chromosome not in coverage:
                coverage[region.chromosome] = {}
            
            for pos in range(region.start, region.end):
                coverage[region.chromosome][pos] = coverage[region.chromosome].get(pos, 0) + 1
        
        print("\n覆盖度统计:")
        for chrom in sorted(coverage.keys()):
            positions = coverage[chrom]
            max_cov = max(positions.values()) if positions else 0
            total_positions = len(positions)
            print(f"  {chrom}: {total_positions} 个有覆盖的位置, 最大覆盖度: {max_cov}")
            
    except ImportError:
        print("注意: 需要安装 gtars Python 包")
        print("安装命令: pip install gtars")
    except Exception as e:
        print(f"API 示例跳过: {e}")
    
    # 清理
    for f in ["fragments.bed", "coverage.wig", "coverage.bw"]:
        if Path(f).exists():
            Path(f).unlink()
    
    print("\n✓ 示例完成!")


if __name__ == "__main__":
    main()
