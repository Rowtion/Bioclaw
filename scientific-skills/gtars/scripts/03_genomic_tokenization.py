#!/usr/bin/env python3
"""
gtars 示例脚本 3: 基因组标记化和 ML 预处理
演示如何使用 TreeTokenizer 进行基因组数据的 ML 预处理
"""

import numpy as np
from pathlib import Path


def create_sample_regions_for_tokenization():
    """创建用于标记化的示例基因组区间"""
    # 模拟增强子/启动子区域
    regions = [
        ("chr1", 1000000, 1001000),
        ("chr1", 2000000, 2000800),
        ("chr1", 3000000, 3001200),
        ("chr2", 500000, 501000),
        ("chr2", 800000, 800600),
        ("chr3", 100000, 100800),
    ]
    return regions


def main():
    print("=" * 60)
    print("gtars 基因组标记化示例")
    print("=" * 60)
    
    # 创建示例训练区域
    print("\n1. 创建示例基因组区间...")
    regions = create_sample_regions_for_tokenization()
    
    with open("training_regions.bed", 'w') as f:
        for chrom, start, end in regions:
            f.write(f"{chrom}\t{start}\t{end}\n")
    
    print(f"✓ 已创建 training_regions.bed ({len(regions)} 个区域)")
    
    # 显示区域
    print("\n训练区域:")
    for i, (chrom, start, end) in enumerate(regions, 1):
        length = end - start
        print(f"  {i}. {chrom}:{start:,}-{end:,} (长度: {length}bp)")
    
    # 示例: 使用 TreeTokenizer (需要 gtars 包)
    print("\n" + "=" * 60)
    print("2. TreeTokenizer 示例")
    print("=" * 60)
    
    try:
        from gtars.tokenizers import TreeTokenizer
        
        # 从 BED 文件创建 tokenizer
        print("\n创建 TreeTokenizer...")
        tokenizer = TreeTokenizer.from_bed_file("training_regions.bed")
        print("✓ Tokenizer 创建成功")
        
        # 标记化查询区域
        query_regions = [
            ("chr1", 1000500, 1000800),   # 在训练区域内
            ("chr1", 1500000, 1500100),   # 不在训练区域内
            ("chr2", 500200, 500500),     # 在训练区域内
        ]
        
        print("\n查询区域标记化:")
        for chrom, start, end in query_regions:
            try:
                token = tokenizer.tokenize(chrom, start, end)
                print(f"  {chrom}:{start:,}-{end:,} -> Token ID: {token}")
            except Exception as e:
                print(f"  {chrom}:{start:,}-{end:,} -> 标记化失败: {e}")
        
        # 批量标记化
        print("\n批量标记化:")
        tokens = []
        for chrom, start, end in query_regions:
            try:
                token = tokenizer.tokenize(chrom, start, end)
                tokens.append(token)
            except:
                tokens.append(-1)  # 未知区域
        
        print(f"  Token 序列: {tokens}")
        
    except ImportError:
        print("\n注意: 需要安装 gtars 包")
        print("安装命令: pip install gtars")
        print("\n模拟标记化过程:")
        
        # 模拟简单的标记化逻辑
        token_map = {}
        for i, (chrom, start, end) in enumerate(regions):
            token_map[(chrom, start, end)] = i
        
        print("\n模拟 Token 映射:")
        for i, (chrom, start, end) in enumerate(regions):
            print(f"  Region {i}: {chrom}:{start:,}-{end:,}")
    
    except Exception as e:
        print(f"Tokenization 示例跳过: {e}")
    
    # 示例 3: 与 geniml 集成
    print("\n" + "=" * 60)
    print("3. ML 预处理流程示例")
    print("=" * 60)
    
    print("""
典型的基因组 ML 预处理流程:

1. 数据准备:
   - 加载基因组区间 (BED 文件)
   - 使用 TreeTokenizer 进行标记化
   
2. 创建训练数据:
   ```python
   from gtars.tokenizers import TreeTokenizer
   
   tokenizer = TreeTokenizer.from_bed_file("training_peaks.bed")
   
   # 标记化所有区域
   tokens = [tokenizer.tokenize(r.chrom, r.start, r.end) 
             for r in regions]
   
   # 转换为 numpy 数组用于 ML
   X = np.array(tokens)
   ```

3. 与深度学习框架集成:
   - 使用 geniml 库进行高级 ML 分析
   - 或直接使用 PyTorch/TensorFlow
   
4. 位置编码:
   - 添加染色体位置信息
   - 创建用于 Transformer 的输入嵌入
""")
    
    # 示例: 简单的位置编码
    print("\n简单位置编码示例:")
    
    # 归一化染色体位置
    chrom_lengths = {
        "chr1": 248956422,
        "chr2": 242193529,
        "chr3": 198295559,
    }
    
    for chrom, start, end in regions[:3]:
        mid = (start + end) / 2
        norm_pos = mid / chrom_lengths.get(chrom, 1)
        print(f"  {chrom}:{start:,}-{end:,} -> 归一化位置: {norm_pos:.6f}")
    
    # 清理
    if Path("training_regions.bed").exists():
        Path("training_regions.bed").unlink()
    
    print("\n✓ 示例完成!")
    print("\n参考: 要与 geniml 集成，请参见 https://github.com/databio/geniml")


if __name__ == "__main__":
    main()
