#!/usr/bin/env python3
"""
PathML 病理图像分析工具
用于全切片图像 (WSI) 分析和机器学习
"""

from typing import Optional, Dict, List, Tuple, Union, Any, Callable
from dataclasses import dataclass
from pathlib import Path
import numpy as np

# 尝试导入 pathml
try:
    from pathml.core import SlideData
    from pathml.preprocessing import Pipeline, StainNormalizationHE, TissueDetectionHE
    PATHML_AVAILABLE = True
except ImportError:
    PATHML_AVAILABLE = False
    print("警告: pathml 未安装。运行: uv pip install pathml")

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


@dataclass
class SlideInfo:
    """切片信息"""
    name: str
    dimensions: Tuple[int, int]
    level_count: int
    level_dimensions: List[Tuple[int, int]]
    level_downsamples: List[float]
    
    def to_dict(self) -> Dict:
        return {
            'name': self.name,
            'dimensions': self.dimensions,
            'level_count': self.level_count,
            'level_dimensions': self.level_dimensions,
            'level_downsamples': self.level_downsamples
        }


class SlideLoader:
    """切片加载器"""
    
    SUPPORTED_FORMATS = [
        '.svs', '.ndpi', '.scn', '.zvi', '.czi', 
        '.tiff', '.tif', '.vms', '.vmu', '.dcm'
    ]
    
    def __init__(self, slide_path: str):
        """
        初始化切片加载器
        
        Args:
            slide_path: 切片文件路径
        """
        if not PATHML_AVAILABLE:
            raise ImportError("pathml 未安装。运行: uv pip install pathml")
        
        self.slide_path = Path(slide_path)
        if not self.slide_path.exists():
            raise FileNotFoundError(f"切片文件不存在: {slide_path}")
        
        self.slide: Optional[SlideData] = None
    
    def load(self) -> SlideData:
        """
        加载切片
        
        Returns:
            SlideData 对象
        """
        self.slide = SlideData.from_slide(str(self.slide_path))
        return self.slide
    
    def get_info(self) -> SlideInfo:
        """
        获取切片信息
        
        Returns:
            切片信息
        """
        if self.slide is None:
            self.load()
        
        # 提取信息
        shape = self.slide.shape
        
        return SlideInfo(
            name=self.slide_path.name,
            dimensions=(shape[1], shape[0]),
            level_count=1,  # 简化为单级别
            level_dimensions=[(shape[1], shape[0])],
            level_downsamples=[1.0]
        )
    
    def get_thumbnail(self, size: Tuple[int, int] = (512, 512)) -> np.ndarray:
        """
        获取缩略图
        
        Args:
            size: 缩略图大小
        
        Returns:
            缩略图数组
        """
        if self.slide is None:
            self.load()
        
        # 获取整张图像并缩放
        region = self.slide.extract_region(location=(0, 0), 
                                           size=self.slide.shape[:2])
        
        # 简单缩放
        from PIL import Image
        img = Image.fromarray(region)
        img.thumbnail(size)
        return np.array(img)


class SlidePreprocessor:
    """切片预处理器"""
    
    def __init__(self, slide: SlideData):
        """
        初始化预处理器
        
        Args:
            slide: SlideData 对象
        """
        if not PATHML_AVAILABLE:
            raise ImportError("pathml 未安装")
        
        self.slide = slide
        self.pipeline: Optional[Pipeline] = None
    
    def create_standard_pipeline(self, stain_method: str = 'macenko') -> Pipeline:
        """
        创建标准预处理管道
        
        Args:
            stain_method: 染色归一化方法 ('macenko', 'vahadane')
        
        Returns:
            Pipeline 对象
        """
        if not PATHML_AVAILABLE:
            raise ImportError("pathml 未安装")
        
        self.pipeline = Pipeline([
            TissueDetectionHE(),
            StainNormalizationHE(target='normalize', 
                               stain_estimation_method=stain_method)
        ])
        
        return self.pipeline
    
    def run(self, **kwargs):
        """
        运行预处理管道
        
        Args:
            **kwargs: 传递给 pipeline.run 的参数
        """
        if self.pipeline is None:
            self.create_standard_pipeline()
        
        self.pipeline.run(self.slide, **kwargs)
    
    def extract_tiles(self, tile_size: int = 256, 
                     stride: int = 256,
                     level: int = 0) -> List[np.ndarray]:
        """
        提取图像块
        
        Args:
            tile_size: 块大小
            stride: 步长
            level: 金字塔级别
        
        Returns:
            图像块列表
        """
        tiles = []
        shape = self.slide.shape
        
        # 计算块位置
        for y in range(0, shape[0] - tile_size + 1, stride):
            for x in range(0, shape[1] - tile_size + 1, stride):
                tile = self.slide.extract_region(
                    location=(y, x),
                    size=(tile_size, tile_size)
                )
                tiles.append(tile)
        
        return tiles


class TileAnalyzer:
    """图像块分析器"""
    
    def __init__(self):
        """初始化分析器"""
        pass
    
    def analyze_intensity(self, tile: np.ndarray) -> Dict:
        """
        分析图像块强度
        
        Args:
            tile: 图像块
        
        Returns:
            强度统计
        """
        if tile.ndim == 3:
            # 多通道图像
            results = {}
            for i, channel_name in enumerate(['R', 'G', 'B']):
                channel = tile[:, :, i]
                results[channel_name] = {
                    'mean': float(np.mean(channel)),
                    'std': float(np.std(channel)),
                    'min': float(np.min(channel)),
                    'max': float(np.max(channel))
                }
            return results
        else:
            # 单通道
            return {
                'mean': float(np.mean(tile)),
                'std': float(np.std(tile)),
                'min': float(np.min(tile)),
                'max': float(np.max(tile))
            }
    
    def detect_tissue(self, tile: np.ndarray, 
                     threshold: float = 0.1) -> bool:
        """
        检测是否包含组织
        
        Args:
            tile: 图像块
            threshold: 组织比例阈值
        
        Returns:
            是否包含组织
        """
        # 简单阈值检测
        if tile.ndim == 3:
            gray = np.mean(tile, axis=2)
        else:
            gray = tile
        
        # 计算非背景像素比例
        tissue_ratio = np.mean(gray < 240)
        return tissue_ratio > threshold
    
    def calculate_histogram(self, tile: np.ndarray, 
                           bins: int = 256) -> Dict:
        """
        计算颜色直方图
        
        Args:
            tile: 图像块
            bins: 直方图 bin 数
        
        Returns:
            直方图数据
        """
        if tile.ndim != 3:
            tile = np.stack([tile] * 3, axis=-1)
        
        histograms = {}
        colors = ['red', 'green', 'blue']
        
        for i, color in enumerate(colors):
            hist, bin_edges = np.histogram(tile[:, :, i], bins=bins, range=(0, 256))
            histograms[color] = {
                'histogram': hist.tolist(),
                'bin_edges': bin_edges.tolist()
            }
        
        return histograms


class SlideVisualizer:
    """切片可视化器"""
    
    def __init__(self, slide: SlideData):
        """
        初始化可视化器
        
        Args:
            slide: SlideData 对象
        """
        self.slide = slide
    
    def show_overview(self, figsize: Tuple[int, int] = (12, 10),
                     save_path: Optional[str] = None):
        """
        显示切片概览
        
        Args:
            figsize: 图形大小
            save_path: 保存路径
        """
        if not MATPLOTLIB_AVAILABLE:
            print("错误: matplotlib 未安装")
            return
        
        loader = SlideLoader(self.slide.path) if hasattr(self.slide, 'path') else None
        if loader:
            thumb = loader.get_thumbnail((1024, 1024))
        else:
            # 直接获取
            shape = self.slide.shape
            thumb = self.slide.extract_region(
                location=(0, 0),
                size=(min(shape[0], 1024), min(shape[1], 1024))
            )
        
        plt.figure(figsize=figsize)
        plt.imshow(thumb)
        plt.title('Slide Overview')
        plt.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"概览图已保存到: {save_path}")
        else:
            plt.show()
    
    def show_tile(self, location: Tuple[int, int], 
                 size: Tuple[int, int] = (256, 256),
                 save_path: Optional[str] = None):
        """
        显示特定区域的图像块
        
        Args:
            location: 位置 (y, x)
            size: 大小 (height, width)
            save_path: 保存路径
        """
        if not MATPLOTLIB_AVAILABLE:
            print("错误: matplotlib 未安装")
            return
        
        tile = self.slide.extract_region(location=location, size=size)
        
        plt.figure(figsize=(8, 8))
        plt.imshow(tile)
        plt.title(f'Tile at {location}, size {size}')
        plt.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
        else:
            plt.show()


class BatchProcessor:
    """批量处理器"""
    
    def __init__(self, output_dir: str):
        """
        初始化批量处理器
        
        Args:
            output_dir: 输出目录
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def process_slides(self, slide_paths: List[str], 
                      processor_func: Callable[[SlideData], Any]) -> Dict[str, Any]:
        """
        批量处理切片
        
        Args:
            slide_paths: 切片路径列表
            processor_func: 处理函数
        
        Returns:
            处理结果字典
        """
        results = {}
        
        for slide_path in slide_paths:
            try:
                loader = SlideLoader(slide_path)
                slide = loader.load()
                result = processor_func(slide)
                results[slide_path] = {
                    'status': 'success',
                    'result': result
                }
            except Exception as e:
                results[slide_path] = {
                    'status': 'error',
                    'error': str(e)
                }
        
        return results


def demo():
    """演示功能"""
    print("=" * 60)
    print("PathML 病理图像分析工具演示")
    print("=" * 60)
    
    if not PATHML_AVAILABLE:
        print("\n错误: pathml 未安装")
        print("请运行: uv pip install pathml")
        print("\n注意: 此工具需要 Java 环境支持")
        return
    
    print("\nPathML 工具可用功能:")
    print()
    print("1. SlideLoader - 切片加载:")
    print("   • load() - 加载 WSI 文件")
    print("   • get_info() - 获取切片信息")
    print("   • get_thumbnail() - 生成缩略图")
    print()
    print("2. SlidePreprocessor - 预处理:")
    print("   • create_standard_pipeline() - 创建标准管道")
    print("   • extract_tiles() - 提取图像块")
    print("   • run() - 执行预处理")
    print()
    print("3. TileAnalyzer - 图像块分析:")
    print("   • analyze_intensity() - 强度分析")
    print("   • detect_tissue() - 组织检测")
    print("   • calculate_histogram() - 颜色直方图")
    print()
    print("4. SlideVisualizer - 可视化:")
    print("   • show_overview() - 显示概览")
    print("   • show_tile() - 显示图像块")
    print()
    print("5. BatchProcessor - 批量处理:")
    print("   • process_slides() - 批量处理多个切片")
    print()
    print("支持的格式:", ", ".join(SlideLoader.SUPPORTED_FORMATS))
    print()
    print("=" * 60)


if __name__ == "__main__":
    demo()
