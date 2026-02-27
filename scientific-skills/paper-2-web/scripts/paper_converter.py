#!/usr/bin/env python3
"""
Paper2All 论文转换工具
用于将学术论文转换为网站、视频和海报
"""

import subprocess
from typing import Optional, Dict, List, Union, Tuple
from dataclasses import dataclass
from pathlib import Path
import json
import os


@dataclass
class ConversionConfig:
    """转换配置"""
    input_dir: str
    output_dir: str
    model_choice: int = 1  # 1: GPT-4, 2: GPT-4.1
    generate_website: bool = True
    generate_poster: bool = True
    generate_video: bool = True
    enable_talking_head: bool = False
    poster_width: int = 48  # 英寸
    poster_height: int = 36  # 英寸
    video_duration: Optional[int] = None
    enable_logo_search: bool = False
    
    def to_dict(self) -> Dict:
        return {
            'input_dir': self.input_dir,
            'output_dir': self.output_dir,
            'model_choice': self.model_choice,
            'components': {
                'website': self.generate_website,
                'poster': self.generate_poster,
                'video': self.generate_video
            },
            'poster_size': {
                'width': self.poster_width,
                'height': self.poster_height
            },
            'video': {
                'duration': self.video_duration,
                'talking_head': self.enable_talking_head
            }
        }


class PaperConverter:
    """论文转换器"""
    
    MODEL_OPTIONS = {
        1: 'GPT-4 (推荐，质量与成本平衡)',
        2: 'GPT-4.1 (最新特性，成本较高)'
    }
    
    def __init__(self, config: ConversionConfig):
        """
        初始化论文转换器
        
        Args:
            config: 转换配置
        """
        self.config = config
        self._validate_config()
    
    def _validate_config(self):
        """验证配置"""
        input_path = Path(self.config.input_dir)
        if not input_path.exists():
            raise ValueError(f"输入目录不存在: {self.config.input_dir}")
        
        # 检查输入文件
        tex_files = list(input_path.glob("*.tex"))
        pdf_files = list(input_path.glob("*.pdf"))
        
        if not tex_files and not pdf_files:
            print(f"警告: 未在 {self.config.input_dir} 中找到 .tex 或 .pdf 文件")
    
    def _build_command(self) -> List[str]:
        """
        构建命令行参数
        
        Returns:
            命令列表
        """
        cmd = [
            'python', 'pipeline_all.py',
            '--input-dir', self.config.input_dir,
            '--output-dir', self.config.output_dir,
            '--model-choice', str(self.config.model_choice)
        ]
        
        if self.config.generate_website:
            cmd.append('--generate-website')
        
        if self.config.generate_poster:
            cmd.append('--generate-poster')
            cmd.extend(['--poster-width-inches', str(self.config.poster_width)])
            cmd.extend(['--poster-height-inches', str(self.config.poster_height)])
        
        if self.config.generate_video:
            cmd.append('--generate-video')
        
        if self.config.enable_talking_head:
            cmd.append('--enable-talking-head')
        
        if self.config.enable_logo_search:
            cmd.append('--enable-logo-search')
        
        return cmd
    
    def convert(self, dry_run: bool = False) -> Dict:
        """
        执行转换
        
        Args:
            dry_run: 仅显示命令，不执行
        
        Returns:
            转换结果
        """
        cmd = self._build_command()
        
        print(f"执行命令: {' '.join(cmd)}")
        
        if dry_run:
            return {'status': 'dry_run', 'command': ' '.join(cmd)}
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1小时超时
            )
            
            return {
                'status': 'success' if result.returncode == 0 else 'error',
                'returncode': result.returncode,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
        except subprocess.TimeoutExpired:
            return {'status': 'timeout', 'error': '转换超时'}
        except Exception as e:
            return {'status': 'error', 'error': str(e)}
    
    def preview_config(self):
        """预览配置"""
        print("=" * 60)
        print("Paper2All 转换配置")
        print("=" * 60)
        print()
        print(f"输入目录: {self.config.input_dir}")
        print(f"输出目录: {self.config.output_dir}")
        print(f"模型选择: {self.MODEL_OPTIONS.get(self.config.model_choice, '未知')}")
        print()
        print("生成组件:")
        print(f"  • 网站: {'是' if self.config.generate_website else '否'}")
        print(f"  • 海报: {'是' if self.config.generate_poster else '否'}")
        if self.config.generate_poster:
            print(f"    尺寸: {self.config.poster_width}\" x {self.config.poster_height}\"")
        print(f"  • 视频: {'是' if self.config.generate_video else '否'}")
        if self.config.generate_video and self.config.enable_talking_head:
            print("    包含虚拟主播")
        print()
        print("=" * 60)


class BatchConverter:
    """批量转换器"""
    
    def __init__(self, base_config: ConversionConfig):
        """
        初始化批量转换器
        
        Args:
            base_config: 基础配置
        """
        self.base_config = base_config
    
    def convert_multiple(self, paper_dirs: List[str], 
                        parallel: bool = False) -> List[Dict]:
        """
        批量转换多个论文
        
        Args:
            paper_dirs: 论文目录列表
            parallel: 是否并行处理
        
        Returns:
            转换结果列表
        """
        results = []
        
        for paper_dir in paper_dirs:
            config = ConversionConfig(
                input_dir=paper_dir,
                output_dir=f"{self.base_config.output_dir}/{Path(paper_dir).name}",
                model_choice=self.base_config.model_choice,
                generate_website=self.base_config.generate_website,
                generate_poster=self.base_config.generate_poster,
                generate_video=self.base_config.generate_video
            )
            
            converter = PaperConverter(config)
            result = converter.convert()
            result['paper'] = paper_dir
            results.append(result)
        
        return results


class PaperOutputChecker:
    """论文输出生成检查器"""
    
    def __init__(self, output_dir: str):
        """
        初始化检查器
        
        Args:
            output_dir: 输出目录
        """
        self.output_dir = Path(output_dir)
    
    def check_outputs(self) -> Dict:
        """
        检查生成的输出文件
        
        Returns:
            检查结果
        """
        results = {
            'website': self._check_website(),
            'poster': self._check_poster(),
            'video': self._check_video()
        }
        return results
    
    def _check_website(self) -> Dict:
        """检查网站输出"""
        website_dir = self.output_dir / 'website'
        if not website_dir.exists():
            return {'exists': False, 'files': []}
        
        files = list(website_dir.glob('*'))
        return {
            'exists': True,
            'files': [f.name for f in files],
            'has_index': (website_dir / 'index.html').exists()
        }
    
    def _check_poster(self) -> Dict:
        """检查海报输出"""
        poster_dir = self.output_dir / 'poster'
        if not poster_dir.exists():
            return {'exists': False, 'files': []}
        
        pdf_files = list(poster_dir.glob('*.pdf'))
        png_files = list(poster_dir.glob('*.png'))
        
        return {
            'exists': True,
            'pdf_files': [f.name for f in pdf_files],
            'png_files': [f.name for f in png_files]
        }
    
    def _check_video(self) -> Dict:
        """检查视频输出"""
        video_dir = self.output_dir / 'video'
        if not video_dir.exists():
            return {'exists': False, 'files': []}
        
        video_files = list(video_dir.glob('*.mp4'))
        
        return {
            'exists': True,
            'video_files': [f.name for f in video_files]
        }
    
    def print_report(self):
        """打印检查报告"""
        results = self.check_outputs()
        
        print("=" * 60)
        print("Paper2All 输出检查报告")
        print("=" * 60)
        print()
        
        for component, status in results.items():
            print(f"{component.upper()}:")
            if status['exists']:
                print(f"  状态: ✓ 已生成")
                if 'files' in status:
                    print(f"  文件数: {len(status['files'])}")
                if status.get('has_index'):
                    print(f"  包含 index.html")
            else:
                print(f"  状态: ✗ 未生成")
            print()
        
        print("=" * 60)


def quick_convert(paper_dir: str, 
                 output_dir: Optional[str] = None,
                 website: bool = True,
                 poster: bool = True,
                 video: bool = False) -> Dict:
    """
    快速转换论文
    
    Args:
        paper_dir: 论文目录
        output_dir: 输出目录 (默认: {paper_dir}_output)
        website: 是否生成网站
        poster: 是否生成海报
        video: 是否生成视频
    
    Returns:
        转换结果
    """
    if output_dir is None:
        output_dir = f"{paper_dir}_output"
    
    config = ConversionConfig(
        input_dir=paper_dir,
        output_dir=output_dir,
        generate_website=website,
        generate_poster=poster,
        generate_video=video
    )
    
    converter = PaperConverter(config)
    return converter.convert()


def demo():
    """演示功能"""
    print("=" * 60)
    print("Paper2All 论文转换工具演示")
    print("=" * 60)
    
    # 1. 配置示例
    print("\n1. 转换配置示例:")
    config = ConversionConfig(
        input_dir="./my_paper",
        output_dir="./output",
        model_choice=1,
        generate_website=True,
        generate_poster=True,
        generate_video=False,
        poster_width=48,
        poster_height=36
    )
    
    converter = PaperConverter(config)
    converter.preview_config()
    
    # 2. 命令预览
    print("\n2. 生成的命令预览:")
    cmd = converter._build_command()
    print(f"   {' '.join(cmd)}")
    
    # 3. 输出检查示例
    print("\n3. 输出检查示例:")
    print("   PaperOutputChecker('./output').print_report()")
    
    # 4. 批量转换示例
    print("\n4. 批量转换示例:")
    print("   paper_dirs = ['./paper1', './paper2', './paper3']")
    print("   batch = BatchConverter(config)")
    print("   results = batch.convert_multiple(paper_dirs)")
    
    print("\n" + "=" * 60)
    print("\n注意: 此工具需要 Paper2All 项目环境")
    print("GitHub: https://github.com/YuhangChen1/Paper2All")
    print("=" * 60)


if __name__ == "__main__":
    demo()
