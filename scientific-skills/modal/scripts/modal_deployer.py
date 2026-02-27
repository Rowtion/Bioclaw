#!/usr/bin/env python3
"""
Modal 云函数部署工具
用于部署和管理 Modal 无服务器应用
"""

import modal
from typing import Optional, List, Dict, Any, Callable
from pathlib import Path
import tempfile
import os


class ModalAppBuilder:
    """Modal 应用构建器"""
    
    def __init__(self, app_name: str, python_version: str = "3.12"):
        """
        初始化 Modal 应用
        
        Args:
            app_name: 应用名称
            python_version: Python 版本
        """
        self.app_name = app_name
        self.python_version = python_version
        self.packages: List[str] = []
        self.image = None
        self.app = None
        self._build_image()
    
    def _build_image(self):
        """构建基础镜像"""
        self.image = modal.Image.debian_slim(
            python_version=self.python_version
        )
        self.app = modal.App(self.app_name, image=self.image)
    
    def add_packages(self, *packages: str) -> 'ModalAppBuilder':
        """
        添加 Python 包
        
        Args:
            *packages: 包名列表
        
        Returns:
            self (链式调用)
        """
        self.packages.extend(packages)
        self.image = self.image.uv_pip_install(*packages)
        self.app = modal.App(self.app_name, image=self.image)
        return self
    
    def add_system_packages(self, *packages: str) -> 'ModalAppBuilder':
        """
        添加系统包
        
        Args:
            *packages: 系统包名
        
        Returns:
            self (链式调用)
        """
        self.image = self.image.apt_install(*packages)
        self.app = modal.App(self.app_name, image=self.image)
        return self
    
    def add_local_code(self, local_path: str, remote_path: str = "/root") -> 'ModalAppBuilder':
        """
        添加本地代码
        
        Args:
            local_path: 本地路径
            remote_path: 远程路径
        
        Returns:
            self (链式调用)
        """
        self.image = self.image.add_local_dir(local_path, remote_path)
        self.app = modal.App(self.app_name, image=self.image)
        return self
    
    def create_function(self, 
                       func: Callable = None,
                       gpu: Optional[str] = None,
                       cpu: float = 0.125,
                       memory: int = 128,
                       timeout: int = 300,
                       secrets: Optional[List] = None) -> Callable:
        """
        创建 Modal 函数装饰器
        
        Args:
            func: 要包装的函数
            gpu: GPU 类型 (T4, L4, A10, A100, H100 等)
            cpu: CPU 核心数
            memory: 内存 (MiB)
            timeout: 超时时间 (秒)
            secrets: Secret 列表
        
        Returns:
            装饰后的函数
        """
        decorator = self.app.function(
            gpu=gpu,
            cpu=cpu,
            memory=memory,
            timeout=timeout,
            secrets=secrets or []
        )
        
        if func is None:
            return decorator
        return decorator(func)
    
    def create_web_endpoint(self,
                           func: Callable = None,
                           method: str = "POST",
                           gpu: Optional[str] = None) -> Callable:
        """
        创建 Web 端点
        
        Args:
            func: 处理函数
            method: HTTP 方法
            gpu: GPU 类型
        
        Returns:
            装饰后的函数
        """
        decorator = self.app.function(gpu=gpu)(
            modal.web_endpoint(method=method)(func) if func else modal.web_endpoint(method=method)
        )
        return decorator


class ModalMLDeployer:
    """Modal ML 模型部署器"""
    
    def __init__(self, model_name: str):
        """
        初始化模型部署器
        
        Args:
            model_name: 模型名称
        """
        self.model_name = model_name
        self.builder = ModalAppBuilder(f"ml-{model_name}")
        self.volume = None
    
    def setup_gpu(self, gpu_type: str = "L4") -> 'ModalMLDeployer':
        """
        配置 GPU
        
        Args:
            gpu_type: GPU 类型
        
        Returns:
            self (链式调用)
        """
        self.gpu_type = gpu_type
        return self
    
    def add_ml_packages(self, framework: str = "pytorch") -> 'ModalMLDeployer':
        """
        添加 ML 框架包
        
        Args:
            framework: 框架名称 (pytorch, tensorflow, transformers)
        
        Returns:
            self (链式调用)
        """
        if framework == "pytorch":
            self.builder.add_packages("torch", "torchvision")
        elif framework == "tensorflow":
            self.builder.add_packages("tensorflow")
        elif framework == "transformers":
            self.builder.add_packages("transformers", "torch", "accelerate")
        return self
    
    def create_volume(self, volume_name: str) -> 'ModalMLDeployer':
        """
        创建持久化存储卷
        
        Args:
            volume_name: 卷名称
        
        Returns:
            self (链式调用)
        """
        self.volume = modal.Volume.from_name(volume_name, create_if_missing=True)
        return self
    
    def deploy(self, inference_func: Callable) -> str:
        """
        部署模型
        
        Args:
            inference_func: 推理函数
        
        Returns:
            部署信息
        """
        # 这里返回部署代码示例
        return f"""
# 使用以下代码部署模型:

{self.builder.app}

@self.builder.app.cls(gpu="{self.gpu_type}")
class {self.model_name.title()}Model:
    @modal.enter()
    def load_model(self):
        # 加载模型代码
        pass
    
    @modal.method()
    def predict(self, input_data):
        # 推理代码
        return inference_func(input_data)

# 部署命令: modal deploy script.py
"""


def create_batch_processor(app_name: str, 
                          process_func: Callable,
                          cpu: float = 2.0,
                          memory: int = 4096) -> modal.App:
    """
    创建批处理处理器
    
    Args:
        app_name: 应用名称
        process_func: 处理函数
        cpu: CPU 核心数
        memory: 内存 (MiB)
    
    Returns:
        Modal App
    """
    app = modal.App(app_name)
    
    @app.function(cpu=cpu, memory=memory)
    def process_item(item):
        return process_func(item)
    
    @app.local_entrypoint()
    def main(items: List[Any]):
        # 并行处理
        results = list(process_item.map(items))
        return results
    
    return app


def create_scheduled_job(app_name: str,
                        job_func: Callable,
                        schedule: str,
                        cpu: float = 1.0,
                        memory: int = 2048) -> modal.App:
    """
    创建定时任务
    
    Args:
        app_name: 应用名称
        job_func: 任务函数
        schedule: cron 表达式
        cpu: CPU 核心数
        memory: 内存 (MiB)
    
    Returns:
        Modal App
    """
    app = modal.App(app_name)
    
    @app.function(
        schedule=modal.Cron(schedule),
        cpu=cpu,
        memory=memory
    )
    def scheduled_task():
        job_func()
        return "Job completed"
    
    return app


def quick_deploy_function(func: Callable,
                         app_name: str = "quick-app",
                         gpu: Optional[str] = None,
                         packages: List[str] = None) -> str:
    """
    快速部署函数
    
    Args:
        func: 要部署的函数
        app_name: 应用名称
        gpu: GPU 类型
        packages: 依赖包列表
    
    Returns:
        部署代码字符串
    """
    pkg_str = ", ".join(f'"{p}"' for p in (packages or []))
    
    code = f'''
import modal

image = modal.Image.debian_slim().uv_pip_install({pkg_str})
app = modal.App("{app_name}", image=image)

@app.function(gpu={f'"{gpu}"' if gpu else None})
def {func.__name__}(*args, **kwargs):
    # 导入依赖
    return func(*args, **kwargs)

@app.local_entrypoint()
def main():
    result = {func.__name__}.remote()
    print(result)

# 运行: modal run script.py
# 部署: modal deploy script.py
'''
    return code


# 预定义配置模板
GPU_PRESETS = {
    "inference_small": {"gpu": "T4", "cpu": 2.0, "memory": 4096},
    "inference_medium": {"gpu": "L4", "cpu": 4.0, "memory": 8192},
    "inference_large": {"gpu": "A100", "cpu": 8.0, "memory": 32768},
    "training_small": {"gpu": "A100", "cpu": 8.0, "memory": 32768},
    "training_large": {"gpu": "H100:4", "cpu": 32.0, "memory": 131072},
}


def get_gpu_preset(preset_name: str) -> Dict[str, Any]:
    """
    获取 GPU 预设配置
    
    Args:
        preset_name: 预设名称
    
    Returns:
        配置字典
    """
    return GPU_PRESETS.get(preset_name, GPU_PRESETS["inference_medium"])


def demo():
    """演示功能"""
    print("=" * 60)
    print("Modal 部署工具示例")
    print("=" * 60)
    
    # 1. 基本应用构建器
    print("\n1. 创建 Modal 应用:")
    builder = ModalAppBuilder("my-app")
    builder.add_packages("pandas", "numpy", "scikit-learn")
    print(f"   应用名称: my-app")
    print(f"   已添加包: pandas, numpy, scikit-learn")
    
    # 2. ML 部署器
    print("\n2. ML 模型部署配置:")
    ml_deployer = ModalMLDeployer("sentiment-classifier")
    ml_deployer.setup_gpu("L4").add_ml_packages("transformers")
    print(f"   模型: sentiment-classifier")
    print(f"   GPU: L4")
    print(f"   框架: transformers")
    
    # 3. GPU 预设
    print("\n3. GPU 预设:")
    for name, config in GPU_PRESETS.items():
        print(f"   {name}: GPU={config['gpu']}, CPU={config['cpu']}, Memory={config['memory']}MiB")
    
    # 4. 生成示例代码
    print("\n4. 快速部署示例代码:")
    
    def example_inference(text: str) -> str:
        return f"Processed: {text}"
    
    code = quick_deploy_function(
        example_inference,
        app_name="inference-app",
        gpu="T4",
        packages=["torch", "transformers"]
    )
    print(code[:500] + "...")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
