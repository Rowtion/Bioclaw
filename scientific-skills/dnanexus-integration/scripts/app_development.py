#!/usr/bin/env python3
"""
DNAnexus 应用开发辅助脚本
功能：
1. 创建应用模板
2. 验证dxapp.json配置
3. 构建和部署应用
4. 测试应用

依赖：dxpy
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Optional


def create_app_template(
    app_name: str,
    output_dir: str,
    language: str = "python",
    inputs: Optional[List[Dict]] = None,
    outputs: Optional[List[Dict]] = None
):
    """
    创建DNAnexus应用模板
    
    Args:
        app_name: 应用名称
        output_dir: 输出目录
        language: 编程语言 (python, bash)
        inputs: 输入定义列表
        outputs: 输出定义列表
    """
    app_dir = Path(output_dir) / app_name
    src_dir = app_dir / "src"
    src_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"创建应用模板: {app_name}")
    
    # 默认输入输出
    if inputs is None:
        inputs = [
            {
                "name": "input_file",
                "label": "Input File",
                "class": "file",
                "optional": False,
                "patterns": ["*.txt", "*.fastq", "*.bam"]
            },
            {
                "name": "param_value",
                "label": "Parameter Value",
                "class": "int",
                "default": 10,
                "optional": True
            }
        ]
    
    if outputs is None:
        outputs = [
            {
                "name": "output_file",
                "label": "Output File",
                "class": "file"
            }
        ]
    
    # 创建dxapp.json
    dxapp = {
        "name": app_name,
        "title": app_name.replace("_", " ").title(),
        "summary": f"Analysis app: {app_name}",
        "description": "Add your description here",
        "version": "0.0.1",
        "inputSpec": inputs,
        "outputSpec": outputs,
        "runSpec": {
            "interpreter": language,
            "file": f"src/{app_name}.py" if language == "python" else f"src/{app_name}.sh",
            "distribution": "Ubuntu",
            "release": "20.04"
        },
        "developers": [],
        "authorizedUsers": []
    }
    
    dxapp_path = app_dir / "dxapp.json"
    with open(dxapp_path, 'w') as f:
        json.dump(dxapp, f, indent=2)
    print(f"  ✓ {dxapp_path}")
    
    # 创建主脚本
    if language == "python":
        script_content = f'''#!/usr/bin/env python3
import dxpy
import subprocess

@dxpy.entry_point('main')
def main(input_file, param_value=10):
    """
    主入口点
    
    Args:
        input_file: 输入文件
        param_value: 参数值
    
    Returns:
        输出文件
    """
    # 下载输入文件
    dxpy.download_dxfile(input_file["$dnanexus_link"], "input.txt")
    
    # 处理逻辑
    # TODO: 添加您的分析代码
    
    # 示例: 简单处理
    with open("input.txt", "r") as f:
        content = f.read()
    
    result = f"Processed with param={{param_value}}\\n{{content}}"
    
    with open("output.txt", "w") as f:
        f.write(result)
    
    # 上传输出
    output_file = dxpy.upload_local_file("output.txt")
    
    return {{
        "output_file": dxpy.dxlink(output_file)
    }}

dxpy.run()
'''
        script_path = src_dir / f"{app_name}.py"
    else:
        script_content = f'''#!/bin/bash
# DNAnexus App: {app_name}

set -e -o pipefail

# 下载输入
dx download "$input_file" -o input.txt

# 处理逻辑
# TODO: 添加您的分析代码

echo "Processing with param=$param_value"
cat input.txt > output.txt

# 上传输出
dx-upload-all-outputs
'''
        script_path = src_dir / f"{app_name}.sh"
    
    with open(script_path, 'w') as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)
    print(f"  ✓ {script_path}")
    
    # 创建Readme
    readme = f"""# {app_name}

## 描述

添加应用描述

## 输入

"""
    for inp in inputs:
        readme += f"- **{inp['name']}**: {inp.get('label', 'N/A')}\n"
    
    readme += "\n## 输出\n\n"
    for out in outputs:
        readme += f"- **{out['name']}**: {out.get('label', 'N/A')}\n"
    
    readme += f"""
## 使用

```bash
# 构建应用
dx build

# 运行应用
dx run {app_name} -iinput_file=file-xxxx
```

## 开发者

- 您的姓名

"""
    
    readme_path = app_dir / "Readme.md"
    with open(readme_path, 'w') as f:
        f.write(readme)
    print(f"  ✓ {readme_path}")
    
    print(f"\n应用模板创建完成: {app_dir}")
    print("\n下一步:")
    print(f"  cd {app_dir}")
    print(f"  # 编辑 src/{app_name}.py 添加您的分析逻辑")
    print(f"  dx build  # 构建应用")


def validate_dxapp(app_dir: str) -> bool:
    """
    验证dxapp.json配置
    
    Args:
        app_dir: 应用目录
    
    Returns:
        是否有效
    """
    dxapp_path = Path(app_dir) / "dxapp.json"
    
    if not dxapp_path.exists():
        print(f"错误: 未找到 {dxapp_path}")
        return False
    
    print(f"验证配置: {dxapp_path}")
    
    with open(dxapp_path) as f:
        try:
            dxapp = json.load(f)
        except json.JSONDecodeError as e:
            print(f"  ✗ JSON解析错误: {e}")
            return False
    
    # 检查必需字段
    required = ["name", "title", "inputSpec", "outputSpec", "runSpec"]
    missing = [f for f in required if f not in dxapp]
    
    if missing:
        print(f"  ✗ 缺少必需字段: {missing}")
        return False
    
    # 检查runSpec
    run_spec = dxapp.get("runSpec", {})
    if "interpreter" not in run_spec:
        print("  ✗ runSpec缺少interpreter字段")
        return False
    
    if "file" not in run_spec:
        print("  ✗ runSpec缺少file字段")
        return False
    
    # 检查脚本文件
    script_file = Path(app_dir) / run_spec["file"]
    if not script_file.exists():
        print(f"  ! 警告: 脚本文件不存在 {script_file}")
    
    print("  ✓ 配置有效")
    return True


def build_app(app_dir: str, as_app: bool = False) -> bool:
    """
    构建DNAnexus应用
    
    Args:
        app_dir: 应用目录
        as_app: 是否构建为可发布应用（而非小程序）
    
    Returns:
        是否成功
    """
    import subprocess
    
    # 先验证
    if not validate_dxapp(app_dir):
        print("配置验证失败，停止构建")
        return False
    
    print(f"构建应用: {app_dir}")
    
    cmd = ["dx", "build"]
    if as_app:
        cmd.append("--app")
    cmd.append(app_dir)
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print("  ✓ 构建成功")
        print(result.stdout)
        return True
    else:
        print("  ✗ 构建失败")
        print(result.stderr)
        return False


def test_app_locally(
    app_dir: str,
    test_inputs: Optional[Dict] = None
):
    """
    本地测试应用（使用dx-docker）
    
    Args:
        app_dir: 应用目录
        test_inputs: 测试输入
    """
    print(f"本地测试应用: {app_dir}")
    print("  注意: 需要安装dx-docker")
    
    # 创建测试输入文件
    if test_inputs:
        test_dir = Path(app_dir) / "test_inputs"
        test_dir.mkdir(exist_ok=True)
        
        for key, value in test_inputs.items():
            test_file = test_dir / f"{key}.json"
            with open(test_file, 'w') as f:
                json.dump(value, f)
        
        print(f"  测试输入已创建: {test_dir}")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="DNAnexus 应用开发工具"
    )
    subparsers = parser.add_subparsers(dest="command", help="可用命令")
    
    # 创建模板
    create_parser = subparsers.add_parser("create", help="创建应用模板")
    create_parser.add_argument("name", help="应用名称")
    create_parser.add_argument("--output", "-o", default=".", help="输出目录")
    create_parser.add_argument("--language", choices=["python", "bash"], default="python")
    
    # 验证配置
    validate_parser = subparsers.add_parser("validate", help="验证dxapp.json")
    validate_parser.add_argument("dir", help="应用目录")
    
    # 构建应用
    build_parser = subparsers.add_parser("build", help="构建应用")
    build_parser.add_argument("dir", help="应用目录")
    build_parser.add_argument("--app", action="store_true", help="构建为可发布应用")
    
    args = parser.parse_args()
    
    if args.command == "create":
        create_app_template(args.name, args.output, args.language)
    elif args.command == "validate":
        validate_dxapp(args.dir)
    elif args.command == "build":
        build_app(args.dir, args.app)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
