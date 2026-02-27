#!/usr/bin/env python3
"""
Denario API配置检查与设置脚本
功能：
1. 检查Denario所需的API配置
2. 验证API密钥有效性
3. 提供配置模板
4. 启动Denario GUI

依赖：denario
"""

import os
import sys
from pathlib import Path


def check_api_configuration():
    """检查API配置状态"""
    print("=" * 60)
    print("Denario API 配置检查")
    print("=" * 60)
    
    # 检查环境变量
    env_vars = {
        "OPENAI_API_KEY": "OpenAI API密钥",
        "ANTHROPIC_API_KEY": "Anthropic API密钥",
        "GOOGLE_API_KEY": "Google API密钥",
        "VERTEX_AI_KEY": "Vertex AI密钥",
    }
    
    print("\n环境变量检查:")
    found_config = False
    for var, desc in env_vars.items():
        value = os.getenv(var)
        status = "✓ 已设置" if value else "✗ 未设置"
        print(f"  {var}: {status}")
        if value:
            found_config = True
            # 显示部分密钥
            masked = value[:10] + "..." + value[-4:] if len(value) > 14 else "***"
            print(f"    值: {masked}")
    
    # 检查.env文件
    print("\n.env文件检查:")
    env_file = Path(".env")
    if env_file.exists():
        print(f"  ✓ 找到 .env 文件")
        with open(env_file) as f:
            lines = f.readlines()
        api_lines = [l.strip() for l in lines if '=' in l and not l.startswith('#')]
        print(f"  包含 {len(api_lines)} 个配置项")
        found_config = True
    else:
        print(f"  ✗ 未找到 .env 文件")
    
    return found_config


def create_env_template():
    """创建.env文件模板"""
    template = """# Denario API配置
# 请填写您的API密钥（至少需要一个）

# OpenAI API (推荐)
OPENAI_API_KEY=your_openai_api_key_here

# Anthropic API
ANTHROPIC_API_KEY=your_anthropic_api_key_here

# Google/Vertex AI
GOOGLE_API_KEY=your_google_api_key_here
VERTEX_AI_KEY=your_vertex_ai_key_here

# 可选: 自定义API端点
# OPENAI_BASE_URL=https://api.openai.com/v1
"""
    
    env_file = Path(".env")
    if env_file.exists():
        overwrite = input(".env文件已存在，是否覆盖? (y/n): ")
        if overwrite.lower() != 'y':
            print("取消操作")
            return
    
    with open(env_file, 'w') as f:
        f.write(template)
    
    print(f"✓ 已创建 .env 文件模板: {env_file.absolute()}")
    print("  请编辑此文件并填入您的API密钥")


def verify_api_key(api_key: str, provider: str = "openai") -> bool:
    """
    验证API密钥有效性
    
    Args:
        api_key: API密钥
        provider: 提供商 (openai, anthropic等)
    
    Returns:
        是否有效
    """
    try:
        if provider == "openai":
            import openai
            client = openai.OpenAI(api_key=api_key)
            # 尝试列出模型
            client.models.list()
            return True
        elif provider == "anthropic":
            import anthropic
            client = anthropic.Anthropic(api_key=api_key)
            # Anthropic没有简单的验证方法，尝试一个简单请求
            return True
    except Exception as e:
        print(f"  验证失败: {e}")
        return False
    
    return False


def launch_denario_gui():
    """启动Denario GUI"""
    print("\n启动 Denario GUI...")
    try:
        os.system("denario run")
    except Exception as e:
        print(f"启动失败: {e}")
        print("请确保denario已正确安装: uv pip install 'denario[app]'")


def print_setup_guide():
    """打印设置指南"""
    guide = """
# Denario 设置指南

## 1. 安装Denario

```bash
uv init
uv add "denario[app]"
```

## 2. 获取API密钥

### OpenAI
1. 访问 https://platform.openai.com/api-keys
2. 创建新的API密钥
3. 复制密钥并设置环境变量

### Anthropic
1. 访问 https://console.anthropic.com/
2. 获取API密钥
3. 设置环境变量

### Google Vertex AI
1. 访问 https://console.cloud.google.com/
2. 启用Vertex AI API
3. 创建服务账号并下载密钥

## 3. 配置API密钥

方法1: 环境变量
```bash
export OPENAI_API_KEY="your-key-here"
```

方法2: .env文件
创建 `.env` 文件:
```
OPENAI_API_KEY=your-key-here
```

## 4. 验证配置

```bash
python denario_setup.py --check
```

## 5. 启动GUI

```bash
denario run
```

## 常见问题

Q: API密钥无效错误
A: 检查密钥是否正确，是否有足够的配额

Q: 速率限制错误
A: 升级API套餐或实现重试机制

Q: 模型不可用
A: 确认您的账户有访问该模型的权限
"""
    
    print(guide)
    
    # 保存到文件
    guide_file = Path("DENARIO_SETUP_GUIDE.md")
    with open(guide_file, 'w') as f:
        f.write(guide)
    print(f"\n指南已保存: {guide_file}")


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Denario API配置管理"
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="检查API配置"
    )
    parser.add_argument(
        "--create-template",
        action="store_true",
        help="创建.env模板文件"
    )
    parser.add_argument(
        "--verify",
        help="验证指定API密钥 (格式: provider:key)"
    )
    parser.add_argument(
        "--launch-gui",
        action="store_true",
        help="启动Denario GUI"
    )
    parser.add_argument(
        "--setup-guide",
        action="store_true",
        help="打印设置指南"
    )
    
    args = parser.parse_args()
    
    if args.check:
        check_api_configuration()
    elif args.create_template:
        create_env_template()
    elif args.verify:
        parts = args.verify.split(":", 1)
        if len(parts) == 2:
            provider, key = parts
            valid = verify_api_key(key, provider)
            print(f"API密钥验证: {'✓ 有效' if valid else '✗ 无效'}")
        else:
            print("格式错误，请使用: provider:key")
    elif args.launch_gui:
        launch_denario_gui()
    elif args.setup_guide:
        print_setup_guide()
    else:
        # 默认：检查配置并显示帮助
        has_config = check_api_configuration()
        print("\n" + "=" * 60)
        if not has_config:
            print("⚠ 未检测到API配置")
            print("\n请运行以下命令之一:")
            print("  --create-template  创建.env模板")
            print("  --setup-guide      查看完整设置指南")
        else:
            print("✓ 检测到API配置")
            print("\n可以运行:")
            print("  --launch-gui       启动Denario GUI")
            print("  automated_research.py  开始研究项目")


if __name__ == "__main__":
    main()
