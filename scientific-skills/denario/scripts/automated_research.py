#!/usr/bin/env python3
"""
Denario 自动化科研流程脚本
功能：
1. 从数据描述生成研究假设
2. 开发研究方法论
3. 执行计算实验
4. 生成学术论文

依赖：denario, pandas, matplotlib, scikit-learn
"""

import os
from pathlib import Path
from typing import Optional, Dict, Any


def setup_research_project(
    project_name: str,
    project_dir: str = "./denario_projects"
) -> str:
    """
    设置研究项目目录
    
    Args:
        project_name: 项目名称
        project_dir: 项目根目录
    
    Returns:
        项目完整路径
    """
    from denario import Denario
    
    project_path = Path(project_dir) / project_name
    project_path.mkdir(parents=True, exist_ok=True)
    
    print(f"项目目录: {project_path}")
    return str(project_path)


def run_full_research_pipeline(
    project_dir: str,
    data_description: str,
    journal: str = "APS"
):
    """
    运行完整的研究流程
    
    Args:
        project_dir: 项目目录
        data_description: 数据描述文本
        journal: 目标期刊 (APS, Nature等)
    """
    try:
        from denario import Denario, Journal
        
        print(f"=" * 60)
        print(f"启动 Denario 自动化研究流程")
        print(f"=" * 60)
        
        # 初始化项目
        print("\n[1/5] 初始化研究项目...")
        den = Denario(project_dir=project_dir)
        den.set_data_description(data_description)
        print(f"  ✓ 项目已初始化: {project_dir}")
        
        # 生成研究假设
        print("\n[2/5] 生成研究假设...")
        try:
            den.get_idea()
            print(f"  ✓ 研究假设已生成")
        except Exception as e:
            print(f"  ! 生成假设时出错: {e}")
            print("  ! 请检查数据描述和API配置")
            return
        
        # 开发方法论
        print("\n[3/5] 开发研究方法论...")
        try:
            den.get_method()
            print(f"  ✓ 方法论已开发")
        except Exception as e:
            print(f"  ! 开发方法论时出错: {e}")
            return
        
        # 执行实验
        print("\n[4/5] 执行计算实验...")
        try:
            den.get_results()
            print(f"  ✓ 实验结果已生成")
        except Exception as e:
            print(f"  ! 执行实验时出错: {e}")
            return
        
        # 生成论文
        print("\n[5/5] 生成学术论文...")
        try:
            journal_enum = getattr(Journal, journal.upper(), Journal.APS)
            den.get_paper(journal=journal_enum)
            print(f"  ✓ 论文已生成 (格式: {journal})")
        except Exception as e:
            print(f"  ! 生成论文时出错: {e}")
            return
        
        print("\n" + "=" * 60)
        print("研究流程完成!")
        print("=" * 60)
        print(f"\n项目文件位置: {project_dir}")
        print("请查看生成的论文和分析结果。")
        
    except ImportError:
        print("错误: 未安装denario。请运行: uv pip install 'denario[app]'")
        raise


def run_custom_research_pipeline(
    project_dir: str,
    data_description: str,
    custom_idea: Optional[str] = None,
    custom_method: Optional[str] = None,
    custom_results: Optional[str] = None,
    journal: str = "APS"
):
    """
    运行自定义研究流程（可指定自定义输入）
    
    Args:
        project_dir: 项目目录
        data_description: 数据描述
        custom_idea: 自定义研究假设文件路径（可选）
        custom_method: 自定义方法论文件路径（可选）
        custom_results: 自定义结果文件路径（可选）
        journal: 目标期刊
    """
    try:
        from denario import Denario, Journal
        
        print(f"=" * 60)
        print(f"启动 Denario 自定义研究流程")
        print(f"=" * 60)
        
        den = Denario(project_dir=project_dir)
        den.set_data_description(data_description)
        
        # 研究假设
        if custom_idea:
            print(f"\n[假设] 使用自定义假设: {custom_idea}")
            den.set_idea(custom_idea)
        else:
            print("\n[假设] 自动生成研究假设...")
            den.get_idea()
        
        # 方法论
        if custom_method:
            print(f"\n[方法] 使用自定义方法论: {custom_method}")
            den.set_method(custom_method)
        else:
            print("\n[方法] 自动生成方法论...")
            den.get_method()
        
        # 结果
        if custom_results:
            print(f"\n[结果] 使用自定义结果: {custom_results}")
            den.set_results(custom_results)
        else:
            print("\n[结果] 自动执行实验...")
            den.get_results()
        
        # 生成论文
        print("\n[论文] 生成学术论文...")
        journal_enum = getattr(Journal, journal.upper(), Journal.APS)
        den.get_paper(journal=journal_enum)
        
        print("\n" + "=" * 60)
        print("自定义研究流程完成!")
        print("=" * 60)
        
    except ImportError:
        print("错误: 未安装denario")
        raise


def main():
    """主函数 - 命令行接口"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Denario 自动化科研流程"
    )
    parser.add_argument(
        "--project",
        required=True,
        help="项目名称"
    )
    parser.add_argument(
        "--data-desc",
        required=True,
        help="数据描述文本或文件路径"
    )
    parser.add_argument(
        "--journal",
        default="APS",
        choices=["APS", "NATURE", "SCIENCE"],
        help="目标期刊格式"
    )
    parser.add_argument(
        "--custom-idea",
        help="自定义研究假设文件(.md)"
    )
    parser.add_argument(
        "--custom-method",
        help="自定义方法论文件(.md)"
    )
    parser.add_argument(
        "--custom-results",
        help="自定义结果文件(.md)"
    )
    parser.add_argument(
        "--output-dir",
        default="./denario_projects",
        help="输出目录"
    )
    
    args = parser.parse_args()
    
    # 读取数据描述
    if Path(args.data_desc).exists():
        with open(args.data_desc, 'r') as f:
            data_description = f.read()
    else:
        data_description = args.data_desc
    
    # 设置项目
    project_dir = setup_research_project(args.project, args.output_dir)
    
    # 运行流程
    if args.custom_idea or args.custom_method or args.custom_results:
        run_custom_research_pipeline(
            project_dir=project_dir,
            data_description=data_description,
            custom_idea=args.custom_idea,
            custom_method=args.custom_method,
            custom_results=args.custom_results,
            journal=args.journal
        )
    else:
        run_full_research_pipeline(
            project_dir=project_dir,
            data_description=data_description,
            journal=args.journal
        )


if __name__ == "__main__":
    main()
