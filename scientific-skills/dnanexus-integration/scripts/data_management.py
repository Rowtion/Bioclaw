#!/usr/bin/env python3
"""
DNAnexus 数据管理脚本
功能：
1. 文件上传/下载
2. 搜索数据对象
3. 批量操作
4. 项目数据管理

依赖：dxpy
"""

import os
import sys
from pathlib import Path
from typing import List, Dict, Optional, Any
import json


def login_to_dnanexus():
    """登录DNAnexus平台"""
    try:
        import dxpy
        # 检查是否已登录
        try:
            dxpy.whoami()
            print("✓ 已登录DNAnexus")
            return True
        except dxpy.exceptions.InvalidAuthentication:
            print("请先登录DNAnexus:")
            print("  dx login")
            return False
    except ImportError:
        print("错误: 未安装dxpy。请运行: uv pip install dxpy")
        return False


def upload_file(
    local_path: str,
    project_id: Optional[str] = None,
    folder: str = "/",
    name: Optional[str] = None,
    properties: Optional[Dict] = None
) -> str:
    """
    上传文件到DNAnexus
    
    Args:
        local_path: 本地文件路径
        project_id: 项目ID (如 project-xxxx)
        folder: 目标文件夹
        name: 远程文件名（默认为本地文件名）
        properties: 文件属性字典
    
    Returns:
        文件ID
    """
    import dxpy
    
    print(f"上传文件: {local_path}")
    
    file_id = dxpy.upload_local_file(
        local_path,
        project=project_id,
        folder=folder,
        name=name,
        properties=properties or {}
    )
    
    print(f"  ✓ 上传完成: {file_id}")
    return file_id


def download_file(
    file_id: str,
    dest_path: str,
    project_id: Optional[str] = None
):
    """
    从DNAnexus下载文件
    
    Args:
        file_id: 文件ID (如 file-xxxx)
        dest_path: 本地保存路径
        project_id: 项目ID
    """
    import dxpy
    
    print(f"下载文件: {file_id}")
    
    dxpy.download_dxfile(
        file_id,
        dest_path,
        project=project_id
    )
    
    print(f"  ✓ 已保存: {dest_path}")


def search_files(
    name_pattern: Optional[str] = None,
    properties: Optional[Dict] = None,
    project_id: Optional[str] = None,
    folder: Optional[str] = None,
    limit: int = 100
) -> List[Dict]:
    """
    搜索DNAnexus中的文件
    
    Args:
        name_pattern: 文件名模式（如 *.fastq.gz）
        properties: 属性筛选
        project_id: 项目ID
        folder: 文件夹路径
        limit: 最大返回数量
    
    Returns:
        文件信息列表
    """
    import dxpy
    
    print(f"搜索文件...")
    
    query = {
        "classname": "file",
        "limit": limit
    }
    
    if name_pattern:
        query["name"] = name_pattern
    if project_id:
        query["project"] = project_id
    if folder:
        query["folder"] = folder
    if properties:
        query["properties"] = properties
    
    results = list(dxpy.find_data_objects(**query))
    
    print(f"  找到 {len(results)} 个文件")
    return results


def batch_upload(
    local_dir: str,
    project_id: str,
    remote_folder: str = "/data",
    pattern: str = "*"
) -> List[str]:
    """
    批量上传文件
    
    Args:
        local_dir: 本地目录
        project_id: 项目ID
        remote_folder: 远程目标文件夹
        pattern: 文件匹配模式
    
    Returns:
        上传的文件ID列表
    """
    import dxpy
    
    local_path = Path(local_dir)
    files = list(local_path.glob(pattern))
    
    print(f"批量上传: {len(files)} 个文件")
    print(f"  从: {local_path}")
    print(f"  到: {project_id}:{remote_folder}")
    
    uploaded_ids = []
    for file_path in files:
        if file_path.is_file():
            try:
                file_id = dxpy.upload_local_file(
                    str(file_path),
                    project=project_id,
                    folder=remote_folder,
                    name=file_path.name
                )
                uploaded_ids.append(file_id)
                print(f"  ✓ {file_path.name}")
            except Exception as e:
                print(f"  ✗ {file_path.name}: {e}")
    
    print(f"\n上传完成: {len(uploaded_ids)}/{len(files)}")
    return uploaded_ids


def batch_download(
    project_id: str,
    remote_folder: str,
    local_dir: str,
    pattern: str = "*"
):
    """
    批量下载文件
    
    Args:
        project_id: 项目ID
        remote_folder: 远程文件夹
        local_dir: 本地保存目录
        pattern: 文件名匹配模式
    """
    import dxpy
    
    # 确保本地目录存在
    Path(local_dir).mkdir(parents=True, exist_ok=True)
    
    # 搜索文件
    files = search_files(
        name_pattern=pattern,
        project_id=project_id,
        folder=remote_folder
    )
    
    print(f"\n批量下载: {len(files)} 个文件")
    
    downloaded = 0
    for file_info in files:
        file_id = file_info["id"]
        file_obj = dxpy.DXFile(file_id)
        file_desc = file_obj.describe()
        filename = file_desc["name"]
        
        dest_path = Path(local_dir) / filename
        
        try:
            dxpy.download_dxfile(file_id, str(dest_path))
            print(f"  ✓ {filename}")
            downloaded += 1
        except Exception as e:
            print(f"  ✗ {filename}: {e}")
    
    print(f"\n下载完成: {downloaded}/{len(files)}")


def create_folder(project_id: str, folder_path: str):
    """
    创建文件夹
    
    Args:
        project_id: 项目ID
        folder_path: 文件夹路径
    """
    import dxpy
    
    dxpy.api.project_new_folder(
        project_id,
        {
            "folder": folder_path,
            "parents": True
        }
    )
    print(f"✓ 创建文件夹: {folder_path}")


def list_project_contents(
    project_id: str,
    folder: str = "/",
    recurse: bool = False
) -> List[Dict]:
    """
    列出项目内容
    
    Args:
        project_id: 项目ID
        folder: 文件夹路径
        recurse: 是否递归列出
    
    Returns:
        对象信息列表
    """
    import dxpy
    
    print(f"列出项目内容: {project_id}:{folder}")
    
    results = list(dxpy.find_data_objects(
        project=project_id,
        folder=folder,
        recurse=recurse
    ))
    
    print(f"  找到 {len(results)} 个对象")
    
    for obj in results[:20]:  # 显示前20个
        obj_id = obj["id"]
        obj_class = obj["describe"]["class"]
        obj_name = obj["describe"].get("name", "N/A")
        print(f"  {obj_class}: {obj_name} ({obj_id})")
    
    if len(results) > 20:
        print(f"  ... 还有 {len(results) - 20} 个对象")
    
    return results


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="DNAnexus 数据管理工具"
    )
    subparsers = parser.add_subparsers(dest="command", help="可用命令")
    
    # 登录命令
    subparsers.add_parser("login", help="检查登录状态")
    
    # 上传命令
    upload_parser = subparsers.add_parser("upload", help="上传文件")
    upload_parser.add_argument("file", help="本地文件路径")
    upload_parser.add_argument("--project", required=True, help="项目ID")
    upload_parser.add_argument("--folder", default="/", help="目标文件夹")
    upload_parser.add_argument("--name", help="远程文件名")
    
    # 下载命令
    download_parser = subparsers.add_parser("download", help="下载文件")
    download_parser.add_argument("file_id", help="文件ID")
    download_parser.add_argument("--dest", required=True, help="本地保存路径")
    download_parser.add_argument("--project", help="项目ID")
    
    # 搜索命令
    search_parser = subparsers.add_parser("search", help="搜索文件")
    search_parser.add_argument("--pattern", help="文件名模式")
    search_parser.add_argument("--project", help="项目ID")
    search_parser.add_argument("--folder", help="文件夹")
    search_parser.add_argument("--limit", type=int, default=100)
    
    # 批量上传
    batch_up_parser = subparsers.add_parser("batch-upload", help="批量上传")
    batch_up_parser.add_argument("--local-dir", required=True, help="本地目录")
    batch_up_parser.add_argument("--project", required=True, help="项目ID")
    batch_up_parser.add_argument("--remote-folder", default="/data")
    batch_up_parser.add_argument("--pattern", default="*")
    
    # 批量下载
    batch_down_parser = subparsers.add_parser("batch-download", help="批量下载")
    batch_down_parser.add_argument("--project", required=True, help="项目ID")
    batch_down_parser.add_argument("--remote-folder", required=True)
    batch_down_parser.add_argument("--local-dir", required=True)
    batch_down_parser.add_argument("--pattern", default="*")
    
    # 列出内容
    list_parser = subparsers.add_parser("list", help="列出项目内容")
    list_parser.add_argument("--project", required=True, help="项目ID")
    list_parser.add_argument("--folder", default="/")
    list_parser.add_argument("--recurse", action="store_true")
    
    args = parser.parse_args()
    
    if args.command == "login":
        login_to_dnanexus()
    elif args.command == "upload":
        if login_to_dnanexus():
            upload_file(
                args.file,
                project_id=args.project,
                folder=args.folder,
                name=args.name
            )
    elif args.command == "download":
        if login_to_dnanexus():
            download_file(
                args.file_id,
                args.dest,
                project_id=args.project
            )
    elif args.command == "search":
        if login_to_dnanexus():
            results = search_files(
                name_pattern=args.pattern,
                project_id=args.project,
                folder=args.folder,
                limit=args.limit
            )
            print("\n搜索结果:")
            for r in results:
                desc = r.get("describe", {})
                print(f"  {desc.get('name', 'N/A')} ({r['id']})")
    elif args.command == "batch-upload":
        if login_to_dnanexus():
            batch_upload(
                args.local_dir,
                args.project,
                args.remote_folder,
                args.pattern
            )
    elif args.command == "batch-download":
        if login_to_dnanexus():
            batch_download(
                args.project,
                args.remote_folder,
                args.local_dir,
                args.pattern
            )
    elif args.command == "list":
        if login_to_dnanexus():
            list_project_contents(
                args.project,
                args.folder,
                args.recurse
            )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
