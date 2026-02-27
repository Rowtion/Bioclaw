#!/usr/bin/env python3
"""
DNAnexus 工作流执行脚本
功能：
1. 运行分析应用/小程序
2. 监控作业状态
3. 创建简单工作流链
4. 批量处理多个样本

依赖：dxpy
"""

import time
from typing import List, Dict, Optional, Any
from pathlib import Path


def run_applet(
    applet_id: str,
    inputs: Dict[str, Any],
    project_id: Optional[str] = None,
    instance_type: Optional[str] = None,
    name: Optional[str] = None
) -> str:
    """
    运行DNAnexus小程序
    
    Args:
        applet_id: 小程序ID (如 applet-xxxx)
        inputs: 输入参数字典
        project_id: 项目ID
        instance_type: 实例类型 (如 mem1_ssd1_v2_x8)
        name: 作业名称
    
    Returns:
        作业ID
    """
    import dxpy
    
    print(f"运行小程序: {applet_id}")
    print(f"  输入: {inputs}")
    
    applet = dxpy.DXApplet(applet_id, project=project_id)
    
    job = applet.run(
        inputs,
        project=project_id,
        instance_type=instance_type,
        name=name
    )
    
    print(f"  ✓ 作业已提交: {job.get_id()}")
    return job.get_id()


def run_app(
    app_name: str,
    inputs: Dict[str, Any],
    project_id: Optional[str] = None,
    name: Optional[str] = None
) -> str:
    """
    运行DNAnexus应用
    
    Args:
        app_name: 应用名称 (如 app-cloud_workstation)
        inputs: 输入参数字典
        project_id: 项目ID
        name: 作业名称
    
    Returns:
        作业ID
    """
    import dxpy
    
    print(f"运行应用: {app_name}")
    
    job = dxpy.run(
        app_name,
        inputs,
        project=project_id,
        name=name
    )
    
    print(f"  ✓ 作业已提交: {job.get_id()}")
    return job.get_id()


def wait_for_job(
    job_id: str,
    poll_interval: int = 30,
    timeout: Optional[int] = None
) -> str:
    """
    等待作业完成
    
    Args:
        job_id: 作业ID
        poll_interval: 轮询间隔（秒）
        timeout: 超时时间（秒）
    
    Returns:
        最终状态
    """
    import dxpy
    
    print(f"等待作业完成: {job_id}")
    
    job = dxpy.DXJob(job_id)
    start_time = time.time()
    
    while True:
        desc = job.describe()
        state = desc["state"]
        
        if state in ["done", "failed", "terminated"]:
            print(f"\n  作业完成，状态: {state}")
            return state
        
        print(f"  状态: {state} (等待 {poll_interval}s...)")
        time.sleep(poll_interval)
        
        if timeout and (time.time() - start_time) > timeout:
            print(f"\n  等待超时")
            return "timeout"


def get_job_outputs(job_id: str) -> Dict:
    """
    获取作业输出
    
    Args:
        job_id: 作业ID
    
    Returns:
        输出字典
    """
    import dxpy
    
    job = dxpy.DXJob(job_id)
    desc = job.describe()
    
    outputs = desc.get("output", {})
    print(f"作业输出 ({job_id}):")
    for key, value in outputs.items():
        print(f"  {key}: {value}")
    
    return outputs


def run_workflow_chain(
    steps: List[Dict],
    project_id: str
) -> List[str]:
    """
    运行工作流链（多步骤分析）
    
    Args:
        steps: 步骤配置列表，每项包含:
            - applet_id: 小程序ID
            - inputs: 输入（可使用前一步的输出引用）
            - name: 步骤名称
        project_id: 项目ID
    
    Returns:
        作业ID列表
    """
    import dxpy
    
    print(f"运行工作流链: {len(steps)} 个步骤")
    
    job_ids = []
    previous_job = None
    
    for i, step in enumerate(steps, 1):
        print(f"\n步骤 {i}/{len(steps)}: {step.get('name', 'unnamed')}")
        
        # 解析输入（处理对前一步输出的引用）
        inputs = step["inputs"].copy()
        
        if previous_job:
            # 替换引用，如 {"$job_output": "output_name"}
            for key, value in inputs.items():
                if isinstance(value, dict) and "$job_output" in value:
                    output_name = value["$job_output"]
                    inputs[key] = previous_job.get_output_ref(output_name)
        
        # 运行
        job_id = run_applet(
            step["applet_id"],
            inputs,
            project_id=project_id,
            name=step.get("name")
        )
        
        job_ids.append(job_id)
        previous_job = dxpy.DXJob(job_id)
    
    print(f"\n工作流链已提交: {len(job_ids)} 个作业")
    return job_ids


def batch_process(
    applet_id: str,
    input_list: List[Dict],
    project_id: str,
    max_parallel: int = 5
) -> List[str]:
    """
    批量处理多个样本
    
    Args:
        applet_id: 小程序ID
        input_list: 输入配置列表（每项是一个输入字典）
        project_id: 项目ID
        max_parallel: 最大并行作业数
    
    Returns:
        作业ID列表
    """
    import dxpy
    
    print(f"批量处理: {len(input_list)} 个样本")
    print(f"  小程序: {applet_id}")
    print(f"  最大并行: {max_parallel}")
    
    jobs = []
    completed = []
    
    for i, inputs in enumerate(input_list):
        print(f"\n提交作业 {i+1}/{len(input_list)}")
        
        # 等待直到有可用槽位
        while len([j for j in jobs if j.describe()["state"] in ["idle", "running", "runnable"]]) >= max_parallel:
            time.sleep(10)
            # 清理已完成的作业
            jobs = [j for j in jobs if j.describe()["state"] not in ["done", "failed", "terminated"]]
        
        job_id = run_applet(
            applet_id,
            inputs,
            project_id=project_id,
            name=f"batch_job_{i+1}"
        )
        jobs.append(dxpy.DXJob(job_id))
        completed.append(job_id)
    
    # 等待所有作业完成
    print(f"\n等待所有 {len(completed)} 个作业完成...")
    for job_id in completed:
        wait_for_job(job_id)
    
    return completed


def list_jobs(
    project_id: str,
    state: Optional[str] = None,
    limit: int = 50
) -> List[Dict]:
    """
    列出项目中的作业
    
    Args:
        project_id: 项目ID
        state: 状态筛选 (done, failed, running, idle)
        limit: 最大返回数量
    
    Returns:
        作业信息列表
    """
    import dxpy
    
    print(f"列出作业: {project_id}")
    
    query = {"project": project_id, "limit": limit}
    if state:
        query["state"] = state
    
    jobs = list(dxpy.find_jobs(**query))
    
    print(f"  找到 {len(jobs)} 个作业")
    
    for job in jobs[:20]:
        desc = job.describe()
        print(f"  {desc.get('name', 'N/A')} - {desc['state']} ({job.get_id()})")
    
    return jobs


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="DNAnexus 工作流执行工具"
    )
    subparsers = parser.add_subparsers(dest="command", help="可用命令")
    
    # 运行小程序
    run_parser = subparsers.add_parser("run", help="运行小程序")
    run_parser.add_argument("applet_id", help="小程序ID")
    run_parser.add_argument("--input", "-i", action="append", help="输入参数 (key=value)")
    run_parser.add_argument("--project", required=True, help="项目ID")
    run_parser.add_argument("--name", help="作业名称")
    run_parser.add_argument("--wait", action="store_true", help="等待完成")
    
    # 等待作业
    wait_parser = subparsers.add_parser("wait", help="等待作业完成")
    wait_parser.add_argument("job_id", help="作业ID")
    wait_parser.add_argument("--interval", type=int, default=30, help="轮询间隔")
    
    # 获取输出
    output_parser = subparsers.add_parser("outputs", help="获取作业输出")
    output_parser.add_argument("job_id", help="作业ID")
    
    # 列出作业
    list_parser = subparsers.add_parser("list", help="列出作业")
    list_parser.add_argument("--project", required=True, help="项目ID")
    list_parser.add_argument("--state", choices=["done", "failed", "running", "idle"])
    list_parser.add_argument("--limit", type=int, default=50)
    
    # 批量处理
    batch_parser = subparsers.add_parser("batch", help="批量处理")
    batch_parser.add_argument("applet_id", help="小程序ID")
    batch_parser.add_argument("--project", required=True, help="项目ID")
    batch_parser.add_argument("--input-file", required=True, help="输入配置JSON文件")
    batch_parser.add_argument("--parallel", type=int, default=5, help="最大并行数")
    
    args = parser.parse_args()
    
    if args.command == "run":
        # 解析输入
        inputs = {}
        if args.input:
            for inp in args.input:
                key, value = inp.split("=", 1)
                # 尝试解析JSON
                try:
                    inputs[key] = json.loads(value)
                except:
                    inputs[key] = value
        
        job_id = run_applet(
            args.applet_id,
            inputs,
            project_id=args.project,
            name=args.name
        )
        
        if args.wait:
            wait_for_job(job_id)
    
    elif args.command == "wait":
        wait_for_job(args.job_id, poll_interval=args.interval)
    
    elif args.command == "outputs":
        get_job_outputs(args.job_id)
    
    elif args.command == "list":
        list_jobs(args.project, args.state, args.limit)
    
    elif args.command == "batch":
        # 读取输入配置
        with open(args.input_file) as f:
            input_list = json.load(f)
        
        batch_process(args.applet_id, input_list, args.project, args.parallel)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    import json
    main()
