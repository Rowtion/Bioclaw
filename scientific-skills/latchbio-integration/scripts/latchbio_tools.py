#!/usr/bin/env python3
"""
LatchBio Platform Tools
Cloud bioinformatics platform for workflow execution and data management.
"""

import subprocess
import os
from typing import List, Dict


def login():
    """Authenticate with LatchBio platform."""
    subprocess.run(["latch", "login"], check=True)
    print("Logged in to LatchBio")


def init_workflow(name: str):
    """Initialize a new workflow."""
    subprocess.run(["latch", "init", name], check=True)
    print(f"Initialized workflow: {name}")


def register_workflow(workflow_dir: str):
    """Register workflow to LatchBio."""
    subprocess.run(["latch", "register", workflow_dir], check=True)
    print("Workflow registered")


def execute_workflow(workflow_id: str, params: Dict = None):
    """Execute a registered workflow."""
    cmd = ["latch", "execute", workflow_id]
    if params:
        for k, v in params.items():
            cmd.extend(["--param", f"{k}={v}"])
    subprocess.run(cmd, check=True)
    print(f"Executed: {workflow_id}")


def list_workflows() -> List[Dict]:
    """List registered workflows."""
    result = subprocess.run(["latch", "workflow", "list"], 
                          capture_output=True, text=True)
    # Parse output
    workflows = []
    for line in result.stdout.split('\n')[2:]:  # Skip headers
        if line.strip():
            parts = line.split()
            if len(parts) >= 2:
                workflows.append({"id": parts[0], "name": parts[1]})
    return workflows


def upload_data(local_path: str, remote_path: str = ""):
    """Upload data to LatchBio."""
    subprocess.run(["latch", "cp", local_path, f"latch://{remote_path}"], check=True)
    print(f"Uploaded: {local_path}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="LatchBio Tools")
    parser.add_argument("command", 
                       choices=["login", "init", "register", "execute", "list", "upload"])
    parser.add_argument("--name", help="Workflow/data name")
    parser.add_argument("--dir", help="Workflow directory")
    parser.add_argument("--workflow-id", help="Workflow ID")
    parser.add_argument("--local-path", help="Local file path")
    parser.add_argument("--remote-path", help="Remote destination")
    
    args = parser.parse_args()
    
    if args.command == "login":
        login()
    elif args.command == "init":
        init_workflow(args.name)
    elif args.command == "register":
        register_workflow(args.dir)
    elif args.command == "execute":
        execute_workflow(args.workflow_id)
    elif args.command == "list":
        workflows = list_workflows()
        for w in workflows:
            print(f"{w['id']}: {w['name']}")
    elif args.command == "upload":
        upload_data(args.local_path, args.remote_path)


if __name__ == "__main__":
    main()
