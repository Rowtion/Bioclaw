#!/usr/bin/env python3
"""
Protocols.io API Tools - Protocol management and integration utilities
Search, create, and manage scientific protocols via protocols.io API
"""

import os
import argparse
import json
import requests
from typing import List, Dict, Optional
from pathlib import Path

BASE_URL = "https://protocols.io/api/v3"


def get_headers(token: Optional[str] = None) -> Dict:
    """Get API headers with authentication."""
    token = token or os.getenv("PROTOCOLSIO_TOKEN")
    if not token:
        raise ValueError("API token required. Set PROTOCOLSIO_TOKEN environment variable.")
    return {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json"
    }


def search_protocols(query: str, token: Optional[str] = None, 
                    filter_type: str = "public", page_size: int = 10) -> List[Dict]:
    """Search for protocols by keyword."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/protocols"
    params = {
        "filter": filter_type,
        "key": query,
        "page_size": page_size,
        "content_format": "html"
    }
    
    try:
        response = requests.get(url, headers=headers, params=params)
        response.raise_for_status()
        data = response.json()
        
        protocols = []
        if "items" in data:
            for item in data["items"]:
                protocols.append({
                    "id": item.get("id"),
                    "title": item.get("title"),
                    "doi": item.get("doi"),
                    "authors": [a.get("name") for a in item.get("authors", [])],
                    "created_at": item.get("created_on"),
                    "url": item.get("url", {}).get("public", "")
                })
        return protocols
    except Exception as e:
        print(f"Error searching protocols: {e}")
        return []


def get_protocol(protocol_id: str, token: Optional[str] = None) -> Dict:
    """Get detailed protocol information."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/protocols/{protocol_id}"
    params = {"content_format": "html"}
    
    try:
        response = requests.get(url, headers=headers, params=params)
        response.raise_for_status()
        data = response.json()
        
        protocol = data.get("protocol", {})
        return {
            "id": protocol.get("id"),
            "title": protocol.get("title"),
            "description": protocol.get("description"),
            "doi": protocol.get("doi"),
            "authors": [a.get("name") for a in protocol.get("authors", [])],
            "steps": len(protocol.get("steps", [])),
            "materials": protocol.get("materials", []),
            "created_at": protocol.get("created_on"),
            "published_at": protocol.get("published_on"),
            "url": protocol.get("url", {}).get("public", "")
        }
    except Exception as e:
        print(f"Error fetching protocol: {e}")
        return {}


def get_protocol_steps(protocol_id: str, token: Optional[str] = None) -> List[Dict]:
    """Get all steps of a protocol."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/protocols/{protocol_id}"
    params = {"content_format": "html"}
    
    try:
        response = requests.get(url, headers=headers, params=params)
        response.raise_for_status()
        data = response.json()
        
        protocol = data.get("protocol", {})
        steps = []
        
        for step in protocol.get("steps", []):
            steps.append({
                "id": step.get("id"),
                "position": step.get("position"),
                "title": step.get("title", ""),
                "description": step.get("description", {}).get("body", ""),
                "components": step.get("components", [])
            })
        
        return steps
    except Exception as e:
        print(f"Error fetching steps: {e}")
        return []


def create_protocol(title: str, description: str, token: Optional[str] = None) -> str:
    """Create a new protocol."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/protocols"
    data = {
        "title": title,
        "description": description,
        "tags": []
    }
    
    try:
        response = requests.post(url, headers=headers, json=data)
        response.raise_for_status()
        result = response.json()
        
        protocol_id = result.get("item", {}).get("id")
        print(f"Created protocol: {protocol_id}")
        return protocol_id
    except Exception as e:
        print(f"Error creating protocol: {e}")
        return ""


def add_step(protocol_id: str, title: str, description: str,
             position: Optional[int] = None, token: Optional[str] = None) -> bool:
    """Add a step to a protocol."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/protocols/{protocol_id}/steps"
    data = {
        "title": title,
        "description": description
    }
    if position is not None:
        data["position"] = position
    
    try:
        response = requests.post(url, headers=headers, json=data)
        response.raise_for_status()
        print(f"Added step to protocol {protocol_id}")
        return True
    except Exception as e:
        print(f"Error adding step: {e}")
        return False


def list_workspaces(token: Optional[str] = None) -> List[Dict]:
    """List user workspaces."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/workspaces"
    
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        data = response.json()
        
        workspaces = []
        if "items" in data:
            for item in data["items"]:
                workspaces.append({
                    "id": item.get("id"),
                    "name": item.get("name"),
                    "description": item.get("description"),
                    "member_count": item.get("members_count"),
                    "protocol_count": item.get("protocols_count")
                })
        return workspaces
    except Exception as e:
        print(f"Error listing workspaces: {e}")
        return []


def get_workspace_protocols(workspace_id: str, token: Optional[str] = None) -> List[Dict]:
    """Get protocols in a workspace."""
    headers = get_headers(token)
    
    url = f"{BASE_URL}/workspaces/{workspace_id}/protocols"
    
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        data = response.json()
        
        protocols = []
        if "items" in data:
            for item in data["items"]:
                protocols.append({
                    "id": item.get("id"),
                    "title": item.get("title"),
                    "doi": item.get("doi"),
                    "created_at": item.get("created_on")
                })
        return protocols
    except Exception as e:
        print(f"Error fetching workspace protocols: {e}")
        return []


def export_protocol(protocol_id: str, output_path: str, 
                   format_type: str = "json", token: Optional[str] = None):
    """Export protocol to file."""
    protocol = get_protocol(protocol_id, token)
    steps = get_protocol_steps(protocol_id, token)
    
    export_data = {
        "protocol": protocol,
        "steps": steps
    }
    
    if format_type == "json":
        with open(output_path, 'w') as f:
            json.dump(export_data, f, indent=2)
    elif format_type == "markdown":
        # Convert to markdown
        md_lines = [
            f"# {protocol.get('title', 'Untitled')}",
            "",
            f"**DOI:** {protocol.get('doi', 'N/A')}",
            f"**Authors:** {', '.join(protocol.get('authors', []))}",
            "",
            "## Description",
            protocol.get('description', ''),
            "",
            "## Steps",
            ""
        ]
        
        for step in steps:
            md_lines.append(f"### Step {step['position']}: {step['title']}")
            md_lines.append("")
            md_lines.append(step['description'])
            md_lines.append("")
        
        with open(output_path, 'w') as f:
            f.write("\n".join(md_lines))
    
    print(f"Protocol exported to {output_path}")


def clone_protocol(protocol_id: str, new_title: Optional[str] = None,
                  token: Optional[str] = None) -> str:
    """Clone an existing protocol."""
    # Get original protocol
    original = get_protocol(protocol_id, token)
    steps = get_protocol_steps(protocol_id, token)
    
    # Create new protocol
    title = new_title or f"Copy of {original.get('title', 'Untitled')}"
    new_id = create_protocol(title, original.get('description', ''), token)
    
    if new_id and steps:
        # Add steps
        for step in steps:
            add_step(new_id, step['title'], step['description'], token=token)
    
    return new_id


def main():
    parser = argparse.ArgumentParser(description="Protocols.io API Tools")
    parser.add_argument("--token", help="API token (or set PROTOCOLSIO_TOKEN env var)")
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Search command
    search_parser = subparsers.add_parser("search", help="Search protocols")
    search_parser.add_argument("query", help="Search query")
    search_parser.add_argument("-n", "--limit", type=int, default=10, help="Max results")
    search_parser.add_argument("-f", "--filter", default="public", choices=["public", "private"])
    
    # Get command
    get_parser = subparsers.add_parser("get", help="Get protocol details")
    get_parser.add_argument("protocol_id", help="Protocol ID")
    get_parser.add_argument("--steps", action="store_true", help="Include steps")
    
    # Create command
    create_parser = subparsers.add_parser("create", help="Create new protocol")
    create_parser.add_argument("-t", "--title", required=True, help="Protocol title")
    create_parser.add_argument("-d", "--description", default="", help="Description")
    
    # Add step command
    step_parser = subparsers.add_parser("add-step", help="Add step to protocol")
    step_parser.add_argument("protocol_id", help="Protocol ID")
    step_parser.add_argument("-t", "--title", required=True, help="Step title")
    step_parser.add_argument("-d", "--description", default="", help="Step description")
    
    # Workspaces command
    ws_parser = subparsers.add_parser("workspaces", help="List workspaces")
    
    # Workspace protocols command
    ws_prot_parser = subparsers.add_parser("workspace-protocols", help="List workspace protocols")
    ws_prot_parser.add_argument("workspace_id", help="Workspace ID")
    
    # Export command
    export_parser = subparsers.add_parser("export", help="Export protocol")
    export_parser.add_argument("protocol_id", help="Protocol ID")
    export_parser.add_argument("-o", "--output", required=True, help="Output file")
    export_parser.add_argument("-f", "--format", choices=["json", "markdown"], default="json")
    
    # Clone command
    clone_parser = subparsers.add_parser("clone", help="Clone a protocol")
    clone_parser.add_argument("protocol_id", help="Protocol ID to clone")
    clone_parser.add_argument("-t", "--title", help="New title")
    
    args = parser.parse_args()
    
    if args.command == "search":
        results = search_protocols(args.query, args.token, args.filter, args.limit)
        print(json.dumps(results, indent=2))
    
    elif args.command == "get":
        protocol = get_protocol(args.protocol_id, args.token)
        if args.steps:
            protocol["steps"] = get_protocol_steps(args.protocol_id, args.token)
        print(json.dumps(protocol, indent=2))
    
    elif args.command == "create":
        protocol_id = create_protocol(args.title, args.description, args.token)
        if protocol_id:
            print(f"Created: {protocol_id}")
    
    elif args.command == "add-step":
        add_step(args.protocol_id, args.title, args.description, token=args.token)
    
    elif args.command == "workspaces":
        workspaces = list_workspaces(args.token)
        print(json.dumps(workspaces, indent=2))
    
    elif args.command == "workspace-protocols":
        protocols = get_workspace_protocols(args.workspace_id, args.token)
        print(json.dumps(protocols, indent=2))
    
    elif args.command == "export":
        export_protocol(args.protocol_id, args.output, args.format, args.token)
    
    elif args.command == "clone":
        new_id = clone_protocol(args.protocol_id, args.title, args.token)
        if new_id:
            print(f"Cloned to: {new_id}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
