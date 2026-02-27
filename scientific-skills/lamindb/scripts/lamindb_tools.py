#!/usr/bin/env python3
"""
LaminDB - Biological Data Management
Manage biological data with integrity, provenance, and validity.
"""

import lamindb as ln
from typing import List, Dict, Optional
import pandas as pd


def init_instance(name: str, storage: str = None):
    """Initialize a LaminDB instance."""
    ln.setup.init(storage=storage or f"./{name}_storage")
    print(f"Initialized instance: {name}")


def register_artifact(filepath: str, description: str = None, 
                     run_id: str = None) -> str:
    """Register a data artifact."""
    artifact = ln.Artifact(filepath, description=description)
    artifact.save()
    print(f"Registered: {artifact.key} (uid: {artifact.uid})")
    return artifact.uid


def query_artifacts(entity_type: str = None, limit: int = 10) -> List[Dict]:
    """Query registered artifacts."""
    artifacts = ln.Artifact.filter().limit(limit).all()
    results = []
    for art in artifacts:
        results.append({
            "uid": art.uid,
            "key": art.key,
            "description": art.description,
            "size": art.size,
            "created_at": art.created_at
        })
    return results


def transform_data(input_uid: str, transform_fn, output_path: str):
    """Transform data with provenance tracking."""
    # Load input
    input_art = ln.Artifact.get(input_uid)
    
    # Create transform record
    transform = ln.Transform(name="data_processing")
    transform.save()
    
    # Run transform
    run = ln.Run(transform=transform)
    
    # Process
    input_art.load()
    result = transform_fn(input_art)
    
    # Save output
    output_art = ln.Artifact(output_path, run=run)
    output_art.save()
    
    print(f"Transform complete. Output: {output_art.uid}")
    return output_art.uid


def main():
    import argparse
    parser = argparse.ArgumentParser(description="LaminDB Tools")
    parser.add_argument("command", choices=["init", "register", "query"])
    parser.add_argument("--name", help="Instance name")
    parser.add_argument("--file", help="File path")
    parser.add_argument("--description", help="Artifact description")
    parser.add_argument("--limit", type=int, default=10)
    
    args = parser.parse_args()
    
    if args.command == "init":
        init_instance(args.name)
    elif args.command == "register":
        register_artifact(args.file, args.description)
    elif args.command == "query":
        results = query_artifacts(limit=args.limit)
        for r in results:
            print(f"{r['uid'][:8]}... | {r['key']} | {r['description']}")


if __name__ == "__main__":
    main()
