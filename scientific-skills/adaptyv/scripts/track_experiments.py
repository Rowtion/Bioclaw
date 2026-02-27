#!/usr/bin/env python3
"""
Experiment Tracker - Monitor experiment status and download results
Track multiple experiments and notify when results are ready.
"""

import os
import json
import time
import argparse
from datetime import datetime
from typing import List, Dict
from adaptyv_client import AdaptyvClient


class ExperimentTracker:
    """Track status of multiple experiments."""
    
    def __init__(self, tracking_file: str = "experiments.json"):
        self.client = AdaptyvClient()
        self.tracking_file = tracking_file
        self.experiments = self._load_tracking()
    
    def _load_tracking(self) -> Dict:
        """Load experiment tracking data."""
        if os.path.exists(self.tracking_file):
            with open(self.tracking_file, 'r') as f:
                return json.load(f)
        return {"experiments": []}
    
    def _save_tracking(self):
        """Save experiment tracking data."""
        with open(self.tracking_file, 'w') as f:
            json.dump(self.experiments, f, indent=2)
    
    def add_experiment(self, experiment_id: str, name: str = ""):
        """Add experiment to tracking."""
        exp_data = {
            "id": experiment_id,
            "name": name or experiment_id,
            "added": datetime.now().isoformat(),
            "status": "submitted",
            "last_checked": None,
            "completed": False
        }
        self.experiments["experiments"].append(exp_data)
        self._save_tracking()
        print(f"Added experiment: {name or experiment_id}")
    
    def check_all(self):
        """Check status of all tracked experiments."""
        print(f"Checking {len(self.experiments['experiments'])} experiments...")
        print("="*60)
        
        for exp in self.experiments["experiments"]:
            if exp.get("completed"):
                continue
            
            try:
                status = self.client.get_experiment_status(exp["id"])
                current_status = status.get("status", "unknown")
                
                exp["last_checked"] = datetime.now().isoformat()
                exp["status"] = current_status
                
                if current_status in ["completed", "failed", "cancelled"]:
                    exp["completed"] = True
                    if current_status == "completed":
                        exp["completed_at"] = datetime.now().isoformat()
                
                status_icon = "✓" if current_status == "completed" else "⏳"
                print(f"{status_icon} {exp['name']}: {current_status}")
                
            except Exception as e:
                print(f"✗ {exp['name']}: Error checking status - {e}")
        
        self._save_tracking()
        print("="*60)
        self._print_summary()
    
    def _print_summary(self):
        """Print summary of all experiments."""
        total = len(self.experiments["experiments"])
        completed = sum(1 for e in self.experiments["experiments"] if e.get("completed"))
        pending = total - completed
        
        print(f"\nSummary:")
        print(f"  Total: {total}")
        print(f"  Completed: {completed}")
        print(f"  Pending: {pending}")
        
        if pending > 0:
            print(f"\n  Typical completion time: ~21 days")
    
    def download_results(self, experiment_id: str, output_dir: str = "results"):
        """Download results for a completed experiment."""
        os.makedirs(output_dir, exist_ok=True)
        
        try:
            results = self.client.get_results(experiment_id)
            
            output_file = os.path.join(output_dir, f"{experiment_id}_results.json")
            with open(output_file, 'w') as f:
                json.dump(results, f, indent=2)
            
            print(f"✓ Results downloaded: {output_file}")
            return results
            
        except Exception as e:
            print(f"✗ Failed to download results: {e}")
            return None
    
    def list_experiments(self):
        """List all tracked experiments."""
        print("\nTracked Experiments:")
        print("="*60)
        for exp in self.experiments["experiments"]:
            status_icon = "✓" if exp.get("completed") else "⏳"
            print(f"{status_icon} {exp['name']} ({exp['id'][:8]}...)")
            print(f"   Status: {exp['status']}")
            print(f"   Added: {exp['added'][:10]}")
            if exp.get("last_checked"):
                print(f"   Last checked: {exp['last_checked'][:10]}")
            print()


def main():
    parser = argparse.ArgumentParser(
        description="Track Adaptyv experiment status"
    )
    parser.add_argument("--tracking-file", default="experiments.json",
                       help="Experiment tracking file")
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Add command
    add_parser = subparsers.add_parser("add", help="Add experiment to tracking")
    add_parser.add_argument("experiment_id", help="Experiment ID")
    add_parser.add_argument("--name", help="Friendly name for experiment")
    
    # Check command
    subparsers.add_parser("check", help="Check status of all experiments")
    
    # List command
    subparsers.add_parser("list", help="List all tracked experiments")
    
    # Download command
    download_parser = subparsers.add_parser("download", help="Download experiment results")
    download_parser.add_argument("experiment_id", help="Experiment ID")
    download_parser.add_argument("--output", default="results", help="Output directory")
    
    args = parser.parse_args()
    
    tracker = ExperimentTracker(args.tracking_file)
    
    if args.command == "add":
        tracker.add_experiment(args.experiment_id, args.name)
    elif args.command == "check":
        tracker.check_all()
    elif args.command == "list":
        tracker.list_experiments()
    elif args.command == "download":
        tracker.download_results(args.experiment_id, args.output)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
