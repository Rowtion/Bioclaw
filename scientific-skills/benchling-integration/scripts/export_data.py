#!/usr/bin/env python3
"""
Export Data from Benchling

This script demonstrates exporting data from Benchling
for external analysis.

Usage:
    export BENCHLING_API_KEY="your_key"
    python export_data.py --output export.csv
"""

import os
import sys
import csv
import argparse


def example_export_sequences():
    """Example: Export sequences to CSV."""
    print("\n" + "="*60)
    print("Exporting Sequences to CSV")
    print("="*60)
    
    example_code = '''
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
import csv

benchling = Benchling(
    url=os.environ['BENCHLING_URL'],
    auth_method=ApiKeyAuth(os.environ['BENCHLING_API_KEY'])
)

# Export sequences
sequences = benchling.dna_sequences.list()

with open('sequences_export.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['ID', 'Name', 'Length', 'Bases'])
    
    for page in sequences:
        for seq in page:
            writer.writerow([
                seq.id,
                seq.name,
                len(seq.bases),
                seq.bases[:50] + '...' if len(seq.bases) > 50 else seq.bases
            ])

print("Exported to: sequences_export.csv")
'''
    
    print(example_code)


def example_export_with_metadata():
    """Example: Export with custom metadata."""
    print("\n" + "="*60)
    print("Exporting with Metadata")
    print("="*60)
    
    example_code = '''
# Export with custom fields
export_data = []

for page in benchling.dna_sequences.list():
    for seq in page:
        if seq.schema_id == "your_target_schema":
            export_data.append({
                'id': seq.id,
                'name': seq.name,
                'bases': seq.bases,
                'length': len(seq.bases),
                'folder': seq.folder_id,
                'created': seq.created_at
            })

# Save to JSON
import json
with open('sequences.json', 'w') as f:
    json.dump(export_data, f, indent=2, default=str)
'''
    
    print(example_code)


def main():
    """Display export examples."""
    parser = argparse.ArgumentParser(description="Export Benchling data")
    parser.add_argument("--output", "-o", default="export.csv", help="Output file")
    args = parser.parse_args()
    
    print("Benchling Data Export Examples")
    print("=" * 60)
    
    example_export_sequences()
    example_export_with_metadata()
    
    print("\n" + "="*60)
    print("Examples completed!")
    print(f"Output would be saved to: {args.output}")
    print("="*60)


if __name__ == "__main__":
    main()
