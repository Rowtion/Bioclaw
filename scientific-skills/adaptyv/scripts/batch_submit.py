#!/usr/bin/env python3
"""
Batch Protein Submission - Submit multiple protein sequences for testing
Reads sequences from FASTA file and submits experiments in batch.
"""

import os
import json
import argparse
from pathlib import Path
from typing import List, Dict
from adaptyv_client import AdaptyvClient


def parse_fasta(filepath: str) -> List[Dict[str, str]]:
    """Parse FASTA file and return list of sequence records."""
    sequences = []
    current_id = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences.append({
                        'id': current_id,
                        'sequence': ''.join(current_seq)
                    })
                current_id = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)
        
        if current_id:
            sequences.append({
                'id': current_id,
                'sequence': ''.join(current_seq)
            })
    
    return sequences


def format_fasta(record: Dict[str, str]) -> str:
    """Format sequence record as FASTA string."""
    return f">{record['id']}\n{record['sequence']}"


def submit_batch(fasta_file: str, 
                 experiment_type: str = "binding",
                 output_file: str = "batch_submissions.json"):
    """Submit all sequences from FASTA file."""
    
    client = AdaptyvClient()
    sequences = parse_fasta(fasta_file)
    
    print(f"Found {len(sequences)} sequences in {fasta_file}")
    
    results = []
    for i, seq_record in enumerate(sequences, 1):
        print(f"\n[{i}/{len(sequences)}] Submitting: {seq_record['id']}")
        
        try:
            fasta_str = format_fasta(seq_record)
            response = client.submit_experiment(
                sequences=fasta_str,
                experiment_type=experiment_type
            )
            
            results.append({
                "input_id": seq_record['id'],
                "experiment_id": response.get("experiment_id"),
                "status": response.get("status"),
                "submitted": True
            })
            print(f"  ✓ Submitted! ID: {response.get('experiment_id')}")
            
        except Exception as e:
            results.append({
                "input_id": seq_record['id'],
                "error": str(e),
                "submitted": False
            })
            print(f"  ✗ Failed: {e}")
    
    # Save results
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    successful = sum(1 for r in results if r.get("submitted"))
    print(f"\n{'='*50}")
    print(f"Batch submission complete!")
    print(f"Successful: {successful}/{len(sequences)}")
    print(f"Results saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Batch submit protein sequences to Adaptyv"
    )
    parser.add_argument("fasta_file", help="Input FASTA file with sequences")
    parser.add_argument("--type", default="binding",
                       choices=["binding", "expression", "thermostability", "enzyme_activity"],
                       help="Experiment type")
    parser.add_argument("--output", default="batch_submissions.json",
                       help="Output JSON file for submission tracking")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.fasta_file):
        print(f"Error: File not found: {args.fasta_file}")
        return
    
    submit_batch(args.fasta_file, args.type, args.output)


if __name__ == "__main__":
    main()
