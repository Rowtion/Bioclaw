#!/usr/bin/env python3
"""
ZINC Database CLI Tool
Command-line interface for common ZINC operations.
"""

import argparse
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from search_utils import search_by_zinc_id, search_by_smiles, get_random_compounds
    from structure_utils import add_tranche_properties, classify_by_lipinski
    from workflow_utils import prepare_docking_library
    import pandas as pd
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure you have installed: pandas, requests")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='ZINC Database Command Line Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Search by ZINC ID
  python zinc_cli.py search-id ZINC000000000001
  
  # Search by SMILES
  python zinc_cli.py search-smiles "c1ccccc1"
  
  # Find analogs (similarity search)
  python zinc_cli.py search-smiles "CC(C)Cc1ccc(cc1)C(C)C(=O)O" --dist 5
  
  # Get random compounds
  python zinc_cli.py random --count 100 --subset lead-like
  
  # Prepare docking library
  python zinc_cli.py prepare-lib input.csv --output docking_library
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Search by ZINC ID
    id_parser = subparsers.add_parser('search-id', help='Search by ZINC ID')
    id_parser.add_argument('zinc_id', help='ZINC ID or comma-separated list')
    id_parser.add_argument('--output', '-o', help='Output CSV file')
    
    # Search by SMILES
    smiles_parser = subparsers.add_parser('search-smiles', help='Search by SMILES')
    smiles_parser.add_argument('smiles', help='SMILES string')
    smiles_parser.add_argument('--dist', type=int, default=0, help='Tanimoto distance')
    smiles_parser.add_argument('--output', '-o', help='Output CSV file')
    
    # Random compounds
    random_parser = subparsers.add_parser('random', help='Get random compounds')
    random_parser.add_argument('--count', type=int, default=100, help='Number of compounds')
    random_parser.add_argument('--subset', choices=['lead-like', 'drug-like', 'fragment'],
                               help='Filter by subset')
    random_parser.add_argument('--output', '-o', help='Output CSV file')
    
    # Prepare library
    lib_parser = subparsers.add_parser('prepare-lib', help='Prepare docking library')
    lib_parser.add_argument('input', help='Input CSV file')
    lib_parser.add_argument('--output', '-o', required=True, help='Output base path')
    lib_parser.add_argument('--format', choices=['csv', 'smi', 'both'], default='both',
                           help='Output format')
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return
    
    try:
        if args.command == 'search-id':
            # Handle multiple IDs
            zinc_ids = [id.strip() for id in args.zinc_id.split(',')]
            print(f"Searching for {len(zinc_ids)} ZINC ID(s)...")
            
            results = search_by_zinc_id(zinc_ids)
            
            # Add properties
            if 'tranche' in results.columns:
                results = add_tranche_properties(results)
                results = classify_by_lipinski(results)
            
            print(f"\nFound {len(results)} compound(s)")
            print(results.to_string())
            
            if args.output:
                results.to_csv(args.output, index=False)
                print(f"\nResults saved to {args.output}")
        
        elif args.command == 'search-smiles':
            print(f"Searching for SMILES: {args.smiles}")
            if args.dist > 0:
                print(f"Using Tanimoto distance: {args.dist}")
            
            results = search_by_smiles(args.smiles, dist=args.dist)
            
            if 'tranche' in results.columns:
                results = add_tranche_properties(results)
                results = classify_by_lipinski(results)
            
            print(f"\nFound {len(results)} compound(s)")
            print(results[['zinc_id', 'smiles', 'tranche']].head(20).to_string())
            
            if args.output:
                results.to_csv(args.output, index=False)
                print(f"\nResults saved to {args.output}")
        
        elif args.command == 'random':
            print(f"Getting {args.count} random compounds...")
            if args.subset:
                print(f"Subset: {args.subset}")
            
            results = get_random_compounds(count=args.count, subset=args.subset)
            
            if 'tranche' in results.columns:
                results = add_tranche_properties(results)
                results = classify_by_lipinski(results)
            
            print(f"\nRetrieved {len(results)} compound(s)")
            print(results[['zinc_id', 'smiles', 'tranche']].head(10).to_string())
            
            if args.output:
                results.to_csv(args.output, index=False)
                print(f"\nResults saved to {args.output}")
        
        elif args.command == 'prepare-lib':
            print(f"Loading compounds from {args.input}...")
            df = pd.read_csv(args.input)
            print(f"Loaded {len(df)} compounds")
            
            formats = []
            if args.format in ['csv', 'both']:
                formats.append('csv')
            if args.format in ['smi', 'both']:
                formats.append('smi')
            
            output_files = prepare_docking_library(df, args.output, formats=formats)
            
            print(f"\nCreated files:")
            for fmt, path in output_files.items():
                print(f"  {fmt}: {path}")
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
