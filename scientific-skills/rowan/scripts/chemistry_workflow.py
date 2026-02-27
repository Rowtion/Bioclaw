#!/usr/bin/env python3
"""
Rowan Quantum Chemistry Workflow Scripts
Utilities for cloud-based computational chemistry workflows.
"""

import argparse
import json
import os
from pathlib import Path


def setup_rowan():
    """Verify Rowan API setup."""
    try:
        import rowan
    except ImportError:
        raise ImportError("rowan-python is required. Install with: uv pip install rowan-python")
    
    api_key = os.environ.get('ROWAN_API_KEY')
    if not api_key:
        print("Warning: ROWAN_API_KEY not set in environment")
        print("Set with: export ROWAN_API_KEY='your_key'")
        return False
    
    try:
        user = rowan.whoami()
        print(f"Connected as: {user.username}")
        print(f"Credits available: {user.credits}")
        return True
    except Exception as e:
        print(f"Connection failed: {e}")
        return False


def predict_pka(smiles: str, output_file: str = None) -> dict:
    """
    Predict pKa for a molecule.
    
    Args:
        smiles: SMILES string
        output_file: Optional output JSON file
        
    Returns:
        pKa prediction results
    """
    try:
        import rowan
        import stjames
    except ImportError:
        raise ImportError("rowan-python and stjames are required")
    
    # Create molecule
    mol = stjames.Molecule.from_smiles(smiles)
    
    # Submit workflow
    workflow = rowan.submit_pka_workflow(
        initial_molecule=mol,
        name=f"pKa: {smiles[:30]}"
    )
    
    print(f"Workflow submitted: {workflow.uuid}")
    print("Waiting for results...")
    
    # Wait for completion
    workflow.wait_for_result()
    workflow.fetch_latest(in_place=True)
    
    results = {
        'smiles': smiles,
        'strongest_acid_pka': workflow.data.get('strongest_acid'),
        'all_pkas': workflow.data.get('all_pkas', []),
        'workflow_id': workflow.uuid,
    }
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Results saved to: {output_file}")
    
    return results


def run_conformer_search(smiles: str, output_dir: str = None) -> dict:
    """
    Run conformer search for a molecule.
    
    Args:
        smiles: SMILES string
        output_dir: Optional output directory
        
    Returns:
        Conformer search results
    """
    try:
        import rowan
        import stjames
    except ImportError:
        raise ImportError("rowan-python and stjames are required")
    
    mol = stjames.Molecule.from_smiles(smiles)
    
    workflow = rowan.submit_conformer_search_workflow(
        initial_molecule=mol,
        name=f"Conformer search: {smiles[:30]}"
    )
    
    print(f"Workflow submitted: {workflow.uuid}")
    workflow.wait_for_result()
    workflow.fetch_latest(in_place=True)
    
    conformers = workflow.data.get('conformers', [])
    
    results = {
        'smiles': smiles,
        'num_conformers': len(conformers),
        'conformers': [
            {
                'energy': c.get('energy'),
                'index': i
            }
            for i, c in enumerate(conformers[:10])  # First 10
        ],
        'workflow_id': workflow.uuid,
    }
    
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, 'conformer_results.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"Results saved to: {output_file}")
    
    return results


def run_geometry_optimization(smiles: str, output_file: str = None) -> dict:
    """
    Run geometry optimization.
    
    Args:
        smiles: SMILES string
        output_file: Optional output file
        
    Returns:
        Optimization results
    """
    try:
        import rowan
        import stjames
    except ImportError:
        raise ImportError("rowan-python and stjames are required")
    
    mol = stjames.Molecule.from_smiles(smiles)
    
    workflow = rowan.submit_basic_calculation_workflow(
        initial_molecule=mol,
        name=f"Optimization: {smiles[:30]}",
        workflow_type="optimization"
    )
    
    print(f"Workflow submitted: {workflow.uuid}")
    workflow.wait_for_result()
    workflow.fetch_latest(in_place=True)
    
    final_mol = workflow.data.get('final_molecule')
    
    results = {
        'smiles': smiles,
        'final_energy': final_mol.energy if final_mol else None,
        'workflow_id': workflow.uuid,
    }
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
    
    return results


def batch_pka_prediction(smiles_list: list, output_file: str = None) -> list:
    """
    Run batch pKa predictions.
    
    Args:
        smiles_list: List of SMILES strings
        output_file: Optional output JSON file
        
    Returns:
        List of results
    """
    try:
        import rowan
        import stjames
    except ImportError:
        raise ImportError("rowan-python and stjames are required")
    
    # Create molecules
    molecules = [stjames.Molecule.from_smiles(smi) for smi in smiles_list]
    
    # Submit batch
    print(f"Submitting {len(molecules)} pKa calculations...")
    results = rowan.batch_pka(molecules)
    
    output = []
    for smiles, result in zip(smiles_list, results):
        output.append({
            'smiles': smiles,
            'strongest_acid_pka': result.strongest_acid if result else None,
        })
    
    if output_file:
        with open(output_file, 'w') as f:
            json.dump(output, f, indent=2)
        print(f"Batch results saved to: {output_file}")
    
    return output


def list_workflows(status: str = None, limit: int = 10) -> list:
    """
    List recent workflows.
    
    Args:
        status: Filter by status (running, completed, failed)
        limit: Maximum number to return
        
    Returns:
        List of workflows
    """
    try:
        import rowan
    except ImportError:
        raise ImportError("rowan-python is required")
    
    workflows = rowan.list_workflows(size=limit, status=status)
    
    results = []
    for wf in workflows:
        results.append({
            'uuid': wf.uuid,
            'name': wf.name,
            'status': wf.status,
            'created': wf.created.isoformat() if hasattr(wf, 'created') else None,
        })
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Rowan Quantum Chemistry Workflows')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Setup command
    setup_parser = subparsers.add_parser('setup', help='Verify Rowan setup')
    
    # pKa command
    pka_parser = subparsers.add_parser('pka', help='Predict pKa')
    pka_parser.add_argument('smiles', help='SMILES string')
    pka_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Conformer command
    conf_parser = subparsers.add_parser('conformer', help='Conformer search')
    conf_parser.add_argument('smiles', help='SMILES string')
    conf_parser.add_argument('-o', '--output', help='Output directory')
    
    # Optimize command
    opt_parser = subparsers.add_parser('optimize', help='Geometry optimization')
    opt_parser.add_argument('smiles', help='SMILES string')
    opt_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Batch pKa prediction')
    batch_parser.add_argument('input', help='Input file with SMILES (one per line)')
    batch_parser.add_argument('-o', '--output', help='Output JSON file')
    
    # List command
    list_parser = subparsers.add_parser('list', help='List workflows')
    list_parser.add_argument('-s', '--status', help='Filter by status')
    list_parser.add_argument('-n', '--limit', type=int, default=10, help='Limit results')
    
    args = parser.parse_args()
    
    if args.command == 'setup':
        setup_rowan()
        
    elif args.command == 'pka':
        results = predict_pka(args.smiles, args.output)
        print(f"\npKa Results:")
        print(f"  SMILES: {results['smiles']}")
        print(f"  Strongest acid pKa: {results['strongest_acid_pka']}")
        
    elif args.command == 'conformer':
        results = run_conformer_search(args.smiles, args.output)
        print(f"\nConformer Results:")
        print(f"  Found {results['num_conformers']} conformers")
        
    elif args.command == 'optimize':
        results = run_geometry_optimization(args.smiles, args.output)
        print(f"\nOptimization Results:")
        print(f"  Final energy: {results['final_energy']}")
        
    elif args.command == 'batch':
        with open(args.input) as f:
            smiles_list = [line.strip() for line in f if line.strip()]
        results = batch_pka_prediction(smiles_list, args.output)
        print(f"\nBatch results for {len(results)} molecules:")
        for r in results:
            print(f"  {r['smiles'][:40]:<40} pKa={r['strongest_acid_pka']}")
            
    elif args.command == 'list':
        workflows = list_workflows(args.status, args.limit)
        print(f"\nRecent Workflows:")
        print(f"{'UUID':<36} {'Status':<12} {'Name'}")
        print("-" * 80)
        for wf in workflows:
            print(f"{wf['uuid']:<36} {wf['status']:<12} {wf['name'][:40]}")
            
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
