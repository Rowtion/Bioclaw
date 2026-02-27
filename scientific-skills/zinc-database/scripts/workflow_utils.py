#!/usr/bin/env python3
"""
ZINC Database Workflow Utilities
Virtual screening workflows and compound library preparation.
"""

import pandas as pd
from typing import List, Optional
from pathlib import Path


def prepare_docking_library(
    compounds_df: pd.DataFrame,
    output_path: str,
    formats: List[str] = ['csv', 'smi'],
    include_properties: bool = True
) -> dict:
    """
    Prepare a compound library for molecular docking.
    
    Args:
        compounds_df: DataFrame with compounds (must have 'zinc_id' and 'smiles' columns)
        output_path: Base path for output files (without extension)
        formats: Output formats ('csv', 'smi', 'sdf')
        include_properties: Include molecular properties in output
    
    Returns:
        Dictionary with output file paths
    """
    output_files = {}
    base_path = Path(output_path)
    
    # Ensure required columns exist
    if 'zinc_id' not in compounds_df.columns or 'smiles' not in compounds_df.columns:
        raise ValueError("DataFrame must have 'zinc_id' and 'smiles' columns")
    
    # CSV format
    if 'csv' in formats:
        csv_path = f"{base_path}.csv"
        if include_properties:
            compounds_df.to_csv(csv_path, index=False)
        else:
            compounds_df[['zinc_id', 'smiles']].to_csv(csv_path, index=False)
        output_files['csv'] = csv_path
    
    # SMILES format (standard .smi file)
    if 'smi' in formats:
        smi_path = f"{base_path}.smi"
        with open(smi_path, 'w') as f:
            for _, row in compounds_df.iterrows():
                f.write(f"{row['smiles']}\t{row['zinc_id']}\n")
        output_files['smi'] = smi_path
    
    return output_files


def filter_for_docking(
    df: pd.DataFrame,
    max_mw: float = 500,
    max_logp: float = 5,
    max_rotatable_bonds: Optional[int] = None,
    min_drug_likeness: bool = True
) -> pd.DataFrame:
    """
    Filter compounds for molecular docking based on drug-like properties.
    
    Args:
        df: DataFrame with molecular properties
        max_mw: Maximum molecular weight
        max_logp: Maximum LogP
        max_rotatable_bonds: Maximum number of rotatable bonds
        min_drug_likeness: Only include drug-like compounds
    
    Returns:
        Filtered DataFrame
    """
    filtered = df.copy()
    
    # Filter by molecular weight
    if 'mw' in filtered.columns:
        filtered = filtered[filtered['mw'] <= max_mw]
    
    # Filter by LogP
    if 'logP' in filtered.columns:
        filtered = filtered[filtered['logP'] <= max_logp]
    
    # Filter by drug-likeness
    if min_drug_likeness and 'drug_like' in filtered.columns:
        filtered = filtered[filtered['drug_like'] == True]
    
    return filtered


def generate_3d_download_urls(
    zinc_ids: List[str],
    file_format: str = 'mol2'
) -> List[str]:
    """
    Generate URLs for downloading 3D structures from ZINC file repository.
    
    Args:
        zinc_ids: List of ZINC IDs
        file_format: File format ('mol2', 'sdf', 'db2.gz')
    
    Returns:
        List of download URLs
    """
    base_url = "https://files.docking.org/zinc22"
    urls = []
    
    for zinc_id in zinc_ids:
        # Extract tranche info from ZINC ID if available
        # URL format depends on tranche structure
        # This is a simplified version - actual URLs require tranche lookup
        url = f"{base_url}/{zinc_id}.{file_format}"
        urls.append(url)
    
    return urls


def create_virtual_screening_list(
    compounds_df: pd.DataFrame,
    target_properties: Optional[dict] = None,
    diversity_threshold: float = 0.7
) -> pd.DataFrame:
    """
    Create a curated list for virtual screening.
    
    Args:
        compounds_df: DataFrame with compounds
        target_properties: Target property ranges (e.g., {'mw': (300, 400)})
        diversity_threshold: Minimum diversity for compound selection
    
    Returns:
        Curated DataFrame for virtual screening
    """
    vs_list = compounds_df.copy()
    
    # Apply target property filters
    if target_properties:
        for prop, (min_val, max_val) in target_properties.items():
            if prop in vs_list.columns:
                vs_list = vs_list[
                    (vs_list[prop] >= min_val) & 
                    (vs_list[prop] <= max_val)
                ]
    
    # Remove duplicates by SMILES
    if 'smiles' in vs_list.columns:
        vs_list = vs_list.drop_duplicates(subset=['smiles'])
    
    return vs_list


def split_library_for_parallel_docking(
    compounds_df: pd.DataFrame,
    n_splits: int = 10,
    output_dir: str = "./split_libraries"
) -> List[str]:
    """
    Split compound library into chunks for parallel docking.
    
    Args:
        compounds_df: DataFrame with compounds
        n_splits: Number of splits
        output_dir: Directory to save split files
    
    Returns:
        List of output file paths
    """
    from pathlib import Path
    import os
    
    os.makedirs(output_dir, exist_ok=True)
    
    chunk_size = len(compounds_df) // n_splits
    output_files = []
    
    for i in range(n_splits):
        start_idx = i * chunk_size
        if i == n_splits - 1:
            end_idx = len(compounds_df)
        else:
            end_idx = (i + 1) * chunk_size
        
        chunk = compounds_df.iloc[start_idx:end_idx]
        output_path = f"{output_dir}/library_chunk_{i+1:03d}.csv"
        chunk.to_csv(output_path, index=False)
        output_files.append(output_path)
    
    return output_files


def analyze_library_diversity(
    compounds_df: pd.DataFrame,
    smiles_col: str = 'smiles'
) -> dict:
    """
    Analyze chemical diversity of a compound library.
    
    Args:
        compounds_df: DataFrame with SMILES
        smiles_col: Column name for SMILES
    
    Returns:
        Dictionary with diversity metrics
    """
    if smiles_col not in compounds_df.columns:
        return {'error': f'Column {smiles_col} not found'}
    
    metrics = {
        'total_compounds': len(compounds_df),
        'unique_smiles': compounds_df[smiles_col].nunique(),
        'duplicates': len(compounds_df) - compounds_df[smiles_col].nunique()
    }
    
    # Subset distribution if available
    if 'subset' in compounds_df.columns:
        metrics['subset_distribution'] = compounds_df['subset'].value_counts().to_dict()
    
    # Property ranges if available
    for prop in ['mw', 'logP', 'h_donors']:
        if prop in compounds_df.columns:
            metrics[f'{prop}_range'] = {
                'min': compounds_df[prop].min(),
                'max': compounds_df[prop].max(),
                'mean': compounds_df[prop].mean()
            }
    
    return metrics


def merge_libraries(
    library_paths: List[str],
    output_path: str,
    remove_duplicates: bool = True
) -> pd.DataFrame:
    """
    Merge multiple compound libraries into one.
    
    Args:
        library_paths: List of library file paths (CSV format)
        output_path: Path for merged output
        remove_duplicates: Remove duplicate compounds by SMILES
    
    Returns:
        Merged DataFrame
    """
    dfs = []
    for path in library_paths:
        df = pd.read_csv(path)
        dfs.append(df)
    
    merged = pd.concat(dfs, ignore_index=True)
    
    if remove_duplicates and 'smiles' in merged.columns:
        merged = merged.drop_duplicates(subset=['smiles'])
    
    merged.to_csv(output_path, index=False)
    return merged


if __name__ == "__main__":
    print("ZINC Database Workflow Utilities")
    print("=" * 40)
    
    # Create sample data
    print("\nCreating sample library...")
    sample_lib = pd.DataFrame({
        'zinc_id': [f'ZINC{i:012d}' for i in range(1, 101)],
        'smiles': ['CCO', 'c1ccccc1', 'CC(C)C', 'CCCC'] * 25,
        'mw': [150, 250, 180, 220] * 25,
        'logP': [1.5, 2.5, 2.0, 2.2] * 25,
        'h_donors': [1, 0, 0, 0] * 25
    })
    
    # Prepare docking library
    print("\nPreparing docking library...")
    output_files = prepare_docking_library(
        sample_lib,
        output_path="example_library",
        formats=['csv', 'smi']
    )
    print(f"Created files: {output_files}")
    
    # Analyze diversity
    print("\nLibrary diversity analysis:")
    metrics = analyze_library_diversity(sample_lib)
    for key, value in metrics.items():
        print(f"  {key}: {value}")
    
    # Split for parallel docking
    print("\nSplitting library for parallel docking...")
    split_files = split_library_for_parallel_docking(
        sample_lib,
        n_splits=5,
        output_dir="./split_libraries"
    )
    print(f"Created {len(split_files)} split files")
