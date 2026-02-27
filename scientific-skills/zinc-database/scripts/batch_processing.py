#!/usr/bin/env python3
"""
ZINC Database Batch Processing
Process large compound lists efficiently.
"""

import pandas as pd
from typing import List, Callable, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import time


def batch_process(
    items: List,
    process_func: Callable,
    batch_size: int = 100,
    max_workers: int = 4,
    delay: float = 0.1
) -> List:
    """
    Process items in batches with parallel execution.
    
    Args:
        items: List of items to process
        process_func: Function to apply to each item
        batch_size: Number of items per batch
        max_workers: Number of parallel workers
        delay: Delay between batches (for rate limiting)
    
    Returns:
        List of results
    """
    results = []
    
    def process_batch(batch):
        batch_results = []
        for item in batch:
            try:
                result = process_func(item)
                batch_results.append(result)
            except Exception as e:
                print(f"Error processing {item}: {e}")
                batch_results.append(None)
        return batch_results
    
    # Split into batches
    batches = [items[i:i + batch_size] for i in range(0, len(items), batch_size)]
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_batch, batch): i 
                   for i, batch in enumerate(batches)}
        
        for future in as_completed(futures):
            batch_idx = futures[future]
            try:
                batch_results = future.result()
                results.extend(batch_results)
                print(f"Completed batch {batch_idx + 1}/{len(batches)}")
            except Exception as e:
                print(f"Batch {batch_idx} failed: {e}")
            
            time.sleep(delay)
    
    return results


def process_zinc_ids_file(
    input_file: str,
    output_file: str,
    process_func: Callable,
    id_col: str = 'zinc_id',
    batch_size: int = 100
):
    """
    Process a file of ZINC IDs and save results.
    
    Args:
        input_file: Input CSV file with ZINC IDs
        output_file: Output CSV file path
        process_func: Function to process each ID
        id_col: Column name for ZINC IDs
        batch_size: Number of IDs per batch
    """
    df = pd.read_csv(input_file)
    
    if id_col not in df.columns:
        raise ValueError(f"Column '{id_col}' not found in input file")
    
    results = batch_process(
        df[id_col].tolist(),
        process_func,
        batch_size=batch_size
    )
    
    df['result'] = results
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")


def deduplicate_compounds(
    df: pd.DataFrame,
    smiles_col: str = 'smiles',
    keep: str = 'first'
) -> pd.DataFrame:
    """
    Remove duplicate compounds from a DataFrame.
    
    Args:
        df: DataFrame with SMILES column
        smiles_col: Column name for SMILES
        keep: Which duplicate to keep ('first', 'last', False)
    
    Returns:
        Deduplicated DataFrame
    """
    if smiles_col not in df.columns:
        return df
    
    original_count = len(df)
    deduped = df.drop_duplicates(subset=[smiles_col], keep=keep)
    removed = original_count - len(deduped)
    
    print(f"Removed {removed} duplicates ({removed/original_count*100:.1f}%)")
    return deduped


def sample_library(
    df: pd.DataFrame,
    n: int,
    method: str = 'random',
    random_state: int = 42
) -> pd.DataFrame:
    """
    Sample compounds from a library.
    
    Args:
        df: DataFrame with compounds
        n: Number of compounds to sample
        method: Sampling method ('random', 'diverse')
        random_state: Random seed
    
    Returns:
        Sampled DataFrame
    """
    if method == 'random':
        return df.sample(n=min(n, len(df)), random_state=random_state)
    elif method == 'diverse':
        # Simple diversity sampling: spread across MW range
        df = df.copy()
        if 'mw' in df.columns:
            df = df.sort_values('mw')
            indices = list(range(0, len(df), len(df) // n))[:n]
            return df.iloc[indices]
        else:
            return df.sample(n=min(n, len(df)), random_state=random_state)
    else:
        raise ValueError(f"Unknown sampling method: {method}")


def merge_results(
    result_files: List[str],
    output_file: str,
    merge_on: str = 'zinc_id'
):
    """
    Merge multiple result files into one.
    
    Args:
        result_files: List of CSV files to merge
        output_file: Output CSV file path
        merge_on: Column to merge on
    """
    dfs = [pd.read_csv(f) for f in result_files]
    
    # Start with first DataFrame
    merged = dfs[0]
    
    # Merge remaining DataFrames
    for df in dfs[1:]:
        merged = merged.merge(df, on=merge_on, how='outer')
    
    merged.to_csv(output_file, index=False)
    print(f"Merged {len(result_files)} files into {output_file}")
    print(f"Total compounds: {len(merged)}")


def validate_compound_batch(
    zinc_ids: List[str],
    require_smiles: bool = True,
    require_properties: bool = False
) -> pd.DataFrame:
    """
    Validate a batch of compounds by checking ZINC database.
    
    Args:
        zinc_ids: List of ZINC IDs to validate
        require_smiles: Require valid SMILES
        require_properties: Require property data
    
    Returns:
        DataFrame with validation results
    """
    # Import here to avoid circular dependency
    from .search_utils import batch_search_zinc_ids
    
    results = batch_search_zinc_ids(zinc_ids)
    
    validation = []
    for zinc_id in zinc_ids:
        match = results[results['zinc_id'] == zinc_id]
        
        is_valid = len(match) > 0
        if require_smiles and is_valid:
            is_valid = pd.notna(match.iloc[0].get('smiles'))
        
        validation.append({
            'zinc_id': zinc_id,
            'found_in_zinc': len(match) > 0,
            'valid': is_valid
        })
    
    return pd.DataFrame(validation)


if __name__ == "__main__":
    print("ZINC Database Batch Processing Utilities")
    print("=" * 40)
    
    # Example: Process function
    def example_process(x):
        time.sleep(0.01)  # Simulate work
        return x * 2
    
    # Batch processing example
    print("\nBatch processing example...")
    items = list(range(100))
    results = batch_process(
        items,
        example_process,
        batch_size=10,
        max_workers=2,
        delay=0.01
    )
    print(f"Processed {len(results)} items")
    print(f"First 5 results: {results[:5]}")
    
    # Sampling example
    print("\nSampling example...")
    sample_df = pd.DataFrame({
        'zinc_id': [f'ZINC{i:012d}' for i in range(1000)],
        'smiles': ['CCO'] * 1000,
        'mw': range(100, 1100)
    })
    
    sampled = sample_library(sample_df, n=10, method='diverse')
    print(f"Sampled {len(sampled)} compounds")
    print(sampled[['zinc_id', 'mw']].head())
