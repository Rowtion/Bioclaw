#!/usr/bin/env python3
"""
Concatenate AnnData Objects

This script demonstrates merging multiple AnnData objects
along observations or variables.

Usage:
    python concatenate_datasets.py
"""

import anndata as ad
import numpy as np
import pandas as pd


def create_batch_data(batch_name, n_cells=50, n_genes=100):
    """Create a batch of sample data."""
    X = np.random.poisson(3, size=(n_cells, n_genes))
    obs = pd.DataFrame({
        'batch': [batch_name] * n_cells,
        'cell_type': np.random.choice(['T', 'B', 'Mono'], n_cells),
    }, index=[f'{batch_name}_cell_{i}' for i in range(n_cells)])
    var = pd.DataFrame({
        'gene_name': [f'gene_{i}' for i in range(n_genes)],
    }, index=[f'gene_{i}' for i in range(n_genes)])
    
    return ad.AnnData(X=X, obs=obs, var=var)


def concatenate_observations():
    """Concatenate along observations (add more cells)."""
    print("\n" + "="*60)
    print("Concatenating Observations (Adding Cells)")
    print("="*60)
    
    # Create multiple batches
    batch1 = create_batch_data('Batch1', n_cells=50)
    batch2 = create_batch_data('Batch2', n_cells=60)
    batch3 = create_batch_data('Batch3', n_cells=40)
    
    print(f"Batch 1: {batch1.shape}")
    print(f"Batch 2: {batch2.shape}")
    print(f"Batch 3: {batch3.shape}")
    
    # Concatenate with inner join (default)
    combined_inner = ad.concat([batch1, batch2, batch3], axis=0, join='inner')
    print(f"\nCombined (inner join): {combined_inner.shape}")
    
    # Concatenate with outer join
    combined_outer = ad.concat([batch1, batch2, batch3], axis=0, join='outer')
    print(f"Combined (outer join): {combined_outer.shape}")
    
    # Concatenate with labels
    combined_labeled = ad.concat(
        [batch1, batch2, batch3],
        axis=0,
        label='dataset',
        keys=['Batch1', 'Batch2', 'Batch3']
    )
    print(f"\nCombined with labels: {combined_labeled.shape}")
    print(f"Dataset categories: {combined_labeled.obs['dataset'].unique()}")


def concatenate_variables():
    """Concatenate along variables (add more genes)."""
    print("\n" + "="*60)
    print("Concatenating Variables (Adding Genes)")
    print("="*60)
    
    # Create data with different genes
    n_cells = 100
    
    rna_data = ad.AnnData(
        X=np.random.poisson(5, size=(n_cells, 1000)),
        var=pd.DataFrame(index=[f'rna_{i}' for i in range(1000)])
    )
    
    protein_data = ad.AnnData(
        X=np.random.uniform(0, 10, size=(n_cells, 50)),
        var=pd.DataFrame(index=[f'prot_{i}' for i in range(50)])
    )
    
    print(f"RNA data: {rna_data.shape}")
    print(f"Protein data: {protein_data.shape}")
    
    # Concatenate along variables (axis=1)
    multimodal = ad.concat([rna_data, protein_data], axis=1)
    print(f"\nMultimodal data: {multimodal.shape}")


def main():
    """Run concatenation examples."""
    print("AnnData Concatenation Examples")
    print("=" * 60)
    
    concatenate_observations()
    concatenate_variables()
    
    print("\n" + "="*60)
    print("Concatenation examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
