#!/usr/bin/env python3
"""
Create and manipulate AnnData objects

This script demonstrates creating AnnData objects from various sources
and basic manipulation operations.

Usage:
    python create_anndata.py
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


def create_from_arrays():
    """Create AnnData from numpy arrays."""
    print("\n" + "="*60)
    print("Creating AnnData from Arrays")
    print("="*60)
    
    # Create data matrix (cells x genes)
    n_cells = 100
    n_genes = 2000
    
    X = np.random.poisson(5, size=(n_cells, n_genes))
    
    # Create observation (cell) metadata
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['T cell', 'B cell', 'Monocyte'], n_cells),
        'sample': np.random.choice(['Sample_A', 'Sample_B'], n_cells),
        'n_genes': np.sum(X > 0, axis=1),
        'n_counts': np.sum(X, axis=1),
    }, index=[f'cell_{i}' for i in range(n_cells)])
    
    # Create variable (gene) metadata
    var = pd.DataFrame({
        'gene_name': [f'Gene_{i}' for i in range(n_genes)],
        'chromosome': np.random.choice(['chr1', 'chr2', 'chrX'], n_genes),
    }, index=[f'ENSG{i:05d}' for i in range(n_genes)])
    
    # Create AnnData object
    adata = ad.AnnData(X=X, obs=obs, var=var)
    
    print(f"Created AnnData object:")
    print(f"  Shape: {adata.shape}")
    print(f"  Observations: {adata.n_obs}")
    print(f"  Variables: {adata.n_vars}")
    print(f"\nObservation columns: {list(adata.obs.columns)}")
    print(f"Variable columns: {list(adata.var.columns)}")
    
    return adata


def create_sparse_anndata():
    """Create sparse AnnData for memory efficiency."""
    print("\n" + "="*60)
    print("Creating Sparse AnnData")
    print("="*60)
    
    # Create sparse matrix (typical for single-cell data)
    n_cells = 1000
    n_genes = 5000
    
    # Create sparse matrix with ~10% non-zero values
    X_dense = np.random.poisson(1, size=(n_cells, n_genes))
    X_dense[X_dense > 5] = 0  # Sparsify
    X_sparse = csr_matrix(X_dense)
    
    adata = ad.AnnData(X=X_sparse)
    
    print(f"Dense size: {X_dense.nbytes / 1024**2:.2f} MB")
    print(f"Sparse size: {X_sparse.data.nbytes / 1024**2:.2f} MB")
    print(f"Compression: {X_dense.nbytes / X_sparse.data.nbytes:.1f}x")
    
    return adata


def add_metadata_layers(adata):
    """Demonstrate adding layers and unstructured data."""
    print("\n" + "="*60)
    print("Adding Layers and Metadata")
    print("="*60)
    
    # Add layers for different data representations
    adata.layers['counts'] = adata.X.copy()
    adata.layers['log1p'] = np.log1p(adata.X)
    adata.layers['scaled'] = (adata.X - adata.X.mean(axis=0)) / adata.X.std(axis=0)
    
    print(f"Added layers: {list(adata.layers.keys())}")
    
    # Add observation-level matrices (e.g., PCA, UMAP)
    n_pcs = 50
    adata.obsm['X_pca'] = np.random.randn(adata.n_obs, n_pcs)
    adata.obsm['X_umap'] = np.random.randn(adata.n_obs, 2)
    
    print(f"Added obsm: {list(adata.obsm.keys())}")
    
    # Add unstructured metadata
    adata.uns['pca'] = {'variance_ratio': np.random.rand(n_pcs)}
    adata.uns['neighbors'] = {'connectivities_key': 'connectivities'}
    
    print(f"Added uns: {list(adata.uns.keys())}")
    
    return adata


def subset_operations(adata):
    """Demonstrate subsetting operations."""
    print("\n" + "="*60)
    print("Subsetting Operations")
    print("="*60)
    
    # Subset by cell type
    if 'cell_type' in adata.obs.columns:
        t_cells = adata[adata.obs['cell_type'] == 'T cell']
        print(f"T cells: {t_cells.n_obs}")
        
        # Multiple conditions
        subset = adata[
            (adata.obs['cell_type'] == 'T cell') & 
            (adata.obs['sample'] == 'Sample_A')
        ]
        print(f"T cells from Sample_A: {subset.n_obs}")
    
    # Subset by gene expression
    high_expr_genes = adata[:, adata.X.sum(axis=0) > np.median(adata.X.sum(axis=0))]
    print(f"High expression genes: {high_expr_genes.n_vars}")
    
    # Index-based subsetting
    first_10_cells = adata[:10, :]
    print(f"First 10 cells: {first_10_cells.n_obs}")
    
    return adata


def main():
    """Run all examples."""
    print("AnnData Creation and Manipulation Examples")
    print("=" * 60)
    
    # Example 1: Create from arrays
    adata = create_from_arrays()
    
    # Example 2: Create sparse
    sparse_adata = create_sparse_anndata()
    
    # Example 3: Add metadata
    adata = add_metadata_layers(adata)
    
    # Example 4: Subsetting
    adata = subset_operations(adata)
    
    print("\n" + "="*60)
    print("Examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
