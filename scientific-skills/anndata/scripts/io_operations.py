#!/usr/bin/env python3
"""
AnnData I/O Operations

This script demonstrates reading and writing AnnData objects
in various formats.

Usage:
    python io_operations.py
"""

import anndata as ad
import numpy as np
import pandas as pd
import os


def create_sample_data():
    """Create sample AnnData for testing."""
    n_cells = 50
    n_genes = 100
    
    X = np.random.poisson(3, size=(n_cells, n_genes))
    obs = pd.DataFrame({
        'cell_type': np.random.choice(['Type_A', 'Type_B'], n_cells),
    }, index=[f'cell_{i}' for i in range(n_cells)])
    var = pd.DataFrame({
        'gene_name': [f'gene_{i}' for i in range(n_genes)],
    }, index=[f'ENSG{i:05d}' for i in range(n_genes)])
    
    return ad.AnnData(X=X, obs=obs, var=var)


def h5ad_operations():
    """Demonstrate h5ad read/write operations."""
    print("\n" + "="*60)
    print("H5AD Read/Write Operations")
    print("="*60)
    
    adata = create_sample_data()
    
    # Write h5ad file
    output_file = 'sample_data.h5ad'
    adata.write_h5ad(output_file)
    print(f"Written: {output_file} ({os.path.getsize(output_file)/1024:.1f} KB)")
    
    # Read h5ad file
    adata_read = ad.read_h5ad(output_file)
    print(f"Read back: {adata_read.shape}")
    
    # Write with compression
    compressed_file = 'sample_data_compressed.h5ad'
    adata.write_h5ad(compressed_file, compression='gzip')
    print(f"Written (compressed): {compressed_file} ({os.path.getsize(compressed_file)/1024:.1f} KB)")
    
    # Cleanup
    os.remove(output_file)
    os.remove(compressed_file)
    print("Cleanup complete")


def backed_mode_example():
    """Demonstrate backed mode for large files."""
    print("\n" + "="*60)
    print("Backed Mode Operations")
    print("="*60)
    
    adata = create_sample_data()
    output_file = 'backed_example.h5ad'
    adata.write_h5ad(output_file)
    
    # Open in backed mode (read-only)
    adata_backed = ad.read_h5ad(output_file, backed='r')
    print(f"Backed mode (read-only): {adata_backed.shape}")
    print(f"Is backed: {adata_backed.isbacked}")
    
    # Access data without loading into memory
    print(f"First 5 cells, first 5 genes shape: {adata_backed[:5, :5].X.shape}")
    
    # Load subset into memory
    subset = adata_backed[:10, :].to_memory()
    print(f"Loaded subset: {subset.shape}, is backed: {subset.isbacked}")
    
    adata_backed.file.close()
    os.remove(output_file)


def convert_csv_example():
    """Demonstrate CSV conversion."""
    print("\n" + "="*60)
    print("CSV Conversion")
    print("="*60)
    
    adata = create_sample_data()
    
    # Export to CSV
    output_dir = 'csv_export'
    adata.write_csvs(output_dir, skip_data=False)
    print(f"Exported to directory: {output_dir}")
    
    # List files
    for f in os.listdir(output_dir):
        print(f"  - {f}")
    
    # Import from CSV
    adata_imported = ad.read_csv(output_dir + '/X.csv').T
    print(f"Imported shape: {adata_imported.shape}")
    
    # Cleanup
    import shutil
    shutil.rmtree(output_dir)
    print("Cleanup complete")


def main():
    """Run all I/O examples."""
    print("AnnData I/O Operations")
    print("=" * 60)
    
    h5ad_operations()
    backed_mode_example()
    convert_csv_example()
    
    print("\n" + "="*60)
    print("I/O examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
