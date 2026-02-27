#!/usr/bin/env python3
"""
scvi-tools: Single-Cell Variational Inference Examples
=======================================================
Deep generative models for single-cell RNA-seq and multi-omics analysis.
"""

import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# Example 1: Basic scVI Analysis (Batch Correction & Latent Representation)
# ==============================================================================

def example_scvi_basic():
    """Basic scVI workflow for single-cell RNA-seq analysis."""
    import scanpy as sc
    import scvi
    
    print("=" * 60)
    print("Example 1: Basic scVI Analysis")
    print("=" * 60)
    
    # Load example dataset (heart cell atlas subsample)
    adata = scvi.data.heart_cell_atlas_subsampled()
    
    print(f"Original data shape: {adata.shape}")
    print(f"Observations: {adata.n_obs}")
    print(f"Variables (genes): {adata.n_vars}")
    print(f"Batch categories: {adata.obs['batch'].unique()}")
    
    # Preprocessing
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.highly_variable_genes(adata, n_top_genes=1200, subset=True)
    
    print(f"\nAfter preprocessing: {adata.shape}")
    
    # Setup anndata for scVI
    scvi.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch"
    )
    
    # Create and train model
    model = scvi.model.SCVI(adata)
    print("\nTraining scVI model...")
    model.train(max_epochs=20, accelerator="cpu")
    
    # Get latent representation
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    
    print(f"Latent representation shape: {latent.shape}")
    
    # Downstream analysis with scanpy
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    
    print(f"\nClusters identified: {adata.obs['leiden'].nunique()}")
    
    return model, adata


# ==============================================================================
# Example 2: Differential Expression Analysis
# ==============================================================================

def example_differential_expression():
    """Differential expression with scVI uncertainty quantification."""
    import scanpy as sc
    import scvi
    
    print("\n" + "=" * 60)
    print("Example 2: Differential Expression Analysis")
    print("=" * 60)
    
    # Load data
    adata = scvi.data.heart_cell_atlas_subsampled()
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.highly_variable_genes(adata, n_top_genes=1200, subset=True)
    
    # Setup and train
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
    model = scvi.model.SCVI(adata)
    model.train(max_epochs=20, accelerator="cpu")
    
    # Get cell type annotations (or use clustering)
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    sc.pp.neighbors(adata, use_rep="X_scVI")
    sc.tl.leiden(adata)
    
    # Check if we have cell types to compare
    if 'cell_type' in adata.obs.columns:
        cell_types = adata.obs['cell_type'].unique()
        if len(cell_types) >= 2:
            print(f"\nComparing cell types: {cell_types[0]} vs {cell_types[1]}")
            
            # Differential expression
            de_results = model.differential_expression(
                groupby="cell_type",
                group1=cell_types[0],
                group2=cell_types[1],
                mode="change"
            )
            
            print(f"\nTop differentially expressed genes:")
            print(de_results.head(10)[['proba_de', 'lfc_mean', 'lfc_std']])
            
            return de_results
    
    print("Using leiden clusters for DE analysis...")
    clusters = adata.obs['leiden'].unique()
    if len(clusters) >= 2:
        de_results = model.differential_expression(
            groupby="leiden",
            group1=clusters[0],
            group2=clusters[1]
        )
        print(f"\nDE results between clusters {clusters[0]} and {clusters[1]}:")
        print(de_results.head())
        return de_results
    
    return None


# ==============================================================================
# Example 3: scANVI - Semi-Supervised Cell Type Annotation
# ==============================================================================

def example_scanvi_annotation():
    """scANVI for semi-supervised cell type annotation."""
    import scanpy as sc
    import scvi
    
    print("\n" + "=" * 60)
    print("Example 3: scANVI - Semi-Supervised Annotation")
    print("=" * 60)
    
    # Load data with labels
    adata = scvi.data.heart_cell_atlas_subsampled()
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.highly_variable_genes(adata, n_top_genes=1200, subset=True)
    
    # Check for cell type annotations
    if 'cell_type' not in adata.obs.columns:
        print("No cell_type column found. Creating synthetic labels...")
        # Create synthetic cell types from batch
        adata.obs['cell_type'] = adata.obs['batch'].astype(str)
    
    # Mask some labels for semi-supervised learning
    cell_types = adata.obs['cell_type'].copy()
    n_cells = len(cell_types)
    unlabeled_idx = np.random.choice(n_cells, size=int(0.3*n_cells), replace=False)
    cell_types.iloc[unlabeled_idx] = "Unknown"
    adata.obs["cell_type_scanvi"] = cell_types
    
    print(f"Labeled cells: {(cell_types != 'Unknown').sum()}")
    print(f"Unlabeled cells: {(cell_types == 'Unknown').sum()}")
    
    # Setup scANVI
    scvi.model.SCANVI.setup_anndata(
        adata,
        layer="counts",
        batch_key="batch",
        labels_key="cell_type_scanvi",
        unlabeled_category="Unknown"
    )
    
    # Create scANVI model
    model = scvi.model.SCANVI(adata)
    print("\nTraining scANVI model...")
    model.train(max_epochs=20, accelerator="cpu")
    
    # Predict labels for unlabeled cells
    predictions = model.predict(adata)
    adata.obs["scanvi_predictions"] = predictions
    
    print(f"\nLabel predictions completed")
    print(f"Predicted cell types: {pd.Series(predictions).unique()}")
    
    return model, adata


# ==============================================================================
# Example 4: TotalVI - CITE-seq Protein + RNA Integration
# ==============================================================================

def example_totalvi_citeseq():
    """TotalVI for CITE-seq protein and RNA joint modeling."""
    import scanpy as sc
    import scvi
    
    print("\n" + "=" * 60)
    print("Example 4: TotalVI - CITE-seq Analysis")
    print("=" * 60)
    
    # Try to load CITE-seq data
    try:
        adata = scvi.data.pbmcs_10x_cite_seq()
        print(f"Loaded CITE-seq data: {adata.shape}")
        print(f"Protein features: {adata.obsm['protein_expression'].shape[1]}")
        
        # Preprocessing
        sc.pp.filter_genes(adata, min_counts=3)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
        
        # Setup TotalVI
        scvi.model.TOTALVI.setup_anndata(
            adata,
            protein_expression_obsm_key="protein_expression"
        )
        
        # Create and train model
        model = scvi.model.TOTALVI(adata)
        print("\nTraining TotalVI model...")
        model.train(max_epochs=20, accelerator="cpu")
        
        # Get joint latent representation
        latent = model.get_latent_representation()
        adata.obsm["X_totalVI"] = latent
        
        print(f"TotalVI latent shape: {latent.shape}")
        
        # Denoised protein expression
        denoised_protein = model.get_protein_foreground_probability()
        print(f"Denoised protein shape: {denoised_protein.shape}")
        
        return model, adata
        
    except Exception as e:
        print(f"CITE-seq data not available: {e}")
        print("Skipping TotalVI example")
        return None, None


# ==============================================================================
# Example 5: Model Saving and Loading
# ==============================================================================

def example_model_save_load():
    """Save and load trained scVI models."""
    import scanpy as sc
    import scvi
    import tempfile
    import os
    
    print("\n" + "=" * 60)
    print("Example 5: Model Save and Load")
    print("=" * 60)
    
    # Load and prepare data
    adata = scvi.data.heart_cell_atlas_subsampled()
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.highly_variable_genes(adata, n_top_genes=1000, subset=True)
    
    # Setup and train
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
    model = scvi.model.SCVI(adata)
    model.train(max_epochs=10, accelerator="cpu")
    
    # Save model
    temp_dir = tempfile.mkdtemp()
    model_path = os.path.join(temp_dir, "scvi_model")
    
    print(f"\nSaving model to: {model_path}")
    model.save(model_path, overwrite=True)
    
    # Load model
    print(f"Loading model from: {model_path}")
    loaded_model = scvi.model.SCVI.load(model_path, adata=adata)
    
    print("Model successfully saved and loaded!")
    
    # Verify loaded model works
    latent = loaded_model.get_latent_representation()
    print(f"Loaded model latent shape: {latent.shape}")
    
    # Cleanup
    import shutil
    shutil.rmtree(temp_dir)
    
    return loaded_model


# ==============================================================================
# Example 6: Solo - Doublet Detection
# ==============================================================================

def example_solo_doublet_detection():
    """SOLO for doublet detection in single-cell data."""
    import scanpy as sc
    import scvi
    from scvi.external import SOLO
    
    print("\n" + "=" * 60)
    print("Example 6: SOLO - Doublet Detection")
    print("=" * 60)
    
    # Load data
    adata = scvi.data.heart_cell_atlas_subsampled()
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.highly_variable_genes(adata, n_top_genes=1200, subset=True)
    
    # First train a SCVI model
    scvi.model.SCVI.setup_anndata(adata, layer="counts")
    vae = scvi.model.SCVI(adata)
    vae.train(max_epochs=20, accelerator="cpu")
    
    # Train SOLO model for doublet detection
    print("\nTraining SOLO doublet classifier...")
    solo = SOLO.from_scvi_model(vae)
    solo.train(max_epochs=20, accelerator="cpu")
    
    # Predict doublets
    doublet_preds = solo.predict(adata)
    adata.obs["solo_doublet"] = doublet_preds
    
    # Get probabilities
    doublet_probs = solo.predict_proba(adata)
    adata.obs["solo_doublet_prob"] = doublet_probs
    
    n_doublets = (doublet_preds == "Doublet").sum()
    print(f"\nPredicted doublets: {n_doublets} ({100*n_doublets/len(adata):.1f}%)")
    
    return solo, adata


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("scvi-tools: Single-Cell Variational Inference Examples")
    print("=" * 70)
    
    try:
        # Run examples
        example_scvi_basic()
        example_differential_expression()
        example_scanvi_annotation()
        example_totalvi_citeseq()
        example_model_save_load()
        example_solo_doublet_detection()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)
        
    except ImportError as e:
        print(f"\nImport Error: {e}")
        print("Please install required packages:")
        print("  pip install scvi-tools scanpy")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
