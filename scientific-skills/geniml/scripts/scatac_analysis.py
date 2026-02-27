"""
Geniml Example: scATAC-seq Analysis
Single-cell ATAC-seq analysis with scEmbed.
"""
import scanpy as sc
from geniml.scembed import ScEmbed
from geniml.io import tokenize_cells


def preprocess_scatac(adata_file: str, universe_file: str, 
                      output_tokens: str):
    """
    Preprocess scATAC-seq data and tokenize cells.
    
    Args:
        adata_file: Path to AnnData file
        universe_file: Universe BED file
        output_tokens: Output token file
    """
    # Tokenize cells
    tokenize_cells(
        adata=adata_file,
        universe_file=universe_file,
        output=output_tokens
    )
    print(f"Cells tokenized: {output_tokens}")


def train_scembed(tokens_file: str, embedding_dim: int = 100, 
                  epochs: int = 100):
    """
    Train scEmbed model.
    
    Args:
        tokens_file: Token file
        embedding_dim: Embedding dimension
        epochs: Training epochs
        
    Returns:
        Trained model
    """
    model = ScEmbed(embedding_dim=embedding_dim)
    model.train(dataset=tokens_file, epochs=epochs)
    
    print("scEmbed model trained")
    return model


def analyze_embeddings(model, adata_file: str, output_file: str):
    """
    Generate embeddings and analyze with scanpy.
    
    Args:
        model: Trained ScEmbed model
        adata_file: AnnData file
        output_file: Output file
    """
    # Load data
    adata = sc.read_h5ad(adata_file)
    
    # Generate embeddings
    embeddings = model.encode(adata)
    adata.obsm['scembed_X'] = embeddings
    
    # Cluster with scanpy
    sc.pp.neighbors(adata, use_rep='scembed_X')
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    
    # Save results
    adata.write(output_file)
    print(f"Results saved: {output_file}")
    
    # Visualize
    sc.pl.umap(adata, color='leiden', save='_clusters.png')


if __name__ == "__main__":
    print("Geniml scATAC-seq Analysis Example")
    # model = train_scembed('tokens.parquet')
    # analyze_embeddings(model, 'scatac_data.h5ad', 'output.h5ad')
