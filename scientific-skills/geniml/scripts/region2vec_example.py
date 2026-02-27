"""
Geniml Example: Region2Vec Training
Train region embeddings on BED files.
"""
from geniml.tokenization import hard_tokenization
from geniml.region2vec import region2vec
from geniml.evaluation import evaluate_embeddings


def tokenize_bed_files(src_folder: str, dst_folder: str, universe_file: str):
    """
    Tokenize BED files using a universe reference.
    
    Args:
        src_folder: Source folder with BED files
        dst_folder: Destination folder for tokens
        universe_file: Universe BED file
    """
    hard_tokenization(
        src_folder=src_folder,
        dst_folder=dst_folder,
        universe_file=universe_file,
        p_value_threshold=1e-9
    )
    print(f"Tokenization complete. Tokens saved to {dst_folder}")


def train_region2vec(token_folder: str, save_dir: str, 
                     embedding_dim: int = 100, num_shufflings: int = 1000):
    """
    Train Region2Vec model on tokenized BED files.
    
    Args:
        token_folder: Folder with tokens
        save_dir: Directory to save model
        embedding_dim: Embedding dimension
        num_shufflings: Number of shufflings
    """
    region2vec(
        token_folder=token_folder,
        save_dir=save_dir,
        num_shufflings=num_shufflings,
        embedding_dim=embedding_dim
    )
    print(f"Model trained and saved to {save_dir}")


def evaluate_model(embeddings_file: str, labels_file: str):
    """
    Evaluate trained embeddings.
    
    Args:
        embeddings_file: Path to embeddings file
        labels_file: Path to labels file
        
    Returns:
        Evaluation metrics
    """
    metrics = evaluate_embeddings(
        embeddings_file=embeddings_file,
        labels_file=labels_file
    )
    
    print("Evaluation Metrics:")
    for metric, value in metrics.items():
        print(f"  {metric}: {value}")
    
    return metrics


if __name__ == "__main__":
    print("Geniml Region2Vec Example")
    # tokenize_bed_files('bed_files/', 'tokens/', 'universe.bed')
    # train_region2vec('tokens/', 'model/')
    # evaluate_model('model/embeddings.npy', 'metadata.csv')
