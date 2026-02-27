"""
ESM Example: Protein Embeddings with ESM C
Generate embeddings for protein sequences using ESM C.
"""
import torch
from esm.models.esmc import ESMC
from esm.sdk.api import ESMProtein


def generate_protein_embeddings(sequences: list):
    """
    Generate embeddings for a list of protein sequences.
    
    Args:
        sequences: List of amino acid sequences
        
    Returns:
        List of embeddings
    """
    print("Loading ESM C model...")
    model = ESMC.from_pretrained("esmc-300m")
    
    if torch.cuda.is_available():
        model = model.to("cuda")
    
    embeddings = []
    
    for seq in sequences:
        protein = ESMProtein(sequence=seq)
        protein_tensor = model.encode(protein)
        
        # Generate embeddings
        output = model.forward(protein_tensor)
        
        # Get pooled representation (mean across sequence)
        pooled = output.mean(dim=1).detach().cpu().numpy()
        embeddings.append(pooled)
        
        print(f"Processed sequence of length {len(seq)}")
    
    return embeddings


def batch_embed_sequences(sequences: list, batch_size: int = 8):
    """
    Process sequences in batches for efficiency.
    
    Args:
        sequences: List of sequences
        batch_size: Batch size
        
    Returns:
        List of embeddings
    """
    model = ESMC.from_pretrained("esmc-300m")
    
    if torch.cuda.is_available():
        model = model.to("cuda")
    
    all_embeddings = []
    
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        proteins = [ESMProtein(sequence=seq) for seq in batch]
        
        # Encode batch
        tensors = [model.encode(p) for p in proteins]
        
        # Process each
        for tensor in tensors:
            output = model.forward(tensor)
            pooled = output.mean(dim=1).detach().cpu().numpy()
            all_embeddings.append(pooled)
        
        print(f"Processed batch {i//batch_size + 1}")
    
    return all_embeddings


if __name__ == "__main__":
    sample_sequences = [
        "MPRTKEINDAGLIVHSP",
        "AGLIVHSPQKTEFLNDGR",
        "KTEFLNDGRMPRTKEIND"
    ]
    print("ESM C Embeddings Example")
    # embeddings = generate_protein_embeddings(sample_sequences)
