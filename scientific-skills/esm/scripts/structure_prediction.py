"""
ESM Example: Structure Prediction
Predict protein structure from sequence using ESM3.
"""
import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig


def predict_structure(sequence: str):
    """
    Predict 3D structure from protein sequence.
    
    Args:
        sequence: Amino acid sequence
        
    Returns:
        Protein object with predicted structure
    """
    print("Loading ESM3 model...")
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-sm-open-v1")
    
    if torch.cuda.is_available():
        model = model.to("cuda")
    
    # Create protein from sequence
    protein = ESMProtein(sequence=sequence)
    
    # Predict structure
    print("Predicting structure...")
    config = GenerationConfig(
        track="structure",
        num_steps=len(sequence)
    )
    protein_with_structure = model.generate(protein, config)
    
    # Access coordinates
    coords = protein_with_structure.coordinates
    print(f"Predicted {coords.shape[0]} atoms")
    
    # Export to PDB format
    pdb_string = protein_with_structure.to_pdb()
    
    return protein_with_structure, pdb_string


def inverse_folding(pdb_file: str):
    """
    Design sequence for a target structure (inverse folding).
    
    Args:
        pdb_file: Path to PDB structure file
        
    Returns:
        Designed protein
    """
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-sm-open-v1")
    
    if torch.cuda.is_available():
        model = model.to("cuda")
    
    # Load structure from PDB
    protein_with_structure = ESMProtein.from_pdb(pdb_file)
    
    # Remove sequence to design new one
    protein_with_structure.sequence = None
    
    # Generate sequence that folds to this structure
    config = GenerationConfig(
        track="sequence",
        num_steps=50,
        temperature=0.7
    )
    designed_protein = model.generate(protein_with_structure, config)
    
    print(f"Designed sequence: {designed_protein.sequence}")
    
    return designed_protein


if __name__ == "__main__":
    print("ESM3 Structure Prediction Example")
    # protein, pdb = predict_structure("MPRTKEINDAGLIVHSP")
