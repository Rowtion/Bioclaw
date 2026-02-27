"""
ESM Example: Protein Sequence Generation with ESM3
Generate novel protein sequences using ESM3.
"""
import torch
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig


def generate_protein_sequence():
    """
    Generate a protein sequence from a partial sequence prompt.
    """
    # Load model (requires GPU and model download)
    print("Loading ESM3 model...")
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-sm-open-v1")
    
    if torch.cuda.is_available():
        model = model.to("cuda")
        print("Using CUDA")
    
    # Create protein with masked positions (underscore = mask)
    protein = ESMProtein(sequence="MPRT___KEND")  # 3 positions to fill
    
    # Generate completion
    print("Generating sequence...")
    config = GenerationConfig(track="sequence", num_steps=8)
    generated = model.generate(protein, config)
    
    print(f"Input:    MPRT___KEND")
    print(f"Output:   {generated.sequence}")
    
    return generated


def design_protein_variants(base_sequence: str, num_variants: int = 5):
    """
    Generate variants of a protein sequence.
    
    Args:
        base_sequence: Starting sequence
        num_variants: Number of variants to generate
    """
    model: ESM3InferenceClient = ESM3.from_pretrained("esm3-sm-open-v1")
    
    if torch.cuda.is_available():
        model = model.to("cuda")
    
    variants = []
    
    for i in range(num_variants):
        # Mask random positions
        masked = list(base_sequence)
        import random
        mask_positions = random.sample(range(len(masked)), min(5, len(masked)//10))
        for pos in mask_positions:
            masked[pos] = '_'
        
        protein = ESMProtein(sequence=''.join(masked))
        config = GenerationConfig(track="sequence", num_steps=8, temperature=0.8)
        generated = model.generate(protein, config)
        
        variants.append({
            'variant_id': i + 1,
            'sequence': generated.sequence,
            'masked_positions': mask_positions
        })
        print(f"Variant {i+1}: {generated.sequence}")
    
    return variants


if __name__ == "__main__":
    print("ESM3 Protein Generation Example")
    print("Note: This requires the ESM package and sufficient GPU memory.")
    # generate_protein_sequence()
