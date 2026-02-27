#!/usr/bin/env python3
"""
ZINC Database RDKit Integration
Advanced cheminformatics analysis using RDKit.
"""

import pandas as pd
from typing import List, Optional, Dict
import warnings

# RDKit imports with error handling
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem, DataStructs, Draw
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcNumRotatableBonds
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    warnings.warn("RDKit not available. Install with: conda install -c conda-forge rdkit")


def smiles_to_mol(smiles: str) -> Optional['Chem.Mol']:
    """
    Convert SMILES to RDKit molecule object.
    
    Args:
        smiles: SMILES string
    
    Returns:
        RDKit Mol object or None if invalid
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except:
        return None


def calculate_descriptors(smiles: str) -> Optional[Dict[str, float]]:
    """
    Calculate molecular descriptors for a SMILES string.
    
    Args:
        smiles: SMILES string
    
    Returns:
        Dictionary of molecular descriptors
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    mol = smiles_to_mol(smiles)
    if mol is None:
        return None
    
    descriptors = {
        'mw': Descriptors.MolWt(mol),
        'logP': Descriptors.MolLogP(mol),
        'tpsa': Descriptors.TPSA(mol),
        'hbd': Descriptors.NumHDonors(mol),
        'hba': Descriptors.NumHAcceptors(mol),
        'rotatable_bonds': CalcNumRotatableBonds(mol),
        'rings': Descriptors.RingCount(mol),
        'aromatic_rings': Descriptors.NumAromaticRings(mol),
        'formal_charge': Chem.GetFormalCharge(mol)
    }
    
    return descriptors


def add_rdkit_descriptors(df: pd.DataFrame, smiles_col: str = 'smiles') -> pd.DataFrame:
    """
    Add RDKit molecular descriptors to a DataFrame.
    
    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column
    
    Returns:
        DataFrame with added descriptor columns
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    df = df.copy()
    
    # Calculate descriptors for each SMILES
    descriptors_list = []
    for smiles in df[smiles_col]:
        desc = calculate_descriptors(smiles)
        descriptors_list.append(desc)
    
    # Convert to DataFrame and merge
    desc_df = pd.DataFrame(descriptors_list)
    
    # Add prefix to descriptor columns
    desc_df.columns = [f'rdkit_{col}' for col in desc_df.columns]
    
    result = pd.concat([df, desc_df], axis=1)
    return result


def calculate_similarity(smiles1: str, smiles2: str, fingerprint_type: str = 'morgan') -> Optional[float]:
    """
    Calculate Tanimoto similarity between two compounds.
    
    Args:
        smiles1: First SMILES string
        smiles2: Second SMILES string
        fingerprint_type: Type of fingerprint ('morgan', 'maccs')
    
    Returns:
        Tanimoto similarity coefficient (0-1)
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    mol1 = smiles_to_mol(smiles1)
    mol2 = smiles_to_mol(smiles2)
    
    if mol1 is None or mol2 is None:
        return None
    
    if fingerprint_type == 'morgan':
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2, nBits=1024)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2, nBits=1024)
    elif fingerprint_type == 'maccs':
        from rdkit.Chem import MACCSkeys
        fp1 = MACCSkeys.GenMACCSKeys(mol1)
        fp2 = MACCSkeys.GenMACCSKeys(mol2)
    else:
        raise ValueError(f"Unknown fingerprint type: {fingerprint_type}")
    
    return DataStructs.TanimotoSimilarity(fp1, fp2)


def find_similar_compounds(
    query_smiles: str,
    candidates_df: pd.DataFrame,
    smiles_col: str = 'smiles',
    threshold: float = 0.7,
    top_n: int = 10
) -> pd.DataFrame:
    """
    Find compounds similar to a query compound.
    
    Args:
        query_smiles: Query SMILES string
        candidates_df: DataFrame of candidate compounds
        smiles_col: Name of SMILES column
        threshold: Minimum similarity threshold
        top_n: Return top N most similar compounds
    
    Returns:
        DataFrame with similar compounds and similarity scores
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    query_mol = smiles_to_mol(query_smiles)
    if query_mol is None:
        raise ValueError("Invalid query SMILES")
    
    query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024)
    
    similarities = []
    for idx, row in candidates_df.iterrows():
        mol = smiles_to_mol(row[smiles_col])
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            sim = DataStructs.TanimotoSimilarity(query_fp, fp)
            similarities.append((idx, sim))
    
    # Filter by threshold and sort
    similarities = [(idx, sim) for idx, sim in similarities if sim >= threshold]
    similarities.sort(key=lambda x: x[1], reverse=True)
    
    # Get top N
    top_indices = [idx for idx, _ in similarities[:top_n]]
    
    result = candidates_df.loc[top_indices].copy()
    result['similarity'] = [sim for _, sim in similarities[:top_n]]
    
    return result


def generate_molecule_image(
    smiles: str,
    output_path: str,
    size: tuple = (300, 300)
) -> bool:
    """
    Generate 2D image of a molecule.
    
    Args:
        smiles: SMILES string
        output_path: Path to save image
        size: Image size (width, height)
    
    Returns:
        True if successful
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    mol = smiles_to_mol(smiles)
    if mol is None:
        return False
    
    img = Draw.MolToImage(mol, size=size)
    img.save(output_path)
    return True


def validate_structures(df: pd.DataFrame, smiles_col: str = 'smiles') -> pd.DataFrame:
    """
    Validate SMILES structures using RDKit.
    
    Args:
        df: DataFrame with SMILES column
        smiles_col: Name of SMILES column
    
    Returns:
        DataFrame with validation results
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    df = df.copy()
    df['is_valid'] = df[smiles_col].apply(lambda s: smiles_to_mol(s) is not None)
    
    return df


def cluster_compounds(
    df: pd.DataFrame,
    smiles_col: str = 'smiles',
    n_clusters: int = 5
) -> pd.DataFrame:
    """
    Cluster compounds by structural similarity.
    
    Args:
        df: DataFrame with SMILES
        smiles_col: Name of SMILES column
        n_clusters: Number of clusters
    
    Returns:
        DataFrame with cluster assignments
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit not available")
    
    from rdkit.ML.Cluster import Butina
    
    # Generate fingerprints
    fps = []
    valid_indices = []
    
    for idx, smiles in enumerate(df[smiles_col]):
        mol = smiles_to_mol(smiles)
        if mol is not None:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(fp)
            valid_indices.append(idx)
    
    # Calculate distance matrix
    dists = []
    for i in range(len(fps)):
        for j in range(i + 1, len(fps)):
            dists.append(1 - DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    
    # Cluster
    cs = Butina.ClusterData(dists, len(fps), 0.7, isDistData=True)
    
    # Assign clusters
    df = df.copy()
    df['cluster'] = -1
    for cluster_id, cluster in enumerate(cs):
        for member_idx in cluster:
            df.loc[valid_indices[member_idx], 'cluster'] = cluster_id
    
    return df


if __name__ == "__main__":
    print("ZINC Database RDKit Integration")
    print("=" * 40)
    
    if not RDKIT_AVAILABLE:
        print("\nRDKit not available. Install with:")
        print("  conda install -c conda-forge rdkit")
        print("\nExample usage (when RDKit is available):")
        print("""
        # Calculate descriptors
        desc = calculate_descriptors('CC(C)Cc1ccc(cc1)C(C)C(=O)O')
        
        # Add to DataFrame
        df_with_desc = add_rdkit_descriptors(df)
        
        # Calculate similarity
        sim = calculate_similarity('c1ccccc1', 'c1ccc(cc1)O')
        
        # Find similar compounds
        similar = find_similar_compounds(query, candidates_df, threshold=0.7)
        """)
    else:
        # Example with RDKit available
        print("\nRDKit is available!")
        
        # Test SMILES
        test_smiles = 'CC(C)Cc1ccc(cc1)C(C)C(=O)O'  # Ibuprofen
        
        print(f"\nCalculating descriptors for: {test_smiles}")
        desc = calculate_descriptors(test_smiles)
        for key, value in desc.items():
            print(f"  {key}: {value}")
        
        # Test similarity
        smiles1 = 'c1ccccc1'  # Benzene
        smiles2 = 'c1ccc(cc1)O'  # Phenol
        sim = calculate_similarity(smiles1, smiles2)
        print(f"\nSimilarity between benzene and phenol: {sim:.3f}")
