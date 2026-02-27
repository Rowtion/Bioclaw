#!/usr/bin/env python3
"""
ZINC Database Structure Utilities
Parse tranche data and analyze molecular properties.
"""

import re
import pandas as pd
from typing import Dict, Optional


def parse_tranche(tranche_str: str) -> Optional[Dict[str, float]]:
    """
    Parse ZINC tranche code to extract molecular properties.
    
    Tranche format: H##P###M###-phase
    - H##: Number of hydrogen bond donors (00-99)
    - P###: LogP × 10 (e.g., P035 = LogP 3.5)
    - M###: Molecular weight in Daltons (e.g., M400 = 400 Da)
    - phase: Reactivity classification
    
    Args:
        tranche_str: ZINC tranche string
    
    Returns:
        Dictionary with parsed properties or None if invalid
    """
    match = re.match(r'H(\d+)P(\d+)M(\d+)-(\d+)', tranche_str)
    if match:
        return {
            'h_donors': int(match.group(1)),
            'logP': int(match.group(2)) / 10.0,
            'mw': int(match.group(3)),
            'phase': int(match.group(4))
        }
    return None


def add_tranche_properties(df: pd.DataFrame, tranche_col: str = 'tranche') -> pd.DataFrame:
    """
    Add parsed tranche properties as new columns to a DataFrame.
    
    Args:
        df: DataFrame with tranche column
        tranche_col: Name of the tranche column
    
    Returns:
        DataFrame with added property columns
    """
    df = df.copy()
    
    # Parse tranche for each row
    properties = df[tranche_col].apply(parse_tranche)
    
    # Extract properties into separate columns
    df['h_donors'] = properties.apply(lambda x: x['h_donors'] if x else None)
    df['logP'] = properties.apply(lambda x: x['logP'] if x else None)
    df['mw'] = properties.apply(lambda x: x['mw'] if x else None)
    df['phase'] = properties.apply(lambda x: x['phase'] if x else None)
    
    return df


def classify_by_lipinski(df: pd.DataFrame, mw_col: str = 'mw', logp_col: str = 'logP',
                         h_donors_col: str = 'h_donors') -> pd.DataFrame:
    """
    Classify compounds by Lipinski's Rule of Five.
    
    Rules:
    - Molecular weight < 500 Da
    - LogP < 5
    - H-bond donors < 5
    - H-bond acceptors < 10 (requires additional data)
    
    Args:
        df: DataFrame with molecular properties
        mw_col: Column name for molecular weight
        logp_col: Column name for LogP
        h_donors_col: Column name for H-bond donors
    
    Returns:
        DataFrame with 'lipinski_violations' and 'drug_like' columns
    """
    df = df.copy()
    
    violations = pd.Series(0, index=df.index)
    
    if mw_col in df.columns:
        violations += (df[mw_col] >= 500).astype(int)
    if logp_col in df.columns:
        violations += (df[logp_col] >= 5).astype(int)
    if h_donors_col in df.columns:
        violations += (df[h_donors_col] >= 5).astype(int)
    
    df['lipinski_violations'] = violations
    df['drug_like'] = violations <= 1  # Pass if 0 or 1 violations
    
    return df


def classify_by_subset(df: pd.DataFrame, mw_col: str = 'mw', logp_col: str = 'logP') -> pd.DataFrame:
    """
    Classify compounds into ZINC subsets based on properties.
    
    Subsets:
    - fragment: MW < 250
    - lead-like: MW 250-350, LogP ≤ 3.5
    - drug-like: MW 350-500, follows Lipinski
    - beyond: MW > 500
    
    Args:
        df: DataFrame with molecular properties
        mw_col: Column name for molecular weight
        logp_col: Column name for LogP
    
    Returns:
        DataFrame with 'subset' column
    """
    df = df.copy()
    
    def get_subset(row):
        mw = row.get(mw_col)
        logp = row.get(logp_col)
        
        if mw is None:
            return 'unknown'
        
        if mw < 250:
            return 'fragment'
        elif 250 <= mw < 350 and (logp is None or logp <= 3.5):
            return 'lead-like'
        elif 350 <= mw <= 500:
            return 'drug-like'
        else:
            return 'beyond'
    
    df['subset'] = df.apply(get_subset, axis=1)
    return df


def filter_by_properties(
    df: pd.DataFrame,
    min_mw: Optional[float] = None,
    max_mw: Optional[float] = None,
    min_logp: Optional[float] = None,
    max_logp: Optional[float] = None,
    max_h_donors: Optional[int] = None
) -> pd.DataFrame:
    """
    Filter compounds by molecular properties.
    
    Args:
        df: DataFrame with parsed properties
        min_mw: Minimum molecular weight
        max_mw: Maximum molecular weight
        min_logp: Minimum LogP
        max_logp: Maximum LogP
        max_h_donors: Maximum H-bond donors
    
    Returns:
        Filtered DataFrame
    """
    filtered = df.copy()
    
    if min_mw is not None and 'mw' in filtered.columns:
        filtered = filtered[filtered['mw'] >= min_mw]
    if max_mw is not None and 'mw' in filtered.columns:
        filtered = filtered[filtered['mw'] <= max_mw]
    if min_logp is not None and 'logP' in filtered.columns:
        filtered = filtered[filtered['logP'] >= min_logp]
    if max_logp is not None and 'logP' in filtered.columns:
        filtered = filtered[filtered['logP'] <= max_logp]
    if max_h_donors is not None and 'h_donors' in filtered.columns:
        filtered = filtered[filtered['h_donors'] <= max_h_donors]
    
    return filtered


def get_property_summary(df: pd.DataFrame) -> pd.DataFrame:
    """
    Get summary statistics of molecular properties.
    
    Args:
        df: DataFrame with parsed properties
    
    Returns:
        DataFrame with summary statistics
    """
    property_cols = ['mw', 'logP', 'h_donors', 'phase']
    available_cols = [c for c in property_cols if c in df.columns]
    
    if not available_cols:
        return pd.DataFrame()
    
    return df[available_cols].describe()


def validate_smiles(smiles: str) -> bool:
    """
    Basic validation of SMILES string.
    
    Args:
        smiles: SMILES string to validate
    
    Returns:
        True if valid, False otherwise
    """
    if not smiles or not isinstance(smiles, str):
        return False
    
    # Basic checks
    if len(smiles) < 1:
        return False
    
    # Check for valid characters
    valid_chars = set('CNOHPScnlbrFI0123456789.@[]()-=#$+/\\%')
    if not all(c in valid_chars for c in smiles.upper()):
        return False
    
    # Check parentheses balance
    if smiles.count('(') != smiles.count(')'):
        return False
    
    # Check brackets balance
    if smiles.count('[') != smiles.count(']'):
        return False
    
    return True


if __name__ == "__main__":
    print("ZINC Database Structure Utilities")
    print("=" * 40)
    
    # Example: Parse tranche
    print("\n1. Parsing tranche string...")
    tranche = "H05P035M400-0"
    props = parse_tranche(tranche)
    print(f"Tranche: {tranche}")
    print(f"  H-bond donors: {props['h_donors']}")
    print(f"  LogP: {props['logP']}")
    print(f"  MW: {props['mw']}")
    print(f"  Phase: {props['phase']}")
    
    # Example: Create sample DataFrame with tranches
    print("\n2. Processing sample compounds...")
    sample_data = pd.DataFrame({
        'zinc_id': ['ZINC001', 'ZINC002', 'ZINC003', 'ZINC004'],
        'smiles': ['CCO', 'c1ccccc1', 'CC(C)C', 'CCCCCCCCCC'],
        'tranche': ['H02P015M150-0', 'H00P025M250-0', 'H01P020M180-0', 'H05P045M450-0']
    })
    
    # Add properties
    sample_data = add_tranche_properties(sample_data)
    print("\nWith parsed properties:")
    print(sample_data[['zinc_id', 'mw', 'logP', 'h_donors']].to_string())
    
    # Classify by Lipinski
    sample_data = classify_by_lipinski(sample_data)
    print("\nLipinski classification:")
    print(sample_data[['zinc_id', 'lipinski_violations', 'drug_like']].to_string())
    
    # Classify by subset
    sample_data = classify_by_subset(sample_data)
    print("\nSubset classification:")
    print(sample_data[['zinc_id', 'mw', 'subset']].to_string())
    
    # Property summary
    print("\nProperty summary:")
    print(get_property_summary(sample_data))
