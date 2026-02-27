#!/usr/bin/env python3
"""
ZINC Database Search Utilities
Search ZINC22 by ZINC ID, SMILES, and supplier codes.
"""

import requests
import pandas as pd
from io import StringIO
from typing import List, Optional, Union
import time


BASE_URL = "https://cartblanche22.docking.org"


def search_by_zinc_id(
    zinc_id: Union[str, List[str]],
    output_fields: str = "zinc_id,smiles,supplier_code,catalogs,tranche"
) -> pd.DataFrame:
    """
    Search ZINC database by ZINC ID(s).
    
    Args:
        zinc_id: Single ZINC ID or list of ZINC IDs
        output_fields: Comma-separated list of fields to return
    
    Returns:
        DataFrame with search results
    """
    if isinstance(zinc_id, list):
        zinc_id = ",".join(zinc_id)
    
    url = f"{BASE_URL}/substances.txt"
    params = {
        "zinc_id": zinc_id,
        "output_fields": output_fields
    }
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    return pd.read_csv(StringIO(response.text), sep='\t')


def search_by_smiles(
    smiles: str,
    dist: int = 0,
    adist: int = 0,
    output_fields: str = "zinc_id,smiles,supplier_code,catalogs,tranche"
) -> pd.DataFrame:
    """
    Search ZINC database by SMILES string.
    
    Args:
        smiles: SMILES string (URL-encoded if necessary)
        dist: Tanimoto distance threshold (0=exact match, higher=more similar)
        adist: Alternative distance parameter for broader searches
        output_fields: Comma-separated list of fields to return
    
    Returns:
        DataFrame with search results
    """
    url = f"{BASE_URL}/smiles.txt"
    params = {
        "smiles": smiles,
        "dist": dist,
        "adist": adist,
        "output_fields": output_fields
    }
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    return pd.read_csv(StringIO(response.text), sep='\t')


def search_by_supplier_code(
    supplier_code: Union[str, List[str]],
    output_fields: str = "zinc_id,smiles,supplier_code,catalogs,tranche"
) -> pd.DataFrame:
    """
    Search ZINC database by supplier catalog code(s).
    
    Args:
        supplier_code: Single supplier code or list of codes
        output_fields: Comma-separated list of fields to return
    
    Returns:
        DataFrame with search results
    """
    if isinstance(supplier_code, list):
        supplier_code = ",".join(supplier_code)
    
    url = f"{BASE_URL}/catitems.txt"
    params = {
        "catitem_id": supplier_code,
        "output_fields": output_fields
    }
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    return pd.read_csv(StringIO(response.text), sep='\t')


def get_random_compounds(
    count: int = 100,
    subset: Optional[str] = None,
    output_fields: str = "zinc_id,smiles,tranche"
) -> pd.DataFrame:
    """
    Get random compounds from ZINC database.
    
    Args:
        count: Number of random compounds to retrieve
        subset: Filter by subset ('lead-like', 'drug-like', 'fragment')
        output_fields: Comma-separated list of fields to return
    
    Returns:
        DataFrame with random compounds
    """
    url = f"{BASE_URL}/substance/random.txt"
    params = {
        "count": count,
        "output_fields": output_fields
    }
    
    if subset:
        params["subset"] = subset
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    return pd.read_csv(StringIO(response.text), sep='\t')


def batch_search_zinc_ids(
    zinc_ids: List[str],
    batch_size: int = 100,
    delay: float = 0.5,
    output_fields: str = "zinc_id,smiles,supplier_code,catalogs,tranche"
) -> pd.DataFrame:
    """
    Batch search multiple ZINC IDs with rate limiting.
    
    Args:
        zinc_ids: List of ZINC IDs
        batch_size: Number of IDs per request
        delay: Delay between requests in seconds
        output_fields: Fields to return
    
    Returns:
        Combined DataFrame with all results
    """
    results = []
    
    for i in range(0, len(zinc_ids), batch_size):
        batch = zinc_ids[i:i + batch_size]
        df = search_by_zinc_id(batch, output_fields)
        results.append(df)
        
        if i + batch_size < len(zinc_ids):
            time.sleep(delay)
    
    return pd.concat(results, ignore_index=True)


def find_analogs(
    smiles: str,
    distance: int = 5,
    output_fields: str = "zinc_id,smiles,supplier_code,catalogs,tranche"
) -> pd.DataFrame:
    """
    Find structural analogs of a compound using similarity search.
    
    Args:
        smiles: Query SMILES string
        distance: Tanimoto distance (1-3 for close analogs, 5-10 for diverse)
        output_fields: Fields to return
    
    Returns:
        DataFrame with analog compounds
    """
    return search_by_smiles(smiles, dist=distance, output_fields=output_fields)


if __name__ == "__main__":
    print("ZINC Database Search Utilities")
    print("=" * 40)
    
    # Example 1: Search by ZINC ID
    print("\n1. Searching by ZINC ID...")
    try:
        result = search_by_zinc_id("ZINC000000000001")
        print(f"Found {len(result)} compound(s)")
        print(result[['zinc_id', 'smiles']].to_string())
    except Exception as e:
        print(f"Error: {e}")
    
    # Example 2: Search by SMILES (exact match)
    print("\n2. Searching by SMILES (benzene)...")
    try:
        result = search_by_smiles("c1ccccc1", dist=0)
        print(f"Found {len(result)} compound(s)")
        print(result[['zinc_id', 'smiles']].head().to_string())
    except Exception as e:
        print(f"Error: {e}")
    
    # Example 3: Get random compounds
    print("\n3. Getting random lead-like compounds...")
    try:
        result = get_random_compounds(count=5, subset="lead-like")
        print(f"Retrieved {len(result)} random compounds")
        print(result[['zinc_id', 'smiles']].to_string())
    except Exception as e:
        print(f"Error: {e}")
