"""
GEO Database Example: Query and Download
Search and download gene expression data from NCBI GEO.
"""
import GEOparse
import pandas as pd


def download_geo_series(geo_id: str, destdir: str = "./geo_data"):
    """
    Download and parse a GEO Series.
    
    Args:
        geo_id: GEO Series ID (e.g., GSE123456)
        destdir: Destination directory
        
    Returns:
        GEOparse object
    """
    print(f"Downloading {geo_id}...")
    gse = GEOparse.get_GEO(geo=geo_id, destdir=destdir)
    
    print(f"\nSeries: {gse.metadata['title'][0]}")
    print(f"Organism: {gse.metadata['organism'][0]}")
    print(f"Samples: {len(gse.gsms)}")
    
    return gse


def extract_expression_data(gse):
    """
    Extract expression matrix from GEO Series.
    
    Args:
        gse: GEOparse series object
        
    Returns:
        Expression DataFrame
    """
    # Extract expression matrix
    if hasattr(gse, 'pivot_samples'):
        expression_df = gse.pivot_samples('VALUE')
        print(f"\nExpression matrix: {expression_df.shape}")
        return expression_df
    else:
        # Build from individual samples
        expression_data = {}
        for gsm_name, gsm in gse.gsms.items():
            if hasattr(gsm, 'table'):
                expression_data[gsm_name] = gsm.table['VALUE']
        
        expression_df = pd.DataFrame(expression_data)
        print(f"Expression matrix: {expression_df.shape}")
        return expression_df


def get_sample_metadata(gse):
    """
    Extract sample metadata.
    
    Args:
        gse: GEOparse series object
        
    Returns:
        Metadata DataFrame
    """
    metadata = []
    
    for gsm_name, gsm in gse.gsms.items():
        meta = {
            'sample_id': gsm_name,
            'title': gsm.metadata.get('title', [''])[0],
            'source': gsm.metadata.get('source_name_ch1', [''])[0],
        }
        
        # Extract characteristics
        chars = gsm.metadata.get('characteristics_ch1', [])
        for char in chars:
            if ':' in char:
                key, value = char.split(':', 1)
                meta[key.strip()] = value.strip()
        
        metadata.append(meta)
    
    df = pd.DataFrame(metadata)
    return df


if __name__ == "__main__":
    print("GEO Database Example")
    # gse = download_geo_series('GSE123456')
    # expr = extract_expression_data(gse)
    # meta = get_sample_metadata(gse)
