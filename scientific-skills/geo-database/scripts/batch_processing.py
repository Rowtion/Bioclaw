"""
GEO Database Example: Batch Processing
Process multiple GEO datasets in batch.
"""
import GEOparse
import pandas as pd


def batch_download(gse_list: list, destdir: str = "./geo_data"):
    """
    Download multiple GEO series.
    
    Args:
        gse_list: List of GSE IDs
        destdir: Destination directory
        
    Returns:
        Summary DataFrame
    """
    results = {}
    
    for gse_id in gse_list:
        try:
            print(f"Processing {gse_id}...")
            gse = GEOparse.get_GEO(geo=gse_id, destdir=destdir)
            
            results[gse_id] = {
                'title': gse.metadata.get('title', ['N/A'])[0],
                'organism': gse.metadata.get('organism', ['N/A'])[0],
                'num_samples': len(gse.gsms),
                'platform': list(gse.gpls.keys())[0] if gse.gpls else 'N/A',
                'status': 'success'
            }
            
            # Save expression data
            if hasattr(gse, 'pivot_samples'):
                expr = gse.pivot_samples('VALUE')
                expr.to_csv(f"{destdir}/{gse_id}_expression.csv")
                results[gse_id]['num_genes'] = len(expr)
        
        except Exception as e:
            print(f"Error: {e}")
            results[gse_id] = {'status': 'error', 'error': str(e)}
    
    df = pd.DataFrame(results).T
    df.to_csv(f"{destdir}/batch_summary.csv")
    
    return df


def meta_analysis(gse_list: list, gene_of_interest: str):
    """
    Meta-analysis of a gene across multiple studies.
    
    Args:
        gse_list: List of GSE IDs
        gene_of_interest: Gene symbol
        
    Returns:
        Meta-analysis results
    """
    results = []
    
    for gse_id in gse_list:
        try:
            gse = GEOparse.get_GEO(geo=gse_id, destdir="./geo_data")
            
            # Get expression data
            if hasattr(gse, 'pivot_samples'):
                expr = gse.pivot_samples('VALUE')
                
                if gene_of_interest in expr.index:
                    values = expr.loc[gene_of_interest]
                    
                    results.append({
                        'study': gse_id,
                        'gene': gene_of_interest,
                        'mean_expression': values.mean(),
                        'std_expression': values.std(),
                        'num_samples': len(values)
                    })
        
        except Exception as e:
            print(f"Error in {gse_id}: {e}")
    
    return pd.DataFrame(results)


if __name__ == "__main__":
    print("GEO Batch Processing Example")
    # summary = batch_download(['GSE100001', 'GSE100002'])
    # meta = meta_analysis(['GSE100001', 'GSE100002'], 'TP53')
