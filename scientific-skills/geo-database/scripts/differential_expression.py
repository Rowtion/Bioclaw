"""
GEO Database Example: Differential Expression
Perform differential expression analysis on GEO data.
"""
import GEOparse
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


def load_and_filter(gse_id: str, control_keyword: str, treatment_keyword: str):
    """
    Load GEO data and filter samples by condition.
    
    Args:
        gse_id: GEO Series ID
        control_keyword: Keyword for control samples
        treatment_keyword: Keyword for treatment samples
        
    Returns:
        Expression data and sample groups
    """
    gse = GEOparse.get_GEO(geo=gse_id, destdir="./geo_data")
    
    # Get expression data
    expr_df = gse.pivot_samples('VALUE')
    
    # Classify samples
    control_samples = []
    treatment_samples = []
    
    for gsm_name, gsm in gse.gsms.items():
        title = gsm.metadata.get('title', [''])[0].lower()
        
        if control_keyword.lower() in title:
            control_samples.append(gsm_name)
        elif treatment_keyword.lower() in title:
            treatment_samples.append(gsm_name)
    
    print(f"Control samples: {len(control_samples)}")
    print(f"Treatment samples: {len(treatment_samples)}")
    
    return expr_df, control_samples, treatment_samples


def differential_expression(expr_df, control_samples, treatment_samples):
    """
    Perform differential expression analysis.
    
    Args:
        expr_df: Expression DataFrame
        control_samples: List of control sample IDs
        treatment_samples: List of treatment sample IDs
        
    Returns:
        Results DataFrame
    """
    results = []
    
    for gene in expr_df.index:
        control_expr = expr_df.loc[gene, control_samples]
        treatment_expr = expr_df.loc[gene, treatment_samples]
        
        # Calculate statistics
        fold_change = treatment_expr.mean() - control_expr.mean()
        t_stat, p_value = stats.ttest_ind(treatment_expr, control_expr)
        
        results.append({
            'gene': gene,
            'log2_fold_change': fold_change,
            'p_value': p_value,
            'control_mean': control_expr.mean(),
            'treatment_mean': treatment_expr.mean()
        })
    
    df = pd.DataFrame(results)
    
    # Multiple testing correction
    _, df['q_value'], _, _ = multipletests(df['p_value'], method='fdr_bh')
    
    # Filter significant genes
    significant = df[(df['q_value'] < 0.05) & (abs(df['log2_fold_change']) > 1)]
    
    print(f"Significant genes: {len(significant)}")
    
    return df, significant


if __name__ == "__main__":
    print("GEO Differential Expression Example")
    # expr, control, treatment = load_and_filter('GSE123456', 'control', 'treatment')
    # results, sig = differential_expression(expr, control, treatment)
