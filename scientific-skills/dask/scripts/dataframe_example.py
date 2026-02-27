"""
Dask DataFrame Example: Parallel CSV Processing
Demonstrates reading and processing multiple CSV files with Dask DataFrames.
"""
import dask.dataframe as dd
import numpy as np


def process_large_csv(file_pattern: str):
    """
    Process multiple CSV files in parallel.
    
    Args:
        file_pattern: Glob pattern for CSV files (e.g., 'data/*.csv')
    """
    # Read multiple CSV files as a single DataFrame
    ddf = dd.read_csv(file_pattern)
    
    print(f"Total partitions: {ddf.npartitions}")
    print(f"Columns: {list(ddf.columns)}")
    print(f"Estimated rows: {len(ddf):,}")
    
    # Lazy operations - nothing computed yet
    filtered = ddf[ddf['value'] > 100]
    grouped = filtered.groupby('category').agg({
        'value': ['mean', 'sum', 'count']
    })
    
    # Trigger computation
    result = grouped.compute()
    print("\nGrouped Results:")
    print(result)
    
    return result


def parallel_etl(input_path: str, output_path: str):
    """
    ETL pipeline: Clean and transform data, then save to Parquet.
    
    Args:
        input_path: Path to input CSV files
        output_path: Path for output Parquet files
    """
    # Read raw data
    ddf = dd.read_csv(input_path)
    
    # Transformations (lazy)
    ddf = ddf[ddf['status'] == 'valid']
    ddf['amount'] = ddf['amount'].astype('float64')
    ddf = ddf.dropna(subset=['important_col'])
    ddf['processed_at'] = dd.to_datetime(ddf['timestamp'])
    
    # Save as Parquet (efficient columnar format)
    ddf.to_parquet(output_path, engine='pyarrow')
    print(f"ETL complete. Output saved to {output_path}")


if __name__ == "__main__":
    # Example usage
    # process_large_csv('data/*.csv')
    # parallel_etl('raw_data/*.csv', 'processed_data/')
    print("Dask DataFrame examples loaded. Uncomment function calls to run.")
