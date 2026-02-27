#!/usr/bin/env python3
"""
Polars Tools - High-performance DataFrame operations
Fast data processing with lazy evaluation and parallel execution
"""

import argparse
import json
import time
from typing import List, Dict, Optional
from pathlib import Path


def read_csv_lazy(file_path: str, n_rows: Optional[int] = None) -> "pl.LazyFrame":
    """Read CSV file lazily for efficient processing."""
    try:
        import polars as pl
        return pl.scan_csv(file_path, n_rows=n_rows)
    except ImportError:
        print("Error: Polars not installed. Run: uv pip install polars")
        return None


def benchmark_csv(file_path: str) -> Dict:
    """Benchmark CSV reading performance."""
    try:
        import polars as pl
        import pandas as pd
    except ImportError:
        print("Error: Required packages not installed")
        return {}
    
    results = {}
    
    # Polars lazy
    start = time.time()
    df_lazy = pl.scan_csv(file_path)
    df_lazy.collect()
    results["polars_lazy"] = time.time() - start
    
    # Polars eager
    start = time.time()
    pl.read_csv(file_path)
    results["polars_eager"] = time.time() - start
    
    # Pandas
    start = time.time()
    pd.read_csv(file_path)
    results["pandas"] = time.time() - start
    
    return results


def profile_dataframe(file_path: str) -> Dict:
    """Generate profile report for dataframe."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return {}
    
    df = pl.read_csv(file_path)
    
    profile = {
        "shape": df.shape,
        "columns": df.columns,
        "dtypes": {col: str(dtype) for col, dtype in zip(df.columns, df.dtypes)},
        "null_counts": df.null_count().to_dicts()[0],
        "numeric_summary": {}
    }
    
    # Summary for numeric columns
    numeric_cols = [c for c in df.columns if df[c].dtype in 
                   [pl.Float32, pl.Float64, pl.Int32, pl.Int64]]
    
    if numeric_cols:
        summary = df[numeric_cols].describe()
        profile["numeric_summary"] = summary.to_dicts()
    
    return profile


def filter_and_transform(input_path: str, output_path: str, 
                        filter_expr: str, columns: List[str]):
    """Apply filter and column selection."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return
    
    lf = pl.scan_csv(input_path)
    
    # Apply filter
    if filter_expr:
        lf = lf.filter(eval(f"pl.{filter_expr}"))
    
    # Select columns
    if columns:
        lf = lf.select(columns)
    
    # Collect and save
    df = lf.collect()
    df.write_csv(output_path)
    
    print(f"Processed {len(df)} rows, saved to {output_path}")


def aggregate_by_group(input_path: str, output_path: str,
                      group_cols: List[str], agg_exprs: Dict[str, str]):
    """Perform groupby aggregation."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return
    
    df = pl.read_csv(input_path)
    
    # Build aggregation expressions
    agg_list = []
    for col, func in agg_exprs.items():
        if func == "sum":
            agg_list.append(pl.col(col).sum().alias(f"{col}_sum"))
        elif func == "mean":
            agg_list.append(pl.col(col).mean().alias(f"{col}_mean"))
        elif func == "count":
            agg_list.append(pl.col(col).count().alias(f"{col}_count"))
        elif func == "min":
            agg_list.append(pl.col(col).min().alias(f"{col}_min"))
        elif func == "max":
            agg_list.append(pl.col(col).max().alias(f"{col}_max"))
    
    result = df.group_by(group_cols).agg(agg_list)
    result.write_csv(output_path)
    
    print(f"Aggregated data saved to {output_path}")
    print(f"Groups: {len(result)}")


def join_dataframes(left_path: str, right_path: str, output_path: str,
                   left_on: str, right_on: Optional[str] = None,
                   how: str = "inner"):
    """Join two dataframes."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return
    
    left = pl.read_csv(left_path)
    right = pl.read_csv(right_path)
    
    right_on = right_on or left_on
    
    result = left.join(right, left_on=left_on, right_on=right_on, how=how)
    result.write_csv(output_path)
    
    print(f"Joined data saved to {output_path}")
    print(f"Result rows: {len(result)}")


def convert_to_parquet(csv_path: str, parquet_path: str):
    """Convert CSV to Parquet for faster access."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return
    
    start = time.time()
    df = pl.read_csv(csv_path)
    df.write_parquet(parquet_path)
    elapsed = time.time() - start
    
    csv_size = Path(csv_path).stat().st_size / (1024 * 1024)
    parquet_size = Path(parquet_path).stat().st_size / (1024 * 1024)
    
    print(f"Converted in {elapsed:.2f}s")
    print(f"CSV size: {csv_size:.2f} MB")
    print(f"Parquet size: {parquet_size:.2f} MB")
    print(f"Compression ratio: {csv_size/parquet_size:.2f}x")


def pivot_table(input_path: str, output_path: str,
               index_col: str, columns_col: str, values_col: str):
    """Create pivot table."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return
    
    df = pl.read_csv(input_path)
    
    pivoted = df.pivot(
        values=values_col,
        index=index_col,
        columns=columns_col
    )
    
    pivoted.write_csv(output_path)
    print(f"Pivot table saved to {output_path}")


def generate_sample_data(output_path: str, n_rows: int = 100000):
    """Generate large sample dataset for testing."""
    try:
        import polars as pl
        import numpy as np
    except ImportError:
        print("Error: Required packages not installed")
        return
    
    np.random.seed(42)
    
    df = pl.DataFrame({
        "id": range(n_rows),
        "category": np.random.choice(["A", "B", "C", "D"], n_rows),
        "region": np.random.choice(["North", "South", "East", "West"], n_rows),
        "value1": np.random.randn(n_rows),
        "value2": np.random.exponential(10, n_rows),
        "value3": np.random.randint(1, 100, n_rows),
        "date": pl.datetime_range(
            start="2024-01-01",
            end="2024-12-31",
            eager=True
        ).sample(n_rows, with_replacement=True)
    })
    
    df.write_csv(output_path)
    print(f"Generated {n_rows} rows, saved to {output_path}")


def query_sql(file_path: str, query: str) -> "pl.DataFrame":
    """Execute SQL query on dataframe."""
    try:
        import polars as pl
    except ImportError:
        print("Error: Polars not installed")
        return None
    
    df = pl.read_csv(file_path)
    
    # Register table and execute query
    result = pl.SQLContext(frame=df).execute(query)
    
    return result


def main():
    parser = argparse.ArgumentParser(description="Polars DataFrame Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Benchmark command
    bench_parser = subparsers.add_parser("benchmark", help="Benchmark CSV reading")
    bench_parser.add_argument("file", help="CSV file to benchmark")
    
    # Profile command
    profile_parser = subparsers.add_parser("profile", help="Profile dataframe")
    profile_parser.add_argument("file", help="CSV file to profile")
    profile_parser.add_argument("-o", "--output", help="Output JSON file")
    
    # Filter command
    filter_parser = subparsers.add_parser("filter", help="Filter and transform data")
    filter_parser.add_argument("input", help="Input CSV")
    filter_parser.add_argument("output", help="Output CSV")
    filter_parser.add_argument("--expr", help="Filter expression (e.g., 'col(\"age\") > 25')")
    filter_parser.add_argument("--cols", help="Comma-separated columns to keep")
    
    # Aggregate command
    agg_parser = subparsers.add_parser("aggregate", help="Groupby aggregation")
    agg_parser.add_argument("input", help="Input CSV")
    agg_parser.add_argument("output", help="Output CSV")
    agg_parser.add_argument("--group", required=True, help="Comma-separated group columns")
    agg_parser.add_argument("--agg", required=True, help="Aggregations: col:func,col2:func2")
    
    # Join command
    join_parser = subparsers.add_parser("join", help="Join dataframes")
    join_parser.add_argument("left", help="Left CSV")
    join_parser.add_argument("right", help="Right CSV")
    join_parser.add_argument("output", help="Output CSV")
    join_parser.add_argument("--left-on", required=True, help="Left join key")
    join_parser.add_argument("--right-on", help="Right join key")
    join_parser.add_argument("--how", default="inner", choices=["inner", "left", "right", "outer"])
    
    # Convert command
    convert_parser = subparsers.add_parser("convert", help="Convert CSV to Parquet")
    convert_parser.add_argument("input", help="Input CSV")
    convert_parser.add_argument("output", help="Output Parquet")
    
    # Pivot command
    pivot_parser = subparsers.add_parser("pivot", help="Create pivot table")
    pivot_parser.add_argument("input", help="Input CSV")
    pivot_parser.add_argument("output", help="Output CSV")
    pivot_parser.add_argument("--index", required=True, help="Index column")
    pivot_parser.add_argument("--columns", required=True, help="Columns column")
    pivot_parser.add_argument("--values", required=True, help="Values column")
    
    # Generate sample command
    gen_parser = subparsers.add_parser("generate", help="Generate sample data")
    gen_parser.add_argument("-o", "--output", default="sample.csv")
    gen_parser.add_argument("-n", type=int, default=100000, help="Number of rows")
    
    # SQL command
    sql_parser = subparsers.add_parser("sql", help="Execute SQL query")
    sql_parser.add_argument("file", help="CSV file")
    sql_parser.add_argument("query", help="SQL query")
    sql_parser.add_argument("-o", "--output", help="Output CSV file")
    
    args = parser.parse_args()
    
    if args.command == "benchmark":
        results = benchmark_csv(args.file)
        print(json.dumps(results, indent=2))
    
    elif args.command == "profile":
        profile = profile_dataframe(args.file)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(profile, f, indent=2, default=str)
            print(f"Profile saved to {args.output}")
        else:
            print(json.dumps(profile, indent=2, default=str))
    
    elif args.command == "filter":
        cols = args.cols.split(",") if args.cols else None
        filter_and_transform(args.input, args.output, args.expr, cols)
    
    elif args.command == "aggregate":
        group_cols = args.group.split(",")
        agg_dict = {}
        for agg in args.agg.split(","):
            col, func = agg.split(":")
            agg_dict[col] = func
        aggregate_by_group(args.input, args.output, group_cols, agg_dict)
    
    elif args.command == "join":
        join_dataframes(args.left, args.right, args.output, 
                       args.left_on, args.right_on, args.how)
    
    elif args.command == "convert":
        convert_to_parquet(args.input, args.output)
    
    elif args.command == "pivot":
        pivot_table(args.input, args.output, args.index, args.columns, args.values)
    
    elif args.command == "generate":
        generate_sample_data(args.output, args.n)
    
    elif args.command == "sql":
        result = query_sql(args.file, args.query)
        if result is not None:
            print(result)
            if args.output:
                result.write_csv(args.output)
                print(f"Results saved to {args.output}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
