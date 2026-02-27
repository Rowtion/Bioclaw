#!/usr/bin/env python3
"""
Vaex Tools
High-performance DataFrame library for large datasets.
"""

from typing import List, Dict
import json


def load_large_dataset(filepath: str) -> object:
    """Load large dataset with Vaex."""
    try:
        import vaex
        if filepath.endswith('.csv'):
            return vaex.read_csv(filepath, convert=True)
        elif filepath.endswith('.hdf5') or filepath.endswith('.hdf'):
            return vaex.open(filepath)
        elif filepath.endswith('.arrow') or filepath.endswith('.feather'):
            return vaex.read_arrow(filepath)
        else:
            print(f"Unsupported format: {filepath}")
            return None
    except ImportError:
        print("Error: vaex not installed. Run: pip install vaex")
        return None


def describe_dataset(df) -> Dict:
    """Get dataset statistics."""
    if df is None:
        return {}
    
    return {
        "shape": df.shape,
        "columns": df.column_names,
        "dtypes": {col: str(dtype) for col, dtype in zip(df.column_names, df.dtypes)},
        "virtual_columns": list(df.virtual_columns.keys()) if hasattr(df, 'virtual_columns') else []
    }


def fast_groupby(df, column: str, agg_column: str = None) -> Dict:
    """Perform fast groupby aggregation."""
    try:
        if agg_column:
            result = df.groupby(column, agg={agg_column: 'mean'})
        else:
            result = df.groupby(column)
        return result.to_dict()
    except Exception as e:
        print(f"Error: {e}")
        return {}


def filter_data(df, expression: str):
    """Filter data using expression."""
    try:
        return df[expression]
    except Exception as e:
        print(f"Error: {e}")
        return None


def create_virtual_column(df, name: str, expression: str):
    """Create virtual column."""
    try:
        df[name] = df[expression]
        print(f"Created virtual column: {name}")
        return df
    except Exception as e:
        print(f"Error: {e}")
        return df


def export_data(df, output_path: str):
    """Export to various formats."""
    try:
        if output_path.endswith('.csv'):
            df.export_csv(output_path)
        elif output_path.endswith('.hdf5'):
            df.export_hdf5(output_path)
        elif output_path.endswith('.arrow'):
            df.export_arrow(output_path)
        else:
            print(f"Unsupported output format: {output_path}")
            return False
        print(f"Exported to: {output_path}")
        return True
    except Exception as e:
        print(f"Error: {e}")
        return False


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Vaex Tools")
    parser.add_argument("command", choices=["load", "describe", "groupby", "filter", "export"])
    parser.add_argument("--file", help="Input file")
    parser.add_argument("--column", help="Column name")
    parser.add_argument("--agg-column", help="Aggregation column")
    parser.add_argument("--expression", help="Filter expression")
    parser.add_argument("--virtual-name", help="Virtual column name")
    parser.add_argument("--virtual-expr", help="Virtual column expression")
    parser.add_argument("--output", help="Output file")
    
    args = parser.parse_args()
    
    if args.command == "load":
        if not args.file:
            print("Error: --file required")
            return
        df = load_large_dataset(args.file)
        if df is not None:
            print(f"Loaded: {df.shape}")
    
    elif args.command == "describe":
        if not args.file:
            print("Error: --file required")
            return
        df = load_large_dataset(args.file)
        info = describe_dataset(df)
        print(json.dumps(info, indent=2))
    
    elif args.command == "filter":
        if not args.file or not args.expression:
            print("Error: --file and --expression required")
            return
        df = load_large_dataset(args.file)
        filtered = filter_data(df, args.expression)
        if filtered is not None:
            print(f"Filtered: {filtered.shape}")
    
    elif args.command == "export":
        if not args.file or not args.output:
            print("Error: --file and --output required")
            return
        df = load_large_dataset(args.file)
        export_data(df, args.output)


if __name__ == "__main__":
    main()
