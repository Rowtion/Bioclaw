#!/usr/bin/env python3
"""
Plotly Tools - Interactive visualization utilities for scientific data
Creates publication-quality interactive charts and dashboards
"""

import argparse
import json
import pandas as pd
import numpy as np
from typing import List, Dict, Optional
from pathlib import Path


def create_scatter_plot(data_path: str, x_col: str, y_col: str, color_col: Optional[str] = None,
                       output: str = "scatter.html", title: str = "Scatter Plot"):
    """Create interactive scatter plot."""
    try:
        import plotly.express as px
    except ImportError:
        print("Error: Plotly not installed. Run: uv pip install plotly")
        return
    
    df = pd.read_csv(data_path)
    
    fig = px.scatter(df, x=x_col, y=y_col, color=color_col,
                     title=title, 
                     hover_data=df.columns.tolist())
    
    fig.update_layout(
        template="plotly_white",
        width=900,
        height=600
    )
    
    fig.write_html(output)
    print(f"Plot saved to {output}")


def create_line_plot(data_path: str, x_col: str, y_col: str, 
                    group_col: Optional[str] = None,
                    output: str = "line.html", title: str = "Line Plot"):
    """Create interactive line plot."""
    try:
        import plotly.express as px
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    df = pd.read_csv(data_path)
    
    fig = px.line(df, x=x_col, y=y_col, color=group_col,
                  title=title, markers=True)
    
    fig.update_layout(
        template="plotly_white",
        width=900,
        height=600,
        hovermode="x unified"
    )
    
    fig.write_html(output)
    print(f"Plot saved to {output}")


def create_heatmap(data_path: str, output: str = "heatmap.html", 
                  title: str = "Heatmap", z_col: Optional[str] = None):
    """Create interactive heatmap."""
    try:
        import plotly.express as px
        import plotly.graph_objects as go
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    df = pd.read_csv(data_path)
    
    # If 3 columns, assume x, y, z for heatmap
    if len(df.columns) == 3 and z_col:
        pivot_df = df.pivot(index=df.columns[0], columns=df.columns[1], values=z_col)
        fig = px.imshow(pivot_df, title=title, aspect="auto")
    else:
        # Correlation heatmap
        corr = df.select_dtypes(include=[np.number]).corr()
        fig = px.imshow(corr, title=title, aspect="auto",
                       color_continuous_scale="RdBu_r", zmid=0)
        fig.update_traces(text=corr.round(2), texttemplate="%{text}")
    
    fig.update_layout(width=800, height=700)
    fig.write_html(output)
    print(f"Heatmap saved to {output}")


def create_bar_chart(data_path: str, x_col: str, y_col: str,
                    output: str = "bar.html", title: str = "Bar Chart"):
    """Create interactive bar chart."""
    try:
        import plotly.express as px
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    df = pd.read_csv(data_path)
    
    fig = px.bar(df, x=x_col, y=y_col, title=title,
                 text_auto=True)
    
    fig.update_layout(
        template="plotly_white",
        width=900,
        height=600
    )
    
    fig.write_html(output)
    print(f"Bar chart saved to {output}")


def create_histogram(data_path: str, col: str, bins: int = 30,
                    output: str = "histogram.html", title: str = "Histogram"):
    """Create interactive histogram."""
    try:
        import plotly.express as px
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    df = pd.read_csv(data_path)
    
    fig = px.histogram(df, x=col, nbins=bins, title=title,
                       marginal="box", opacity=0.7)
    
    fig.update_layout(
        template="plotly_white",
        width=900,
        height=600,
        bargap=0.1
    )
    
    fig.write_html(output)
    print(f"Histogram saved to {output}")


def create_3d_surface(output: str = "surface.html", title: str = "3D Surface"):
    """Create 3D surface plot example."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    # Create sample surface data
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.sqrt(X**2 + Y**2))
    
    fig = go.Figure(data=[go.Surface(x=x, y=y, z=Z)])
    
    fig.update_layout(
        title=title,
        width=900,
        height=700,
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z"
        )
    )
    
    fig.write_html(output)
    print(f"3D surface saved to {output}")


def create_subplots(data_path: str, columns: List[str], 
                   output: str = "subplots.html", title: str = "Dashboard"):
    """Create subplot dashboard."""
    try:
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    df = pd.read_csv(data_path)
    
    n_cols = min(3, len(columns))
    n_rows = (len(columns) + n_cols - 1) // n_cols
    
    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=columns,
        specs=[[{"type": "scatter"} for _ in range(n_cols)] for _ in range(n_rows)]
    )
    
    for i, col in enumerate(columns):
        row = i // n_cols + 1
        col_idx = i % n_cols + 1
        
        fig.add_trace(
            go.Histogram(x=df[col], name=col),
            row=row, col=col_idx
        )
    
    fig.update_layout(height=300 * n_rows, width=1000, title_text=title)
    fig.write_html(output)
    print(f"Dashboard saved to {output}")


def export_static(input_html: str, output_image: str, width: int = 1200, height: int = 800):
    """Export HTML plot to static image."""
    try:
        import plotly.io as pio
    except ImportError:
        print("Error: kaleido required for image export. Run: uv pip install kaleido")
        return
    
    try:
        from PIL import Image
        img_bytes = pio.to_image(input_html, format="png", width=width, height=height)
        with open(output_image, 'wb') as f:
            f.write(img_bytes)
        print(f"Image exported to {output_image}")
    except Exception as e:
        print(f"Error exporting image: {e}")


def create_timeseries(data_path: str, date_col: str, value_col: str,
                     output: str = "timeseries.html", title: str = "Time Series"):
    """Create interactive time series plot."""
    try:
        import plotly.express as px
    except ImportError:
        print("Error: Plotly not installed")
        return
    
    df = pd.read_csv(data_path)
    df[date_col] = pd.to_datetime(df[date_col])
    
    fig = px.line(df, x=date_col, y=value_col, title=title)
    
    fig.update_xaxes(rangeslider_visible=True)
    fig.update_layout(
        template="plotly_white",
        width=1000,
        height=600
    )
    
    fig.write_html(output)
    print(f"Time series saved to {output}")


def generate_sample_data(output: str = "sample_data.csv", n: int = 100):
    """Generate sample data for testing."""
    np.random.seed(42)
    
    data = {
        "x": np.random.randn(n),
        "y": np.random.randn(n) * 2 + 5,
        "category": np.random.choice(["A", "B", "C"], n),
        "value": np.random.exponential(10, n),
        "date": pd.date_range(start="2024-01-01", periods=n, freq="D")
    }
    
    df = pd.DataFrame(data)
    df.to_csv(output, index=False)
    print(f"Sample data saved to {output}")
    return output


def main():
    parser = argparse.ArgumentParser(description="Plotly Visualization Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Scatter command
    scatter_parser = subparsers.add_parser("scatter", help="Create scatter plot")
    scatter_parser.add_argument("data", help="CSV data file")
    scatter_parser.add_argument("-x", required=True, help="X column")
    scatter_parser.add_argument("-y", required=True, help="Y column")
    scatter_parser.add_argument("-c", "--color", help="Color column")
    scatter_parser.add_argument("-o", "--output", default="scatter.html")
    scatter_parser.add_argument("-t", "--title", default="Scatter Plot")
    
    # Line command
    line_parser = subparsers.add_parser("line", help="Create line plot")
    line_parser.add_argument("data", help="CSV data file")
    line_parser.add_argument("-x", required=True, help="X column")
    line_parser.add_argument("-y", required=True, help="Y column")
    line_parser.add_argument("-g", "--group", help="Group column")
    line_parser.add_argument("-o", "--output", default="line.html")
    line_parser.add_argument("-t", "--title", default="Line Plot")
    
    # Heatmap command
    heatmap_parser = subparsers.add_parser("heatmap", help="Create heatmap")
    heatmap_parser.add_argument("data", help="CSV data file")
    heatmap_parser.add_argument("-z", help="Value column for pivot")
    heatmap_parser.add_argument("-o", "--output", default="heatmap.html")
    heatmap_parser.add_argument("-t", "--title", default="Heatmap")
    
    # Bar command
    bar_parser = subparsers.add_parser("bar", help="Create bar chart")
    bar_parser.add_argument("data", help="CSV data file")
    bar_parser.add_argument("-x", required=True, help="X column")
    bar_parser.add_argument("-y", required=True, help="Y column")
    bar_parser.add_argument("-o", "--output", default="bar.html")
    bar_parser.add_argument("-t", "--title", default="Bar Chart")
    
    # Histogram command
    hist_parser = subparsers.add_parser("histogram", help="Create histogram")
    hist_parser.add_argument("data", help="CSV data file")
    hist_parser.add_argument("-c", "--column", required=True, help="Column name")
    hist_parser.add_argument("-b", "--bins", type=int, default=30)
    hist_parser.add_argument("-o", "--output", default="histogram.html")
    hist_parser.add_argument("-t", "--title", default="Histogram")
    
    # 3D surface command
    surface_parser = subparsers.add_parser("surface", help="Create 3D surface")
    surface_parser.add_argument("-o", "--output", default="surface.html")
    surface_parser.add_argument("-t", "--title", default="3D Surface")
    
    # Timeseries command
    ts_parser = subparsers.add_parser("timeseries", help="Create time series")
    ts_parser.add_argument("data", help="CSV data file")
    ts_parser.add_argument("-d", "--date", required=True, help="Date column")
    ts_parser.add_argument("-v", "--value", required=True, help="Value column")
    ts_parser.add_argument("-o", "--output", default="timeseries.html")
    ts_parser.add_argument("-t", "--title", default="Time Series")
    
    # Dashboard command
    dash_parser = subparsers.add_parser("dashboard", help="Create dashboard with subplots")
    dash_parser.add_argument("data", help="CSV data file")
    dash_parser.add_argument("-c", "--columns", required=True, help="Comma-separated column names")
    dash_parser.add_argument("-o", "--output", default="dashboard.html")
    dash_parser.add_argument("-t", "--title", default="Dashboard")
    
    # Generate sample data command
    sample_parser = subparsers.add_parser("generate-sample", help="Generate sample data")
    sample_parser.add_argument("-o", "--output", default="sample_data.csv")
    sample_parser.add_argument("-n", type=int, default=100, help="Number of rows")
    
    args = parser.parse_args()
    
    if args.command == "scatter":
        create_scatter_plot(args.data, args.x, args.y, args.color, args.output, args.title)
    elif args.command == "line":
        create_line_plot(args.data, args.x, args.y, args.group, args.output, args.title)
    elif args.command == "heatmap":
        create_heatmap(args.data, args.output, args.title, args.z)
    elif args.command == "bar":
        create_bar_chart(args.data, args.x, args.y, args.output, args.title)
    elif args.command == "histogram":
        create_histogram(args.data, args.column, args.bins, args.output, args.title)
    elif args.command == "surface":
        create_3d_surface(args.output, args.title)
    elif args.command == "timeseries":
        create_timeseries(args.data, args.date, args.value, args.output, args.title)
    elif args.command == "dashboard":
        cols = args.columns.split(",")
        create_subplots(args.data, cols, args.output, args.title)
    elif args.command == "generate-sample":
        generate_sample_data(args.output, args.n)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
