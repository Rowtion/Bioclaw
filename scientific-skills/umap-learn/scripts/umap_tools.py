#!/usr/bin/env python3
"""
UMAP-learn Tools
Dimensionality reduction and visualization for high-dimensional data.
"""

import numpy as np
from typing import List, Dict, Tuple
import json


def reduce_dimensions(data: np.ndarray, n_components: int = 2, 
                     n_neighbors: int = 15, min_dist: float = 0.1) -> np.ndarray:
    """Reduce dimensions using UMAP."""
    try:
        import umap
        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=42
        )
        embedding = reducer.fit_transform(data)
        return embedding
    except ImportError:
        print("Error: umap-learn not installed. Run: pip install umap-learn")
        return None


def visualize_embedding(embedding: np.ndarray, labels: List = None, 
                       save_path: str = None):
    """Visualize UMAP embedding."""
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 8))
        if labels is not None:
            scatter = plt.scatter(embedding[:, 0], embedding[:, 1], 
                                c=labels, cmap='Spectral', s=5)
            plt.colorbar(scatter)
        else:
            plt.scatter(embedding[:, 0], embedding[:, 1], s=5)
        
        plt.title('UMAP Projection')
        plt.xlabel('UMAP 1')
        plt.ylabel('UMAP 2')
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Saved to: {save_path}")
        else:
            plt.show()
    except ImportError:
        print("Error: matplotlib not installed")


def compare_parameters(data: np.ndarray, param_grid: Dict) -> Dict:
    """Compare different UMAP parameters."""
    results = {}
    
    for n_neighbors in param_grid.get('n_neighbors', [5, 15, 30]):
        for min_dist in param_grid.get('min_dist', [0.1, 0.5]):
            key = f"nn{n_neighbors}_md{min_dist}"
            embedding = reduce_dimensions(data, n_neighbors=n_neighbors, min_dist=min_dist)
            if embedding is not None:
                results[key] = {
                    "shape": embedding.shape,
                    "n_neighbors": n_neighbors,
                    "min_dist": min_dist
                }
    
    return results


def load_and_reduce(filepath: str, **umap_params) -> np.ndarray:
    """Load data from file and reduce dimensions."""
    try:
        import pandas as pd
        data = pd.read_csv(filepath).values
        return reduce_dimensions(data, **umap_params)
    except Exception as e:
        print(f"Error loading file: {e}")
        return None


def main():
    import argparse
    parser = argparse.ArgumentParser(description="UMAP Tools")
    parser.add_argument("command", choices=["reduce", "visualize", "compare"])
    parser.add_argument("--file", help="Input CSV file")
    parser.add_argument("--output", help="Output file for embedding")
    parser.add_argument("--n-components", type=int, default=2)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.1)
    parser.add_argument("--labels", help="Labels file (optional)")
    parser.add_argument("--plot", help="Save plot to file")
    
    args = parser.parse_args()
    
    if args.command == "reduce":
        if not args.file:
            print("Error: --file required")
            return
        embedding = load_and_reduce(
            args.file, 
            n_components=args.n_components,
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist
        )
        if embedding is not None:
            np.savetxt(args.output or "umap_embedding.csv", embedding, delimiter=',')
            print(f"Embedding shape: {embedding.shape}")
    
    elif args.command == "visualize":
        if not args.file:
            print("Error: --file required")
            return
        embedding = np.loadtxt(args.file, delimiter=',')
        labels = None
        if args.labels:
            labels = np.loadtxt(args.labels)
        visualize_embedding(embedding, labels, args.plot)
    
    elif args.command == "compare":
        if not args.file:
            print("Error: --file required")
            return
        import pandas as pd
        data = pd.read_csv(args.file).values
        results = compare_parameters(data, {
            'n_neighbors': [5, 15, 30],
            'min_dist': [0.1, 0.5]
        })
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
