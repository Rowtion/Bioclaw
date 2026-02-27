#!/usr/bin/env python3
"""
MatchMS - Mass Spectrometry Data Processing
Process and compare mass spectrometry data for metabolomics.
"""

import numpy as np
from typing import List, Tuple, Dict
import json


def load_mgf(filepath: str) -> List[Dict]:
    """Load spectra from MGF file."""
    try:
        from matchms.importing import load_from_mgf
        return list(load_from_mgf(filepath))
    except ImportError:
        print("Error: matchms not installed. Run: pip install matchms")
        return []


def normalize_intensities(spectrum):
    """Normalize spectrum intensities."""
    from matchms.filtering import normalize_intensities as norm
    return norm(spectrum)


def calculate_similarity(spec1, spec2, method: str = "cosine") -> float:
    """Calculate similarity between two spectra."""
    from matchms.similarity import CosineGreedy, ModifiedCosine
    
    if method == "cosine":
        similarity = CosineGreedy(tolerance=0.005)
    elif method == "modified_cosine":
        similarity = ModifiedCosine(tolerance=0.005)
    else:
        raise ValueError(f"Unknown method: {method}")
    
    score = similarity.pair(spec1, spec2)
    return score[0]


def compare_spectra_library(query_spec, library_path: str, top_n: int = 10) -> List[Tuple]:
    """Compare query spectrum against a spectral library."""
    library = load_mgf(library_path)
    scores = []
    
    for lib_spec in library:
        score = calculate_similarity(query_spec, lib_spec)
        scores.append((lib_spec.metadata.get("compound_name", "Unknown"), score))
    
    scores.sort(key=lambda x: x[1], reverse=True)
    return scores[:top_n]


def process_spectra_batch(input_path: str, output_path: str):
    """Process batch of spectra."""
    spectra = load_mgf(input_path)
    processed = []
    
    for spec in spectra:
        # Normalize
        spec = normalize_intensities(spec)
        # Remove noise
        from matchms.filtering import select_by_relative_intensity
        spec = select_by_relative_intensity(spec, intensity_from=0.01)
        processed.append(spec)
    
    # Save processed
    from matchms.exporting import save_as_mgf
    save_as_mgf(processed, output_path)
    print(f"Processed {len(processed)} spectra")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="MatchMS Spectra Tools")
    parser.add_argument("command", choices=["load", "compare", "process"])
    parser.add_argument("--input", help="Input MGF file")
    parser.add_argument("--library", help="Library MGF file")
    parser.add_argument("--output", help="Output file")
    parser.add_argument("--query-idx", type=int, default=0, help="Query spectrum index")
    parser.add_argument("--top-n", type=int, default=10)
    
    args = parser.parse_args()
    
    if args.command == "load":
        spectra = load_mgf(args.input)
        print(f"Loaded {len(spectra)} spectra")
        for i, s in enumerate(spectra[:5]):
            print(f"  {i}: {s.metadata.get('compound_name', 'N/A')}")
    
    elif args.command == "compare":
        spectra = load_mgf(args.input)
        if args.query_idx >= len(spectra):
            print("Invalid query index")
            return
        results = compare_spectra_library(spectra[args.query_idx], args.library, args.top_n)
        print("Top matches:")
        for name, score in results:
            print(f"  {score:.3f} | {name}")
    
    elif args.command == "process":
        process_spectra_batch(args.input, args.output)


if __name__ == "__main__":
    main()
