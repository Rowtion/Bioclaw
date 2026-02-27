#!/usr/bin/env python3
"""
Scikit-Bio Tools
Bioinformatics utilities for sequence analysis, phylogenetics, and diversity.
"""

from typing import List, Dict, Tuple
import json


def load_sequences(filepath: str, format: str = "fasta") -> List[Dict]:
    """Load biological sequences from file."""
    try:
        from skbio import read
        sequences = []
        for seq in read(filepath, format=format):
            sequences.append({
                "id": seq.metadata.get("id", "unknown"),
                "sequence": str(seq),
                "length": len(seq)
            })
        return sequences
    except ImportError:
        print("Error: scikit-bio not installed. Run: pip install scikit-bio")
        return []


def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of DNA sequence."""
    try:
        from skbio import DNA
        dna = DNA(sequence)
        return dna.gc_content()
    except:
        # Fallback calculation
        gc = sequence.upper().count('G') + sequence.upper().count('C')
        return gc / len(sequence) if sequence else 0


def align_sequences(sequences: List[str]) -> Dict:
    """Perform multiple sequence alignment."""
    try:
        from skbio.alignment import global_pairwise_align_nucleotide
        from skbio import DNA
        
        if len(sequences) < 2:
            return {"error": "Need at least 2 sequences"}
        
        # Pairwise alignment of first two sequences
        seq1, seq2 = DNA(sequences[0]), DNA(sequences[1])
        alignment, score, _ = global_pairwise_align_nucleotide(seq1, seq2)
        
        return {
            "alignment": str(alignment),
            "score": float(score),
            "identity": alignment[0].match_frequency(alignment[1])
        }
    except ImportError:
        print("Error: scikit-bio not installed")
        return {}


def calculate_alpha_diversity(otu_table: List[List[int]], metric: str = "shannon") -> List[float]:
    """Calculate alpha diversity for samples."""
    try:
        from skbio.diversity import alpha_diversity
        import numpy as np
        
        otu_array = np.array(otu_table)
        ids = [f"sample_{i}" for i in range(len(otu_table))]
        
        result = alpha_diversity(metric, otu_array, ids)
        return result.tolist()
    except ImportError:
        print("Error: scikit-bio not installed")
        return []


def build_phylogenetic_tree(sequences: List[Dict]) -> str:
    """Build phylogenetic tree from sequences."""
    try:
        from skbio.tree import nj
        from skbio import DistanceMatrix
        from skbio.sequence import DistanceMatrix as SeqDistMatrix
        
        # Calculate distance matrix
        ids = [s["id"] for s in sequences]
        # Simplified - in practice would use proper distance calculation
        
        return "Tree construction requires distance matrix input"
    except ImportError:
        print("Error: scikit-bio not installed")
        return ""


def analyze_sequence_stats(sequences: List[Dict]) -> Dict:
    """Analyze statistics of sequence collection."""
    if not sequences:
        return {}
    
    lengths = [s["length"] for s in sequences]
    
    return {
        "num_sequences": len(sequences),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "mean_length": sum(lengths) / len(lengths),
        "total_bases": sum(lengths)
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scikit-Bio Tools")
    parser.add_argument("command", 
                       choices=["load", "gc", "align", "diversity", "stats"])
    parser.add_argument("--file", help="Input file")
    parser.add_argument("--sequence", help="DNA sequence")
    parser.add_argument("--format", default="fasta", help="File format")
    parser.add_argument("--metric", default="shannon", help="Diversity metric")
    
    args = parser.parse_args()
    
    if args.command == "load":
        if not args.file:
            print("Error: --file required")
            return
        sequences = load_sequences(args.file, args.format)
        print(f"Loaded {len(sequences)} sequences")
        for s in sequences[:5]:
            print(f"  {s['id']}: {s['length']} bp")
    
    elif args.command == "gc":
        if not args.sequence:
            print("Error: --sequence required")
            return
        gc = calculate_gc_content(args.sequence)
        print(f"GC Content: {gc*100:.2f}%")
    
    elif args.command == "stats":
        if not args.file:
            print("Error: --file required")
            return
        sequences = load_sequences(args.file, args.format)
        stats = analyze_sequence_stats(sequences)
        print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
