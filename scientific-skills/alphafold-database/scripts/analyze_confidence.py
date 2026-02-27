#!/usr/bin/env python3
"""
Analyze AlphaFold Confidence Metrics

This script analyzes pLDDT scores and Predicted Aligned Error (PAE)
matrices from AlphaFold predictions.

Usage:
    python analyze_confidence.py --confidence AF-P00520-F1_confidence.json
"""

import argparse
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def load_confidence_scores(confidence_file):
    """
    Load confidence scores from AlphaFold JSON file.
    
    Args:
        confidence_file: Path to confidence JSON file
        
    Returns:
        Dictionary with confidence data
    """
    with open(confidence_file, 'r') as f:
        data = json.load(f)
    return data


def load_pae_matrix(pae_file):
    """
    Load PAE matrix from AlphaFold JSON file.
    
    Args:
        pae_file: Path to PAE JSON file
        
    Returns:
        Numpy array of PAE values
    """
    with open(pae_file, 'r') as f:
        data = json.load(f)
    
    # PAE data is stored as a list of lists
    pae_matrix = np.array(data.get('predicted_aligned_error', data.get('distance', [])))
    return pae_matrix


def analyze_plddt(confidence_data):
    """
    Analyze per-residue pLDDT (predicted Local Distance Difference Test) scores.
    
    Args:
        confidence_data: Loaded confidence JSON data
        
    Returns:
        Dictionary with analysis results
    """
    plddt_scores = np.array(confidence_data.get('confidenceScore', []))
    
    if len(plddt_scores) == 0:
        print("Warning: No pLDDT scores found in confidence file")
        return None
    
    # Confidence categories
    very_high = plddt_scores > 90
    high = (plddt_scores >= 70) & (plddt_scores <= 90)
    low = (plddt_scores >= 50) & (plddt_scores < 70)
    very_low = plddt_scores < 50
    
    analysis = {
        'scores': plddt_scores,
        'mean': np.mean(plddt_scores),
        'median': np.median(plddt_scores),
        'min': np.min(plddt_scores),
        'max': np.max(plddt_scores),
        'std': np.std(plddt_scores),
        'very_high_confidence': np.sum(very_high),
        'high_confidence': np.sum(high),
        'low_confidence': np.sum(low),
        'very_low_confidence': np.sum(very_low),
        'very_high_pct': np.sum(very_high) / len(plddt_scores) * 100,
        'high_pct': np.sum(high) / len(plddt_scores) * 100,
        'low_pct': np.sum(low) / len(plddt_scores) * 100,
        'very_low_pct': np.sum(very_low) / len(plddt_scores) * 100,
    }
    
    return analysis


def analyze_pae(pae_matrix):
    """
    Analyze Predicted Aligned Error matrix.
    
    Args:
        pae_matrix: 2D numpy array of PAE values
        
    Returns:
        Dictionary with analysis results
    """
    if pae_matrix.size == 0:
        print("Warning: Empty PAE matrix")
        return None
    
    analysis = {
        'matrix': pae_matrix,
        'mean': np.mean(pae_matrix),
        'median': np.median(pae_matrix),
        'min': np.min(pae_matrix),
        'max': np.max(pae_matrix),
        'std': np.std(pae_matrix),
        'q25': np.percentile(pae_matrix, 25),
        'q75': np.percentile(pae_matrix, 75),
    }
    
    # Calculate fraction with low PAE (<5 Å - confident positioning)
    low_pae_fraction = np.sum(pae_matrix < 5) / pae_matrix.size
    analysis['low_pae_fraction'] = low_pae_fraction
    
    return analysis


def plot_plddt_analysis(analysis, output_file='plddt_analysis.png'):
    """
    Create visualization of pLDDT scores.
    
    Args:
        analysis: pLDDT analysis results
        output_file: Path to save plot
    """
    scores = analysis['scores']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: pLDDT scores along sequence
    ax = axes[0, 0]
    ax.plot(scores, linewidth=1.5, color='blue')
    ax.axhline(90, color='green', linestyle='--', alpha=0.7, label='Very High (>90)')
    ax.axhline(70, color='yellow', linestyle='--', alpha=0.7, label='High (70-90)')
    ax.axhline(50, color='orange', linestyle='--', alpha=0.7, label='Low (50-70)')
    ax.axhline(50, color='red', linestyle='--', alpha=0.7, label='Very Low (<50)')
    ax.fill_between(range(len(scores)), 90, 100, alpha=0.2, color='green')
    ax.fill_between(range(len(scores)), 70, 90, alpha=0.2, color='yellow')
    ax.fill_between(range(len(scores)), 50, 70, alpha=0.2, color='orange')
    ax.fill_between(range(len(scores)), 0, 50, alpha=0.2, color='red')
    ax.set_xlabel('Residue')
    ax.set_ylabel('pLDDT Score')
    ax.set_title('Per-Residue pLDDT Scores')
    ax.set_ylim(0, 100)
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Confidence distribution histogram
    ax = axes[0, 1]
    ax.hist(scores, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(analysis['mean'], color='red', linestyle='--', linewidth=2, label=f"Mean: {analysis['mean']:.1f}")
    ax.axvline(analysis['median'], color='orange', linestyle='--', linewidth=2, label=f"Median: {analysis['median']:.1f}")
    ax.set_xlabel('pLDDT Score')
    ax.set_ylabel('Frequency')
    ax.set_title('pLDDT Score Distribution')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Confidence categories pie chart
    ax = axes[1, 0]
    categories = ['Very High\n(>90)', 'High\n(70-90)', 'Low\n(50-70)', 'Very Low\n(<50)']
    values = [
        analysis['very_high_confidence'],
        analysis['high_confidence'],
        analysis['low_confidence'],
        analysis['very_low_confidence']
    ]
    colors = ['green', 'yellow', 'orange', 'red']
    explode = (0.05, 0, 0, 0.05)
    
    ax.pie(values, labels=categories, colors=colors, autopct='%1.1f%%',
           startangle=90, explode=explode, shadow=True)
    ax.set_title('Confidence Category Distribution')
    
    # Plot 4: Summary statistics table
    ax = axes[1, 1]
    ax.axis('off')
    
    stats_text = f"""
    pLDDT Statistics
    
    Total Residues: {len(scores)}
    
    Mean Score: {analysis['mean']:.2f}
    Median Score: {analysis['median']:.2f}
    Min Score: {analysis['min']:.2f}
    Max Score: {analysis['max']:.2f}
    Std Dev: {analysis['std']:.2f}
    
    Confidence Distribution:
      Very High (>90): {analysis['very_high_pct']:.1f}%
      High (70-90): {analysis['high_pct']:.1f}%
      Low (50-70): {analysis['low_pct']:.1f}%
      Very Low (<50): {analysis['very_low_pct']:.1f}%
    
    Interpretation:
    """
    
    if analysis['very_high_pct'] + analysis['high_pct'] > 80:
        stats_text += "    High overall confidence - suitable for\n    detailed structural analysis"
    elif analysis['very_high_pct'] + analysis['high_pct'] > 50:
        stats_text += "    Moderate confidence - reliable backbone\n    structure expected"
    else:
        stats_text += "    Low confidence - use with caution,\n    may contain disordered regions"
    
    ax.text(0.1, 0.9, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"pLDDT analysis plot saved to: {output_file}")
    plt.close()


def plot_pae_analysis(analysis, output_file='pae_analysis.png'):
    """
    Create visualization of PAE matrix.
    
    Args:
        analysis: PAE analysis results
        output_file: Path to save plot
    """
    pae_matrix = analysis['matrix']
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Plot 1: PAE heatmap
    ax = axes[0]
    im = ax.imshow(pae_matrix, cmap='viridis_r', vmin=0, vmax=30)
    ax.set_xlabel('Residue i')
    ax.set_ylabel('Residue j')
    ax.set_title('Predicted Aligned Error (PAE) Matrix')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('PAE (Å)')
    
    # Plot 2: PAE statistics
    ax = axes[1]
    ax.axis('off')
    
    stats_text = f"""
    PAE Statistics
    
    Matrix Shape: {pae_matrix.shape}
    
    Mean PAE: {analysis['mean']:.2f} Å
    Median PAE: {analysis['median']:.2f} Å
    Min PAE: {analysis['min']:.2f} Å
    Max PAE: {analysis['max']:.2f} Å
    Std Dev: {analysis['std']:.2f} Å
    
    Percentiles:
      25th: {analysis['q25']:.2f} Å
      75th: {analysis['q75']:.2f} Å
    
    Low PAE fraction (<5Å): {analysis['low_pae_fraction']:.1%}
    
    Interpretation:
      Low PAE (<5 Å): Confident relative 
        positioning of domains
      Medium PAE (5-15 Å): Moderate confidence
      High PAE (>15 Å): Uncertain relative 
        positions, domains may be mobile
    """
    
    ax.text(0.1, 0.9, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"PAE analysis plot saved to: {output_file}")
    plt.close()


def print_summary(plddt_analysis, pae_analysis):
    """
    Print summary of confidence analysis.
    
    Args:
        plddt_analysis: pLDDT analysis results
        pae_analysis: PAE analysis results
    """
    print("\n" + "="*60)
    print("Confidence Analysis Summary")
    print("="*60)
    
    if plddt_analysis:
        print("\n📊 pLDDT (Local Confidence):")
        print(f"   Mean: {plddt_analysis['mean']:.2f} (range: 0-100)")
        print(f"   Distribution:")
        print(f"     - Very High (>90): {plddt_analysis['very_high_pct']:.1f}%")
        print(f"     - High (70-90): {plddt_analysis['high_pct']:.1f}%")
        print(f"     - Low (50-70): {plddt_analysis['low_pct']:.1f}%")
        print(f"     - Very Low (<50): {plddt_analysis['very_low_pct']:.1f}%")
    
    if pae_analysis:
        print("\n📐 PAE (Domain Positioning Confidence):")
        print(f"   Mean: {pae_analysis['mean']:.2f} Å")
        print(f"   Low PAE fraction (<5Å): {pae_analysis['low_pae_fraction']:.1%}")
    
    print("\n" + "="*60)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Analyze AlphaFold confidence metrics"
    )
    parser.add_argument(
        "--confidence", "-c",
        help="Path to confidence JSON file"
    )
    parser.add_argument(
        "--pae", "-p",
        help="Path to PAE JSON file"
    )
    parser.add_argument(
        "--output-prefix", "-o",
        default="analysis",
        help="Prefix for output plot files"
    )
    
    args = parser.parse_args()
    
    if not args.confidence and not args.pae:
        parser.error("Please provide at least --confidence or --pae file")
    
    print("="*60)
    print("AlphaFold Confidence Analysis")
    print("="*60)
    
    plddt_analysis = None
    pae_analysis = None
    
    # Analyze pLDDT scores
    if args.confidence:
        print(f"\n📁 Loading confidence file: {args.confidence}")
        try:
            confidence_data = load_confidence_scores(args.confidence)
            plddt_analysis = analyze_plddt(confidence_data)
            
            if plddt_analysis:
                output_file = f"{args.output_prefix}_plddt.png"
                plot_plddt_analysis(plddt_analysis, output_file)
        except Exception as e:
            print(f"Error analyzing confidence file: {e}")
    
    # Analyze PAE matrix
    if args.pae:
        print(f"\n📁 Loading PAE file: {args.pae}")
        try:
            pae_matrix = load_pae_matrix(args.pae)
            pae_analysis = analyze_pae(pae_matrix)
            
            if pae_analysis:
                output_file = f"{args.output_prefix}_pae.png"
                plot_pae_analysis(pae_analysis, output_file)
        except Exception as e:
            print(f"Error analyzing PAE file: {e}")
    
    # Print summary
    print_summary(plddt_analysis, pae_analysis)
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()
