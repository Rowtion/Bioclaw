#!/usr/bin/env python3
"""
Sequence Analysis with Biopython

This script demonstrates DNA and protein sequence analysis
including GC content, translation, and motif finding.

Usage:
    python sequence_analysis.py
"""

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np


def dna_sequence_analysis():
    """Analyze DNA sequences."""
    print("\n" + "="*60)
    print("DNA Sequence Analysis")
    print("="*60)
    
    # Create DNA sequence
    dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    
    print(f"\nOriginal sequence: {dna_seq}")
    print(f"Length: {len(dna_seq)} bp")
    
    # GC content
    gc = gc_fraction(dna_seq)
    print(f"GC content: {gc:.2%}")
    
    # Complement and reverse complement
    print(f"\nComplement: {dna_seq.complement()}")
    print(f"Reverse complement: {dna_seq.reverse_complement()}")
    
    # Transcription
    rna_seq = dna_seq.transcribe()
    print(f"\nTranscribed (RNA): {rna_seq}")
    
    # Translation
    protein_seq = dna_seq.translate()
    print(f"Translated (Protein): {protein_seq}")
    
    # Molecular weight
    mw_dna = molecular_weight(dna_seq, seq_type='DNA')
    print(f"\nMolecular weight: {mw_dna:.2f} Da")


def protein_sequence_analysis():
    """Analyze protein sequences."""
    print("\n" + "="*60)
    print("Protein Sequence Analysis")
    print("="*60)
    
    # Protein sequence
    protein_seq = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    
    print(f"\nSequence: {protein_seq[:50]}...")
    print(f"Length: {len(protein_seq)} aa")
    
    # Analyze protein properties
    analysis = ProteinAnalysis(protein_seq)
    
    print(f"\nMolecular weight: {analysis.molecular_weight():.2f} Da")
    print(f"Isoelectric point: {analysis.isoelectric_point():.2f}")
    print(f"Aromaticity: {analysis.aromaticity():.3f}")
    print(f"Instability index: {analysis.instability_index():.2f}")
    
    # Amino acid composition
    aa_counts = analysis.count_amino_acids()
    print(f"\nTop 5 amino acids:")
    sorted_aa = sorted(aa_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    for aa, count in sorted_aa:
        print(f"  {aa}: {count}")


def find_orfs():
    """Find open reading frames."""
    print("\n" + "="*60)
    print("Finding Open Reading Frames (ORFs)")
    print("="*60)
    
    dna_seq = Seq(
        "ATGAAATAGATGCGCTAGATGGCCCGGGTAA"
    )
    
    print(f"Sequence: {dna_seq}")
    print(f"Length: {len(dna_seq)} bp")
    
    # Find all ORFs starting with ATG
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    
    orfs = []
    for i in range(len(dna_seq) - 2):
        if str(dna_seq[i:i+3]) == start_codon:
            # Look for stop codon
            for j in range(i + 3, len(dna_seq) - 2, 3):
                codon = str(dna_seq[j:j+3])
                if codon in stop_codons:
                    orf = dna_seq[i:j+3]
                    orfs.append((i, j+3, len(orf), str(orf.translate())))
                    break
    
    print(f"\nFound {len(orfs)} ORF(s):")
    for start, end, length, protein in orfs:
        print(f"  Position {start}-{end}: {length} bp -> {protein}")


def main():
    """Run sequence analysis examples."""
    print("Biopython Sequence Analysis")
    print("=" * 60)
    
    dna_sequence_analysis()
    protein_sequence_analysis()
    find_orfs()
    
    print("\n" + "="*60)
    print("Sequence analysis completed!")
    print("="*60)


if __name__ == "__main__":
    main()
