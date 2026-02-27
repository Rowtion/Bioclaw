#!/usr/bin/env python3
"""
BLAST Search with Biopython

This script demonstrates running BLAST searches and
parsing results using Biopython.

Usage:
    python blast_search.py
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


def blast_online():
    """Run online BLAST search."""
    print("\n" + "="*60)
    print("Running Online BLAST Search")
    print("="*60)
    
    # Short sequence for demo (full searches can take time)
    sequence = "ATCGATCGATCGATCGATCGATCG"
    
    print(f"\nQuery sequence: {sequence}")
    print("Running BLASTN against nt database...")
    print("(This may take a moment...)")
    
    try:
        # Run BLAST
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
        
        # Parse results
        blast_record = NCBIXML.read(result_handle)
        
        print(f"\nQuery: {blast_record.query}")
        print(f"Version: {blast_record.version}")
        
        # Display top hits
        print("\nTop 5 hits:")
        for i, alignment in enumerate(blast_record.alignments[:5], 1):
            print(f"\n{i}. {alignment.title[:60]}...")
            
            # Get first HSP (High Scoring Pair)
            if alignment.hsps:
                hsp = alignment.hsps[0]
                print(f"   E-value: {hsp.expect:.2e}")
                print(f"   Score: {hsp.score}")
                print(f"   Identities: {hsp.identities}/{hsp.align_length}")
        
        result_handle.close()
        
    except Exception as e:
        print(f"\nError: {e}")
        print("Note: BLAST searches require internet connection.")


def parse_blast_xml():
    """Parse BLAST XML output."""
    print("\n" + "="*60)
    print("Parsing BLAST XML Results")
    print("="*60)
    
    example_xml = """<?xml version="1.0"?>
    <!-- Example BLAST XML structure -->
    """
    
    print("\nExample code for parsing BLAST XML:")
    code = '''
from Bio.Blast import NCBIXML

# Parse BLAST XML file
with open("blast_results.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)
    
    # Iterate through alignments
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.01:  # Filter by E-value
                print(f"Sequence: {alignment.title}")
                print(f"E-value: {hsp.expect}")
                print(f"Alignment: {hsp.query}")
'''
    print(code)


def filter_blast_results():
    """Example: Filter BLAST results."""
    print("\n" + "="*60)
    print("Filtering BLAST Results")
    print("="*60)
    
    code = '''
# Filter hits by E-value threshold
e_value_threshold = 0.001

significant_hits = []
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < e_value_threshold:
            significant_hits.append({
                'title': alignment.title,
                'e_value': hsp.expect,
                'score': hsp.score,
                'identities': hsp.identities
            })

# Sort by E-value
significant_hits.sort(key=lambda x: x['e_value'])

print(f"Found {len(significant_hits)} significant hits")
'''
    print(code)


def main():
    """Run BLAST examples."""
    print("Biopython BLAST Search Examples")
    print("=" * 60)
    
    # Note: BLAST searches can be slow
    print("\nNote: Online BLAST searches can take several minutes.")
    print("The demo uses a short sequence, but may still take time.")
    
    response = input("\nRun online BLAST search? (y/n): ")
    if response.lower() == 'y':
        blast_online()
    else:
        print("\nSkipping online BLAST.")
    
    parse_blast_xml()
    filter_blast_results()
    
    print("\n" + "="*60)
    print("BLAST examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
