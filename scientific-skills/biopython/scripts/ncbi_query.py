#!/usr/bin/env python3
"""
Query NCBI Databases

This script demonstrates querying NCBI databases using
Biopython's Entrez module.

Usage:
    python ncbi_query.py
"""

from Bio import Entrez, SeqIO
import time

# Set your email (required by NCBI)
Entrez.email = "your.email@example.com"
# Optional: Set API key for higher rate limits
# Entrez.api_key = "your_api_key"


def search_pubmed():
    """Search PubMed for publications."""
    print("\n" + "="*60)
    print("Searching PubMed")
    print("="*60)
    
    # Search query
    query = "CRISPR[Title/Abstract] AND 2023[PDAT]"
    
    handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
    results = Entrez.read(handle)
    handle.close()
    
    print(f"\nQuery: {query}")
    print(f"Total results: {results['Count']}")
    print(f"Retrieved IDs: {len(results['IdList'])}")
    
    # Fetch summaries
    if results['IdList']:
        id_list = ",".join(results['IdList'][:5])
        handle = Entrez.esummary(db="pubmed", id=id_list)
        summaries = Entrez.read(handle)
        handle.close()
        
        print("\nTop 5 publications:")
        for i, summary in enumerate(summaries[:5], 1):
            print(f"\n{i}. {summary.get('Title', 'N/A')}")
            print(f"   Authors: {', '.join(summary.get('AuthorList', [])[:3])}...")
            print(f"   Journal: {summary.get('Source', 'N/A')}")


def fetch_sequence():
    """Fetch a sequence from GenBank."""
    print("\n" + "="*60)
    print("Fetching Sequence from GenBank")
    print("="*60)
    
    # Example: Insulin gene
    accession = "NM_000207"
    
    print(f"\nFetching: {accession}")
    
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    print(f"\nRecord ID: {record.id}")
    print(f"Name: {record.name}")
    print(f"Description: {record.description}")
    print(f"Length: {len(record.seq)} bp")
    print(f"Sequence: {record.seq[:50]}...")
    
    # Show features
    print(f"\nFeatures ({len(record.features)}):")
    for feature in record.features[:5]:
        print(f"  {feature.type}: {feature.location}")


def search_protein():
    """Search protein database."""
    print("\n" + "="*60)
    print("Searching Protein Database")
    print("="*60)
    
    query = "insulin[Protein Name] AND mammals[Organism]"
    
    handle = Entrez.esearch(db="protein", term=query, retmax=5)
    results = Entrez.read(handle)
    handle.close()
    
    print(f"\nQuery: {query}")
    print(f"Found: {results['Count']} proteins")
    
    if results['IdList']:
        print("\nProtein IDs:")
        for pid in results['IdList']:
            print(f"  - {pid}")


def main():
    """Run NCBI query examples."""
    print("Biopython NCBI Database Queries")
    print("=" * 60)
    print("\nIMPORTANT: Set your email address in the script!")
    print("Edit this file and set: Entrez.email = 'your@email.com'")
    print()
    
    # Check if email is set
    if Entrez.email == "your.email@example.com":
        print("WARNING: Please set your email address before running queries.")
        print("NCBI requires a valid email address for API access.")
        return
    
    try:
        search_pubmed()
        time.sleep(1)  # Be nice to NCBI servers
        
        fetch_sequence()
        time.sleep(1)
        
        search_protein()
        
    except Exception as e:
        print(f"\nError: {e}")
        print("Note: Internet connection required for NCBI queries.")
    
    print("\n" + "="*60)
    print("NCBI query examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
