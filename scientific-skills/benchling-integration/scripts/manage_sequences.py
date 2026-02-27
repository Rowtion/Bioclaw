#!/usr/bin/env python3
"""
Benchling Sequence Management

This script demonstrates creating and updating DNA/protein
sequences in Benchling using the Benchling SDK.

Usage:
    export BENCHLING_API_KEY="your_key"
    export BENCHLING_URL="https://your-tenant.benchling.com"
    python manage_sequences.py
"""

import os
import sys

# Check for credentials
if not os.environ.get('BENCHLING_API_KEY') or not os.environ.get('BENCHLING_URL'):
    print("Error: Please set BENCHLING_API_KEY and BENCHLING_URL environment variables")
    print("\nExample:")
    print("  export BENCHLING_API_KEY='your_api_key'")
    print("  export BENCHLING_URL='https://your-tenant.benchling.com'")
    sys.exit(1)


def example_create_dna():
    """Example: Create a DNA sequence in Benchling."""
    print("\n" + "="*60)
    print("Creating DNA Sequence")
    print("="*60)
    
    example_code = '''
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.models import DnaSequenceCreate

# Initialize Benchling client
benchling = Benchling(
    url=os.environ['BENCHLING_URL'],
    auth_method=ApiKeyAuth(os.environ['BENCHLING_API_KEY'])
)

# Create DNA sequence
sequence = benchling.dna_sequences.create(
    DnaSequenceCreate(
        name="GFP Reporter Construct",
        bases="ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA",
        is_circular=True,
        folder_id="your_folder_id",  # Replace with actual folder ID
    )
)

print(f"Created sequence: {sequence.id}")
print(f"Name: {sequence.name}")
'''
    
    print(example_code)
    print("\nNote: Replace 'your_folder_id' with an actual folder ID from your Benchling tenant")


def example_list_sequences():
    """Example: List DNA sequences."""
    print("\n" + "="*60)
    print("Listing DNA Sequences")
    print("="*60)
    
    example_code = '''
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth

benchling = Benchling(
    url=os.environ['BENCHLING_URL'],
    auth_method=ApiKeyAuth(os.environ['BENCHLING_API_KEY'])
)

# List sequences with pagination
sequences = benchling.dna_sequences.list()

print("DNA Sequences:")
for page in sequences:
    for seq in page:
        print(f"  - {seq.name} ({seq.id})")

# Get total count
total = sequences.estimated_count()
print(f"\nTotal sequences: {total}")
'''
    
    print(example_code)


def example_update_sequence():
    """Example: Update a sequence."""
    print("\n" + "="*60)
    print("Updating Sequence")
    print("="*60)
    
    example_code = '''
from benchling_sdk.models import DnaSequenceUpdate

# Update sequence
updated = benchling.dna_sequences.update(
    sequence_id="seq_abc123",  # Replace with actual sequence ID
    dna_sequence=DnaSequenceUpdate(
        name="Updated Sequence Name",
        fields=benchling.models.fields({
            "gene_name": "GFP",
            "organism": "Aequorea victoria"
        })
    )
)

print(f"Updated: {updated.name}")
'''
    
    print(example_code)


def main():
    """Display Benchling SDK examples."""
    print("Benchling Sequence Management Examples")
    print("=" * 60)
    print("\nThese are code examples showing how to use the Benchling SDK.")
    print("Set your credentials as environment variables to run them.")
    
    example_create_dna()
    example_list_sequences()
    example_update_sequence()
    
    print("\n" + "="*60)
    print("Examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
