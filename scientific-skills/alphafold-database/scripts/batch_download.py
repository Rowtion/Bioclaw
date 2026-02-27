#!/usr/bin/env python3
"""
Batch Download AlphaFold Structures

This script downloads AlphaFold predictions for multiple proteins
from a list of UniProt IDs.

Usage:
    python batch_download.py --input uniprot_ids.txt --output ./batch_structures
"""

import argparse
import csv
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests


ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction"
ALPHAFOLD_FILES_URL = "https://alphafold.ebi.ac.uk/files"


def read_uniprot_list(input_file):
    """
    Read list of UniProt IDs from file.
    
    Args:
        input_file: Path to file containing UniProt IDs (one per line)
        
    Returns:
        List of UniProt IDs
    """
    uniprot_ids = []
    
    with open(input_file, 'r') as f:
        for line in f:
            # Strip whitespace and comments
            line = line.strip()
            if line and not line.startswith('#'):
                # Handle CSV format
                if ',' in line:
                    parts = line.split(',')
                    uniprot_ids.append(parts[0].strip())
                else:
                    uniprot_ids.append(line)
    
    return uniprot_ids


def check_structure_exists(uniprot_id):
    """
    Check if AlphaFold prediction exists for a UniProt ID.
    
    Args:
        uniprot_id: UniProt accession
        
    Returns:
        (exists, metadata) tuple
    """
    url = f"{ALPHAFOLD_API_URL}/{uniprot_id}"
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            return True, data
        else:
            return False, None
    except requests.exceptions.RequestException:
        return False, None


def download_single_file(url, output_path, max_retries=3):
    """
    Download a single file with retry logic.
    
    Args:
        url: URL to download
        output_path: Path to save file
        max_retries: Maximum number of retry attempts
        
    Returns:
        True if successful
    """
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(response.content)
                return True
            else:
                time.sleep(1)
        except requests.exceptions.RequestException:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)  # Exponential backoff
            continue
    
    return False


def download_protein(uniprot_id, output_dir, version="v4", format="cif"):
    """
    Download all files for a single protein.
    
    Args:
        uniprot_id: UniProt accession
        output_dir: Output directory
        version: AlphaFold version
        format: File format
        
    Returns:
        Dictionary with download results
    """
    result = {
        'uniprot_id': uniprot_id,
        'success': False,
        'files': [],
        'error': None
    }
    
    # Check if prediction exists
    exists, metadata = check_structure_exists(uniprot_id)
    
    if not exists:
        result['error'] = "No prediction found"
        return result
    
    # Create output directory for this protein
    protein_dir = os.path.join(output_dir, uniprot_id)
    os.makedirs(protein_dir, exist_ok=True)
    
    # Download each prediction
    for pred in metadata:
        alphafold_id = pred.get('entryId')
        
        # Download structure
        structure_url = f"{ALPHAFOLD_FILES_URL}/{alphafold_id}-model_{version}.{format}"
        structure_path = os.path.join(protein_dir, f"{alphafold_id}.{format}")
        
        if download_single_file(structure_url, structure_path):
            result['files'].append(structure_path)
        
        # Download confidence
        confidence_url = f"{ALPHAFOLD_FILES_URL}/{alphafold_id}-confidence_{version}.json"
        confidence_path = os.path.join(protein_dir, f"{alphafold_id}_confidence.json")
        
        if download_single_file(confidence_url, confidence_path):
            result['files'].append(confidence_path)
    
    result['success'] = len(result['files']) > 0
    return result


def batch_download(uniprot_ids, output_dir, version="v4", format="cif", 
                   max_workers=5, delay=0.5):
    """
    Download structures for multiple proteins in parallel.
    
    Args:
        uniprot_ids: List of UniProt IDs
        output_dir: Output directory
        version: AlphaFold version
        format: File format
        max_workers: Number of parallel downloads
        delay: Delay between requests
        
    Returns:
        List of download results
    """
    os.makedirs(output_dir, exist_ok=True)
    results = []
    
    print(f"Starting batch download of {len(uniprot_ids)} proteins...")
    print(f"Output directory: {output_dir}")
    print(f"Parallel workers: {max_workers}")
    print("-" * 60)
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_id = {
            executor.submit(download_protein, uid, output_dir, version, format): uid
            for uid in uniprot_ids
        }
        
        # Process results as they complete
        for i, future in enumerate(as_completed(future_to_id)):
            uniprot_id = future_to_id[future]
            
            try:
                result = future.result()
                results.append(result)
                
                # Print progress
                status = "✓" if result['success'] else "✗"
                error_msg = f" ({result['error']})" if result['error'] else ""
                print(f"[{i+1}/{len(uniprot_ids)}] {status} {uniprot_id}{error_msg}")
                
            except Exception as e:
                results.append({
                    'uniprot_id': uniprot_id,
                    'success': False,
                    'error': str(e)
                })
                print(f"[{i+1}/{len(uniprot_ids)}] ✗ {uniprot_id} (Error: {e})")
            
            # Rate limiting
            time.sleep(delay)
    
    return results


def save_summary(results, output_file):
    """
    Save download summary to CSV file.
    
    Args:
        results: List of download results
        output_file: Path to output CSV
    """
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['UniProt_ID', 'Success', 'Files_Downloaded', 'Error'])
        
        for result in results:
            writer.writerow([
                result['uniprot_id'],
                result['success'],
                len(result['files']),
                result['error'] or ''
            ])
    
    print(f"\nSummary saved to: {output_file}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Batch download AlphaFold structures"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="File containing UniProt IDs (one per line)"
    )
    parser.add_argument(
        "--output", "-o",
        default="./batch_structures",
        help="Output directory (default: ./batch_structures)"
    )
    parser.add_argument(
        "--version", "-v",
        default="v4",
        help="AlphaFold version"
    )
    parser.add_argument(
        "--format", "-f",
        default="cif",
        choices=["cif", "pdb"],
        help="Structure file format"
    )
    parser.add_argument(
        "--workers", "-w",
        type=int,
        default=5,
        help="Number of parallel downloads (default: 5)"
    )
    parser.add_argument(
        "--delay", "-d",
        type=float,
        default=0.5,
        help="Delay between requests in seconds (default: 0.5)"
    )
    
    args = parser.parse_args()
    
    # Read UniProt IDs
    print("="*60)
    print("AlphaFold Batch Download")
    print("="*60)
    
    uniprot_ids = read_uniprot_list(args.input)
    print(f"\nLoaded {len(uniprot_ids)} UniProt ID(s) from {args.input}")
    
    if len(uniprot_ids) == 0:
        print("Error: No UniProt IDs found in input file")
        sys.exit(1)
    
    # Download structures
    results = batch_download(
        uniprot_ids,
        args.output,
        version=args.version,
        format=args.format,
        max_workers=args.workers,
        delay=args.delay
    )
    
    # Print summary
    successful = sum(1 for r in results if r['success'])
    failed = len(results) - successful
    
    print("\n" + "="*60)
    print("Download Summary")
    print("="*60)
    print(f"Total: {len(results)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    
    # Save summary CSV
    summary_file = os.path.join(args.output, "download_summary.csv")
    save_summary(results, summary_file)
    
    print("\nBatch download complete!")


if __name__ == "__main__":
    main()
