#!/usr/bin/env python3
"""
Download AlphaFold Protein Structures

This script downloads protein structures from the AlphaFold database
by UniProt ID, including PDB/mmCIF files and confidence metrics.

Usage:
    python download_structure.py --uniprot P00520 --output ./structures
"""

import argparse
import os
import requests
import sys
from pathlib import Path


# AlphaFold DB API base URL
ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction"
ALPHAFOLD_FILES_URL = "https://alphafold.ebi.ac.uk/files"


def get_prediction_metadata(uniprot_id):
    """
    Get prediction metadata for a UniProt accession.
    
    Args:
        uniprot_id: UniProt accession (e.g., 'P00520')
        
    Returns:
        List of prediction dictionaries
    """
    url = f"{ALPHAFOLD_API_URL}/{uniprot_id}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error fetching metadata for {uniprot_id}: {e}")
        return None


def download_file(url, output_path, description="file"):
    """
    Download a file from URL with progress indication.
    
    Args:
        url: URL to download from
        output_path: Path to save file
        description: Description for progress messages
        
    Returns:
        True if successful, False otherwise
    """
    try:
        print(f"  Downloading {description}...")
        response = requests.get(url, timeout=60, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_path, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\r    Progress: {percent:.1f}%", end='', flush=True)
        
        print(f"\n  Saved to: {output_path}")
        return True
        
    except requests.exceptions.RequestException as e:
        print(f"\n  Error downloading {description}: {e}")
        return False


def download_structure(alphafold_id, output_dir, version="v4", format="cif"):
    """
    Download structure file for an AlphaFold prediction.
    
    Args:
        alphafold_id: AlphaFold identifier (e.g., 'AF-P00520-F1')
        output_dir: Directory to save file
        version: Database version
        format: File format ('cif' or 'pdb')
        
    Returns:
        Path to downloaded file or None
    """
    url = f"{ALPHAFOLD_FILES_URL}/{alphafold_id}-model_{version}.{format}"
    output_path = os.path.join(output_dir, f"{alphafold_id}.{format}")
    
    if download_file(url, output_path, f"structure ({format})"):
        return output_path
    return None


def download_confidence(alphafold_id, output_dir, version="v4"):
    """
    Download confidence scores JSON.
    
    Args:
        alphafold_id: AlphaFold identifier
        output_dir: Directory to save file
        version: Database version
        
    Returns:
        Path to downloaded file or None
    """
    url = f"{ALPHAFOLD_FILES_URL}/{alphafold_id}-confidence_{version}.json"
    output_path = os.path.join(output_dir, f"{alphafold_id}_confidence.json")
    
    if download_file(url, output_path, "confidence scores"):
        return output_path
    return None


def download_pae(alphafold_id, output_dir, version="v4"):
    """
    Download Predicted Aligned Error (PAE) matrix.
    
    Args:
        alphafold_id: AlphaFold identifier
        output_dir: Directory to save file
        version: Database version
        
    Returns:
        Path to downloaded file or None
    """
    url = f"{ALPHAFOLD_FILES_URL}/{alphafold_id}-predicted_aligned_error_{version}.json"
    output_path = os.path.join(output_dir, f"{alphafold_id}_pae.json")
    
    if download_file(url, output_path, "PAE matrix"):
        return output_path
    return None


def download_all(uniprot_id, output_dir, version="v4", format="cif", 
                 include_confidence=True, include_pae=True):
    """
    Download all available files for a UniProt ID.
    
    Args:
        uniprot_id: UniProt accession
        output_dir: Output directory
        version: Database version
        format: Structure file format
        include_confidence: Download confidence scores
        include_pae: Download PAE matrix
        
    Returns:
        Dictionary with paths to downloaded files
    """
    print(f"\nFetching metadata for {uniprot_id}...")
    predictions = get_prediction_metadata(uniprot_id)
    
    if not predictions:
        print(f"No predictions found for {uniprot_id}")
        return None
    
    print(f"Found {len(predictions)} prediction(s)")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    results = {
        'uniprot_id': uniprot_id,
        'predictions': []
    }
    
    for pred in predictions:
        alphafold_id = pred.get('entryId')
        print(f"\nProcessing: {alphafold_id}")
        print(f"  UniProt ID: {pred.get('uniprotAccession')}")
        print(f"  Protein: {pred.get('gene') or 'N/A'}")
        print(f"  Organism: {pred.get('organismScientificName')}")
        print(f"  Confidence: {pred.get('confidenceScore', 'N/A')}")
        
        pred_result = {
            'alphafold_id': alphafold_id,
            'structure': None,
            'confidence': None,
            'pae': None
        }
        
        # Download structure
        structure_path = download_structure(alphafold_id, output_dir, version, format)
        pred_result['structure'] = structure_path
        
        # Download confidence scores
        if include_confidence:
            confidence_path = download_confidence(alphafold_id, output_dir, version)
            pred_result['confidence'] = confidence_path
        
        # Download PAE
        if include_pae:
            pae_path = download_pae(alphafold_id, output_dir, version)
            pred_result['pae'] = pae_path
        
        results['predictions'].append(pred_result)
    
    return results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Download AlphaFold protein structures"
    )
    parser.add_argument(
        "--uniprot", "-u",
        required=True,
        help="UniProt accession ID (e.g., P00520)"
    )
    parser.add_argument(
        "--output", "-o",
        default="./structures",
        help="Output directory (default: ./structures)"
    )
    parser.add_argument(
        "--version", "-v",
        default="v4",
        help="AlphaFold version (default: v4)"
    )
    parser.add_argument(
        "--format", "-f",
        default="cif",
        choices=["cif", "pdb"],
        help="Structure file format (default: cif)"
    )
    parser.add_argument(
        "--no-confidence",
        action="store_true",
        help="Skip downloading confidence scores"
    )
    parser.add_argument(
        "--no-pae",
        action="store_true",
        help="Skip downloading PAE matrix"
    )
    
    args = parser.parse_args()
    
    print("="*60)
    print("AlphaFold Structure Downloader")
    print("="*60)
    
    results = download_all(
        args.uniprot,
        args.output,
        version=args.version,
        format=args.format,
        include_confidence=not args.no_confidence,
        include_pae=not args.no_pae
    )
    
    if results:
        print("\n" + "="*60)
        print("Download Summary")
        print("="*60)
        for pred in results['predictions']:
            print(f"\n{pred['alphafold_id']}:")
            print(f"  Structure: {pred['structure'] or 'Failed'}")
            if not args.no_confidence:
                print(f"  Confidence: {pred['confidence'] or 'Failed'}")
            if not args.no_pae:
                print(f"  PAE: {pred['pae'] or 'Failed'}")
        print("\nDownload complete!")
    else:
        print("\nDownload failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
