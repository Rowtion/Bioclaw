#!/usr/bin/env python3
"""
PDB Database Tools - Protein structure search and download utilities
Uses RCSB PDB API for accessing 3D protein/nucleic acid structures
"""

import os
import sys
import argparse
import json
import requests
from typing import List, Dict, Optional
from pathlib import Path

BASE_URL = "https://data.rcsb.org/rest/v1"
DOWNLOAD_URL = "https://files.rcsb.org/download"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def search_structures(query: str, max_results: int = 10) -> List[Dict]:
    """Search PDB structures by text query."""
    search_payload = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {"value": query}
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_results
            }
        }
    }
    
    try:
        response = requests.post(
            SEARCH_URL,
            json=search_payload,
            headers={"Content-Type": "application/json"}
        )
        response.raise_for_status()
        data = response.json()
        
        results = []
        if "result_set" in data:
            for entry in data["result_set"]:
                pdb_id = entry.get("identifier", "").upper()
                results.append({
                    "pdb_id": pdb_id,
                    "score": entry.get("score", 0)
                })
        return results
    except Exception as e:
        print(f"Error searching structures: {e}")
        return []


def get_entry_info(pdb_id: str) -> Dict:
    """Get detailed information about a PDB entry."""
    pdb_id = pdb_id.upper()
    try:
        # Get entry info
        entry_url = f"{BASE_URL}/core/entry/{pdb_id}"
        response = requests.get(entry_url)
        response.raise_for_status()
        entry_data = response.json()
        
        # Get structural info
        struct_url = f"{BASE_URL}/core/pubmed/{pdb_id}"
        pubmed_data = {}
        try:
            pm_response = requests.get(struct_url)
            if pm_response.status_code == 200:
                pubmed_data = pm_response.json()
        except:
            pass
        
        return {
            "pdb_id": pdb_id,
            "title": entry_data.get("struct", {}).get("title", "N/A"),
            "deposition_date": entry_data.get("rcsb_accession_info", {}).get("deposit_date", "N/A"),
            "method": entry_data.get("exptl", [{}])[0].get("method", "N/A") if entry_data.get("exptl") else "N/A",
            "resolution": entry_data.get("rcsb_entry_info", {}).get("resolution_combined", "N/A"),
            "organism": entry_data.get("rcsb_entity_source_organism", [{}])[0].get("scientific_name", "N/A") if entry_data.get("rcsb_entity_source_organism") else "N/A",
            "pubmed": pubmed_data
        }
    except Exception as e:
        print(f"Error fetching entry info for {pdb_id}: {e}")
        return {"pdb_id": pdb_id, "error": str(e)}


def download_structure(pdb_id: str, format_type: str = "cif", output_dir: str = ".") -> str:
    """Download structure file in specified format."""
    pdb_id = pdb_id.upper()
    format_extensions = {
        "pdb": ".pdb",
        "cif": ".cif",
        "mmCIF": ".cif",
        "bcif": ".bcif"
    }
    
    ext = format_extensions.get(format_type, ".cif")
    url = f"{DOWNLOAD_URL}/{pdb_id}{ext}"
    
    output_path = Path(output_dir) / f"{pdb_id}{ext}"
    
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        print(f"Downloaded: {output_path}")
        return str(output_path)
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        return ""


def batch_download(pdb_ids: List[str], format_type: str = "cif", output_dir: str = "."):
    """Download multiple structures."""
    os.makedirs(output_dir, exist_ok=True)
    results = []
    
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.strip().upper()
        if len(pdb_id) == 4:
            path = download_structure(pdb_id, format_type, output_dir)
            results.append({"pdb_id": pdb_id, "path": path})
        else:
            print(f"Skipping invalid PDB ID: {pdb_id}")
    
    return results


def get_sequence(pdb_id: str, entity_id: str = "1") -> str:
    """Get protein sequence for a PDB entity."""
    pdb_id = pdb_id.upper()
    try:
        url = f"{BASE_URL}/core/polymer_entity/{pdb_id}_{entity_id}"
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()
        
        sequence = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code", "")
        return sequence
    except Exception as e:
        print(f"Error fetching sequence: {e}")
        return ""


def search_by_sequence(sequence: str, identity_cutoff: float = 0.9, evalue: float = 0.1) -> List[Dict]:
    """Search PDB by sequence similarity using BLAST."""
    from rcsbapi.search import SequenceQuery
    
    try:
        query = SequenceQuery(
            value=sequence,
            evalue_cutoff=evalue,
            identity_cutoff=identity_cutoff
        )
        results = list(query())
        return [{"pdb_id": r.upper(), "method": "sequence_similarity"} for r in results]
    except Exception as e:
        print(f"Sequence search requires rcsb-api package: {e}")
        return []


def main():
    parser = argparse.ArgumentParser(description="PDB Database Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Search command
    search_parser = subparsers.add_parser("search", help="Search PDB structures")
    search_parser.add_argument("query", help="Search query")
    search_parser.add_argument("-n", "--max-results", type=int, default=10, help="Maximum results")
    
    # Info command
    info_parser = subparsers.add_parser("info", help="Get structure information")
    info_parser.add_argument("pdb_id", help="PDB ID (e.g., 4HHB)")
    
    # Download command
    download_parser = subparsers.add_parser("download", help="Download structure file")
    download_parser.add_argument("pdb_id", help="PDB ID")
    download_parser.add_argument("-f", "--format", choices=["pdb", "cif", "bcif"], default="cif")
    download_parser.add_argument("-o", "--output", default=".", help="Output directory")
    
    # Batch download command
    batch_parser = subparsers.add_parser("batch-download", help="Download multiple structures")
    batch_parser.add_argument("pdb_ids", help="Comma-separated PDB IDs")
    batch_parser.add_argument("-f", "--format", choices=["pdb", "cif", "bcif"], default="cif")
    batch_parser.add_argument("-o", "--output", default="pdb_structures", help="Output directory")
    
    # Sequence command
    seq_parser = subparsers.add_parser("sequence", help="Get protein sequence")
    seq_parser.add_argument("pdb_id", help="PDB ID")
    seq_parser.add_argument("-e", "--entity", default="1", help="Entity ID")
    
    args = parser.parse_args()
    
    if args.command == "search":
        print(f"Searching for: {args.query}")
        results = search_structures(args.query, args.max_results)
        print(f"\nFound {len(results)} structures:")
        for r in results:
            print(f"  {r['pdb_id']} (score: {r['score']:.3f})")
    
    elif args.command == "info":
        info = get_entry_info(args.pdb_id)
        print(json.dumps(info, indent=2))
    
    elif args.command == "download":
        download_structure(args.pdb_id, args.format, args.output)
    
    elif args.command == "batch-download":
        pdb_ids = [p.strip() for p in args.pdb_ids.split(",")]
        batch_download(pdb_ids, args.format, args.output)
    
    elif args.command == "sequence":
        seq = get_sequence(args.pdb_id, args.entity)
        if seq:
            print(f">{args.pdb_id}_{args.entity}")
            print(seq)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
