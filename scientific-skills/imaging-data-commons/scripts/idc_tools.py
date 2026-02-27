#!/usr/bin/env python3
"""
Imaging Data Commons (IDC) Tools
Query and download public cancer imaging data from NCI IDC.
"""

import os
import json
from typing import List, Dict, Optional


def query_collections() -> List[Dict]:
    """Query available IDC collections."""
    try:
        from idc_index import IDCClient
        client = IDCClient()
        
        query = """
        SELECT DISTINCT collection_id, 
               COUNT(DISTINCT PatientID) as patient_count,
               COUNT(DISTINCT StudyInstanceUID) as study_count
        FROM index 
        GROUP BY collection_id
        ORDER BY patient_count DESC
        """
        return client.sql_query(query).to_dict('records')
    except ImportError:
        print("Error: idc-index not installed. Run: pip install idc-index")
        return []


def search_studies(collection_id: str = None, modality: str = None, 
                   body_part: str = None, limit: int = 10) -> List[Dict]:
    """Search imaging studies by criteria."""
    try:
        from idc_index import IDCClient
        client = IDCClient()
        
        conditions = []
        if collection_id:
            conditions.append(f"collection_id = '{collection_id}'")
        if modality:
            conditions.append(f"Modality = '{modality}'")
        if body_part:
            conditions.append(f"BodyPartExamined LIKE '%{body_part}%'")
        
        where_clause = " AND ".join(conditions) if conditions else "1=1"
        
        query = f"""
        SELECT collection_id, PatientID, StudyInstanceUID, 
               SeriesInstanceUID, Modality, BodyPartExamined,
               series_size_MB
        FROM index 
        WHERE {where_clause}
        LIMIT {limit}
        """
        return client.sql_query(query).to_dict('records')
    except ImportError:
        print("Error: idc-index not installed")
        return []


def download_series(series_uid: str, download_dir: str = "./idc_downloads"):
    """Download a specific imaging series."""
    try:
        from idc_index import IDCClient
        client = IDCClient()
        
        os.makedirs(download_dir, exist_ok=True)
        
        client.download_from_selection(
            seriesInstanceUID=series_uid,
            downloadDir=download_dir
        )
        print(f"Downloaded to: {download_dir}")
        return True
    except Exception as e:
        print(f"Download failed: {e}")
        return False


def get_viewer_url(series_uid: str) -> str:
    """Get OHIF viewer URL for a series."""
    try:
        from idc_index import IDCClient
        client = IDCClient()
        return client.get_viewer_URL(seriesInstanceUID=series_uid)
    except Exception as e:
        print(f"Error: {e}")
        return ""


def main():
    import argparse
    parser = argparse.ArgumentParser(description="IDC Imaging Data Commons Tools")
    parser.add_argument("command", choices=["collections", "search", "download", "viewer"])
    parser.add_argument("--collection", help="Filter by collection")
    parser.add_argument("--modality", help="Filter by modality (CT, MR, PT, etc)")
    parser.add_argument("--body-part", help="Filter by body part")
    parser.add_argument("--series-uid", help="Series instance UID")
    parser.add_argument("--output", default="./idc_downloads", help="Download directory")
    parser.add_argument("--limit", type=int, default=10, help="Result limit")
    
    args = parser.parse_args()
    
    if args.command == "collections":
        cols = query_collections()
        for c in cols[:args.limit]:
            print(f"{c['collection_id']}: {c['patient_count']} patients, {c['study_count']} studies")
    
    elif args.command == "search":
        results = search_studies(args.collection, args.modality, args.body_part, args.limit)
        for r in results:
            print(f"{r['SeriesInstanceUID'][:20]}... | {r['Modality']} | {r['BodyPartExamined']} | {r['series_size_MB']:.1f} MB")
    
    elif args.command == "download":
        if not args.series_uid:
            print("Error: --series-uid required")
            return
        download_series(args.series_uid, args.output)
    
    elif args.command == "viewer":
        if not args.series_uid:
            print("Error: --series-uid required")
            return
        url = get_viewer_url(args.series_uid)
        print(f"Viewer URL: {url}")


if __name__ == "__main__":
    main()
