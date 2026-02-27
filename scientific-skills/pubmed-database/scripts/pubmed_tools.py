#!/usr/bin/env python3
"""
PubMed Database Tools - E-utilities API access for literature search
Advanced query construction, batch processing, and citation management
"""

import os
import time
import argparse
import json
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional
from urllib.parse import quote_plus
import requests

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def esearch(query: str, db: str = "pubmed", retmax: int = 20,
           retstart: int = 0, use_history: bool = False,
           api_key: Optional[str] = None) -> Dict:
    """Search PubMed and retrieve PMIDs."""
    params = {
        "db": db,
        "term": query,
        "retmax": retmax,
        "retstart": retstart,
        "retmode": "json",
        "usehistory": "y" if use_history else "n"
    }
    
    if api_key:
        params["api_key"] = api_key
    
    url = f"{BASE_URL}/esearch.fcgi"
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error in esearch: {e}")
        return {}


def efetch(pmids: List[str], db: str = "pubmed", rettype: str = "abstract",
          retmode: str = "text", api_key: Optional[str] = None) -> str:
    """Fetch records by PMID."""
    params = {
        "db": db,
        "id": ",".join(pmids),
        "rettype": rettype,
        "retmode": retmode
    }
    
    if api_key:
        params["api_key"] = api_key
    
    url = f"{BASE_URL}/efetch.fcgi"
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.text
    except Exception as e:
        print(f"Error in efetch: {e}")
        return ""


def esummary(pmids: List[str], db: str = "pubmed", api_key: Optional[str] = None) -> Dict:
    """Get document summaries."""
    params = {
        "db": db,
        "id": ",".join(pmids),
        "retmode": "json"
    }
    
    if api_key:
        params["api_key"] = api_key
    
    url = f"{BASE_URL}/esummary.fcgi"
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error in esummary: {e}")
        return {}


def elink(pmids: List[str], dbfrom: str = "pubmed", db: str = "pubmed",
         cmd: str = "neighbor", api_key: Optional[str] = None) -> Dict:
    """Find related articles."""
    params = {
        "dbfrom": dbfrom,
        "db": db,
        "id": ",".join(pmids),
        "cmd": cmd,
        "retmode": "json"
    }
    
    if api_key:
        params["api_key"] = api_key
    
    url = f"{BASE_URL}/elink.fcgi"
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error in elink: {e}")
        return {}


def search_and_fetch(query: str, max_results: int = 10, 
                    api_key: Optional[str] = None) -> List[Dict]:
    """Search and fetch article details in one operation."""
    # Search
    search_result = esearch(query, retmax=max_results, api_key=api_key)
    
    pmids = search_result.get("esearchresult", {}).get("idlist", [])
    if not pmids:
        return []
    
    # Get summaries
    summary_result = esummary(pmids, api_key=api_key)
    
    articles = []
    result_data = summary_result.get("result", {})
    
    for pmid in pmids:
        if pmid in result_data:
            item = result_data[pmid]
            articles.append({
                "pmid": pmid,
                "title": item.get("title", ""),
                "authors": [a.get("name", "") for a in item.get("authors", [])],
                "journal": item.get("fulljournalname", ""),
                "year": item.get("pubdate", "")[:4] if item.get("pubdate") else "",
                "doi": item.get("elocationid", "").replace("doi: ", "") if item.get("elocationid") else ""
            })
    
    return articles


def build_query(keywords: List[str], author: Optional[str] = None,
               journal: Optional[str] = None, year_from: Optional[int] = None,
               year_to: Optional[int] = None, article_type: Optional[str] = None) -> str:
    """Build a complex PubMed query."""
    parts = []
    
    # Keywords
    if keywords:
        keyword_query = " AND ".join([f'"{k}"[tiab]' for k in keywords])
        parts.append(f"({keyword_query})")
    
    # Author
    if author:
        parts.append(f"{author}[au]")
    
    # Journal
    if journal:
        parts.append(f"{journal}[ta]")
    
    # Date range
    if year_from and year_to:
        parts.append(f"{year_from}:{year_to}[dp]")
    elif year_from:
        parts.append(f"{year_from}[dp]")
    
    # Article type
    if article_type:
        parts.append(f"{article_type}[pt]")
    
    return " AND ".join(parts)


def fetch_abstracts(pmids: List[str], api_key: Optional[str] = None) -> Dict[str, str]:
    """Fetch abstracts for given PMIDs."""
    abstracts = {}
    
    # Process in batches of 100
    for i in range(0, len(pmids), 100):
        batch = pmids[i:i+100]
        text = efetch(batch, rettype="abstract", retmode="text", api_key=api_key)
        
        # Parse abstracts (simple parsing)
        current_pmid = None
        current_abstract = []
        
        for line in text.split("\n"):
            if line.startswith("PMID:"):
                if current_pmid:
                    abstracts[current_pmid] = "\n".join(current_abstract).strip()
                current_pmid = line.replace("PMID:", "").strip()
                current_abstract = []
            elif line.strip():
                current_abstract.append(line.strip())
        
        if current_pmid:
            abstracts[current_pmid] = "\n".join(current_abstract).strip()
        
        # Rate limiting
        time.sleep(0.34 if api_key else 0.34)
    
    return abstracts


def export_to_bibtex(articles: List[Dict], output_path: str):
    """Export articles to BibTeX format."""
    with open(output_path, 'w') as f:
        for article in articles:
            pmid = article.get("pmid", "")
            title = article.get("title", "").replace("{", "").replace("}", "")
            journal = article.get("journal", "")
            year = article.get("year", "")
            
            # Create citation key
            first_author = article.get("authors", ["Unknown"])[0].split()[0] if article.get("authors") else "Unknown"
            cite_key = f"{first_author}{year}{pmid}"
            
            f.write(f"@article{{{cite_key},\n")
            f.write(f"  title = {{{title}}},\n")
            f.write(f"  author = {{{' and '.join(article.get('authors', []))}}},\n")
            f.write(f"  journal = {{{journal}}},\n")
            f.write(f"  year = {{{year}}},\n")
            f.write(f"  pmid = {{{pmid}}},\n")
            if article.get("doi"):
                f.write(f"  doi = {{{article['doi']}}},\n")
            f.write("}\n\n")
    
    print(f"BibTeX exported to {output_path}")


def export_to_csv(articles: List[Dict], output_path: str):
    """Export articles to CSV."""
    import csv
    
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["pmid", "title", "authors", "journal", "year", "doi"])
        writer.writeheader()
        
        for article in articles:
            row = article.copy()
            row["authors"] = "; ".join(row.get("authors", []))
            writer.writerow(row)
    
    print(f"CSV exported to {output_path}")


def get_related_articles(pmid: str, api_key: Optional[str] = None) -> List[str]:
    """Get related articles for a PMID."""
    result = elink([pmid], cmd="neighbor", api_key=api_key)
    
    related = []
    linksets = result.get("linksets", [])
    for linkset in linksets:
        for linksetdb in linkset.get("linksetdbs", []):
            if linksetdb.get("linkname") == "pubmed_pubmed":
                related.extend(linksetdb.get("links", []))
    
    return related[:10]  # Return top 10


def systematic_review_search(queries: List[str], api_key: Optional[str] = None) -> Dict:
    """Execute multiple search strategies for systematic review."""
    results = {}
    
    for i, query in enumerate(queries, 1):
        print(f"Executing search {i}/{len(queries)}...")
        search_result = esearch(query, retmax=1000, api_key=api_key)
        
        pmids = search_result.get("esearchresult", {}).get("idlist", [])
        count = search_result.get("esearchresult", {}).get("count", 0)
        
        results[f"query_{i}"] = {
            "query": query,
            "count": count,
            "pmids": pmids
        }
        
        time.sleep(0.34 if api_key else 1.0)
    
    return results


def main():
    parser = argparse.ArgumentParser(description="PubMed Database Tools")
    parser.add_argument("--api-key", help="NCBI API key (or set NCBI_API_KEY env var)")
    
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Search command
    search_parser = subparsers.add_parser("search", help="Search PubMed")
    search_parser.add_argument("query", help="Search query")
    search_parser.add_argument("-n", "--max-results", type=int, default=10)
    search_parser.add_argument("-o", "--output", help="Output JSON file")
    
    # Fetch command
    fetch_parser = subparsers.add_parser("fetch", help="Fetch abstracts")
    fetch_parser.add_argument("pmids", help="Comma-separated PMIDs")
    fetch_parser.add_argument("-o", "--output", help="Output file")
    
    # Build query command
    build_parser = subparsers.add_parser("build-query", help="Build complex query")
    build_parser.add_argument("-k", "--keywords", help="Comma-separated keywords")
    build_parser.add_argument("-a", "--author", help="Author name")
    build_parser.add_argument("-j", "--journal", help="Journal name")
    build_parser.add_argument("--year-from", type=int, help="Start year")
    build_parser.add_argument("--year-to", type=int, help="End year")
    build_parser.add_argument("--type", help="Article type")
    
    # Export command
    export_parser = subparsers.add_parser("export", help="Export to BibTeX or CSV")
    export_parser.add_argument("pmids", help="Comma-separated PMIDs")
    export_parser.add_argument("-f", "--format", choices=["bibtex", "csv"], default="bibtex")
    export_parser.add_argument("-o", "--output", required=True, help="Output file")
    
    # Related command
    related_parser = subparsers.add_parser("related", help="Get related articles")
    related_parser.add_argument("pmid", help="PMID")
    
    args = parser.parse_args()
    
    api_key = args.api_key or os.getenv("NCBI_API_KEY")
    
    if args.command == "search":
        articles = search_and_fetch(args.query, args.max_results, api_key)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(articles, f, indent=2)
            print(f"Results saved to {args.output}")
        else:
            print(json.dumps(articles, indent=2))
    
    elif args.command == "fetch":
        pmids = [p.strip() for p in args.pmids.split(",")]
        abstracts = fetch_abstracts(pmids, api_key)
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(abstracts, f, indent=2)
            print(f"Abstracts saved to {args.output}")
        else:
            for pmid, abstract in abstracts.items():
                print(f"\n{'='*60}")
                print(f"PMID: {pmid}")
                print(f"{'='*60}")
                print(abstract)
    
    elif args.command == "build-query":
        keywords = args.keywords.split(",") if args.keywords else []
        query = build_query(keywords, args.author, args.journal,
                          args.year_from, args.year_to, args.type)
        print(query)
    
    elif args.command == "export":
        pmids = [p.strip() for p in args.pmids.split(",")]
        articles = []
        
        # Get article details
        for i in range(0, len(pmids), 100):
            batch = pmids[i:i+100]
            summary = esummary(batch, api_key=api_key)
            result_data = summary.get("result", {})
            
            for pmid in batch:
                if pmid in result_data:
                    item = result_data[pmid]
                    articles.append({
                        "pmid": pmid,
                        "title": item.get("title", ""),
                        "authors": [a.get("name", "") for a in item.get("authors", [])],
                        "journal": item.get("fulljournalname", ""),
                        "year": item.get("pubdate", "")[:4] if item.get("pubdate") else "",
                        "doi": item.get("elocationid", "").replace("doi: ", "") if item.get("elocationid") else ""
                    })
        
        if args.format == "bibtex":
            export_to_bibtex(articles, args.output)
        else:
            export_to_csv(articles, args.output)
    
    elif args.command == "related":
        related = get_related_articles(args.pmid, api_key)
        print(f"Related articles for PMID {args.pmid}:")
        for pmid in related:
            print(f"  {pmid}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
