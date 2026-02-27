"""
Document Analysis Example: Citation Extraction
Extract and analyze citations from academic documents.
"""
import re
from collections import Counter
from typing import List, Dict


def extract_citations(text: str) -> List[Dict]:
    """
    Extract citations from text using common patterns.
    
    Args:
        text: Document text
        
    Returns:
        List of citation dictionaries
    """
    citations = []
    
    # Pattern for Author-Year citations (Smith, 2020; Smith et al., 2019)
    author_year_pattern = r'([A-Z][a-z]+(?:\s+et\s+al\.)?\s*,\s*\d{4})'
    matches = re.findall(author_year_pattern, text)
    
    for match in matches:
        citations.append({
            'type': 'author-year',
            'citation': match,
            'year': int(re.search(r'\d{4}', match).group())
        })
    
    # Pattern for numbered citations [1], [1, 2], [1-3]
    numbered_pattern = r'\[(\d+(?:[-,]\s*\d+)*)\]'
    matches = re.findall(numbered_pattern, text)
    
    for match in matches:
        citations.append({
            'type': 'numbered',
            'citation': match
        })
    
    return citations


def analyze_citation_patterns(citations: List[Dict]) -> Dict:
    """
    Analyze citation patterns.
    
    Args:
        citations: List of citation dictionaries
        
    Returns:
        Analysis results
    """
    years = [c['year'] for c in citations if 'year' in c]
    
    analysis = {
        'total_citations': len(citations),
        'author_year_citations': len([c for c in citations if c['type'] == 'author-year']),
        'numbered_citations': len([c for c in citations if c['type'] == 'numbered']),
    }
    
    if years:
        analysis['year_range'] = (min(years), max(years))
        analysis['most_common_years'] = Counter(years).most_common(5)
    
    return analysis


def generate_citation_report(text: str) -> str:
    """
    Generate a citation analysis report.
    
    Args:
        text: Document text
        
    Returns:
        Formatted report
    """
    citations = extract_citations(text)
    analysis = analyze_citation_patterns(citations)
    
    report = []
    report.append("Citation Analysis Report")
    report.append("=" * 50)
    report.append(f"Total citations found: {analysis['total_citations']}")
    report.append(f"Author-year citations: {analysis['author_year_citations']}")
    report.append(f"Numbered citations: {analysis['numbered_citations']}")
    
    if 'year_range' in analysis:
        report.append(f"\nYear range: {analysis['year_range'][0]} - {analysis['year_range'][1]}")
        report.append("\nMost cited years:")
        for year, count in analysis['most_common_years']:
            report.append(f"  {year}: {count} citations")
    
    return '\n'.join(report)


if __name__ == "__main__":
    # Example usage
    sample_text = """
    Recent studies (Smith et al., 2019; Jones, 2020) have shown...
    This was confirmed by subsequent research [1, 2, 3].
    According to Wang et al., 2021, these findings suggest...
    """
    
    print(generate_citation_report(sample_text))
