"""
Document Conversion Example: Format Conversion Utilities
Convert between different document formats.
"""
import json
import csv
from pathlib import Path
from typing import List, Dict


def text_to_jsonl(input_file: str, output_file: str, metadata: Dict = None):
    """
    Convert text file to JSONL format with optional metadata.
    
    Args:
        input_file: Input text file
        output_file: Output JSONL file
        metadata: Optional metadata to include
    """
    metadata = metadata or {}
    
    with open(input_file, 'r', encoding='utf-8') as infile, \
         open(output_file, 'w', encoding='utf-8') as outfile:
        
        for i, line in enumerate(infile, 1):
            record = {
                'id': i,
                'text': line.strip(),
                **metadata
            }
            outfile.write(json.dumps(record) + '\n')
    
    print(f"Converted to {output_file}")


def structured_to_csv(data: List[Dict], output_file: str):
    """
    Save structured data to CSV.
    
    Args:
        data: List of dictionaries
        output_file: Output CSV file
    """
    if not data:
        return
    
    fieldnames = data[0].keys()
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)
    
    print(f"Saved {len(data)} records to {output_file}")


def merge_documents(file_list: List[str], output_file: str, separator: str = "\n\n"):
    """
    Merge multiple documents into one.
    
    Args:
        file_list: List of input files
        output_file: Output file
        separator: Separator between documents
    """
    with open(output_file, 'w', encoding='utf-8') as outfile:
        for i, file_path in enumerate(file_list):
            if i > 0:
                outfile.write(separator)
            
            with open(file_path, 'r', encoding='utf-8') as infile:
                outfile.write(infile.read())
    
    print(f"Merged {len(file_list)} documents into {output_file}")


def split_document(input_file: str, output_dir: str, chunks: int = 5):
    """
    Split a document into chunks.
    
    Args:
        input_file: Input file
        output_dir: Output directory
        chunks: Number of chunks
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    chunk_size = len(content) // chunks
    
    for i in range(chunks):
        start = i * chunk_size
        end = start + chunk_size if i < chunks - 1 else len(content)
        
        chunk_file = Path(output_dir) / f"chunk_{i+1}.txt"
        with open(chunk_file, 'w', encoding='utf-8') as f:
            f.write(content[start:end])
    
    print(f"Split into {chunks} chunks in {output_dir}")


if __name__ == "__main__":
    print("Document conversion utilities loaded.")
