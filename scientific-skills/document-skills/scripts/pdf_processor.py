"""
Document Processing Example: PDF Text Extraction
Extract and process text from PDF documents.
"""
import PyPDF2
import re
from pathlib import Path


def extract_text_from_pdf(pdf_path: str) -> str:
    """
    Extract text from a PDF file.
    
    Args:
        pdf_path: Path to PDF file
        
    Returns:
        Extracted text
    """
    text = ""
    with open(pdf_path, 'rb') as file:
        pdf_reader = PyPDF2.PdfReader(file)
        for page in pdf_reader.pages:
            text += page.extract_text() + "\n"
    return text


def extract_sections(text: str, section_headers: list) -> dict:
    """
    Extract specific sections from document text.
    
    Args:
        text: Full document text
        section_headers: List of section header patterns
        
    Returns:
        Dictionary of section names to content
    """
    sections = {}
    current_section = "introduction"
    current_content = []
    
    lines = text.split('\n')
    
    for line in lines:
        line = line.strip()
        # Check if line is a section header
        for header in section_headers:
            if re.match(rf'^{header}\s*[:\.]?', line, re.IGNORECASE):
                sections[current_section] = '\n'.join(current_content)
                current_section = header.lower()
                current_content = []
                break
        else:
            current_content.append(line)
    
    sections[current_section] = '\n'.join(current_content)
    return sections


def batch_process_pdfs(directory: str, output_dir: str):
    """
    Process all PDFs in a directory.
    
    Args:
        directory: Input directory containing PDFs
        output_dir: Output directory for extracted text
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    for pdf_file in Path(directory).glob('*.pdf'):
        print(f"Processing: {pdf_file.name}")
        text = extract_text_from_pdf(str(pdf_file))
        
        output_file = Path(output_dir) / f"{pdf_file.stem}.txt"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(text)
        
        print(f"  Extracted {len(text)} characters")


if __name__ == "__main__":
    # Example usage
    # text = extract_text_from_pdf('document.pdf')
    # sections = extract_sections(text, ['Abstract', 'Introduction', 'Methods', 'Results'])
    print("Document processing utilities loaded.")
