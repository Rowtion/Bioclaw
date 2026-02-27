#!/usr/bin/env python3
"""
PPTX Posters Tools - HTML/CSS poster creation and export utilities
Creates research posters that can be exported to PDF or PPTX
"""

import argparse
import json
from pathlib import Path
from typing import List, Dict, Optional


HTML_TEMPLATE = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        :root {{
            --primary-color: {primary_color};
            --secondary-color: {secondary_color};
            --accent-color: {accent_color};
            --text-color: #1a202c;
            --bg-color: #ffffff;
        }}
        
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', 'Arial', sans-serif;
            background: #f0f0f0;
            padding: 20px;
        }}
        
        .poster {{
            width: {width}px;
            height: {height}px;
            margin: 0 auto;
            background: var(--bg-color);
            box-shadow: 0 0 50px rgba(0,0,0,0.2);
            position: relative;
            overflow: hidden;
        }}
        
        /* Header */
        .header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: white;
            padding: 40px 60px;
            text-align: center;
        }}
        
        .poster-title {{
            font-size: {title_size}pt;
            font-weight: bold;
            margin-bottom: 20px;
            line-height: 1.2;
        }}
        
        .authors {{
            font-size: {author_size}pt;
            margin-bottom: 10px;
        }}
        
        .affiliations {{
            font-size: {affil_size}pt;
            opacity: 0.9;
        }}
        
        /* Content */
        .content {{
            display: grid;
            grid-template-columns: repeat({columns}, 1fr);
            gap: 30px;
            padding: 40px;
            height: calc(100% - 200px);
        }}
        
        .column {{
            display: flex;
            flex-direction: column;
            gap: 30px;
        }}
        
        .block {{
            background: #f8fafc;
            border-radius: 12px;
            padding: 24px;
            border-left: 5px solid var(--accent-color);
        }}
        
        .block-title {{
            font-size: {section_size}pt;
            color: var(--primary-color);
            font-weight: bold;
            margin-bottom: 16px;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .block-content {{
            font-size: {body_size}pt;
            line-height: 1.6;
            color: var(--text-color);
        }}
        
        .block-content ul {{
            margin-left: 24px;
            margin-top: 12px;
        }}
        
        .block-content li {{
            margin-bottom: 8px;
        }}
        
        .block-image {{
            width: 100%;
            border-radius: 8px;
            margin: 16px 0;
        }}
        
        /* Footer */
        .footer {{
            background: var(--primary-color);
            color: white;
            padding: 20px 40px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            font-size: 16pt;
        }}
        
        .footer-left {{
            opacity: 0.8;
        }}
        
        .footer-right {{
            text-align: right;
        }}
        
        @media print {{
            body {{
                background: white;
                padding: 0;
            }}
            .poster {{
                box-shadow: none;
            }}
        }}
    </style>
</head>
<body>
    <div class="poster">
        <header class="header">
            <h1 class="poster-title">{title}</h1>
            <div class="authors">{authors}</div>
            <div class="affiliations">{affiliations}</div>
        </header>
        
        <div class="content">
            {columns_html}
        </div>
        
        <footer class="footer">
            <div class="footer-left">
                {conference} | {date}
            </div>
            <div class="footer-right">
                {contact}
            </div>
        </footer>
    </div>
</body>
</html>
'''


def generate_poster_html(config: Dict) -> str:
    """Generate poster HTML from configuration."""
    # Generate columns HTML
    columns_html = ""
    for col_idx, column in enumerate(config.get("columns", [])):
        blocks_html = ""
        for block in column.get("blocks", []):
            image_html = ""
            if block.get("image"):
                image_html = f'<img src="{block["image"]}" class="block-image" alt="{block["title"]}">'
            
            content_html = block.get("content", "")
            if isinstance(content_html, list):
                content_html = "<ul>" + "".join([f"<li>{item}</li>" for item in content_html]) + "</ul>"
            else:
                content_html = f"<p>{content_html}</p>"
            
            blocks_html += f'''
                <div class="block">
                    <h2 class="block-title">{block.get("title", "")}</h2>
                    <div class="block-content">
                        {image_html}
                        {content_html}
                    </div>
                </div>
            '''
        
        columns_html += f'<div class="column">{blocks_html}</div>'
    
    # Fill template
    html = HTML_TEMPLATE.format(
        title=config.get("title", "Research Poster"),
        authors=config.get("authors", ""),
        affiliations=config.get("affiliations", ""),
        conference=config.get("conference", ""),
        date=config.get("date", ""),
        contact=config.get("contact", ""),
        columns_html=columns_html,
        columns=len(config.get("columns", 3)),
        width=config.get("width", 2592),
        height=config.get("height", 3456),
        title_size=config.get("title_size", 72),
        author_size=config.get("author_size", 36),
        affil_size=config.get("affil_size", 28),
        section_size=config.get("section_size", 32),
        body_size=config.get("body_size", 22),
        primary_color=config.get("primary_color", "#1a365d"),
        secondary_color=config.get("secondary_color", "#2b6cb0"),
        accent_color=config.get("accent_color", "#3182ce")
    )
    
    return html


def create_poster_config(title: str, output_path: str):
    """Create a poster configuration template."""
    config = {
        "title": title,
        "authors": "Author 1, Author 2, Author 3*",
        "affiliations": "Department, Institution, City, Country",
        "conference": "Conference Name 2025",
        "date": "January 2025",
        "contact": "email@institution.edu | @handle",
        "width": 2592,
        "height": 3456,
        "title_size": 72,
        "author_size": 36,
        "affil_size": 28,
        "section_size": 32,
        "body_size": 22,
        "primary_color": "#1a365d",
        "secondary_color": "#2b6cb0",
        "accent_color": "#3182ce",
        "columns": [
            {
                "blocks": [
                    {
                        "title": "Introduction",
                        "content": [
                            "Background and motivation",
                            "Research question",
                            "Study objectives"
                        ]
                    },
                    {
                        "title": "Methods",
                        "content": [
                            "Study design",
                            "Data collection",
                            "Analysis approach"
                        ],
                        "image": "figures/methods.png"
                    }
                ]
            },
            {
                "blocks": [
                    {
                        "title": "Results",
                        "content": "Key findings from the study",
                        "image": "figures/results.png"
                    },
                    {
                        "title": "Discussion",
                        "content": [
                            "Interpretation of results",
                            "Comparison with prior work",
                            "Implications"
                        ]
                    }
                ]
            },
            {
                "blocks": [
                    {
                        "title": "Conclusions",
                        "content": [
                            "Main conclusion 1",
                            "Main conclusion 2",
                            "Future directions"
                        ]
                    },
                    {
                        "title": "Acknowledgments",
                        "content": "Funding sources and acknowledgments"
                    },
                    {
                        "title": "References",
                        "content": "1. Author et al. (2024). Journal Name."
                    }
                ]
            }
        ]
    }
    
    with open(output_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"Poster config template created: {output_path}")


def build_poster(config_path: str, output_path: str):
    """Build poster HTML from config file."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    html = generate_poster_html(config)
    
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"Poster HTML generated: {output_path}")


def validate_poster_config(config_path: str) -> bool:
    """Validate poster configuration."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    required_fields = ["title", "authors", "columns"]
    errors = []
    
    for field in required_fields:
        if field not in config:
            errors.append(f"Missing required field: {field}")
    
    if "columns" in config:
        if len(config["columns"]) < 1 or len(config["columns"]) > 4:
            errors.append("Number of columns should be 1-4")
    
    if errors:
        print("Validation errors:")
        for e in errors:
            print(f"  - {e}")
        return False
    else:
        print("Configuration is valid")
        return True


def export_pdf(input_html: str, output_pdf: str):
    """Export HTML poster to PDF using Chrome headless."""
    import subprocess
    
    cmd = [
        "google-chrome",
        "--headless",
        "--disable-gpu",
        f"--print-to-pdf={output_pdf}",
        "--print-to-pdf-no-header",
        "--no-margins",
        input_html
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"PDF exported: {output_pdf}")
            return True
        else:
            print(f"Export failed: {result.stderr}")
            return False
    except FileNotFoundError:
        print("Chrome not found. Install Chrome or use browser print to PDF.")
        return False


def create_simple_poster(title: str, sections: List[str], output: str):
    """Create a simple poster with minimal configuration."""
    config = {
        "title": title,
        "authors": "Your Name*",
        "affiliations": "Your Institution",
        "conference": "Conference 2025",
        "date": "2025",
        "contact": "your.email@institution.edu",
        "columns": []
    }
    
    # Distribute sections across 3 columns
    sections_per_col = (len(sections) + 2) // 3
    for i in range(3):
        col_sections = sections[i * sections_per_col:(i + 1) * sections_per_col]
        column = {"blocks": []}
        for section in col_sections:
            column["blocks"].append({
                "title": section,
                "content": f"Content for {section}"
            })
        config["columns"].append(column)
    
    html = generate_poster_html(config)
    
    with open(output, 'w') as f:
        f.write(html)
    
    print(f"Simple poster created: {output}")


def main():
    parser = argparse.ArgumentParser(description="PPTX Posters Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Create template command
    template_parser = subparsers.add_parser("create-template", help="Create poster config template")
    template_parser.add_argument("-t", "--title", required=True, help="Poster title")
    template_parser.add_argument("-o", "--output", default="poster_config.json", help="Output JSON file")
    
    # Build command
    build_parser = subparsers.add_parser("build", help="Build poster HTML from config")
    build_parser.add_argument("config", help="Config JSON file")
    build_parser.add_argument("-o", "--output", default="poster.html", help="Output HTML file")
    
    # Validate command
    validate_parser = subparsers.add_parser("validate", help="Validate poster config")
    validate_parser.add_argument("config", help="Config JSON file")
    
    # Export PDF command
    pdf_parser = subparsers.add_parser("export-pdf", help="Export to PDF (requires Chrome)")
    pdf_parser.add_argument("html", help="Input HTML file")
    pdf_parser.add_argument("-o", "--output", default="poster.pdf", help="Output PDF file")
    
    # Simple poster command
    simple_parser = subparsers.add_parser("simple", help="Create simple poster")
    simple_parser.add_argument("-t", "--title", required=True, help="Poster title")
    simple_parser.add_argument("-s", "--sections", required=True, help="Comma-separated section names")
    simple_parser.add_argument("-o", "--output", default="poster.html", help="Output HTML file")
    
    args = parser.parse_args()
    
    if args.command == "create-template":
        create_poster_config(args.title, args.output)
    
    elif args.command == "build":
        build_poster(args.config, args.output)
    
    elif args.command == "validate":
        validate_poster_config(args.config)
    
    elif args.command == "export-pdf":
        export_pdf(args.html, args.output)
    
    elif args.command == "simple":
        sections = [s.strip() for s in args.sections.split(",")]
        create_simple_poster(args.title, sections, args.output)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
