#!/usr/bin/env python3
"""
Scientific Writing Tools
Assist with academic writing, structure, and style.
"""

import re
from typing import List, Dict, Tuple
import json


def analyze_sentence_structure(text: str) -> Dict:
    """Analyze sentence structure for readability."""
    sentences = re.split(r'[.!?]+', text)
    sentences = [s.strip() for s in sentences if s.strip()]
    
    lengths = [len(s.split()) for s in sentences]
    avg_length = sum(lengths) / len(lengths) if lengths else 0
    
    long_sentences = sum(1 for l in lengths if l > 25)
    short_sentences = sum(1 for l in lengths if l < 10)
    
    return {
        "total_sentences": len(sentences),
        "average_length": round(avg_length, 1),
        "long_sentences": long_sentences,
        "short_sentences": short_sentences,
        "readability": "good" if 15 <= avg_length <= 20 else "needs review"
    }


def check_passive_voice(text: str) -> List[Dict]:
    """Identify passive voice constructions."""
    passive_patterns = [
        r'\b\w+\s+(was|were|is|are|been|be|being)\s+\w+ed\b',
        r'\b\w+\s+(has|have|had)\s+been\s+\w+ed\b'
    ]
    
    findings = []
    for pattern in passive_patterns:
        matches = re.finditer(pattern, text, re.IGNORECASE)
        for match in matches:
            findings.append({
                "text": match.group(),
                "position": match.start(),
                "suggestion": "Consider using active voice"
            })
    
    return findings


def generate_section_outline(section_type: str, topic: str) -> List[str]:
    """Generate outline for paper section."""
    outlines = {
        "introduction": [
            f"Background on {topic}",
            f"Current state of research in {topic}",
            f"Knowledge gaps in {topic}",
            "Research objectives",
            "Significance of the study"
        ],
        "methods": [
            "Study design",
            "Participants/Materials",
            "Data collection procedures",
            "Analysis methods",
            "Ethical considerations"
        ],
        "results": [
            "Participant characteristics",
            "Primary outcomes",
            "Secondary outcomes",
            "Statistical significance",
            "Effect sizes"
        ],
        "discussion": [
            "Summary of key findings",
            "Comparison with prior research",
            "Theoretical implications",
            "Practical applications",
            "Limitations",
            "Future directions",
            "Conclusion"
        ]
    }
    
    return outlines.get(section_type, ["Introduction", "Methods", "Results", "Discussion"])


def suggest_transition_words(current_section: str, next_section: str) -> List[str]:
    """Suggest transition phrases between sections."""
    transitions = {
        ("introduction", "methods"): [
            "To address these objectives, we...",
            "The following methods were employed...",
            "We conducted..."
        ],
        ("methods", "results"): [
            "The analysis revealed...",
            "Our findings indicate...",
            "We observed..."
        ],
        ("results", "discussion"): [
            "These results suggest...",
            "The implications of these findings...",
            "Several explanations may account for..."
        ]
    }
    
    key = (current_section.lower(), next_section.lower())
    return transitions.get(key, ["Furthermore...", "In addition...", "Moreover..."])


def estimate_reading_time(text: str) -> Dict:
    """Estimate reading time for text."""
    words = len(text.split())
    minutes = words / 200  # Average reading speed
    
    return {
        "word_count": words,
        "reading_time_minutes": round(minutes, 1),
        "difficulty": "complex" if words > 5000 else "moderate" if words > 2000 else "brief"
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scientific Writing Tools")
    parser.add_argument("command", 
                       choices=["analyze", "passive", "outline", "transition", "reading_time"])
    parser.add_argument("--text", help="Text to analyze")
    parser.add_argument("--section", help="Section type")
    parser.add_argument("--topic", help="Paper topic")
    parser.add_argument("--current", help="Current section")
    parser.add_argument("--next", help="Next section")
    
    args = parser.parse_args()
    
    if args.command == "analyze":
        if not args.text:
            print("Error: --text required")
            return
        result = analyze_sentence_structure(args.text)
        print(json.dumps(result, indent=2))
    
    elif args.command == "passive":
        if not args.text:
            print("Error: --text required")
            return
        findings = check_passive_voice(args.text)
        print(json.dumps(findings, indent=2))
    
    elif args.command == "outline":
        if not args.section or not args.topic:
            print("Error: --section and --topic required")
            return
        outline = generate_section_outline(args.section, args.topic)
        for i, item in enumerate(outline, 1):
            print(f"{i}. {item}")
    
    elif args.command == "transition":
        if not args.current or not args.next:
            print("Error: --current and --next required")
            return
        transitions = suggest_transition_words(args.current, args.next)
        for t in transitions:
            print(f"- {t}")
    
    elif args.command == "reading_time":
        if not args.text:
            print("Error: --text required")
            return
        result = estimate_reading_time(args.text)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
