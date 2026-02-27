#!/usr/bin/env python3
"""
Peer Review Tools - Structured manuscript and grant review utilities
Provides checklist-based evaluation for scientific peer review
"""

import json
import argparse
from datetime import datetime
from typing import Dict, List
from dataclasses import dataclass, asdict
from pathlib import Path


@dataclass
class ReviewSection:
    """A section of the peer review."""
    title: str
    comments: str
    severity: str = "minor"  # critical, major, minor, suggestion


@dataclass
class PeerReview:
    """Structured peer review document."""
    manuscript_title: str
    authors: str
    review_date: str
    reviewer_name: str
    
    # Summary
    summary: str = ""
    recommendation: str = ""  # accept, minor_revision, major_revision, reject
    
    # Sections
    sections: List[ReviewSection] = None
    
    # Checklist results
    methodology_checklist: Dict = None
    statistics_checklist: Dict = None
    ethics_checklist: Dict = None
    
    def __post_init__(self):
        if self.sections is None:
            self.sections = []
        if self.methodology_checklist is None:
            self.methodology_checklist = {}
        if self.statistics_checklist is None:
            self.statistics_checklist = {}
        if self.ethics_checklist is None:
            self.ethics_checklist = {}


class PeerReviewTemplate:
    """Template for structured peer review."""
    
    METHODOLOGY_CHECKLIST = {
        "reproducibility": "Methods sufficiently detailed for replication",
        "appropriateness": "Methods appropriate for research questions",
        "controls": "Adequate controls included",
        "validation": "Validation approaches documented",
        "replicates": "Biological/technical replicates specified",
        "sample_size": "Sample size justified with power analysis",
        "randomization": "Randomization procedures described",
        "blinding": "Blinding procedures described"
    }
    
    STATISTICS_CHECKLIST = {
        "assumptions": "Statistical assumptions checked and reported",
        "effect_sizes": "Effect sizes reported alongside p-values",
        "multiple_testing": "Multiple testing correction applied",
        "confidence_intervals": "Confidence intervals provided",
        "missing_data": "Missing data handling described",
        "software": "Statistical software and versions specified"
    }
    
    ETHICS_CHECKLIST = {
        "approval": "IRB/IACUC approval documented",
        "consent": "Informed consent described (if human subjects)",
        "privacy": "Patient privacy protections described",
        "conflicts": "Conflicts of interest disclosed",
        "funding": "Funding sources disclosed"
    }
    
    @classmethod
    def create_review(cls, title: str, authors: str, reviewer: str) -> PeerReview:
        """Create a new peer review with initialized checklists."""
        review = PeerReview(
            manuscript_title=title,
            authors=authors,
            review_date=datetime.now().strftime("%Y-%m-%d"),
            reviewer_name=reviewer,
            methodology_checklist={k: None for k in cls.METHODOLOGY_CHECKLIST.keys()},
            statistics_checklist={k: None for k in cls.STATISTICS_CHECKLIST.keys()},
            ethics_checklist={k: None for k in cls.ETHICS_CHECKLIST.keys()}
        )
        return review


def generate_review_report(review: PeerReview) -> str:
    """Generate formatted peer review report."""
    lines = []
    lines.append("=" * 70)
    lines.append("PEER REVIEW REPORT")
    lines.append("=" * 70)
    lines.append(f"Manuscript: {review.manuscript_title}")
    lines.append(f"Authors: {review.authors}")
    lines.append(f"Reviewer: {review.reviewer_name}")
    lines.append(f"Date: {review.review_date}")
    lines.append("")
    
    # Summary
    lines.append("-" * 70)
    lines.append("SUMMARY")
    lines.append("-" * 70)
    lines.append(review.summary)
    lines.append("")
    lines.append(f"RECOMMENDATION: {review.recommendation.upper()}")
    lines.append("")
    
    # Major comments
    major_comments = [s for s in review.sections if s.severity in ["critical", "major"]]
    if major_comments:
        lines.append("-" * 70)
        lines.append("MAJOR COMMENTS")
        lines.append("-" * 70)
        for i, comment in enumerate(major_comments, 1):
            lines.append(f"{i}. [{comment.severity.upper()}] {comment.title}")
            lines.append(f"   {comment.comments}")
            lines.append("")
    
    # Minor comments
    minor_comments = [s for s in review.sections if s.severity in ["minor", "suggestion"]]
    if minor_comments:
        lines.append("-" * 70)
        lines.append("MINOR COMMENTS")
        lines.append("-" * 70)
        for i, comment in enumerate(minor_comments, 1):
            lines.append(f"{i}. [{comment.severity.upper()}] {comment.title}")
            lines.append(f"   {comment.comments}")
            lines.append("")
    
    # Checklists
    lines.append("-" * 70)
    lines.append("METHODLOGY CHECKLIST")
    lines.append("-" * 70)
    for key, value in review.methodology_checklist.items():
        status = "✓" if value else "✗" if value is False else "?"
        description = PeerReviewTemplate.METHODOLOGY_CHECKLIST[key]
        lines.append(f"[{status}] {description}")
    lines.append("")
    
    lines.append("-" * 70)
    lines.append("STATISTICS CHECKLIST")
    lines.append("-" * 70)
    for key, value in review.statistics_checklist.items():
        status = "✓" if value else "✗" if value is False else "?"
        description = PeerReviewTemplate.STATISTICS_CHECKLIST[key]
        lines.append(f"[{status}] {description}")
    lines.append("")
    
    lines.append("-" * 70)
    lines.append("ETHICS CHECKLIST")
    lines.append("-" * 70)
    for key, value in review.ethics_checklist.items():
        status = "✓" if value else "✗" if value is False else "?"
        description = PeerReviewTemplate.ETHICS_CHECKLIST[key]
        lines.append(f"[{status}] {description}")
    lines.append("")
    
    return "\n".join(lines)


def create_review_template(output_path: str, title: str = "", authors: str = "", reviewer: str = ""):
    """Create a new review template file."""
    review = PeerReviewTemplate.create_review(title, authors, reviewer)
    
    # Add template sections
    review.sections = [
        ReviewSection("Abstract and Title", "[Add comments]", "minor"),
        ReviewSection("Introduction", "[Add comments]", "minor"),
        ReviewSection("Methods", "[Add comments]", "major"),
        ReviewSection("Results", "[Add comments]", "major"),
        ReviewSection("Discussion", "[Add comments]", "major"),
        ReviewSection("Figures and Tables", "[Add comments]", "minor"),
        ReviewSection("References", "[Add comments]", "minor"),
    ]
    
    data = {
        "manuscript_title": review.manuscript_title,
        "authors": review.authors,
        "review_date": review.review_date,
        "reviewer_name": review.reviewer_name,
        "summary": "[Provide brief summary of manuscript and overall assessment]",
        "recommendation": "[accept/minor_revision/major_revision/reject]",
        "sections": [asdict(s) for s in review.sections],
        "methodology_checklist": review.methodology_checklist,
        "statistics_checklist": review.statistics_checklist,
        "ethics_checklist": review.ethics_checklist
    }
    
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)
    
    print(f"Review template created: {output_path}")


def load_review(json_path: str) -> PeerReview:
    """Load a review from JSON file."""
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    review = PeerReview(
        manuscript_title=data.get("manuscript_title", ""),
        authors=data.get("authors", ""),
        review_date=data.get("review_date", ""),
        reviewer_name=data.get("reviewer_name", ""),
        summary=data.get("summary", ""),
        recommendation=data.get("recommendation", ""),
        sections=[ReviewSection(**s) for s in data.get("sections", [])],
        methodology_checklist=data.get("methodology_checklist", {}),
        statistics_checklist=data.get("statistics_checklist", {}),
        ethics_checklist=data.get("ethics_checklist", {})
    )
    return review


def evaluate_methodology(description: str) -> Dict[str, bool]:
    """Auto-evaluate methodology section (basic heuristic)."""
    keywords = {
        "reproducibility": ["protocol", "procedure", "detailed", "described"],
        "controls": ["control", "blank", "negative control", "positive control"],
        "replicates": ["replicate", "triplicate", "n="],
        "sample_size": ["n=", "sample size", "power", "statistical power"],
        "randomization": ["random", "randomized", "randomly"],
        "blinding": ["blind", "blinded", "double-blind", "observer"]
    }
    
    description_lower = description.lower()
    results = {}
    for criterion, terms in keywords.items():
        results[criterion] = any(term in description_lower for term in terms)
    
    return results


def main():
    parser = argparse.ArgumentParser(description="Peer Review Tools")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Create template command
    create_parser = subparsers.add_parser("create-template", help="Create review template")
    create_parser.add_argument("-o", "--output", required=True, help="Output JSON file")
    create_parser.add_argument("-t", "--title", default="", help="Manuscript title")
    create_parser.add_argument("-a", "--authors", default="", help="Authors")
    create_parser.add_argument("-r", "--reviewer", default="", help="Reviewer name")
    
    # Generate report command
    report_parser = subparsers.add_parser("generate-report", help="Generate review report from JSON")
    report_parser.add_argument("json_file", help="Review JSON file")
    report_parser.add_argument("-o", "--output", help="Output text file (default: stdout)")
    
    # Evaluate methodology command
    eval_parser = subparsers.add_parser("evaluate-methodology", help="Auto-evaluate methodology")
    eval_parser.add_argument("text_file", help="Methods section text file")
    
    args = parser.parse_args()
    
    if args.command == "create-template":
        create_review_template(args.output, args.title, args.authors, args.reviewer)
    
    elif args.command == "generate-report":
        review = load_review(args.json_file)
        report = generate_review_report(review)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(report)
            print(f"Report saved: {args.output}")
        else:
            print(report)
    
    elif args.command == "evaluate-methodology":
        with open(args.text_file, 'r') as f:
            text = f.read()
        results = evaluate_methodology(text)
        
        print("Methodology Evaluation:")
        for criterion, passed in results.items():
            status = "✓" if passed else "✗"
            print(f"  [{status}] {criterion}")
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
