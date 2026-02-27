#!/usr/bin/env python3
"""
Scientific Critical Thinking Tools
Evaluate arguments, identify biases, and assess evidence quality.
"""

from typing import List, Dict, Tuple
import json


def evaluate_argument_strength(premises: List[str], conclusion: str) -> Dict:
    """Evaluate the strength of a scientific argument."""
    score = 0
    feedback = []
    
    # Check number of premises
    if len(premises) >= 3:
        score += 2
        feedback.append("Multiple supporting premises")
    elif len(premises) >= 1:
        score += 1
        feedback.append("Limited supporting premises")
    else:
        feedback.append("No premises provided")
    
    # Check conclusion specificity
    if len(conclusion.split()) >= 5:
        score += 1
        feedback.append("Specific conclusion")
    
    return {
        "score": min(score, 5),
        "max_score": 5,
        "feedback": feedback,
        "strength": "strong" if score >= 4 else "moderate" if score >= 2 else "weak"
    }


def identify_logical_fallacies(text: str) -> List[Dict]:
    """Identify potential logical fallacies in text."""
    fallacies = []
    
    fallacy_patterns = {
        "appeal_to_authority": ["expert says", "according to", "authority claims"],
        "correlation_causation": ["causes", "leads to", "results in"],
        "hasty_generalization": ["always", "never", "everyone", "nobody"],
        "confirmation_bias": ["obviously", "clearly", "undoubtedly"]
    }
    
    text_lower = text.lower()
    for fallacy_type, patterns in fallacy_patterns.items():
        for pattern in patterns:
            if pattern in text_lower:
                fallacies.append({
                    "type": fallacy_type,
                    "pattern": pattern,
                    "suggestion": f"Review usage of '{pattern}' - potential {fallacy_type}"
                })
    
    return fallacies


def assess_evidence_quality(evidence_type: str, sample_size: int = None,
                           peer_reviewed: bool = False) -> Dict:
    """Assess the quality of scientific evidence."""
    quality_score = 0
    factors = []
    
    # Evidence type weighting
    type_scores = {
        "systematic_review": 5,
        "randomized_controlled_trial": 4,
        "cohort_study": 3,
        "case_control": 2,
        "case_report": 1,
        "expert_opinion": 1,
        "anecdotal": 0
    }
    
    quality_score += type_scores.get(evidence_type, 0)
    
    # Sample size
    if sample_size:
        if sample_size >= 1000:
            quality_score += 2
            factors.append("Large sample size")
        elif sample_size >= 100:
            quality_score += 1
            factors.append("Moderate sample size")
    
    # Peer review
    if peer_reviewed:
        quality_score += 1
        factors.append("Peer reviewed")
    
    return {
        "quality_score": min(quality_score, 8),
        "max_score": 8,
        "level": "high" if quality_score >= 6 else "moderate" if quality_score >= 4 else "low",
        "factors": factors
    }


def critique_paper_abstract(abstract: str) -> Dict:
    """Generate critical questions for a paper abstract."""
    questions = []
    
    sections = ["background", "methods", "results", "conclusion"]
    for section in sections:
        if section not in abstract.lower():
            questions.append(f"What is the {section}?")
    
    questions.extend([
        "What are the limitations of this study?",
        "Are the results generalizable?",
        "What alternative explanations exist?",
        "Is the sample representative?",
        "Were appropriate statistical methods used?"
    ])
    
    return {
        "missing_sections": [s for s in sections if s not in abstract.lower()],
        "critical_questions": questions
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scientific Critical Thinking Tools")
    parser.add_argument("command", 
                       choices=["argument", "fallacies", "evidence", "critique"])
    parser.add_argument("--text", help="Text to analyze")
    parser.add_argument("--premises", nargs="+", help="Argument premises")
    parser.add_argument("--conclusion", help="Argument conclusion")
    parser.add_argument("--evidence-type", help="Type of evidence")
    parser.add_argument("--sample-size", type=int, help="Sample size")
    parser.add_argument("--peer-reviewed", action="store_true")
    
    args = parser.parse_args()
    
    if args.command == "argument":
        if not args.premises or not args.conclusion:
            print("Error: --premises and --conclusion required")
            return
        result = evaluate_argument_strength(args.premises, args.conclusion)
        print(json.dumps(result, indent=2))
    
    elif args.command == "fallacies":
        if not args.text:
            print("Error: --text required")
            return
        fallacies = identify_logical_fallacies(args.text)
        print(json.dumps(fallacies, indent=2))
    
    elif args.command == "evidence":
        result = assess_evidence_quality(
            args.evidence_type or "case_report",
            args.sample_size,
            args.peer_reviewed
        )
        print(json.dumps(result, indent=2))
    
    elif args.command == "critique":
        if not args.text:
            print("Error: --text required")
            return
        critique = critique_paper_abstract(args.text)
        print(json.dumps(critique, indent=2))


if __name__ == "__main__":
    main()
