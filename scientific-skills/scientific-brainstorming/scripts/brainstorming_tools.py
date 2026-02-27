#!/usr/bin/env python3
"""
Scientific Brainstorming Tools
Generate research ideas and hypotheses through structured brainstorming.
"""

import json
import random
from typing import List, Dict


BRAINSTORMING_TECHNIQUES = {
    "scamper": {
        "prompts": ["Substitute", "Combine", "Adapt", "Modify", "Put to other uses", "Eliminate", "Reverse"],
        "description": "Systematic innovation technique"
    },
    "six_hats": {
        "prompts": ["White (facts)", "Red (emotions)", "Black (cautions)", "Yellow (benefits)", "Green (creativity)", "Blue (process)"],
        "description": "Six Thinking Hats method"
    },
    "mind_mapping": {
        "prompts": ["Central concept", "Main branches", "Sub-branches", "Connections", "Keywords"],
        "description": "Visual brainstorming"
    }
}


def generate_research_ideas(domain: str, num_ideas: int = 5) -> List[str]:
    """Generate research ideas for a given domain."""
    templates = [
        f"How can {domain} be integrated with machine learning?",
        f"What are the unexplored applications of {domain} in healthcare?",
        f"How can {domain} address climate change challenges?",
        f"What ethical considerations arise from advances in {domain}?",
        f"How can {domain} benefit from interdisciplinary collaboration?",
        f"What novel methodologies can advance {domain}?",
        f"How can {domain} improve data accessibility?",
        f"What are the long-term implications of {domain} development?"
    ]
    return random.sample(templates, min(num_ideas, len(templates)))


def apply_scamper(topic: str) -> Dict[str, str]:
    """Apply SCAMPER technique to a research topic."""
    prompts = BRAINSTORMING_TECHNIQUES["scamper"]["prompts"]
    results = {}
    for prompt in prompts:
        if prompt == "Substitute":
            results[prompt] = f"What components of {topic} can be replaced?"
        elif prompt == "Combine":
            results[prompt] = f"What can be combined with {topic}?"
        elif prompt == "Adapt":
            results[prompt] = f"How can {topic} be adapted for new uses?"
        elif prompt == "Modify":
            results[prompt] = f"How can {topic} be modified or magnified?"
        elif prompt == "Put to other uses":
            results[prompt] = f"What other uses can {topic} serve?"
        elif prompt == "Eliminate":
            results[prompt] = f"What can be removed from {topic}?"
        elif prompt == "Reverse":
            results[prompt] = f"How can {topic} be reversed or rearranged?"
    return results


def generate_collaboration_opportunities(field1: str, field2: str) -> List[str]:
    """Generate interdisciplinary collaboration ideas."""
    ideas = [
        f"Apply {field1} methodologies to solve {field2} problems",
        f"Use {field2} data to validate {field1} models",
        f"Develop hybrid tools combining {field1} and {field2}",
        f"Create educational programs merging {field1} and {field2}",
        f"Establish joint research centers for {field1}-{field2} studies"
    ]
    return ideas


def brainstorm_session(topic: str, technique: str = "scamper") -> Dict:
    """Run a structured brainstorming session."""
    session = {
        "topic": topic,
        "technique": technique,
        "ideas": generate_research_ideas(topic, 5),
        "scamper_analysis": apply_scamper(topic),
        "collaborations": generate_collaboration_opportunities(topic, "AI")
    }
    return session


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Scientific Brainstorming Tools")
    parser.add_argument("command", choices=["ideas", "scamper", "collaborate", "session"])
    parser.add_argument("--topic", required=True, help="Research topic")
    parser.add_argument("--field2", help="Second field for collaboration")
    parser.add_argument("--num", type=int, default=5, help="Number of ideas")
    parser.add_argument("--output", help="Output JSON file")
    
    args = parser.parse_args()
    
    if args.command == "ideas":
        ideas = generate_research_ideas(args.topic, args.num)
        for i, idea in enumerate(ideas, 1):
            print(f"{i}. {idea}")
    
    elif args.command == "scamper":
        results = apply_scamper(args.topic)
        for prompt, question in results.items():
            print(f"{prompt}: {question}")
    
    elif args.command == "collaborate":
        if not args.field2:
            print("Error: --field2 required for collaboration")
            return
        ideas = generate_collaboration_opportunities(args.topic, args.field2)
        for idea in ideas:
            print(f"- {idea}")
    
    elif args.command == "session":
        session = brainstorm_session(args.topic)
        print(json.dumps(session, indent=2))
        if args.output:
            with open(args.output, 'w') as f:
                json.dump(session, f, indent=2)


if __name__ == "__main__":
    main()
