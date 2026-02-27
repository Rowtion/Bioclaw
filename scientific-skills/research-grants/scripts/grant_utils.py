#!/usr/bin/env python3
"""
Research Grant Utilities
Utilities for managing grant proposal workflows.
"""

import argparse
import json
import os
from datetime import datetime, timedelta
from pathlib import Path


class GrantDatabase:
    """Simple JSON-based grant opportunity database."""
    
    def __init__(self, db_file: str = "grants.json"):
        self.db_file = db_file
        self.grants = self._load()
    
    def _load(self) -> dict:
        """Load database from file."""
        if os.path.exists(self.db_file):
            with open(self.db_file) as f:
                return json.load(f)
        return {'grants': [], 'applications': []}
    
    def _save(self):
        """Save database to file."""
        with open(self.db_file, 'w') as f:
            json.dump(self.grants, f, indent=2)
    
    def add_grant(self, agency: str, program: str, deadline: str,
                  amount: str, description: str, url: str = None) -> str:
        """Add a new grant opportunity."""
        grant_id = f"grant_{len(self.grants['grants']):04d}"
        grant = {
            'id': grant_id,
            'agency': agency,
            'program': program,
            'deadline': deadline,
            'amount': amount,
            'description': description,
            'url': url,
            'added': datetime.now().isoformat(),
            'status': 'open'
        }
        self.grants['grants'].append(grant)
        self._save()
        return grant_id
    
    def list_grants(self, agency: str = None, status: str = None) -> list:
        """List grant opportunities with optional filtering."""
        grants = self.grants['grants']
        
        if agency:
            grants = [g for g in grants if g['agency'].lower() == agency.lower()]
        if status:
            grants = [g for g in grants if g['status'] == status]
        
        return grants
    
    def get_upcoming_deadlines(self, days: int = 30) -> list:
        """Get grants with deadlines within specified days."""
        upcoming = []
        now = datetime.now()
        
        for grant in self.grants['grants']:
            if grant['status'] != 'open':
                continue
            try:
                deadline = datetime.fromisoformat(grant['deadline'])
                if now <= deadline <= now + timedelta(days=days):
                    upcoming.append(grant)
            except ValueError:
                continue
        
        return sorted(upcoming, key=lambda x: x['deadline'])
    
    def add_application(self, grant_id: str, title: str, pi: str) -> str:
        """Add a new application."""
        app_id = f"app_{len(self.grants['applications']):04d}"
        application = {
            'id': app_id,
            'grant_id': grant_id,
            'title': title,
            'pi': pi,
            'status': 'draft',
            'created': datetime.now().isoformat(),
            'updated': datetime.now().isoformat(),
            'milestones': []
        }
        self.grants['applications'].append(application)
        self._save()
        return app_id


def calculate_budget(personnel: list, equipment: list = None, 
                     travel: float = 0, other: float = 0,
                     indirect_rate: float = 0.5) -> dict:
    """
    Calculate research budget with indirect costs.
    
    Args:
        personnel: List of dicts with 'name', 'role', 'months', 'rate'
        equipment: List of equipment costs
        travel: Travel budget
        other: Other direct costs
        indirect_rate: Indirect cost rate (e.g., 0.5 for 50%)
        
    Returns:
        Budget breakdown
    """
    # Calculate personnel costs
    personnel_total = 0
    for person in personnel:
        cost = person.get('months', 0) * person.get('rate', 0)
        personnel_total += cost
    
    # Calculate equipment costs
    equipment_total = sum(equipment) if equipment else 0
    
    # Direct costs
    direct_costs = personnel_total + equipment_total + travel + other
    
    # Indirect costs
    indirect_costs = direct_costs * indirect_rate
    
    # Total
    total_costs = direct_costs + indirect_costs
    
    return {
        'personnel': personnel_total,
        'equipment': equipment_total,
        'travel': travel,
        'other_direct': other,
        'total_direct': direct_costs,
        'indirect': indirect_costs,
        'total': total_costs,
    }


def generate_nsf_project_summary(title: str, overview: str, 
                                  intellectual_merit: str,
                                  broader_impacts: str) -> str:
    """
    Generate NSF Project Summary document.
    
    Args:
        title: Project title
        overview: Brief overview
        intellectual_merit: Intellectual merit statement
        broader_impacts: Broader impacts statement
        
    Returns:
        Formatted project summary
    """
    summary = f"""# NSF Project Summary

## Title
{title}

## Overview
{overview}

## Intellectual Merit
{intellectual_merit}

## Broader Impacts
{broader_impacts}

---
Generated: {datetime.now().strftime('%Y-%m-%d')}
"""
    return summary


def generate_nih_specific_aims(title: str, significance: str,
                                innovation: str, aims: list) -> str:
    """
    Generate NIH Specific Aims page.
    
    Args:
        title: Project title
        significance: Significance statement
        innovation: Innovation statement
        aims: List of specific aims dicts with 'title' and 'approach'
        
    Returns:
        Formatted Specific Aims document
    """
    aims_text = "\n\n".join([
        f"**Aim {i+1}**: {aim['title']}\n\n{aim.get('approach', '')}"
        for i, aim in enumerate(aims)
    ])
    
    document = f"""# NIH Specific Aims

## Title
{title}

## Significance
{significance}

## Innovation
{innovation}

## Specific Aims

{aims_text}

---
Generated: {datetime.now().strftime('%Y-%m-%d')}
Page limit: 1 page
"""
    return document


def main():
    parser = argparse.ArgumentParser(description='Research Grant Utilities')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Database commands
    db_parser = subparsers.add_parser('db', help='Database operations')
    db_sub = db_parser.add_subparsers(dest='db_cmd', help='DB commands')
    
    # Add grant
    add_grant = db_sub.add_parser('add', help='Add grant opportunity')
    add_grant.add_argument('-a', '--agency', required=True, help='Funding agency')
    add_grant.add_argument('-p', '--program', required=True, help='Program name')
    add_grant.add_argument('-d', '--deadline', required=True, help='Deadline (ISO format)')
    add_grant.add_argument('-m', '--amount', help='Award amount')
    add_grant.add_argument('--desc', help='Description')
    add_grant.add_argument('-f', '--file', default='grants.json', help='Database file')
    
    # List grants
    list_grants = db_sub.add_parser('list', help='List grants')
    list_grants.add_argument('-a', '--agency', help='Filter by agency')
    list_grants.add_argument('-f', '--file', default='grants.json', help='Database file')
    
    # Deadlines
    deadlines = db_sub.add_parser('deadlines', help='Show upcoming deadlines')
    deadlines.add_argument('-d', '--days', type=int, default=30, help='Days ahead')
    deadlines.add_argument('-f', '--file', default='grants.json', help='Database file')
    
    # Budget calculator
    budget_parser = subparsers.add_parser('budget', help='Calculate budget')
    budget_parser.add_argument('-p', '--personnel', type=float, default=100000,
                               help='Personnel costs')
    budget_parser.add_argument('-e', '--equipment', type=float, default=50000,
                               help='Equipment costs')
    budget_parser.add_argument('-t', '--travel', type=float, default=5000,
                               help='Travel costs')
    budget_parser.add_argument('-o', '--other', type=float, default=10000,
                               help='Other direct costs')
    budget_parser.add_argument('-r', '--indirect-rate', type=float, default=0.5,
                               help='Indirect rate (e.g., 0.5 for 50%)')
    
    # NSF Project Summary
    nsf_parser = subparsers.add_parser('nsf-summary', help='Generate NSF Project Summary')
    nsf_parser.add_argument('-t', '--title', required=True, help='Project title')
    nsf_parser.add_argument('-o', '--overview', required=True, help='Overview text')
    nsf_parser.add_argument('-i', '--intellectual-merit', required=True, 
                            help='Intellectual merit')
    nsf_parser.add_argument('-b', '--broader-impacts', required=True,
                            help='Broader impacts')
    nsf_parser.add_argument('--output', help='Output file')
    
    # NIH Specific Aims
    nih_parser = subparsers.add_parser('nih-aims', help='Generate NIH Specific Aims')
    nih_parser.add_argument('-t', '--title', required=True, help='Project title')
    nih_parser.add_argument('-s', '--significance', required=True, help='Significance')
    nih_parser.add_argument('-i', '--innovation', required=True, help='Innovation')
    
    args = parser.parse_args()
    
    if args.command == 'db':
        db = GrantDatabase(args.file if hasattr(args, 'file') else 'grants.json')
        
        if args.db_cmd == 'add':
            grant_id = db.add_grant(args.agency, args.program, args.deadline,
                                    args.amount, args.desc or '')
            print(f"Grant added: {grant_id}")
            
        elif args.db_cmd == 'list':
            grants = db.list_grants(agency=args.agency)
            print(f"\n{'Agency':<15} {'Program':<30} {'Deadline':<20} {'Status'}")
            print("-" * 80)
            for g in grants:
                print(f"{g['agency']:<15} {g['program']:<30} "
                      f"{g['deadline']:<20} {g['status']}")
                      
        elif args.db_cmd == 'deadlines':
            upcoming = db.get_upcoming_deadlines(args.days)
            print(f"\nUpcoming Deadlines (next {args.days} days):")
            print(f"{'Deadline':<20} {'Agency':<15} {'Program'}")
            print("-" * 60)
            for g in upcoming:
                deadline = datetime.fromisoformat(g['deadline']).strftime('%Y-%m-%d')
                print(f"{deadline:<20} {g['agency']:<15} {g['program']}")
                
    elif args.command == 'budget':
        budget = calculate_budget(
            personnel=[{'months': 12, 'rate': args.personnel / 12}],
            equipment=[args.equipment],
            travel=args.travel,
            other=args.other,
            indirect_rate=args.indirect_rate
        )
        print("\nBudget Summary:")
        print(f"  Personnel:     ${budget['personnel']:,.2f}")
        print(f"  Equipment:     ${budget['equipment']:,.2f}")
        print(f"  Travel:        ${budget['travel']:,.2f}")
        print(f"  Other Direct:  ${budget['other_direct']:,.2f}")
        print(f"  Total Direct:  ${budget['total_direct']:,.2f}")
        print(f"  Indirect:      ${budget['indirect']:,.2f}")
        print(f"  Total:         ${budget['total']:,.2f}")
        
    elif args.command == 'nsf-summary':
        summary = generate_nsf_project_summary(
            args.title, args.overview,
            args.intellectual_merit, args.broader_impacts
        )
        if args.output:
            with open(args.output, 'w') as f:
                f.write(summary)
            print(f"Project summary saved to: {args.output}")
        else:
            print(summary)
            
    elif args.command == 'nih-aims':
        print("Note: Specific aims content should be provided via arguments")
        print("Use --help for usage information")
        
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
