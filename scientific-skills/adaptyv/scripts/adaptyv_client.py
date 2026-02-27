#!/usr/bin/env python3
"""
Adaptyv API Client - Basic Usage Examples
Submit protein sequences for experimental testing via Adaptyv cloud laboratory.
"""

import os
import requests
import json
from typing import Optional, Dict, Any
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Configuration
API_KEY = os.getenv("ADAPTYV_API_KEY")
BASE_URL = os.getenv("ADAPTYV_BASE_URL", 
    "https://kq5jp7qj7wdqklhsxmovkzn4l40obksv.lambda-url.eu-central-1.on.aws")


class AdaptyvClient:
    """Client for interacting with Adaptyv protein testing API."""
    
    def __init__(self, api_key: Optional[str] = None, base_url: Optional[str] = None):
        self.api_key = api_key or API_KEY
        self.base_url = base_url or BASE_URL
        
        if not self.api_key:
            raise ValueError("API key required. Set ADAPTYV_API_KEY environment variable.")
    
    def _headers(self) -> Dict[str, str]:
        return {
            "Authorization": f"Bearer {self.api_key}",
            "Content-Type": "application/json"
        }
    
    def submit_experiment(self, 
                         sequences: str,
                         experiment_type: str = "binding",
                         webhook_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Submit protein sequences for experimental testing.
        
        Args:
            sequences: Protein sequences in FASTA format
            experiment_type: Type of experiment (binding, expression, thermostability, enzyme_activity)
            webhook_url: Optional URL to receive status callbacks
            
        Returns:
            API response with experiment ID and status
        """
        url = f"{self.base_url}/experiments"
        
        payload = {
            "sequences": sequences,
            "experiment_type": experiment_type
        }
        
        if webhook_url:
            payload["webhook_url"] = webhook_url
        
        response = requests.post(url, headers=self._headers(), json=payload)
        response.raise_for_status()
        return response.json()
    
    def get_experiment_status(self, experiment_id: str) -> Dict[str, Any]:
        """Get status of a submitted experiment."""
        url = f"{self.base_url}/experiments/{experiment_id}"
        response = requests.get(url, headers=self._headers())
        response.raise_for_status()
        return response.json()
    
    def list_experiments(self, limit: int = 10) -> Dict[str, Any]:
        """List recent experiments."""
        url = f"{self.base_url}/experiments?limit={limit}"
        response = requests.get(url, headers=self._headers())
        response.raise_for_status()
        return response.json()
    
    def get_results(self, experiment_id: str) -> Dict[str, Any]:
        """Download results for a completed experiment."""
        url = f"{self.base_url}/experiments/{experiment_id}/results"
        response = requests.get(url, headers=self._headers())
        response.raise_for_status()
        return response.json()


def main():
    """Example usage of Adaptyv client."""
    
    # Example protein sequence (single-chain nanobody)
    test_sequence = """>nanobody_example
QVQLQESGGGLVQAGGSLRLSCAASGFTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAADKGYCSGGSCEYWGQGTQVTVSS"""
    
    try:
        # Initialize client
        client = AdaptyvClient()
        
        # Submit experiment
        print("Submitting protein for binding assay...")
        result = client.submit_experiment(
            sequences=test_sequence,
            experiment_type="binding"
        )
        
        experiment_id = result.get("experiment_id")
        print(f"Experiment submitted! ID: {experiment_id}")
        print(f"Status: {result.get('status')}")
        
        # Check status (in real usage, you'd poll this over time)
        print("\nChecking experiment status...")
        status = client.get_experiment_status(experiment_id)
        print(f"Current status: {status.get('status')}")
        print(f"Estimated completion: ~21 days")
        
    except Exception as e:
        print(f"Error: {e}")
        print("\nNote: This is an alpha/beta API. Contact support@adaptyvbio.com for access.")


if __name__ == "__main__":
    main()
