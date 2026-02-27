#!/usr/bin/env python3
"""
Benchling Inventory Operations

This script demonstrates managing inventory items
in Benchling.

Usage:
    export BENCHLING_API_KEY="your_key"
    python inventory_operations.py
"""

import os
import sys


def example_create_container():
    """Example: Create a container."""
    print("\n" + "="*60)
    print("Creating Container")
    print("="*60)
    
    example_code = '''
from benchling_sdk.benchling import Benchling
from benchling_sdk.auth.api_key_auth import ApiKeyAuth
from benchling_sdk.models import ContainerCreate

benchling = Benchling(
    url=os.environ['BENCHLING_URL'],
    auth_method=ApiKeyAuth(os.environ['BENCHLING_API_KEY'])
)

# Create container
container = benchling.containers.create(
    ContainerCreate(
        name="Sample Tube 001",
        schema_id="your_container_schema_id",
        fields=benchling.models.fields({
            "concentration": "100 ng/µL",
            "volume": "50 µL"
        })
    )
)

print(f"Created container: {container.id}")
'''
    
    print(example_code)


def example_transfer_item():
    """Example: Transfer item between locations."""
    print("\n" + "="*60)
    print("Transferring Inventory Item")
    print("="*60)
    
    example_code = '''
# Transfer container to new location
transfer = benchling.containers.transfer(
    container_id="cont_abc123",
    destination_id="box_xyz789"
)

print(f"Transferred to: {transfer.destination_storage_id}")
'''
    
    print(example_code)


def example_list_inventory():
    """Example: List inventory items."""
    print("\n" + "="*60)
    print("Listing Inventory")
    print("="*60)
    
    example_code = '''
# List containers in a box
containers = benchling.containers.list(
    parent_storage_id="box_abc123"
)

for page in containers:
    for container in page:
        print(f"  - {container.name}: {container.barcode}")

# List all boxes
boxes = benchling.boxes.list()
for page in boxes:
    for box in page:
        print(f"Box: {box.name}")
'''
    
    print(example_code)


def main():
    """Display inventory examples."""
    print("Benchling Inventory Operations Examples")
    print("=" * 60)
    
    example_create_container()
    example_transfer_item()
    example_list_inventory()
    
    print("\n" + "="*60)
    print("Examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
