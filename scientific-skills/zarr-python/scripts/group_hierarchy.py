#!/usr/bin/env python3
"""
Zarr Group and Hierarchy Utilities
Organize multiple arrays hierarchically like directories.
"""

import zarr
from typing import Optional, Dict, Any
import json


def create_hierarchy(store_path: str, structure: Dict[str, Any]) -> zarr.Group:
    """
    Create a hierarchical structure of groups and arrays.
    
    Args:
        store_path: Path to store
        structure: Dictionary defining the hierarchy
            {
                'group_name': {
                    'subgroup': {
                        'array_name': {'shape': (100, 100), 'chunks': (10, 10), 'dtype': 'f4'}
                    }
                }
            }
    
    Returns:
        Root Zarr Group
    """
    root = zarr.group(store=store_path)
    
    def _create_recursive(parent, name: str, content):
        if isinstance(content, dict):
            if 'shape' in content:
                # It's an array definition
                parent.create_array(
                    name=name,
                    shape=content['shape'],
                    chunks=content.get('chunks', content['shape']),
                    dtype=content.get('dtype', 'f4')
                )
            else:
                # It's a group
                group = parent.create_group(name)
                for child_name, child_content in content.items():
                    _create_recursive(group, child_name, child_content)
    
    for name, content in structure.items():
        _create_recursive(root, name, content)
    
    return root


def create_climate_dataset(store_path: str) -> zarr.Group:
    """
    Create a sample climate dataset structure with groups.
    
    Args:
        store_path: Path to store
    
    Returns:
        Root Zarr Group with climate data structure
    """
    structure = {
        'temperature': {
            't2m': {'shape': (365, 720, 1440), 'chunks': (1, 720, 1440), 'dtype': 'f4'},
            'tmin': {'shape': (365, 720, 1440), 'chunks': (1, 720, 1440), 'dtype': 'f4'},
            'tmax': {'shape': (365, 720, 1440), 'chunks': (1, 720, 1440), 'dtype': 'f4'}
        },
        'precipitation': {
            'prcp': {'shape': (365, 720, 1440), 'chunks': (1, 720, 1440), 'dtype': 'f4'}
        },
        'humidity': {
            'rh': {'shape': (365, 720, 1440), 'chunks': (1, 720, 1440), 'dtype': 'f4'}
        }
    }
    
    root = create_hierarchy(store_path, structure)
    
    # Add metadata
    root.attrs['project'] = 'Climate Analysis'
    root.attrs['institution'] = 'Research Institute'
    root.attrs['created'] = '2024-01-15'
    
    return root


def print_tree(group: zarr.Group, max_depth: int = 10):
    """
    Print the hierarchy tree of a Zarr group.
    
    Args:
        group: Zarr Group to print
        max_depth: Maximum depth to display
    """
    print(group.tree(level=max_depth))


def copy_group_structure(source_path: str, dest_path: str, mode: str = 'w'):
    """
    Copy the structure (groups and array shapes) from one Zarr store to another.
    Does not copy data, only creates empty arrays with same structure.
    
    Args:
        source_path: Source Zarr store path
        dest_path: Destination Zarr store path
        mode: Mode for opening destination ('w' or 'a')
    """
    source = zarr.open_group(source_path, mode='r')
    dest = zarr.open_group(dest_path, mode=mode)
    
    def _copy_recursive(src_grp, dest_grp):
        # Copy attributes
        for key, value in src_grp.attrs.items():
            dest_grp.attrs[key] = value
        
        # Copy arrays
        for name in src_grp.array_keys():
            arr = src_grp[name]
            dest_grp.create_array(
                name=name,
                shape=arr.shape,
                chunks=arr.chunks,
                dtype=arr.dtype
            )
            # Copy array attributes
            for key, value in arr.attrs.items():
                dest_grp[name].attrs[key] = value
        
        # Recurse into subgroups
        for name in src_grp.group_keys():
            src_sub = src_grp[name]
            dest_sub = dest_grp.create_group(name)
            _copy_recursive(src_sub, dest_sub)
    
    _copy_recursive(source, dest)
    print(f"Structure copied from {source_path} to {dest_path}")


def add_attributes(obj, attributes: Dict[str, Any]):
    """
    Add JSON-serializable attributes to a Zarr array or group.
    
    Args:
        obj: Zarr Array or Group
        attributes: Dictionary of attributes to add
    """
    for key, value in attributes.items():
        obj.attrs[key] = value


def get_attributes(obj) -> Dict[str, Any]:
    """
    Get all attributes from a Zarr array or group.
    
    Args:
        obj: Zarr Array or Group
    
    Returns:
        Dictionary of attributes
    """
    return dict(obj.attrs)


def search_arrays(group: zarr.Group, pattern: str = None) -> list:
    """
    Search for arrays in a group hierarchy.
    
    Args:
        group: Zarr Group to search
        pattern: Optional pattern to match array names
    
    Returns:
        List of array paths
    """
    arrays = []
    
    def _search(grp, path=""):
        for name in grp.array_keys():
            full_path = f"{path}/{name}" if path else name
            if pattern is None or pattern in name:
                arrays.append(full_path)
        
        for name in grp.group_keys():
            subgrp = grp[name]
            subpath = f"{path}/{name}" if path else name
            _search(subgrp, subpath)
    
    _search(group)
    return arrays


def merge_groups(source_paths: list, dest_path: str, group_names: list = None):
    """
    Merge multiple Zarr groups into a single destination.
    
    Args:
        source_paths: List of source Zarr store paths
        dest_path: Destination Zarr store path
        group_names: Optional names for each source group in destination
    """
    if group_names is None:
        group_names = [f"group_{i}" for i in range(len(source_paths))]
    
    dest = zarr.open_group(dest_path, mode='w')
    
    for src_path, grp_name in zip(source_paths, group_names):
        src = zarr.open_group(src_path, mode='r')
        dest_grp = dest.create_group(grp_name)
        
        def _copy_recursive(src_grp, dest_grp):
            for key, value in src_grp.attrs.items():
                dest_grp.attrs[key] = value
            
            for name in src_grp.array_keys():
                arr = src_grp[name]
                dest_arr = dest_grp.create_array(
                    name=name,
                    shape=arr.shape,
                    chunks=arr.chunks,
                    dtype=arr.dtype
                )
                dest_arr[:] = arr[:]  # Copy data
                
                for k, v in arr.attrs.items():
                    dest_arr.attrs[k] = v
            
            for name in src_grp.group_keys():
                src_sub = src_grp[name]
                dest_sub = dest_grp.create_group(name)
                _copy_recursive(src_sub, dest_sub)
        
        _copy_recursive(src, dest_grp)
    
    print(f"Merged {len(source_paths)} groups into {dest_path}")


if __name__ == "__main__":
    print("Zarr Group and Hierarchy Utilities")
    print("=" * 40)
    
    # Create example climate dataset
    root = create_climate_dataset("example_climate.zarr")
    print("\nCreated climate dataset structure:")
    print_tree(root)
    
    # Add attributes
    temp = root['temperature/t2m']
    add_attributes(temp, {
        'description': '2-meter temperature',
        'units': 'K',
        'long_name': 'Temperature at 2 meters'
    })
    
    print("\nArray attributes:")
    print(get_attributes(temp))
    
    print("\nAll arrays in hierarchy:")
    print(search_arrays(root))
