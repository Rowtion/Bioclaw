#!/usr/bin/env python3
"""
Zarr-Dask Integration Utilities
Parallel computing with chunked arrays.
"""

import zarr
import numpy as np
from typing import Optional, Callable, Any


def zarr_to_dask(zarr_path: str, store_type: str = 'local'):
    """
    Load a Zarr array as a Dask array for parallel computation.
    
    Args:
        zarr_path: Path to Zarr array
        store_type: 'local', 's3', or 'gcs'
    
    Returns:
        Dask array
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required. Install with: uv pip install dask")
    
    if store_type == 'local':
        return da.from_zarr(zarr_path)
    elif store_type == 's3':
        try:
            import s3fs
        except ImportError:
            raise ImportError("s3fs required for S3 access")
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=zarr_path, s3=s3)
        return da.from_zarr(store)
    elif store_type == 'gcs':
        try:
            import gcsfs
        except ImportError:
            raise ImportError("gcsfs required for GCS access")
        gcs = gcsfs.GCSFileSystem()
        store = gcsfs.GCSMap(root=zarr_path, gcs=gcs)
        return da.from_zarr(store)
    else:
        raise ValueError(f"Unknown store type: {store_type}")


def dask_to_zarr(dask_array, zarr_path: str, store_type: str = 'local', **kwargs):
    """
    Save a Dask array to Zarr storage.
    
    Args:
        dask_array: Dask array to save
        zarr_path: Path to save Zarr array
        store_type: 'local', 's3', or 'gcs'
        **kwargs: Additional arguments for to_zarr
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required")
    
    if store_type == 'local':
        da.to_zarr(dask_array, zarr_path, **kwargs)
    elif store_type == 's3':
        try:
            import s3fs
        except ImportError:
            raise ImportError("s3fs required")
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=zarr_path, s3=s3)
        da.to_zarr(dask_array, store, **kwargs)
    elif store_type == 'gcs':
        try:
            import gcsfs
        except ImportError:
            raise ImportError("gcsfs required")
        gcs = gcsfs.GCSFileSystem()
        store = gcsfs.GCSMap(root=zarr_path, gcs=gcs)
        da.to_zarr(dask_array, store, **kwargs)


def parallel_compute(zarr_path: str, operation: str, axis: Optional[int] = None, **kwargs):
    """
    Perform parallel computation on a Zarr array using Dask.
    
    Args:
        zarr_path: Path to Zarr array
        operation: Operation to perform ('mean', 'sum', 'std', 'max', 'min', 'dot')
        axis: Axis for reduction operations
        **kwargs: Additional arguments
    
    Returns:
        Computed result
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required")
    
    darr = da.from_zarr(zarr_path)
    
    if operation == 'mean':
        result = darr.mean(axis=axis)
    elif operation == 'sum':
        result = darr.sum(axis=axis)
    elif operation == 'std':
        result = darr.std(axis=axis)
    elif operation == 'max':
        result = darr.max(axis=axis)
    elif operation == 'min':
        result = darr.min(axis=axis)
    elif operation == 'dot':
        other = kwargs.get('other')
        if other is None:
            raise ValueError("'other' required for dot operation")
        result = darr.dot(other)
    else:
        raise ValueError(f"Unknown operation: {operation}")
    
    return result.compute()


def parallel_apply(zarr_path: str, func: Callable, output_path: str, dtype: Optional[str] = None):
    """
    Apply a function to each chunk of a Zarr array in parallel.
    
    Args:
        zarr_path: Path to input Zarr array
        func: Function to apply to each chunk
        output_path: Path to save output Zarr array
        dtype: Output dtype (defaults to input dtype)
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required")
    
    darr = da.from_zarr(zarr_path)
    result = darr.map_blocks(func, dtype=dtype or darr.dtype)
    da.to_zarr(result, output_path)


def create_large_random_array(output_path: str, shape: tuple, chunks: tuple, dtype: str = 'f4'):
    """
    Create a large random array using Dask and save to Zarr.
    Useful for creating test datasets larger than memory.
    
    Args:
        output_path: Path to save Zarr array
        shape: Total shape of array
        chunks: Chunk shape
        dtype: Data type
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required")
    
    darr = da.random.random(shape, chunks=chunks).astype(dtype)
    da.to_zarr(darr, output_path)
    print(f"Created random array: shape={shape}, chunks={chunks}")


def rechunk_array(input_path: str, output_path: str, new_chunks: tuple):
    """
    Rechunk a Zarr array using Dask.
    
    Args:
        input_path: Path to input Zarr array
        output_path: Path to save rechunked array
        new_chunks: New chunk shape
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required")
    
    darr = da.from_zarr(input_path)
    rechunked = darr.rechunk(new_chunks)
    da.to_zarr(rechunked, output_path)
    print(f"Rechunked array from {darr.chunks} to {new_chunks}")


def parallel_matrix_multiply(a_path: str, b_path: str, output_path: str):
    """
    Perform parallel matrix multiplication on two Zarr arrays.
    
    Args:
        a_path: Path to first Zarr array
        b_path: Path to second Zarr array
        output_path: Path to save result
    """
    try:
        import dask.array as da
    except ImportError:
        raise ImportError("dask required")
    
    a = da.from_zarr(a_path)
    b = da.from_zarr(b_path)
    result = a.dot(b)
    da.to_zarr(result, output_path)
    print(f"Matrix multiplication result saved to {output_path}")


if __name__ == "__main__":
    print("Zarr-Dask Integration Utilities")
    print("=" * 40)
    
    # Create a sample array
    print("\nCreating sample Zarr array...")
    z = zarr.open_array("example_dask.zarr", mode='w', shape=(10000, 10000),
                        chunks=(1000, 1000), dtype='f4')
    z[:] = np.random.random((10000, 10000))
    
    # Load as Dask array
    print("Loading as Dask array...")
    darr = zarr_to_dask("example_dask.zarr")
    print(f"Dask array shape: {darr.shape}, chunks: {darr.chunks}")
    
    # Parallel computation
    print("\nComputing mean...")
    mean_val = parallel_compute("example_dask.zarr", "mean")
    print(f"Mean value: {mean_val}")
    
    print("\nComputing sum along axis 0...")
    sum_result = parallel_compute("example_dask.zarr", "sum", axis=0)
    print(f"Sum result shape: {sum_result.shape}")
