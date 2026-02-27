#!/usr/bin/env python3
"""
Zarr Array Creation Utilities
Chunked N-dimensional array creation with optimal settings.
"""

import zarr
import numpy as np
from zarr.codecs.blosc import BloscCodec
from zarr.codecs import GzipCodec, ZstdCodec, BytesCodec
from typing import Optional, List, Union


def create_optimal_array(
    store_path: str,
    shape: tuple,
    dtype: str = "f4",
    target_chunk_mb: float = 1.0,
    compression: str = "zstd",
    fill_value: Optional[Union[int, float]] = None
) -> zarr.Array:
    """
    Create a Zarr array with optimal chunking for the given shape.
    
    Args:
        store_path: Path to store the array
        shape: Shape of the array
        dtype: Data type (default: f4 for float32)
        target_chunk_mb: Target chunk size in MB (default: 1.0)
        compression: Compression codec ('zstd', 'lz4', 'gzip', 'none')
        fill_value: Fill value for empty array (None for no fill)
    
    Returns:
        Zarr Array instance
    """
    # Calculate element size in bytes
    dtype_size = np.dtype(dtype).itemsize
    target_elements = int((target_chunk_mb * 1024 * 1024) / dtype_size)
    
    # Calculate chunk shape
    ndim = len(shape)
    chunk_shape = calculate_chunk_shape(shape, target_elements, ndim)
    
    # Configure compression
    codecs = _get_compression_codecs(compression)
    
    # Create array
    if fill_value is not None:
        z = zarr.create_array(
            store=store_path,
            shape=shape,
            chunks=chunk_shape,
            dtype=dtype,
            codecs=codecs,
            fill_value=fill_value
        )
    else:
        z = zarr.create_array(
            store=store_path,
            shape=shape,
            chunks=chunk_shape,
            dtype=dtype,
            codecs=codecs
        )
    
    return z


def calculate_chunk_shape(shape: tuple, target_elements: int, ndim: int) -> tuple:
    """Calculate balanced chunk shape for N-dimensional array."""
    if ndim == 1:
        return (min(shape[0], target_elements),)
    elif ndim == 2:
        side = int(np.sqrt(target_elements))
        return (
            min(shape[0], side),
            min(shape[1], side)
        )
    elif ndim == 3:
        side = int(target_elements ** (1/3))
        return (
            min(shape[0], side),
            min(shape[1], side),
            min(shape[2], side)
        )
    else:
        # For higher dimensions, use more balanced approach
        side = int(target_elements ** (1/ndim))
        return tuple(min(s, side) for s in shape)


def _get_compression_codecs(compression: str) -> List:
    """Get compression codec based on name."""
    if compression == "zstd":
        return [BloscCodec(cname='zstd', clevel=5, shuffle='shuffle')]
    elif compression == "lz4":
        return [BloscCodec(cname='lz4', clevel=1)]
    elif compression == "gzip":
        return [GzipCodec(level=6)]
    elif compression == "none":
        return [BytesCodec()]
    else:
        return [BloscCodec(cname='zstd', clevel=5, shuffle='shuffle')]


def create_time_series_array(
    store_path: str,
    spatial_shape: tuple,
    dtype: str = "f4",
    time_chunk_size: int = 1
) -> zarr.Array:
    """
    Create array optimized for time series data.
    Time dimension should have small chunks for efficient appending.
    
    Args:
        store_path: Path to store the array
        spatial_shape: Shape of spatial dimensions (e.g., (720, 1440) for lat/lon)
        dtype: Data type
        time_chunk_size: Number of time steps per chunk (default: 1)
    
    Returns:
        Zarr Array with shape (0, *spatial_shape) for appending
    """
    shape = (0,) + spatial_shape
    chunks = (time_chunk_size,) + spatial_shape
    
    z = zarr.open_array(
        store=store_path,
        mode='a',
        shape=shape,
        chunks=chunks,
        dtype=dtype
    )
    return z


def create_row_optimized_array(
    store_path: str,
    shape: tuple,
    dtype: str = "f4"
) -> zarr.Array:
    """
    Create array optimized for row-wise access.
    Chunks span full columns for efficient row reads.
    """
    rows, cols = shape
    # Each chunk spans all columns
    chunk_shape = (max(1, int(1048576 / (cols * np.dtype(dtype).itemsize))), cols)
    
    z = zarr.create_array(
        store=store_path,
        shape=shape,
        chunks=chunk_shape,
        dtype=dtype
    )
    return z


def create_column_optimized_array(
    store_path: str,
    shape: tuple,
    dtype: str = "f4"
) -> zarr.Array:
    """
    Create array optimized for column-wise access.
    Chunks span full rows for efficient column reads.
    """
    rows, cols = shape
    # Each chunk spans all rows
    chunk_shape = (rows, max(1, int(1048576 / (rows * np.dtype(dtype).itemsize))))
    
    z = zarr.create_array(
        store=store_path,
        shape=shape,
        chunks=chunk_shape,
        dtype=dtype
    )
    return z


if __name__ == "__main__":
    # Example usage
    print("Creating sample arrays...")
    
    # 1. General purpose array
    z1 = create_optimal_array(
        "example_general.zarr",
        shape=(10000, 10000),
        dtype="f4",
        target_chunk_mb=1.0,
        compression="zstd"
    )
    print(f"Created general array: {z1.shape}, chunks={z1.chunks}")
    
    # 2. Time series array
    z2 = create_time_series_array(
        "example_timeseries.zarr",
        spatial_shape=(720, 1440),
        dtype="f4",
        time_chunk_size=1
    )
    print(f"Created time series array: {z2.shape}, chunks={z2.chunks}")
    
    # 3. Row-optimized array
    z3 = create_row_optimized_array(
        "example_row_optimized.zarr",
        shape=(10000, 10000),
        dtype="f4"
    )
    print(f"Created row-optimized array: {z3.shape}, chunks={z3.chunks}")
    
    print("\nDone! Arrays created in current directory.")
