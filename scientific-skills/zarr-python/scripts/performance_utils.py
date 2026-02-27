#!/usr/bin/env python3
"""
Zarr Performance Utilities
Compression, sharding, and performance optimization.
"""

import zarr
import numpy as np
from zarr.codecs.blosc import BloscCodec
from zarr.codecs import GzipCodec, ZstdCodec, BytesCodec
from typing import List, Optional


def get_compression_info(zarr_path: str) -> dict:
    """
    Get compression information for a Zarr array.
    
    Args:
        zarr_path: Path to Zarr array
    
    Returns:
        Dictionary with compression statistics
    """
    z = zarr.open_array(zarr_path, mode='r')
    
    return {
        'shape': z.shape,
        'chunks': z.chunks,
        'dtype': str(z.dtype),
        'uncompressed_bytes': z.nbytes,
        'compressed_bytes': z.nbytes_stored,
        'compression_ratio': z.nbytes / z.nbytes_stored if z.nbytes_stored > 0 else 0,
        'codecs': [str(c) for c in z.codecs] if hasattr(z, 'codecs') else 'default'
    }


def print_compression_report(zarr_path: str):
    """
    Print a detailed compression report for a Zarr array.
    
    Args:
        zarr_path: Path to Zarr array
    """
    info = get_compression_info(zarr_path)
    
    print(f"\nCompression Report: {zarr_path}")
    print("=" * 50)
    print(f"Shape:              {info['shape']}")
    print(f"Chunks:             {info['chunks']}")
    print(f"Dtype:              {info['dtype']}")
    print(f"Uncompressed size:  {info['uncompressed_bytes'] / 1e6:.2f} MB")
    print(f"Compressed size:    {info['compressed_bytes'] / 1e6:.2f} MB")
    print(f"Compression ratio:  {info['compression_ratio']:.2f}x")
    print(f"Codecs:             {info['codecs']}")


def create_compressed_array(
    store_path: str,
    shape: tuple,
    chunks: tuple,
    dtype: str = "f4",
    compression: str = "zstd",
    compression_level: int = 5,
    shuffle: bool = True
) -> zarr.Array:
    """
    Create a Zarr array with specified compression settings.
    
    Args:
        store_path: Path to store
        shape: Array shape
        chunks: Chunk shape
        dtype: Data type
        compression: 'zstd', 'lz4', 'gzip', 'zlib', 'blosclz', 'snappy', 'none'
        compression_level: Compression level (varies by codec)
        shuffle: Enable shuffle filter for numeric data
    
    Returns:
        Zarr Array instance
    """
    shuffle_val = 'shuffle' if shuffle else 'noshuffle'
    
    if compression == "none":
        codecs = [BytesCodec()]
    elif compression == "gzip":
        codecs = [GzipCodec(level=compression_level)]
    elif compression == "zstd":
        codecs = [BloscCodec(cname='zstd', clevel=compression_level, shuffle=shuffle_val)]
    elif compression in ['lz4', 'lz4hc', 'blosclz', 'zlib', 'snappy']:
        codecs = [BloscCodec(cname=compression, clevel=compression_level, shuffle=shuffle_val)]
    else:
        codecs = [BloscCodec(cname='zstd', clevel=compression_level, shuffle=shuffle_val)]
    
    z = zarr.create_array(
        store=store_path,
        shape=shape,
        chunks=chunks,
        dtype=dtype,
        codecs=codecs
    )
    return z


def create_sharded_array(
    store_path: str,
    shape: tuple,
    chunks: tuple,
    shards: tuple,
    dtype: str = "f4",
    compression: str = "zstd"
) -> zarr.Array:
    """
    Create a Zarr array with sharding for large-scale storage.
    Sharding groups many small chunks into larger storage objects.
    
    Args:
        store_path: Path to store
        shape: Array shape
        chunks: Chunk shape (small for access patterns)
        shards: Shard shape (larger, groups chunks together)
        dtype: Data type
        compression: Compression codec
    
    Returns:
        Zarr Array instance
    """
    shuffle_val = 'shuffle'
    
    if compression == "zstd":
        codecs = [BloscCodec(cname='zstd', clevel=5, shuffle=shuffle_val)]
    elif compression == "lz4":
        codecs = [BloscCodec(cname='lz4', clevel=1)]
    else:
        codecs = [BloscCodec(cname='zstd', clevel=5, shuffle=shuffle_val)]
    
    z = zarr.create_array(
        store=store_path,
        shape=shape,
        chunks=chunks,
        shards=shards,
        dtype=dtype,
        codecs=codecs
    )
    return z


def benchmark_compression(
    data: np.ndarray,
    chunks: tuple,
    codecs_to_test: List[str] = None
) -> dict:
    """
    Benchmark different compression codecs on sample data.
    
    Args:
        data: Sample data array
        chunks: Chunk shape
        codecs_to_test: List of codecs to test
    
    Returns:
        Dictionary with benchmark results
    """
    if codecs_to_test is None:
        codecs_to_test = ['zstd', 'lz4', 'gzip', 'none']
    
    results = {}
    
    for codec in codecs_to_test:
        import tempfile
        import time
        
        with tempfile.TemporaryDirectory() as tmpdir:
            path = f"{tmpdir}/test.zarr"
            
            # Create array with codec
            z = create_compressed_array(
                path,
                shape=data.shape,
                chunks=chunks,
                dtype=str(data.dtype),
                compression=codec
            )
            
            # Time write
            start = time.time()
            z[:] = data
            write_time = time.time() - start
            
            # Get compression info
            info = get_compression_info(path)
            
            results[codec] = {
                'write_time': write_time,
                'compression_ratio': info['compression_ratio'],
                'compressed_mb': info['compressed_bytes'] / 1e6,
                'uncompressed_mb': info['uncompressed_bytes'] / 1e6
            }
    
    return results


def print_benchmark_results(results: dict):
    """Print compression benchmark results in a formatted table."""
    print("\nCompression Benchmark Results")
    print("=" * 70)
    print(f"{'Codec':<10} {'Ratio':<10} {'Compressed':<12} {'Write Time':<12}")
    print("-" * 70)
    
    for codec, data in results.items():
        print(f"{codec:<10} {data['compression_ratio']:<10.2f} "
              f"{data['compressed_mb']:<12.2f} {data['write_time']:<12.3f}s")


def optimize_chunks_for_access(
    shape: tuple,
    dtype: str,
    access_pattern: str = "balanced",
    target_chunk_mb: float = 1.0
) -> tuple:
    """
    Calculate optimal chunk shape based on access pattern.
    
    Args:
        shape: Array shape
        dtype: Data type
        access_pattern: 'row', 'column', 'balanced', or 'time_series'
        target_chunk_mb: Target chunk size in MB
    
    Returns:
        Optimal chunk shape tuple
    """
    dtype_size = np.dtype(dtype).itemsize
    target_elements = int((target_chunk_mb * 1024 * 1024) / dtype_size)
    ndim = len(shape)
    
    if access_pattern == "row" and ndim >= 2:
        # Optimize for row-wise access: chunks span columns
        rows = max(1, min(shape[0], target_elements // shape[1]))
        return (rows,) + shape[1:]
    
    elif access_pattern == "column" and ndim >= 2:
        # Optimize for column-wise access: chunks span rows
        cols = max(1, min(shape[1], target_elements // shape[0]))
        return (shape[0], cols) + shape[2:]
    
    elif access_pattern == "time_series" and ndim >= 3:
        # Time series: small time dimension, full spatial
        return (1,) + shape[1:]
    
    else:  # balanced
        side = int(target_elements ** (1/ndim))
        return tuple(min(s, side) for s in shape)


def consolidate_store_metadata(store_path: str):
    """
    Consolidate metadata for a Zarr store.
    Reduces I/O operations for hierarchical stores.
    
    Args:
        store_path: Path to Zarr store
    """
    zarr.consolidate_metadata(store_path)
    print(f"Metadata consolidated for {store_path}")


def open_consolidated(store_path: str):
    """
    Open a Zarr store with consolidated metadata.
    Faster for hierarchical stores, especially on cloud storage.
    
    Args:
        store_path: Path to Zarr store
    
    Returns:
        Zarr Group or Array with consolidated metadata
    """
    return zarr.open_consolidated(store_path)


if __name__ == "__main__":
    print("Zarr Performance Utilities")
    print("=" * 40)
    
    # Create test data
    print("\nCreating test data...")
    test_data = np.random.random((1000, 1000)).astype('f4')
    
    # Benchmark compression
    print("\nBenchmarking compression codecs...")
    results = benchmark_compression(test_data, chunks=(100, 100))
    print_benchmark_results(results)
    
    # Create arrays with different optimizations
    print("\n\nCreating optimized arrays...")
    
    # 1. Standard compression
    z1 = create_compressed_array(
        "example_compressed.zarr",
        shape=(5000, 5000),
        chunks=(500, 500),
        compression="zstd",
        compression_level=5
    )
    z1[:] = np.random.random((5000, 5000))
    print_compression_report("example_compressed.zarr")
    
    # 2. Sharded array
    z2 = create_sharded_array(
        "example_sharded.zarr",
        shape=(10000, 10000),
        chunks=(100, 100),
        shards=(1000, 1000),
        compression="zstd"
    )
    print("\nCreated sharded array: 100x100 chunks grouped into 1000x1000 shards")
