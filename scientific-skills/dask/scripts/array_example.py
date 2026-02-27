"""
Dask Array Example: Large-scale Numerical Computing
Demonstrates working with arrays larger than memory using Dask Arrays.
"""
import dask.array as da
import numpy as np


def create_large_random_array():
    """Create a large random array with chunked processing."""
    # Create 100,000 x 100,000 array (10 billion elements)
    # Chunked into 10,000 x 10,000 blocks
    x = da.random.random((100000, 100000), chunks=(10000, 10000))
    
    print(f"Array shape: {x.shape}")
    print(f"Chunk size: {x.chunks}")
    print(f"Total size: {x.nbytes / 1e9:.2f} GB")
    
    # Operations are lazy
    y = x + 100
    z = y.mean(axis=0)
    
    # Compute result
    result = z.compute()
    print(f"Result shape: {result.shape}")
    return result


def parallel_matrix_operations():
    """Demonstrate parallel matrix operations."""
    # Create two large matrices
    a = da.random.random((50000, 50000), chunks=(5000, 5000))
    b = da.random.random((50000, 50000), chunks=(5000, 5000))
    
    # Matrix multiplication (parallelized)
    c = da.dot(a, b)
    
    # Reduction operations
    row_sums = c.sum(axis=1)
    col_means = c.mean(axis=0)
    
    # Compute all results
    print("Computing matrix operations...")
    result_sums, result_means = da.compute(row_sums, col_means)
    
    print(f"Row sums shape: {result_sums.shape}")
    print(f"Column means shape: {result_means.shape}")


def process_image_stack():
    """Process a stack of images in parallel."""
    # Simulate 1000 images of 1024x1024
    images = da.random.random((1000, 1024, 1024), chunks=(10, 1024, 1024))
    
    # Apply filters lazily
    normalized = (images - images.mean(axis=0)) / images.std(axis=0)
    
    # Compute statistics per image
    per_image_mean = normalized.mean(axis=(1, 2))
    per_image_std = normalized.std(axis=(1, 2))
    
    return da.compute(per_image_mean, per_image_std)


if __name__ == "__main__":
    # Example: Create and process large array
    print("Creating large random array...")
    result = create_large_random_array()
    print(f"Mean of means: {result.mean():.4f}")
