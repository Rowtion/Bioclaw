#!/usr/bin/env python3
"""
Zarr Cloud Storage Utilities
S3 and GCS integration for cloud-native workflows.
"""

import zarr
from typing import Optional, Dict, Any


def create_s3_array(
    bucket_path: str,
    shape: tuple,
    chunks: tuple,
    dtype: str = "f4",
    aws_access_key: Optional[str] = None,
    aws_secret_key: Optional[str] = None,
    region: str = "us-east-1",
    endpoint_url: Optional[str] = None,
    anon: bool = False
) -> zarr.Array:
    """
    Create a Zarr array on S3 storage.
    
    Args:
        bucket_path: S3 path (e.g., 'my-bucket/data/array.zarr')
        shape: Array shape
        chunks: Chunk shape
        dtype: Data type
        aws_access_key: AWS access key (optional)
        aws_secret_key: AWS secret key (optional)
        region: AWS region
        endpoint_url: Custom S3 endpoint (e.g., MinIO)
        anon: Use anonymous access
    
    Returns:
        Zarr Array instance
    """
    try:
        import s3fs
    except ImportError:
        raise ImportError("s3fs required. Install with: uv pip install s3fs")
    
    # Configure S3 filesystem
    s3_kwargs = {"anon": anon}
    if aws_access_key and aws_secret_key:
        s3_kwargs["key"] = aws_access_key
        s3_kwargs["secret"] = aws_secret_key
    if endpoint_url:
        s3_kwargs["client_kwargs"] = {"endpoint_url": endpoint_url}
    
    s3 = s3fs.S3FileSystem(**s3_kwargs)
    store = s3fs.S3Map(root=bucket_path, s3=s3)
    
    z = zarr.open_array(
        store=store,
        mode='w',
        shape=shape,
        chunks=chunks,
        dtype=dtype
    )
    return z


def open_s3_array(
    bucket_path: str,
    mode: str = 'r',
    aws_access_key: Optional[str] = None,
    aws_secret_key: Optional[str] = None,
    region: str = "us-east-1",
    endpoint_url: Optional[str] = None,
    anon: bool = False
) -> zarr.Array:
    """
    Open an existing Zarr array on S3.
    
    Args:
        bucket_path: S3 path
        mode: Access mode ('r', 'r+', 'a')
        aws_access_key: AWS access key
        aws_secret_key: AWS secret key
        region: AWS region
        endpoint_url: Custom S3 endpoint
        anon: Use anonymous access
    
    Returns:
        Zarr Array instance
    """
    try:
        import s3fs
    except ImportError:
        raise ImportError("s3fs required. Install with: uv pip install s3fs")
    
    s3_kwargs = {"anon": anon}
    if aws_access_key and aws_secret_key:
        s3_kwargs["key"] = aws_access_key
        s3_kwargs["secret"] = aws_secret_key
    if endpoint_url:
        s3_kwargs["client_kwargs"] = {"endpoint_url": endpoint_url}
    
    s3 = s3fs.S3FileSystem(**s3_kwargs)
    store = s3fs.S3Map(root=bucket_path, s3=s3)
    
    z = zarr.open_array(store=store, mode=mode)
    return z


def create_gcs_array(
    bucket_path: str,
    shape: tuple,
    chunks: tuple,
    dtype: str = "f4",
    project: Optional[str] = None,
    token: Optional[str] = None
) -> zarr.Array:
    """
    Create a Zarr array on Google Cloud Storage.
    
    Args:
        bucket_path: GCS path (e.g., 'my-bucket/data/array.zarr')
        shape: Array shape
        chunks: Chunk shape
        dtype: Data type
        project: GCP project ID
        token: Path to service account JSON or 'anon' for anonymous
    
    Returns:
        Zarr Array instance
    """
    try:
        import gcsfs
    except ImportError:
        raise ImportError("gcsfs required. Install with: uv pip install gcsfs")
    
    gcs = gcsfs.GCSFileSystem(project=project, token=token)
    store = gcsfs.GCSMap(root=bucket_path, gcs=gcs)
    
    z = zarr.open_array(
        store=store,
        mode='w',
        shape=shape,
        chunks=chunks,
        dtype=dtype
    )
    return z


def open_gcs_array(
    bucket_path: str,
    mode: str = 'r',
    project: Optional[str] = None,
    token: Optional[str] = None
) -> zarr.Array:
    """
    Open an existing Zarr array on GCS.
    
    Args:
        bucket_path: GCS path
        mode: Access mode
        project: GCP project ID
        token: Path to service account JSON or 'anon'
    
    Returns:
        Zarr Array instance
    """
    try:
        import gcsfs
    except ImportError:
        raise ImportError("gcsfs required. Install with: uv pip install gcsfs")
    
    gcs = gcsfs.GCSFileSystem(project=project, token=token)
    store = gcsfs.GCSMap(root=bucket_path, gcs=gcs)
    
    z = zarr.open_array(store=store, mode=mode)
    return z


def consolidate_cloud_metadata(store_path: str, backend: str = "s3", **kwargs):
    """
    Consolidate metadata for cloud storage to reduce latency.
    
    Args:
        store_path: Cloud storage path
        backend: 's3' or 'gcs'
        **kwargs: Backend-specific authentication args
    """
    if backend == "s3":
        try:
            import s3fs
        except ImportError:
            raise ImportError("s3fs required")
        s3 = s3fs.S3FileSystem(**kwargs)
        store = s3fs.S3Map(root=store_path, s3=s3)
    elif backend == "gcs":
        try:
            import gcsfs
        except ImportError:
            raise ImportError("gcsfs required")
        gcs = gcsfs.GCSFileSystem(**kwargs)
        store = gcsfs.GCSMap(root=store_path, gcs=gcs)
    else:
        raise ValueError(f"Unknown backend: {backend}")
    
    zarr.consolidate_metadata(store)
    print(f"Metadata consolidated for {store_path}")


def open_consolidated_cloud(store_path: str, backend: str = "s3", **kwargs):
    """
    Open Zarr store with consolidated metadata (faster for cloud).
    
    Args:
        store_path: Cloud storage path
        backend: 's3' or 'gcs'
        **kwargs: Backend-specific authentication args
    
    Returns:
        Zarr Group or Array
    """
    if backend == "s3":
        try:
            import s3fs
        except ImportError:
            raise ImportError("s3fs required")
        s3 = s3fs.S3FileSystem(**kwargs)
        store = s3fs.S3Map(root=store_path, s3=s3)
    elif backend == "gcs":
        try:
            import gcsfs
        except ImportError:
            raise ImportError("gcsfs required")
        gcs = gcsfs.GCSFileSystem(**kwargs)
        store = gcsfs.GCSMap(root=store_path, gcs=gcs)
    else:
        raise ValueError(f"Unknown backend: {backend}")
    
    return zarr.open_consolidated(store)


if __name__ == "__main__":
    print("Zarr Cloud Storage Utilities")
    print("=" * 40)
    print("\nExample usage:")
    print("""
    # S3 Example
    z = create_s3_array(
        bucket_path='my-bucket/data/array.zarr',
        shape=(10000, 10000),
        chunks=(1000, 1000),
        dtype='f4',
        region='us-west-2'
    )
    
    # GCS Example
    z = create_gcs_array(
        bucket_path='my-bucket/data/array.zarr',
        shape=(10000, 10000),
        chunks=(1000, 1000),
        project='my-project'
    )
    """)
