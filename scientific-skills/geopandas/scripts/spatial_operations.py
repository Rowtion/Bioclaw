"""
GeoPandas Example: Spatial Data Operations
Read, process, and analyze geospatial vector data.
"""
import geopandas as gpd
import matplotlib.pyplot as plt


def read_spatial_data(file_path: str):
    """
    Read spatial data from various formats.
    
    Args:
        file_path: Path to spatial file
        
    Returns:
        GeoDataFrame
    """
    gdf = gpd.read_file(file_path)
    
    print(f"DataFrame shape: {gdf.shape}")
    print(f"Columns: {list(gdf.columns)}")
    print(f"CRS: {gdf.crs}")
    print(f"Geometry types: {gdf.geometry.geom_type.unique()}")
    
    return gdf


def reproject_data(gdf: gpd.GeoDataFrame, target_crs: str):
    """
    Reproject data to a different coordinate system.
    
    Args:
        gdf: GeoDataFrame
        target_crs: Target CRS (e.g., 'EPSG:3857')
        
    Returns:
        Reprojected GeoDataFrame
    """
    gdf_reprojected = gdf.to_crs(target_crs)
    print(f"Reprojected from {gdf.crs} to {gdf_reprojected.crs}")
    return gdf_reprojected


def calculate_geometric_properties(gdf: gpd.GeoDataFrame):
    """
    Calculate area, length, and other geometric properties.
    
    Args:
        gdf: GeoDataFrame
        
    Returns:
        GeoDataFrame with new columns
    """
    # Ensure projected CRS for accurate measurements
    if gdf.crs and gdf.crs.is_geographic:
        gdf = gdf.to_crs('EPSG:3857')
    
    gdf['area'] = gdf.geometry.area
    gdf['perimeter'] = gdf.geometry.length
    gdf['centroid_x'] = gdf.geometry.centroid.x
    gdf['centroid_y'] = gdf.geometry.centroid.y
    
    print("\nGeometric Properties:")
    print(f"Total area: {gdf['area'].sum():,.0f} sq units")
    print(f"Average area: {gdf['area'].mean():,.0f} sq units")
    
    return gdf


def buffer_geometries(gdf: gpd.GeoDataFrame, distance: float):
    """
    Create buffer around geometries.
    
    Args:
        gdf: GeoDataFrame
        distance: Buffer distance
        
    Returns:
        Buffered GeoDataFrame
    """
    gdf_buffered = gdf.copy()
    gdf_buffered['geometry'] = gdf.geometry.buffer(distance)
    
    print(f"Created {distance} unit buffer around {len(gdf)} features")
    return gdf_buffered


if __name__ == "__main__":
    print("GeoPandas Spatial Operations Example")
    # gdf = read_spatial_data('data.shp')
    # gdf_projected = reproject_data(gdf, 'EPSG:3857')
    # gdf_props = calculate_geometric_properties(gdf_projected)
