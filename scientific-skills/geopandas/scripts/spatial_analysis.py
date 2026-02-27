"""
GeoPandas Example: Spatial Analysis
Perform spatial joins and overlay operations.
"""
import geopandas as gpd


def spatial_join_points_to_polygons(points_gdf: gpd.GeoDataFrame,
                                     polygons_gdf: gpd.GeoDataFrame):
    """
    Join points to polygons they fall within.
    
    Args:
        points_gdf: Point GeoDataFrame
        polygons_gdf: Polygon GeoDataFrame
        
    Returns:
        Joined GeoDataFrame
    """
    # Ensure same CRS
    if points_gdf.crs != polygons_gdf.crs:
        points_gdf = points_gdf.to_crs(polygons_gdf.crs)
    
    joined = gpd.sjoin(
        points_gdf,
        polygons_gdf,
        predicate='within',
        how='left'
    )
    
    print(f"Joined {len(joined)} points to polygons")
    return joined


def spatial_overlay(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame,
                    how: str = 'intersection'):
    """
    Perform spatial overlay operation.
    
    Args:
        gdf1: First GeoDataFrame
        gdf2: Second GeoDataFrame
        how: Overlay type ('intersection', 'union', 'difference', 'symmetric_difference')
        
    Returns:
        Overlay result
    """
    # Ensure same CRS
    if gdf1.crs != gdf2.crs:
        gdf2 = gdf2.to_crs(gdf1.crs)
    
    result = gpd.overlay(gdf1, gdf2, how=how)
    print(f"Overlay ({how}): {len(result)} features")
    
    return result


def dissolve_by_attribute(gdf: gpd.GeoDataFrame, by: str, aggfunc: str = 'sum'):
    """
    Dissolve geometries by attribute.
    
    Args:
        gdf: GeoDataFrame
        by: Column to dissolve by
        aggfunc: Aggregation function
        
    Returns:
        Dissolved GeoDataFrame
    """
    dissolved = gdf.dissolve(by=by, aggfunc=aggfunc)
    print(f"Dissolved {len(gdf)} features into {len(dissolved)} groups")
    return dissolved


def nearest_neighbor_join(gdf1: gpd.GeoDataFrame, gdf2: gpd.GeoDataFrame,
                          max_distance: float = None):
    """
    Find nearest neighbors between two datasets.
    
    Args:
        gdf1: First GeoDataFrame (left)
        gdf2: Second GeoDataFrame (right)
        max_distance: Maximum distance to consider
        
    Returns:
        Joined GeoDataFrame
    """
    # Ensure same CRS
    if gdf1.crs != gdf2.crs:
        gdf2 = gdf2.to_crs(gdf1.crs)
    
    joined = gpd.sjoin_nearest(
        gdf1,
        gdf2,
        max_distance=max_distance,
        how='left'
    )
    
    print(f"Found nearest neighbors for {len(joined)} features")
    return joined


if __name__ == "__main__":
    print("GeoPandas Spatial Analysis Example")
    # points = gpd.read_file('points.shp')
    # polygons = gpd.read_file('polygons.shp')
    # joined = spatial_join_points_to_polygons(points, polygons)
