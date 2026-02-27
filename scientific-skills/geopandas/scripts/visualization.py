"""
GeoPandas Example: Visualization
Create maps and visualizations.
"""
import geopandas as gpd
import matplotlib.pyplot as plt


def create_choropleth_map(gdf: gpd.GeoDataFrame, column: str, 
                          cmap: str = 'YlOrRd', output: str = 'map.png'):
    """
    Create a choropleth map.
    
    Args:
        gdf: GeoDataFrame
        column: Column to visualize
        cmap: Colormap
        output: Output file
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    gdf.plot(
        column=column,
        cmap=cmap,
        legend=True,
        legend_kwds={'label': column},
        ax=ax
    )
    
    ax.set_title(f'{column} by Region')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Map saved: {output}")


def create_multilayer_map(gdf_list: list, colors: list, 
                          labels: list, output: str = 'multilayer.png'):
    """
    Create a multi-layer map.
    
    Args:
        gdf_list: List of GeoDataFrames
        colors: List of colors
        labels: List of labels
        output: Output file
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for gdf, color, label in zip(gdf_list, colors, labels):
        gdf.plot(color=color, ax=ax, label=label, alpha=0.6)
    
    ax.legend()
    ax.set_title('Multi-layer Map')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    print(f"Multi-layer map saved: {output}")


def create_interactive_map(gdf: gpd.GeoDataFrame, column: str = None):
    """
    Create an interactive folium map.
    
    Args:
        gdf: GeoDataFrame
        column: Column to visualize
        
    Returns:
        Folium map object
    """
    try:
        import folium
        
        # Reproject to WGS84 for folium
        if gdf.crs != 'EPSG:4326':
            gdf = gdf.to_crs('EPSG:4326')
        
        # Create map centered on data
        center = [gdf.geometry.centroid.y.mean(), gdf.geometry.centroid.x.mean()]
        m = folium.Map(location=center, zoom_start=10)
        
        # Add GeoJSON
        folium.GeoJson(gdf).add_to(m)
        
        return m
    except ImportError:
        print("folium not installed. Install with: pip install folium")
        return None


def plot_with_basemap(gdf: gpd.GeoDataFrame, output: str = 'basemap.png'):
    """
    Create map with basemap using contextily.
    
    Args:
        gdf: GeoDataFrame (should be in EPSG:3857 or will be reprojected)
        output: Output file
    """
    try:
        import contextily as ctx
        
        # Reproject to Web Mercator for basemap
        if gdf.crs != 'EPSG:3857':
            gdf = gdf.to_crs('EPSG:3857')
        
        fig, ax = plt.subplots(figsize=(12, 8))
        gdf.plot(ax=ax, alpha=0.5, color='red')
        ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)
        ax.set_title('Map with Basemap')
        ax.axis('off')
        
        plt.tight_layout()
        plt.savefig(output, dpi=300, bbox_inches='tight')
        print(f"Basemap saved: {output}")
    
    except ImportError:
        print("contextily not installed. Install with: pip install contextily")


if __name__ == "__main__":
    print("GeoPandas Visualization Example")
    # gdf = gpd.read_file('data.shp')
    # create_choropleth_map(gdf, 'population')
