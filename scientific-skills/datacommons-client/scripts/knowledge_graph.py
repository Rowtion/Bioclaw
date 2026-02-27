"""
Data Commons Client Example: Knowledge Graph Exploration
Explore entity relationships and properties in Data Commons.
"""
from datacommons_client import DataCommonsClient


def explore_entity_properties():
    """Discover properties of a geographic entity."""
    client = DataCommonsClient()
    
    # Get property labels for California
    labels = client.node.fetch_property_labels(
        node_dcids=["geoId/06"],
        out=True
    )
    
    print("California Properties:")
    print(labels.to_dict())


def navigate_hierarchy():
    """Navigate geographic hierarchy."""
    client = DataCommonsClient()
    
    # Get children of USA (states)
    children = client.node.fetch_place_children(
        node_dcids=["country/USA"]
    )
    
    print("\nUS States (sample):")
    result = children.to_dict()
    for dcid, data in list(result.items())[:5]:
        print(f"  {dcid}")


def get_entity_names():
    """Get human-readable names for DCIDs."""
    client = DataCommonsClient()
    
    names = client.node.fetch_entity_names(
        node_dcids=["geoId/06", "geoId/48", "country/USA"]
    )
    
    print("\nEntity Names:")
    for dcid, name in names.to_dict().items():
        print(f"  {dcid}: {name}")


def discover_available_variables():
    """Find what statistics are available for a place."""
    client = DataCommonsClient()
    
    variables = client.observation.fetch_available_statistical_variables(
        entity_dcids=["geoId/06"]
    )
    
    print("\nAvailable Variables for California (sample):")
    vars_df = variables.to_dataframe()
    print(vars_df.head(10))


if __name__ == "__main__":
    explore_entity_properties()
    navigate_hierarchy()
    get_entity_names()
    discover_available_variables()
