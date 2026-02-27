"""
Data Commons Client Example: Population Statistics Query
Query demographic and economic data from Data Commons.
"""
import os
from datacommons_client import DataCommonsClient


def query_population_data():
    """Query population statistics for US states."""
    client = DataCommonsClient()
    
    # Get population for California and Texas
    response = client.observation.fetch(
        variable_dcids=["Count_Person"],
        entity_dcids=["geoId/06", "geoId/48"],  # CA, TX
        date="latest"
    )
    
    # Convert to DataFrame
    df = response.to_observations_as_records()
    print("Latest Population Data:")
    print(df)
    return df


def query_time_series():
    """Query unemployment rate time series."""
    client = DataCommonsClient()
    
    response = client.observation.fetch(
        variable_dcids=["UnemploymentRate_Person"],
        entity_dcids=["country/USA"],
        date="all"
    )
    
    df = response.to_observations_as_records()
    print("\nUnemployment Rate (All Years):")
    print(df.head(10))
    return df


def query_by_hierarchy():
    """Query all counties in a state."""
    client = DataCommonsClient()
    
    # Query all counties in California with median income
    response = client.observation.fetch(
        variable_dcids=["Median_Income_Household"],
        entity_expression="geoId/06<-containedInPlace+{typeOf:County}",
        date="2020"
    )
    
    df = response.to_observations_as_records()
    print("\nMedian Household Income by County (CA, 2020):")
    print(df.head())
    return df


def resolve_place_names():
    """Convert place names to DCIDs."""
    client = DataCommonsClient()
    
    response = client.resolve.fetch_dcids_by_name(
        names=["California", "Texas", "New York"],
        entity_type="State"
    )
    
    print("\nResolved Place Names:")
    for name, result in response.to_dict().items():
        if result["candidates"]:
            print(f"{name}: {result['candidates'][0]['dcid']}")


if __name__ == "__main__":
    # Set API key if available
    # os.environ["DC_API_KEY"] = "your_key_here"
    
    query_population_data()
    query_time_series()
    query_by_hierarchy()
    resolve_place_names()
