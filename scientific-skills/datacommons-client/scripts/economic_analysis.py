"""
Data Commons Client Example: Economic Analysis
Analyze economic indicators across different regions.
"""
from datacommons_client import DataCommonsClient
import pandas as pd


def compare_states_economy():
    """Compare economic indicators across US states."""
    client = DataCommonsClient()
    
    variables = [
        "Median_Income_Household",
        "UnemploymentRate_Person",
        "Count_Person"
    ]
    
    states = ["geoId/06", "geoId/48", "geoId/36", "geoId/12"]  # CA, TX, NY, FL
    
    response = client.observation.fetch(
        variable_dcids=variables,
        entity_dcids=states,
        date="latest"
    )
    
    df = response.to_observations_as_records()
    
    # Pivot for easier comparison
    pivot = df.pivot_table(
        values='value',
        index='entity',
        columns='variable'
    )
    
    print("Economic Comparison:")
    print(pivot)
    return pivot


def analyze_income_distribution():
    """Analyze income distribution across counties."""
    client = DataCommonsClient()
    
    response = client.observation.fetch(
        variable_dcids=["Median_Income_Household"],
        entity_expression="country/USA<-containedInPlace+{typeOf:County}",
        date="2021"
    )
    
    df = response.to_observations_as_records()
    
    print("\nIncome Distribution Stats:")
    print(df['value'].describe())
    
    return df


if __name__ == "__main__":
    compare_states_economy()
    analyze_income_distribution()
