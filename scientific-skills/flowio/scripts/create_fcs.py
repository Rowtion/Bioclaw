"""
FlowIO Example: Creating FCS Files
Create new FCS files from NumPy arrays.
"""
from flowio import create_fcs
import numpy as np


def create_simple_fcs(output_path: str):
    """
    Create a simple FCS file with synthetic data.
    
    Args:
        output_path: Output file path
    """
    # Create synthetic event data
    np.random.seed(42)
    num_events = 10000
    num_channels = 6
    
    # Generate data with different distributions
    events = np.zeros((num_events, num_channels))
    
    # FSC-A (Forward Scatter) - cell size
    events[:, 0] = np.random.lognormal(4, 0.3, num_events) * 100
    
    # SSC-A (Side Scatter) - granularity
    events[:, 1] = np.random.lognormal(3.5, 0.4, num_events) * 100
    
    # FL1-A (FITC) - fluorescence
    events[:, 2] = np.random.exponential(50, num_events)
    
    # FL2-A (PE) - fluorescence
    events[:, 3] = np.random.exponential(40, num_events)
    
    # FL3-A (PerCP) - fluorescence
    events[:, 4] = np.random.exponential(30, num_events)
    
    # Time
    events[:, 5] = np.linspace(0, 600, num_events)
    
    # Channel names
    channel_names = ['FSC-A', 'SSC-A', 'FL1-A', 'FL2-A', 'FL3-A', 'Time']
    descriptive_names = [
        'Forward Scatter',
        'Side Scatter',
        'FITC',
        'PE',
        'PerCP',
        'Time'
    ]
    
    # Add metadata
    metadata = {
        '$SRC': 'Python FlowIO',
        '$DATE': '19-OCT-2025',
        '$CYT': 'Synthetic Instrument',
        '$INST': 'Laboratory A',
        '$OP': 'Researcher'
    }
    
    # Create FCS file
    create_fcs(
        output_path,
        events,
        channel_names,
        opt_channel_names=descriptive_names,
        metadata=metadata
    )
    
    print(f"Created FCS file: {output_path}")
    print(f"  Events: {num_events:,}")
    print(f"  Channels: {num_channels}")


def create_from_existing_template(template_path: str, output_path: str):
    """
    Create a new FCS file using an existing file as template.
    
    Args:
        template_path: Path to template FCS file
        output_path: Output file path
    """
    from flowio import FlowData
    
    # Read template
    template = FlowData(template_path)
    
    # Create new data with same structure
    events = template.as_array(preprocess=False)
    
    # Modify data (e.g., add noise, scale)
    modified_events = events * 1.1 + np.random.normal(0, 10, events.shape)
    
    # Create new file with template metadata
    create_fcs(
        output_path,
        modified_events,
        template.pnn_labels,
        opt_channel_names=template.pns_labels,
        metadata={**template.text, '$SRC': 'Modified data'}
    )
    
    print(f"Created modified FCS file: {output_path}")


if __name__ == "__main__":
    print("FlowIO FCS Creator Example")
    # create_simple_fcs('output.fcs')
