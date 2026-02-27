"""
FlowIO Example: Reading FCS Files
Read and parse Flow Cytometry Standard (FCS) files.
"""
from flowio import FlowData
import numpy as np


def read_fcs_file(file_path: str):
    """
    Read and display information from an FCS file.
    
    Args:
        file_path: Path to FCS file
    """
    flow_data = FlowData(file_path)
    
    print(f"File: {flow_data.name}")
    print(f"FCS Version: {flow_data.version}")
    print(f"Events: {flow_data.event_count:,}")
    print(f"Channels: {flow_data.channel_count}")
    
    print("\nChannel Information:")
    for i, (pnn, pns) in enumerate(zip(flow_data.pnn_labels, flow_data.pns_labels)):
        ch_type = "scatter" if i in flow_data.scatter_indices else \
                  "fluoro" if i in flow_data.fluoro_indices else \
                  "time" if i == flow_data.time_index else "other"
        print(f"  [{i}] {pnn:10s} | {pns:30s} | {ch_type}")
    
    # Get event data as NumPy array
    events = flow_data.as_array()
    print(f"\nEvent data shape: {events.shape}")
    print(f"Data type: {events.dtype}")
    
    return flow_data, events


def extract_metadata(flow_data: FlowData):
    """
    Extract metadata from FCS file.
    
    Args:
        flow_data: FlowData object
        
    Returns:
        Dictionary of metadata
    """
    metadata = {
        'filename': flow_data.name,
        'fcs_version': flow_data.version,
        'event_count': flow_data.event_count,
        'channel_count': flow_data.channel_count,
        'date': flow_data.text.get('$DATE', 'Unknown'),
        'instrument': flow_data.text.get('$CYT', 'Unknown'),
        'source': flow_data.text.get('$SRC', 'Unknown'),
        'pnn_labels': flow_data.pnn_labels,
        'pns_labels': flow_data.pns_labels
    }
    
    print("\nMetadata Summary:")
    for key, value in metadata.items():
        if isinstance(value, list):
            print(f"  {key}: [{len(value)} items]")
        else:
            print(f"  {key}: {value}")
    
    return metadata


def analyze_channel_statistics(flow_data: FlowData):
    """
    Calculate statistics for each channel.
    
    Args:
        flow_data: FlowData object
        
    Returns:
        Statistics dictionary
    """
    events = flow_data.as_array()
    
    stats = []
    for i, name in enumerate(flow_data.pnn_labels):
        channel_data = events[:, i]
        stats.append({
            'channel': name,
            'mean': float(np.mean(channel_data)),
            'median': float(np.median(channel_data)),
            'std': float(np.std(channel_data)),
            'min': float(np.min(channel_data)),
            'max': float(np.max(channel_data))
        })
    
    print("\nChannel Statistics:")
    for s in stats[:5]:  # Show first 5
        print(f"  {s['channel']}: mean={s['mean']:.2f}, std={s['std']:.2f}")
    
    return stats


if __name__ == "__main__":
    print("FlowIO FCS Reader Example")
    # flow_data, events = read_fcs_file('sample.fcs')
    # metadata = extract_metadata(flow_data)
    # stats = analyze_channel_statistics(flow_data)
