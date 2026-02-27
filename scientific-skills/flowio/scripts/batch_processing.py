"""
FlowIO Example: Batch Processing
Batch process multiple FCS files.
"""
from pathlib import Path
import pandas as pd
from flowio import FlowData


def batch_extract_metadata(directory: str, pattern: str = "*.fcs"):
    """
    Extract metadata from all FCS files in a directory.
    
    Args:
        directory: Input directory
        pattern: File pattern to match
        
    Returns:
        DataFrame with metadata
    """
    summaries = []
    
    for fcs_path in Path(directory).glob(pattern):
        try:
            flow = FlowData(str(fcs_path), only_text=True)
            summaries.append({
                'filename': fcs_path.name,
                'version': flow.version,
                'events': flow.event_count,
                'channels': flow.channel_count,
                'date': flow.text.get('$DATE', 'N/A'),
                'instrument': flow.text.get('$CYT', 'N/A')
            })
            print(f"Processed: {fcs_path.name}")
        except Exception as e:
            print(f"Error processing {fcs_path.name}: {e}")
    
    df = pd.DataFrame(summaries)
    print(f"\nProcessed {len(df)} files")
    return df


def batch_convert_to_csv(input_dir: str, output_dir: str):
    """
    Convert all FCS files to CSV format.
    
    Args:
        input_dir: Input directory with FCS files
        output_dir: Output directory for CSV files
    """
    Path(output_dir).mkdir(exist_ok=True)
    
    for fcs_path in Path(input_dir).glob("*.fcs"):
        try:
            flow = FlowData(str(fcs_path))
            
            # Convert to DataFrame
            df = pd.DataFrame(
                flow.as_array(),
                columns=flow.pnn_labels
            )
            
            # Save to CSV
            output_file = Path(output_dir) / f"{fcs_path.stem}.csv"
            df.to_csv(output_file, index=False)
            
            print(f"Converted: {fcs_path.name} -> {output_file.name}")
        except Exception as e:
            print(f"Error converting {fcs_path.name}: {e}")


def filter_and_export(input_file: str, output_file: str, 
                      channel_idx: int = 0, threshold: float = 100):
    """
    Filter events based on channel threshold and export.
    
    Args:
        input_file: Input FCS file
        output_file: Output FCS file
        channel_idx: Channel index to filter on
        threshold: Minimum value threshold
    """
    flow = FlowData(input_file)
    events = flow.as_array(preprocess=False)
    
    # Apply filter
    mask = events[:, channel_idx] > threshold
    filtered_events = events[mask]
    
    print(f"Original events: {len(events):,}")
    print(f"Filtered events: {len(filtered_events):,}")
    
    # Create filtered FCS
    from flowio import create_fcs
    create_fcs(
        output_file,
        filtered_events,
        flow.pnn_labels,
        opt_channel_names=flow.pns_labels,
        metadata=flow.text
    )
    
    print(f"Saved to: {output_file}")


if __name__ == "__main__":
    print("FlowIO Batch Processing Example")
    # df = batch_extract_metadata('data/')
    # batch_convert_to_csv('fcs_files/', 'csv_files/')
