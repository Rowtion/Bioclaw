#!/usr/bin/env python3
"""
PyOpenMS Utility Scripts
Utilities for mass spectrometry data analysis and proteomics workflows.
"""

import argparse
import os
from pathlib import Path


def analyze_mzml_file(input_file: str, output_summary: str = None) -> dict:
    """
    Analyze an mzML file and extract basic statistics.
    
    Args:
        input_file: Path to mzML file
        output_summary: Optional path to save summary report
        
    Returns:
        Dictionary containing file statistics
    """
    try:
        import pyopenms as ms
    except ImportError:
        raise ImportError("pyopenms is required. Install with: uv pip install pyopenms")
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"File not found: {input_file}")
    
    # Load experiment
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)
    
    # Collect statistics
    stats = {
        'num_spectra': exp.getNrSpectra(),
        'num_chromatograms': exp.getNrChromatograms(),
        'ms_levels': [],
        'rt_range': (float('inf'), float('-inf')),
        'mz_ranges': [],
    }
    
    ms_levels = set()
    min_rt, max_rt = float('inf'), float('-inf')
    
    for spectrum in exp:
        ms_level = spectrum.getMSLevel()
        ms_levels.add(ms_level)
        rt = spectrum.getRT()
        min_rt = min(min_rt, rt)
        max_rt = max(max_rt, rt)
        
        mz, intensity = spectrum.get_peaks()
        if len(mz) > 0:
            stats['mz_ranges'].append({
                'ms_level': ms_level,
                'min_mz': float(min(mz)),
                'max_mz': float(max(mz)),
                'num_peaks': len(mz)
            })
    
    stats['ms_levels'] = sorted(list(ms_levels))
    stats['rt_range'] = (min_rt, max_rt)
    
    # Save summary if requested
    if output_summary:
        with open(output_summary, 'w') as f:
            f.write("# mzML File Analysis Summary\n\n")
            f.write(f"**File**: {input_file}\n\n")
            f.write(f"**Number of Spectra**: {stats['num_spectra']}\n")
            f.write(f"**Number of Chromatograms**: {stats['num_chromatograms']}\n")
            f.write(f"**MS Levels**: {', '.join(map(str, stats['ms_levels']))}\n")
            f.write(f"**RT Range**: {min_rt:.2f} - {max_rt:.2f}\n\n")
            f.write("## m/z Ranges by Spectrum\n\n")
            for i, mz_range in enumerate(stats['mz_ranges'][:10]):  # First 10
                f.write(f"Spectrum {i}: MS{mz_range['ms_level']}, ")
                f.write(f"m/z {mz_range['min_mz']:.2f}-{mz_range['max_mz']:.2f}, ")
                f.write(f"{mz_range['num_peaks']} peaks\n")
    
    return stats


def convert_mzml_to_featurexml(input_file: str, output_file: str, params: dict = None):
    """
    Convert mzML to featureXML format using feature detection.
    
    Args:
        input_file: Path to input mzML file
        output_file: Path to output featureXML file
        params: Optional feature finding parameters
    """
    try:
        import pyopenms as ms
    except ImportError:
        raise ImportError("pyopenms is required. Install with: uv pip install pyopenms")
    
    # Load experiment
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)
    
    # Create feature finder
    ff = ms.FeatureFinder()
    ff_param = ms.FeatureFinder().getDefaultParameters()
    
    if params:
        for key, value in params.items():
            if key in ff_param.keys():
                ff_param.setValue(key, value)
    
    features = ms.FeatureMap()
    
    # Run feature detection
    ff.run("centroided", exp, features, ff_param, ms.FeatureMap())
    
    # Save features
    ms.FeatureXMLFile().store(output_file, features)
    
    print(f"Features detected: {features.size()}")
    print(f"Saved to: {output_file}")
    
    return features.size()


def extract_ms2_spectra(input_file: str, output_file: str) -> int:
    """
    Extract MS2 spectra from mzML file.
    
    Args:
        input_file: Path to input mzML file
        output_file: Path to output mzML file
        
    Returns:
        Number of MS2 spectra extracted
    """
    try:
        import pyopenms as ms
    except ImportError:
        raise ImportError("pyopenms is required. Install with: uv pip install pyopenms")
    
    # Load experiment
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)
    
    # Filter MS2 spectra
    ms2_exp = ms.MSExperiment()
    for spectrum in exp:
        if spectrum.getMSLevel() == 2:
            ms2_exp.addSpectrum(spectrum)
    
    # Save
    ms.MzMLFile().store(output_file, ms2_exp)
    
    num_ms2 = ms2_exp.getNrSpectra()
    print(f"Extracted {num_ms2} MS2 spectra to {output_file}")
    
    return num_ms2


def calculate_tic(input_file: str, output_file: str = None):
    """
    Calculate Total Ion Chromatogram (TIC) from mzML file.
    
    Args:
        input_file: Path to input mzML file
        output_file: Optional path to save TIC data as CSV
    """
    try:
        import pyopenms as ms
    except ImportError:
        raise ImportError("pyopenms is required. Install with: uv pip install pyopenms")
    
    # Load experiment
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)
    
    # Calculate TIC
    tic_data = []
    for spectrum in exp:
        rt = spectrum.getRT()
        mz, intensity = spectrum.get_peaks()
        tic = sum(intensity) if len(intensity) > 0 else 0
        tic_data.append((rt, tic))
    
    # Save to CSV if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write("Retention_Time,TIC\n")
            for rt, tic in tic_data:
                f.write(f"{rt},{tic}\n")
        print(f"TIC data saved to: {output_file}")
    
    return tic_data


def filter_by_mz_range(input_file: str, output_file: str, min_mz: float, max_mz: float):
    """
    Filter spectra by m/z range.
    
    Args:
        input_file: Path to input mzML file
        output_file: Path to output mzML file
        min_mz: Minimum m/z value
        max_mz: Maximum m/z value
    """
    try:
        import pyopenms as ms
    except ImportError:
        raise ImportError("pyopenms is required. Install with: uv pip install pyopenms")
    
    # Load experiment
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)
    
    # Filter spectra
    filtered = ms.MSExperiment()
    for spectrum in exp:
        mz, intensity = spectrum.get_peaks()
        mask = (mz >= min_mz) & (mz <= max_mz)
        
        if mask.any():
            new_spec = ms.MSSpectrum()
            new_spec.setRT(spectrum.getRT())
            new_spec.setMSLevel(spectrum.getMSLevel())
            new_spec.set_peaks((mz[mask], intensity[mask]))
            filtered.addSpectrum(new_spec)
    
    # Save
    ms.MzMLFile().store(output_file, filtered)
    print(f"Filtered spectra: {filtered.getNrSpectra()}")
    print(f"Saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='PyOpenMS Utilities')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze mzML file')
    analyze_parser.add_argument('input', help='Input mzML file')
    analyze_parser.add_argument('-o', '--output', help='Output summary file')
    
    # Extract MS2 command
    ms2_parser = subparsers.add_parser('extract-ms2', help='Extract MS2 spectra')
    ms2_parser.add_argument('input', help='Input mzML file')
    ms2_parser.add_argument('output', help='Output mzML file')
    
    # TIC command
    tic_parser = subparsers.add_parser('tic', help='Calculate Total Ion Chromatogram')
    tic_parser.add_argument('input', help='Input mzML file')
    tic_parser.add_argument('-o', '--output', help='Output CSV file')
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter by m/z range')
    filter_parser.add_argument('input', help='Input mzML file')
    filter_parser.add_argument('output', help='Output mzML file')
    filter_parser.add_argument('--min-mz', type=float, required=True, help='Minimum m/z')
    filter_parser.add_argument('--max-mz', type=float, required=True, help='Maximum m/z')
    
    args = parser.parse_args()
    
    if args.command == 'analyze':
        stats = analyze_mzml_file(args.input, args.output)
        print(f"\nAnalysis complete:")
        print(f"  Spectra: {stats['num_spectra']}")
        print(f"  MS Levels: {stats['ms_levels']}")
        print(f"  RT Range: {stats['rt_range'][0]:.2f} - {stats['rt_range'][1]:.2f}")
        
    elif args.command == 'extract-ms2':
        extract_ms2_spectra(args.input, args.output)
        
    elif args.command == 'tic':
        calculate_tic(args.input, args.output)
        
    elif args.command == 'filter':
        filter_by_mz_range(args.input, args.output, args.min_mz, args.max_mz)
        
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
