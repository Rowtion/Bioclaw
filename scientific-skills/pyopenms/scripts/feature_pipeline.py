#!/usr/bin/env python3
"""
PyOpenMS Feature Detection Pipeline
Complete workflow for LC-MS/MS feature detection and quantification.
"""

import argparse
import os
from pathlib import Path


def run_feature_detection_pipeline(
    input_file: str,
    output_dir: str,
    method: str = "centroided"
) -> dict:
    """
    Run complete feature detection pipeline.
    
    Args:
        input_file: Path to input mzML file
        output_dir: Output directory for results
        method: Feature detection method (centroided or profile)
        
    Returns:
        Dictionary with pipeline results
    """
    try:
        import pyopenms as ms
        import pandas as pd
    except ImportError:
        raise ImportError("pyopenms and pandas are required")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print(f"Loading {input_file}...")
    exp = ms.MSExperiment()
    ms.MzMLFile().load(input_file, exp)
    
    # Run feature detection
    print("Running feature detection...")
    ff = ms.FeatureFinder()
    param = ms.FeatureFinder().getDefaultParameters()
    features = ms.FeatureMap()
    
    ff.run(method, exp, features, param, ms.FeatureMap())
    
    # Save features
    feature_xml = os.path.join(output_dir, "features.featureXML")
    ms.FeatureXMLFile().store(feature_xml, features)
    print(f"Features saved to: {feature_xml}")
    
    # Export to DataFrame
    df = features.get_df()
    csv_path = os.path.join(output_dir, "features.csv")
    df.to_csv(csv_path, index=False)
    print(f"Features exported to: {csv_path}")
    
    # Generate summary
    summary = {
        'total_features': features.size(),
        'rt_range': (float(df['RT'].min()), float(df['RT'].max())) if 'RT' in df.columns else None,
        'mz_range': (float(df['MZ'].min()), float(df['MZ'].max())) if 'MZ' in df.columns else None,
        'feature_xml': feature_xml,
        'feature_csv': csv_path,
    }
    
    # Save summary
    summary_path = os.path.join(output_dir, "summary.txt")
    with open(summary_path, 'w') as f:
        f.write("# Feature Detection Summary\n\n")
        f.write(f"Input file: {input_file}\n")
        f.write(f"Total features detected: {summary['total_features']}\n")
        if summary['rt_range']:
            f.write(f"RT range: {summary['rt_range'][0]:.2f} - {summary['rt_range'][1]:.2f}\n")
        if summary['mz_range']:
            f.write(f"m/z range: {summary['mz_range'][0]:.4f} - {summary['mz_range'][1]:.4f}\n")
    
    return summary


def align_features(feature_files: list, output_file: str):
    """
    Align features across multiple runs.
    
    Args:
        feature_files: List of featureXML file paths
        output_file: Output consensusXML file path
    """
    try:
        import pyopenms as ms
    except ImportError:
        raise ImportError("pyopenms is required")
    
    # Load feature maps
    maps = []
    for f in feature_files:
        fm = ms.FeatureMap()
        ms.FeatureXMLFile().load(f, fm)
        maps.append(fm)
    
    # Align
    aligner = ms.MapAlignmentAlgorithmPoseClustering()
    
    # Link features
    linker = ms.FeatureGroupingAlgorithmQT()
    consensus = ms.ConsensusMap()
    linker.group(maps, consensus)
    
    # Save
    ms.ConsensusXMLFile().store(output_file, consensus)
    print(f"Consensus map saved to: {output_file}")
    print(f"Number of consensus features: {consensus.size()}")


def main():
    parser = argparse.ArgumentParser(description='PyOpenMS Feature Detection Pipeline')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Detect command
    detect_parser = subparsers.add_parser('detect', help='Detect features')
    detect_parser.add_argument('input', help='Input mzML file')
    detect_parser.add_argument('-o', '--output', required=True, help='Output directory')
    detect_parser.add_argument('--method', default='centroided', 
                               choices=['centroided', 'profile'],
                               help='Detection method')
    
    # Align command
    align_parser = subparsers.add_parser('align', help='Align features across runs')
    align_parser.add_argument('features', nargs='+', help='Input featureXML files')
    align_parser.add_argument('-o', '--output', required=True, help='Output consensusXML file')
    
    args = parser.parse_args()
    
    if args.command == 'detect':
        results = run_feature_detection_pipeline(args.input, args.output, args.method)
        print("\nPipeline complete!")
        print(f"  Features detected: {results['total_features']}")
        
    elif args.command == 'align':
        align_features(args.features, args.output)
        
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
