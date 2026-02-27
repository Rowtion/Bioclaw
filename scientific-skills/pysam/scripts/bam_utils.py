#!/usr/bin/env python3
"""
PySAM BAM/CRAM Utilities
Utilities for working with genomic alignment files.
"""

import argparse
import os
from pathlib import Path
from typing import Optional


def get_bam_stats(bam_file: str) -> dict:
    """
    Get basic statistics from a BAM file.
    
    Args:
        bam_file: Path to BAM/CRAM file
        
    Returns:
        Dictionary with BAM statistics
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required. Install with: uv pip install pysam")
    
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        stats = {
            'filename': bam_file,
            'nreferences': samfile.nreferences,
            'references': list(samfile.references),
            'lengths': list(samfile.lengths),
        }
        
        # Count reads (may be slow for large files)
        total_reads = 0
        mapped_reads = 0
        unmapped_reads = 0
        
        for read in samfile:
            total_reads += 1
            if read.is_unmapped:
                unmapped_reads += 1
            else:
                mapped_reads += 1
                
        stats['total_reads'] = total_reads
        stats['mapped_reads'] = mapped_reads
        stats['unmapped_reads'] = unmapped_reads
        stats['mapping_rate'] = mapped_reads / total_reads if total_reads > 0 else 0
        
    return stats


def extract_region(bam_file: str, chrom: str, start: int, end: int, output_file: str):
    """
    Extract reads from a specific genomic region.
    
    Args:
        bam_file: Input BAM/CRAM file
        chrom: Chromosome name
        start: Start position (0-based)
        end: End position
        output_file: Output BAM file
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        with pysam.AlignmentFile(output_file, "wb", template=infile) as outfile:
            for read in infile.fetch(chrom, start, end):
                outfile.write(read)
    
    # Index output file
    pysam.index(output_file)
    print(f"Extracted region {chrom}:{start}-{end} to {output_file}")


def calculate_coverage(bam_file: str, chrom: str, start: int, end: int, output_file: str = None):
    """
    Calculate per-base coverage for a region.
    
    Args:
        bam_file: Input BAM file
        chrom: Chromosome name
        start: Start position
        end: End position
        output_file: Optional output file for coverage data
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        coverage = []
        for pileupcolumn in samfile.pileup(chrom, start, end):
            pos = pileupcolumn.reference_pos
            if start <= pos < end:
                coverage.append((pos, pileupcolumn.nsegments))
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write("Position,Coverage\n")
            for pos, cov in coverage:
                f.write(f"{pos},{cov}\n")
        print(f"Coverage data saved to: {output_file}")
    
    return coverage


def filter_bam(bam_file: str, output_file: str, min_quality: int = 20, 
               remove_duplicates: bool = False):
    """
    Filter BAM file by mapping quality and other criteria.
    
    Args:
        bam_file: Input BAM file
        output_file: Output BAM file
        min_quality: Minimum mapping quality
        remove_duplicates: Remove duplicate reads
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.AlignmentFile(bam_file, "rb") as infile:
        with pysam.AlignmentFile(output_file, "wb", template=infile) as outfile:
            for read in infile:
                if read.mapping_quality < min_quality:
                    continue
                if remove_duplicates and read.is_duplicate:
                    continue
                outfile.write(read)
    
    # Index output file
    pysam.index(output_file)
    print(f"Filtered BAM saved to: {output_file}")


def bam_to_fastq(bam_file: str, output_fastq: str, output_fastq2: str = None):
    """
    Convert BAM to FASTQ format.
    
    Args:
        bam_file: Input BAM file
        output_fastq: Output FASTQ file (R1 for paired-end)
        output_fastq2: Optional output FASTQ file for R2 (paired-end)
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_fastq, 'w') as fq1:
            if output_fastq2:
                fq2 = open(output_fastq2, 'w')
            
            for read in bam:
                if read.is_unmapped:
                    continue
                    
                # Convert to FASTQ
                fastq_entry = f"@{read.query_name}\n{read.query_sequence}\n+\n"
                
                # Get quality scores
                if read.query_qualities:
                    qual = ''.join([chr(q + 33) for q in read.query_qualities])
                else:
                    qual = 'I' * len(read.query_sequence)
                
                fastq_entry += qual + "\n"
                
                # Write to appropriate file
                if output_fastq2 and read.is_read2:
                    fq2.write(fastq_entry)
                else:
                    fq1.write(fastq_entry)
            
            if output_fastq2:
                fq2.close()
    
    print(f"Converted {bam_file} to FASTQ")


def count_reads_per_chromosome(bam_file: str, output_file: str = None):
    """
    Count reads per chromosome.
    
    Args:
        bam_file: Input BAM file
        output_file: Optional output file for results
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        counts = {}
        for read in samfile:
            if not read.is_unmapped:
                ref_name = samfile.get_reference_name(read.reference_id)
                counts[ref_name] = counts.get(ref_name, 0) + 1
    
    # Sort by count (descending)
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write("Chromosome,Read_Count\n")
            for chrom, count in sorted_counts:
                f.write(f"{chrom},{count}\n")
        print(f"Read counts saved to: {output_file}")
    
    return sorted_counts


def main():
    parser = argparse.ArgumentParser(description='PySAM BAM Utilities')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Get BAM statistics')
    stats_parser.add_argument('input', help='Input BAM file')
    
    # Extract command
    extract_parser = subparsers.add_parser('extract', help='Extract region')
    extract_parser.add_argument('input', help='Input BAM file')
    extract_parser.add_argument('chrom', help='Chromosome')
    extract_parser.add_argument('start', type=int, help='Start position (0-based)')
    extract_parser.add_argument('end', type=int, help='End position')
    extract_parser.add_argument('-o', '--output', required=True, help='Output BAM file')
    
    # Coverage command
    cov_parser = subparsers.add_parser('coverage', help='Calculate coverage')
    cov_parser.add_argument('input', help='Input BAM file')
    cov_parser.add_argument('chrom', help='Chromosome')
    cov_parser.add_argument('start', type=int, help='Start position')
    cov_parser.add_argument('end', type=int, help='End position')
    cov_parser.add_argument('-o', '--output', help='Output CSV file')
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter BAM')
    filter_parser.add_argument('input', help='Input BAM file')
    filter_parser.add_argument('-o', '--output', required=True, help='Output BAM file')
    filter_parser.add_argument('-q', '--min-quality', type=int, default=20,
                               help='Minimum mapping quality')
    filter_parser.add_argument('--remove-duplicates', action='store_true',
                               help='Remove duplicate reads')
    
    # Convert command
    conv_parser = subparsers.add_parser('convert', help='Convert to FASTQ')
    conv_parser.add_argument('input', help='Input BAM file')
    conv_parser.add_argument('-o', '--output', required=True, help='Output FASTQ file')
    conv_parser.add_argument('-o2', '--output2', help='Output FASTQ file for R2 (paired)')
    
    # Count command
    count_parser = subparsers.add_parser('count', help='Count reads per chromosome')
    count_parser.add_argument('input', help='Input BAM file')
    count_parser.add_argument('-o', '--output', help='Output CSV file')
    
    args = parser.parse_args()
    
    if args.command == 'stats':
        stats = get_bam_stats(args.input)
        print(f"\nBAM Statistics:")
        print(f"  Total reads: {stats['total_reads']:,}")
        print(f"  Mapped reads: {stats['mapped_reads']:,}")
        print(f"  Unmapped reads: {stats['unmapped_reads']:,}")
        print(f"  Mapping rate: {stats['mapping_rate']:.2%}")
        print(f"  References: {stats['nreferences']}")
        
    elif args.command == 'extract':
        extract_region(args.input, args.chrom, args.start, args.end, args.output)
        
    elif args.command == 'coverage':
        calculate_coverage(args.input, args.chrom, args.start, args.end, args.output)
        
    elif args.command == 'filter':
        filter_bam(args.input, args.output, args.min_quality, args.remove_duplicates)
        
    elif args.command == 'convert':
        bam_to_fastq(args.input, args.output, args.output2)
        
    elif args.command == 'count':
        counts = count_reads_per_chromosome(args.input, args.output)
        print("\nTop 10 chromosomes by read count:")
        for chrom, count in counts[:10]:
            print(f"  {chrom}: {count:,}")
            
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
