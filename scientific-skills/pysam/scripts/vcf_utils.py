#!/usr/bin/env python3
"""
PySAM VCF/BCF Utilities
Utilities for working with variant call files.
"""

import argparse
import os
from pathlib import Path


def get_vcf_stats(vcf_file: str) -> dict:
    """
    Get statistics from VCF file.
    
    Args:
        vcf_file: Path to VCF/BCF file
        
    Returns:
        Dictionary with VCF statistics
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required. Install with: uv pip install pysam")
    
    with pysam.VariantFile(vcf_file) as vcf:
        # Get contigs from header
        contigs = list(vcf.header.contigs)
        samples = list(vcf.header.samples)
        
        # Count variants
        variant_count = 0
        snp_count = 0
        indel_count = 0
        
        for variant in vcf:
            variant_count += 1
            ref = variant.ref
            alts = variant.alts
            
            # Classify variant
            is_indel = False
            for alt in alts:
                if len(ref) != len(alt):
                    is_indel = True
                    break
            
            if is_indel:
                indel_count += 1
            else:
                snp_count += 1
    
    return {
        'filename': vcf_file,
        'contigs': contigs,
        'samples': samples,
        'variant_count': variant_count,
        'snp_count': snp_count,
        'indel_count': indel_count,
    }


def filter_vcf(input_file: str, output_file: str, 
               min_quality: float = 30.0,
               min_depth: int = 10,
               chrom: str = None):
    """
    Filter VCF file by various criteria.
    
    Args:
        input_file: Input VCF/BCF file
        output_file: Output VCF file
        min_quality: Minimum QUAL score
        min_depth: Minimum read depth
        chrom: Only keep variants on this chromosome
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.VariantFile(input_file) as vcf_in:
        with pysam.VariantFile(output_file, 'w', header=vcf_in.header) as vcf_out:
            for variant in vcf_in:
                # Filter by chromosome
                if chrom and variant.chrom != chrom:
                    continue
                
                # Filter by quality
                if variant.qual is not None and variant.qual < min_quality:
                    continue
                
                # Filter by depth (check INFO field)
                if 'DP' in variant.info:
                    if variant.info['DP'] < min_depth:
                        continue
                
                vcf_out.write(variant)
    
    print(f"Filtered VCF saved to: {output_file}")


def extract_sample_genotypes(vcf_file: str, output_file: str, 
                             samples: list = None, region: str = None):
    """
    Extract genotype information for samples.
    
    Args:
        vcf_file: Input VCF file
        output_file: Output CSV file
        samples: List of samples to extract (None = all)
        region: Genomic region (chrom:start-end)
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.VariantFile(vcf_file) as vcf:
        # Get samples
        if samples is None:
            samples = list(vcf.header.samples)
        
        # Parse region
        chrom, start, end = None, None, None
        if region:
            chrom, pos = region.split(':')
            start, end = map(int, pos.split('-'))
        
        with open(output_file, 'w') as out:
            # Write header
            header = ['CHROM', 'POS', 'ID', 'REF', 'ALT'] + samples
            out.write(','.join(header) + '\n')
            
            # Iterate variants
            for variant in vcf.fetch(chrom, start, end) if chrom else vcf:
                row = [
                    variant.chrom,
                    str(variant.pos),
                    str(variant.id) if variant.id else '.',
                    variant.ref,
                    ','.join(variant.alts) if variant.alts else '.'
                ]
                
                # Get genotypes for each sample
                for sample in samples:
                    if sample in variant.samples:
                        gt = variant.samples[sample].get('GT')
                        if gt:
                            gt_str = '/'.join(['.' if x is None else str(x) for x in gt])
                        else:
                            gt_str = './.'
                    else:
                        gt_str = './.'
                    row.append(gt_str)
                
                out.write(','.join(row) + '\n')
    
    print(f"Genotypes extracted to: {output_file}")


def vcf_to_bed(vcf_file: str, output_file: str, variant_type: str = 'all'):
    """
    Convert VCF to BED format.
    
    Args:
        vcf_file: Input VCF file
        output_file: Output BED file
        variant_type: Filter by type (snp, indel, all)
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    with pysam.VariantFile(vcf_file) as vcf:
        with open(output_file, 'w') as out:
            for variant in vcf:
                # Determine variant type
                ref = variant.ref
                alts = variant.alts
                
                is_indel = False
                for alt in alts:
                    if len(ref) != len(alt):
                        is_indel = True
                        break
                
                # Filter by type
                if variant_type == 'snp' and is_indel:
                    continue
                if variant_type == 'indel' and not is_indel:
                    continue
                
                # Calculate end position
                max_alt_len = max(len(alt) for alt in alts)
                end_pos = variant.pos + max(len(ref), max_alt_len) - 1
                
                # Write BED entry
                name = f"{variant.ref}>{','.join(variant.alts)}"
                out.write(f"{variant.chrom}\t{variant.pos-1}\t{end_pos}\t{name}\n")
    
    print(f"BED file saved to: {output_file}")


def annotate_variants(vcf_file: str, output_file: str, 
                      gene_bed: str = None, dbsnp: str = None):
    """
    Annotate variants with gene information and known variants.
    
    Args:
        vcf_file: Input VCF file
        output_file: Output annotated VCF
        gene_bed: BED file with gene coordinates
        dbsnp: VCF file with known variants
    """
    try:
        import pysam
    except ImportError:
        raise ImportError("pysam is required")
    
    # Load gene annotations if provided
    genes = []
    if gene_bed and os.path.exists(gene_bed):
        with open(gene_bed) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    genes.append({
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3]
                    })
    
    with pysam.VariantFile(vcf_file) as vcf_in:
        # Add new INFO fields
        header = vcf_in.header.copy()
        header.add_line('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">')
        header.add_line('##INFO=<ID=IMPACT,Number=1,Type=String,Description="Variant impact">')
        
        with pysam.VariantFile(output_file, 'w', header=header) as vcf_out:
            for variant in vcf_in:
                # Find overlapping genes
                overlapping_genes = []
                for gene in genes:
                    if (variant.chrom == gene['chrom'] and
                        variant.pos >= gene['start'] and
                        variant.pos <= gene['end']):
                        overlapping_genes.append(gene['name'])
                
                # Add annotation
                if overlapping_genes:
                    variant.info['GENE'] = ','.join(overlapping_genes)
                
                vcf_out.write(variant)
    
    print(f"Annotated VCF saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='PySAM VCF Utilities')
    subparsers = parser.add_subparsers(dest='command', help='Commands')
    
    # Stats command
    stats_parser = subparsers.add_parser('stats', help='Get VCF statistics')
    stats_parser.add_argument('input', help='Input VCF file')
    
    # Filter command
    filter_parser = subparsers.add_parser('filter', help='Filter VCF')
    filter_parser.add_argument('input', help='Input VCF file')
    filter_parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    filter_parser.add_argument('-q', '--min-quality', type=float, default=30.0,
                               help='Minimum QUAL score')
    filter_parser.add_argument('-d', '--min-depth', type=int, default=10,
                               help='Minimum depth')
    filter_parser.add_argument('-c', '--chrom', help='Chromosome to keep')
    
    # Genotypes command
    gt_parser = subparsers.add_parser('genotypes', help='Extract genotypes')
    gt_parser.add_argument('input', help='Input VCF file')
    gt_parser.add_argument('-o', '--output', required=True, help='Output CSV file')
    gt_parser.add_argument('-s', '--samples', nargs='+', help='Samples to extract')
    gt_parser.add_argument('-r', '--region', help='Region (chrom:start-end)')
    
    # Convert command
    conv_parser = subparsers.add_parser('convert', help='Convert to BED')
    conv_parser.add_argument('input', help='Input VCF file')
    conv_parser.add_argument('-o', '--output', required=True, help='Output BED file')
    conv_parser.add_argument('-t', '--type', choices=['snp', 'indel', 'all'], 
                             default='all', help='Variant type')
    
    # Annotate command
    ann_parser = subparsers.add_parser('annotate', help='Annotate variants')
    ann_parser.add_argument('input', help='Input VCF file')
    ann_parser.add_argument('-o', '--output', required=True, help='Output VCF file')
    ann_parser.add_argument('-g', '--genes', help='Gene BED file')
    ann_parser.add_argument('-d', '--dbsnp', help='dbSNP VCF file')
    
    args = parser.parse_args()
    
    if args.command == 'stats':
        stats = get_vcf_stats(args.input)
        print(f"\nVCF Statistics:")
        print(f"  Variants: {stats['variant_count']:,}")
        print(f"  SNPs: {stats['snp_count']:,}")
        print(f"  Indels: {stats['indel_count']:,}")
        print(f"  Samples: {', '.join(stats['samples'])}")
        print(f"  Contigs: {len(stats['contigs'])}")
        
    elif args.command == 'filter':
        filter_vcf(args.input, args.output, args.min_quality, args.min_depth, args.chrom)
        
    elif args.command == 'genotypes':
        extract_sample_genotypes(args.input, args.output, args.samples, args.region)
        
    elif args.command == 'convert':
        vcf_to_bed(args.input, args.output, args.type)
        
    elif args.command == 'annotate':
        annotate_variants(args.input, args.output, args.genes, args.dbsnp)
        
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
