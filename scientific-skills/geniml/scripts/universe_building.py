"""
Geniml Example: Universe Building
Build consensus peak sets from BED file collections.
"""
import subprocess
from geniml.universe import UniverseBuilder


def generate_coverage(bed_files: list, chrom_sizes: str, output_dir: str):
    """
    Generate coverage tracks from BED files.
    
    Args:
        bed_files: List of BED file paths
        chrom_sizes: Chromosome sizes file
        output_dir: Output directory
    """
    # Combine BED files
    combined = f"{output_dir}/combined.bed"
    with open(combined, 'w') as outfile:
        for bed in bed_files:
            with open(bed) as infile:
                outfile.write(infile.read())
    
    # Generate coverage with uniwig
    subprocess.run([
        'uniwig', '-m', '25',
        combined, chrom_sizes, f"{output_dir}/coverage/"
    ])
    
    print(f"Coverage tracks generated in {output_dir}/coverage/")


def build_universe_cc(coverage_folder: str, output_file: str,
                      cutoff: int = 5, merge: int = 100):
    """
    Build universe using Coverage Cutoff method.
    
    Args:
        coverage_folder: Folder with coverage tracks
        output_file: Output BED file
        cutoff: Coverage cutoff
        merge: Merge distance
    """
    import subprocess
    subprocess.run([
        'geniml', 'universe', 'build', 'cc',
        '--coverage-folder', coverage_folder,
        '--output-file', output_file,
        '--cutoff', str(cutoff),
        '--merge', str(merge),
        '--filter-size', '50'
    ])
    
    print(f"Universe built: {output_file}")


def evaluate_universe(universe_file: str, coverage_folder: str, 
                      bed_folder: str):
    """
    Evaluate universe quality.
    
    Args:
        universe_file: Universe BED file
        coverage_folder: Coverage folder
        bed_folder: Folder with original BED files
    """
    import subprocess
    subprocess.run([
        'geniml', 'universe', 'evaluate',
        '--universe', universe_file,
        '--coverage-folder', coverage_folder,
        '--bed-folder', bed_folder
    ])


if __name__ == "__main__":
    print("Geniml Universe Building Example")
    # build_universe_cc('coverage/', 'universe.bed')
