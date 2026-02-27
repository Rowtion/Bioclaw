#!/usr/bin/env python3
"""
FITS File Processing

This script demonstrates reading, writing, and analyzing
FITS (Flexible Image Transport System) files using astropy.

Usage:
    python fits_processing.py
"""

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt


def create_sample_fits():
    """Create a sample FITS file for demonstration."""
    print("\n" + "="*60)
    print("Creating Sample FITS File")
    print("="*60)
    
    # Create sample image data
    data = np.random.normal(1000, 50, size=(100, 100))
    
    # Add some features
    y, x = np.indices((100, 100))
    data += 5000 * np.exp(-((x-50)**2 + (y-50)**2) / 100)
    
    # Create HDU
    hdu = fits.PrimaryHDU(data)
    
    # Add header keywords
    hdu.header['OBJECT'] = 'SAMPLE'
    hdu.header['EXPTIME'] = 300.0
    hdu.header['FILTER'] = 'V'
    hdu.header['DATE-OBS'] = '2024-01-15T12:00:00'
    
    # Create HDU list and write
    hdul = fits.HDUList([hdu])
    output_file = 'sample_image.fits'
    hdul.writeto(output_file, overwrite=True)
    hdul.close()
    
    print(f"Created: {output_file}")
    return output_file


def read_fits_file(filename):
    """Read and display FITS file information."""
    print("\n" + "="*60)
    print(f"Reading FITS File: {filename}")
    print("="*60)
    
    with fits.open(filename) as hdul:
        # Display file structure
        print("\nFile Structure:")
        hdul.info()
        
        # Access primary HDU
        primary = hdul[0]
        
        print(f"\nData shape: {primary.data.shape}")
        print(f"Data type: {primary.data.dtype}")
        
        print("\nHeader Keywords:")
        for key in ['OBJECT', 'EXPTIME', 'FILTER', 'DATE-OBS']:
            if key in primary.header:
                print(f"  {key} = {primary.header[key]}")
        
        # Calculate statistics
        data = primary.data
        print(f"\nStatistics:")
        print(f"  Mean: {np.mean(data):.2f}")
        print(f"  Std: {np.std(data):.2f}")
        print(f"  Min: {np.min(data):.2f}")
        print(f"  Max: {np.max(data):.2f}")


def modify_fits_header(filename):
    """Modify FITS header."""
    print("\n" + "="*60)
    print(f"Modifying FITS Header: {filename}")
    print("="*60)
    
    with fits.open(filename, mode='update') as hdul:
        header = hdul[0].header
        
        # Add new keywords
        header['HISTORY'] = 'Modified by fits_processing.py'
        header['COMMENT'] = 'This is a sample FITS file'
        header['PROCESS'] = 'Demo'
        
        print("Added header keywords:")
        print(f"  PROCESS = {header['PROCESS']}")
        
        # Changes are saved automatically when file closes
    
    print("Header updated successfully")


def visualize_fits(filename):
    """Visualize FITS image data."""
    print("\n" + "="*60)
    print(f"Visualizing FITS Image: {filename}")
    print("="*60)
    
    with fits.open(filename) as hdul:
        data = hdul[0].data
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Original image
        im1 = axes[0].imshow(data, origin='lower', cmap='gray')
        axes[0].set_title('FITS Image')
        axes[0].set_xlabel('X pixel')
        axes[0].set_ylabel('Y pixel')
        plt.colorbar(im1, ax=axes[0])
        
        # Histogram
        axes[1].hist(data.flatten(), bins=50, color='steelblue', edgecolor='black')
        axes[1].set_title('Pixel Value Distribution')
        axes[1].set_xlabel('Pixel Value')
        axes[1].set_ylabel('Frequency')
        axes[1].set_yscale('log')
        
        plt.tight_layout()
        plt.savefig('fits_visualization.png', dpi=150)
        print("Saved visualization to: fits_visualization.png")
        plt.close()


def main():
    """Run FITS processing examples."""
    print("Astropy FITS Processing Examples")
    print("=" * 60)
    
    # Create sample file
    filename = create_sample_fits()
    
    # Read file
    read_fits_file(filename)
    
    # Modify header
    modify_fits_header(filename)
    
    # Visualize
    visualize_fits(filename)
    
    # Cleanup
    import os
    os.remove(filename)
    print(f"\nCleaned up: {filename}")
    
    print("\n" + "="*60)
    print("FITS processing examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
