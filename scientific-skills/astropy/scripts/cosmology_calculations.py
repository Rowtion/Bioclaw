#!/usr/bin/env python3
"""
Cosmological Calculations

This script demonstrates cosmological distance and time
calculations using astropy.cosmology.

Usage:
    python cosmology_calculations.py
"""

import numpy as np
import astropy.units as u
from astropy.cosmology import Planck18, WMAP9, z_at_value
import matplotlib.pyplot as plt


def basic_distances():
    """Calculate basic cosmological distances."""
    print("\n" + "="*60)
    print("Cosmological Distances")
    print("="*60)
    
    z = 1.0
    
    # Different distance measures
    print(f"\nAt redshift z = {z}:")
    print(f"  Luminosity distance: {Planck18.luminosity_distance(z):.2f}")
    print(f"  Angular diameter distance: {Planck18.angular_diameter_distance(z):.2f}")
    print(f"  Comoving distance: {Planck18.comoving_distance(z):.2f}")
    
    # Distance modulus
    dist_mod = Planck18.distmod(z)
    print(f"  Distance modulus: {dist_mod:.2f}")


def age_calculations():
    """Calculate ages and lookback times."""
    print("\n" + "="*60)
    print("Age and Lookback Time Calculations")
    print("="*60)
    
    # Age of universe today
    age_now = Planck18.age(0)
    print(f"\nAge of universe now: {age_now:.2f}")
    
    # Age at different redshifts
    redshifts = [0, 1, 2, 5, 10]
    print("\nAge of universe at different redshifts:")
    for z in redshifts:
        age = Planck18.age(z)
        print(f"  z={z:2d}: {age:.2f}")
    
    # Lookback time
    print("\nLookback time:")
    for z in redshifts[1:]:
        t_lb = Planck18.lookback_time(z)
        print(f"  z={z:2d}: {t_lb:.2f}")


def compare_cosmologies():
    """Compare different cosmological models."""
    print("\n" + "="*60)
    print("Comparing Cosmological Models")
    print("="*60)
    
    redshifts = np.linspace(0, 5, 100)
    
    # Calculate distances for different models
    planck_dist = [Planck18.luminosity_distance(z).value for z in redshifts]
    wmap_dist = [WMAP9.luminosity_distance(z).value for z in redshifts]
    
    print(f"\nAt z=2:")
    print(f"  Planck18: {Planck18.luminosity_distance(2):.2f}")
    print(f"  WMAP9: {WMAP9.luminosity_distance(2):.2f}")
    
    print(f"\nAt z=5:")
    print(f"  Planck18: {Planck18.luminosity_distance(5):.2f}")
    print(f"  WMAP9: {WMAP9.luminosity_distance(5):.2f}")
    
    # Plot comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(redshifts, planck_dist, label='Planck18', linewidth=2)
    ax.plot(redshifts, wmap_dist, label='WMAP9', linewidth=2, linestyle='--')
    ax.set_xlabel('Redshift (z)')
    ax.set_ylabel('Luminosity Distance (Mpc)')
    ax.set_title('Cosmological Model Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('cosmology_comparison.png', dpi=150)
    print("\nSaved plot to: cosmology_comparison.png")
    plt.close()


def hubble_parameter():
    """Calculate Hubble parameter evolution."""
    print("\n" + "="*60)
    print("Hubble Parameter Evolution")
    print("="*60)
    
    print(f"\nHubble constant today: {Planck18.H0:.2f}")
    
    redshifts = [0, 1, 2, 5, 10]
    print("\nHubble parameter at different redshifts:")
    for z in redshifts:
        H = Planck18.H(z)
        print(f"  z={z:2d}: H = {H:.2f}")


def main():
    """Run cosmology examples."""
    print("Astropy Cosmology Calculations")
    print("=" * 60)
    
    basic_distances()
    age_calculations()
    compare_cosmologies()
    hubble_parameter()
    
    print("\n" + "="*60)
    print("Cosmology examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
