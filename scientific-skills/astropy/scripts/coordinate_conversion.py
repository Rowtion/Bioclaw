#!/usr/bin/env python3
"""
Astronomical Coordinate Conversion

This script demonstrates converting between celestial coordinate
systems using astropy.

Usage:
    python coordinate_conversion.py
"""

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time


def basic_coordinate_creation():
    """Create coordinates in different formats."""
    print("\n" + "="*60)
    print("Creating SkyCoord Objects")
    print("="*60)
    
    # Method 1: Using degrees
    coord1 = SkyCoord(ra=10.684*u.degree, dec=41.269*u.degree, frame='icrs')
    print(f"Degrees: RA={coord1.ra.degree:.3f}°, Dec={coord1.dec.degree:.3f}°")
    
    # Method 2: Using sexagesimal strings
    coord2 = SkyCoord('00h42m44.3s', '+41d16m09s', frame='icrs')
    print(f"Sexagesimal: {coord2.to_string('hmsdms')}")
    
    # Method 3: From name (requires internet)
    try:
        coord3 = SkyCoord.from_name('M31')
        print(f"M31: RA={coord3.ra.degree:.3f}°, Dec={coord3.dec.degree:.3f}°")
    except Exception as e:
        print(f"Could not resolve name (internet required): {e}")


def coordinate_transformations():
    """Transform between coordinate systems."""
    print("\n" + "="*60)
    print("Coordinate Transformations")
    print("="*60)
    
    # Create coordinate in ICRS
    coord_icrs = SkyCoord(ra=266.41683*u.degree, dec=-29.00781*u.degree, frame='icrs')
    print(f"ICRS: RA={coord_icrs.ra.degree:.5f}°, Dec={coord_icrs.dec.degree:.5f}°")
    
    # Transform to Galactic
    coord_gal = coord_icrs.galactic
    print(f"Galactic: l={coord_gal.l.degree:.5f}°, b={coord_gal.b.degree:.5f}°")
    
    # Transform to FK5
    coord_fk5 = coord_icrs.fk5
    print(f"FK5: RA={coord_fk5.ra.degree:.5f}°, Dec={coord_fk5.dec.degree:.5f}°")


def altaz_conversion():
    """Convert to AltAz coordinates for observation planning."""
    print("\n" + "="*60)
    print("AltAz Coordinate Conversion")
    print("="*60)
    
    # Target coordinate
    target = SkyCoord(ra=83.82208*u.degree, dec=-5.39111*u.degree, frame='icrs')
    print(f"Target: RA={target.ra.degree:.5f}°, Dec={target.dec.degree:.5f}°")
    
    # Observer location (e.g., Mauna Kea)
    location = EarthLocation(
        lat=19.8207*u.degree,
        lon=-155.4681*u.degree,
        height=4205*u.m
    )
    print(f"Observer: Mauna Kea")
    
    # Observation time
    obstime = Time('2024-06-15 08:00:00')
    print(f"Time: {obstime.iso}")
    
    # Convert to AltAz
    altaz_frame = AltAz(obstime=obstime, location=location)
    target_altaz = target.transform_to(altaz_frame)
    
    print(f"\nAltAz Coordinates:")
    print(f"  Altitude: {target_altaz.alt.degree:.2f}°")
    print(f"  Azimuth: {target_altaz.az.degree:.2f}°")
    print(f"  Airmass: {target_altaz.secz:.3f}")


def angular_separations():
    """Calculate angular separations between coordinates."""
    print("\n" + "="*60)
    print("Angular Separations")
    print("="*60)
    
    # Two coordinates
    coord1 = SkyCoord(ra=10*u.degree, dec=20*u.degree, frame='icrs')
    coord2 = SkyCoord(ra=11*u.degree, dec=21*u.degree, frame='icrs')
    
    # Separation
    sep = coord1.separation(coord2)
    print(f"Separation: {sep.degree:.5f}°")
    print(f"            {sep.arcminute:.3f} arcmin")
    print(f"            {sep.arcsecond:.3f} arcsec")
    
    # Position angle
    pa = coord1.position_angle(coord2)
    print(f"Position angle: {pa.degree:.2f}°")


def main():
    """Run coordinate examples."""
    print("Astropy Coordinate Conversion Examples")
    print("=" * 60)
    
    basic_coordinate_creation()
    coordinate_transformations()
    altaz_conversion()
    angular_separations()
    
    print("\n" + "="*60)
    print("Coordinate examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
