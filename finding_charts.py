#!/usr/bin/env python
"""
Finding Chart Generator for Pan-STARRS and Other Survey Images

This module creates astronomical finding charts from FITS images, highlighting
objects of interest and reference stars. It supports various features including
slit overlay, star labeling, and customizable chart appearance.

Author: E. Banados
Date: 2013-2025 
Version: 2.0 (Refactored for clarity and maintainability and more general use)
"""

from __future__ import division, print_function

import sys
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional, Tuple, List, Dict, Any

from astropy.io import ascii, fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy import wcs
from astropy.nddata.utils import Cutout2D
from astropy.visualization import ImageNormalize, ZScaleInterval

import aplpy


# Constants
ARCMIN_IN_DEG = 1.0 / 60.0
DEFAULT_RADIUS_ARCSEC = 4.0
DEFAULT_RADIUS_DEG = DEFAULT_RADIUS_ARCSEC / 3600.0

# Label position offsets in degrees
LABEL_POSITIONS = {
    "left": [10, 0],
    "right": [-10, 0],
    "top": [0, 5],
    "bottom": [0, -5],
}

EXAMPLES = '''
EXAMPLES:

1. Basic finding chart creation:
   fchart_ps1.py ../finals/idrops_allp1.stamps --suffix p1_

2. Highlighting QSO and Reference Stars:
   fchart_ps1.py PSO037.97065-28.83891_qso2chart.info \\
                 --stars_info PSO037.97065-28.83891_stars2chart.info \\
                 -o PSO037.97065-28.83891_fchart.png

3. Highlighting only the QSO:
   fchart_ps1.py PSO071.45075-02.33329_qso2chart.info \\
                 -o PSO071.45075-02.33329_fchart.png

4. Show slit overlay:
   fchart_ps1.py PSO002+25_P2_y_qso.info \\
                 --stars PSO002+25_P2_y_stars.info \\
                 --show_slit 1.3,120,330.23
'''


class FindingChartGenerator:
    """Main class for generating astronomical finding charts."""
    
    def __init__(self, label_position: str = 'left'):
        """
        Initialize the finding chart generator.
        
        Parameters
        ----------
        label_position : str
            Default position for star labels ('left', 'right', 'top', 'bottom')
        """
        self.label_position = label_position
        self.radius = DEFAULT_RADIUS_DEG
    
    def add_stars_info(self, 
                       stars_file: str, 
                       fig: aplpy.FITSFigure, 
                       ra: float, 
                       dec: float, 
                       show_default_slit: bool = False) -> None:
        """
        Add reference star information to the finding chart.
        
        Parameters
        ----------
        stars_file : str
            Path to file containing star information (ra, dec, label columns)
        fig : aplpy.FITSFigure
            The figure object to add stars to
        ra : float
            Right ascension of the target object (degrees)
        dec : float
            Declination of the target object (degrees)
        show_default_slit : bool
            If True, show a default slit at the PA of the star and target
        """
        stars = ascii.read(stars_file)
        nstars = len(stars)
        ypos = 0.055 * nstars
        
        ralabel, declabel = LABEL_POSITIONS[self.label_position]
        coord_qso = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
        
        for star_ra, star_dec, label in zip(stars['ra'], stars['dec'], stars['label']):
            self._process_single_star(
                fig, star_ra, star_dec, label, coord_qso, 
                ralabel, declabel, ypos, show_default_slit
            )
            ypos -= 0.095
    
    def _process_single_star(self, 
                            fig: aplpy.FITSFigure,
                            star_ra: float, 
                            star_dec: float, 
                            label: str,
                            coord_qso: SkyCoord,
                            ralabel: float,
                            declabel: float,
                            ypos: float,
                            show_default_slit: bool) -> None:
        """Process and display a single reference star."""
        print(f"PA and Distance calculated from star {label} to object")
        
        # Create star coordinates
        coord = SkyCoord(ra=star_ra, dec=star_dec, unit=(u.degree, u.degree))
        
        # Format coordinate strings
        rastring = coord.ra.to_string(precision=3, sep=":", unit=u.hour)
        decstring = coord.dec.to_string(precision=3, sep=":", unit=u.degree, alwayssign=True)
        
        # Calculate offsets
        ra_offset, dec_offset = coord.spherical_offsets_to(coord_qso)
        raoff = self._format_offset(ra_offset.to("arcsec").value, 'E', 'W')
        decoff = self._format_offset(dec_offset.to("arcsec").value, 'N', 'S')
        
        # Calculate position angle and separation
        PA = coord.position_angle(coord_qso).degree
        ang_sep_arcsec = coord.separation(coord_qso).arcsec
        
        # Adjust label position based on PA
        ralabel, declabel = self._adjust_label_position(PA, raoff, ralabel, declabel)
        
        # Add visual elements
        self._add_star_visual_elements(
            fig, star_ra, star_dec, label, ralabel, declabel, ang_sep_arcsec
        )
        
        # Add text information
        self._add_star_text_info(
            fig, label, rastring, decstring, raoff, decoff, PA, ang_sep_arcsec, ypos
        )
        
        # Show default slit if requested
        if show_default_slit:
            self.plot_slit(fig, coord_qso.ra.degree, coord_qso.dec.degree, (5.0, 150, PA))
    
    def _format_offset(self, offset: float, pos_suffix: str, neg_suffix: str) -> str:
        """Format offset value with appropriate direction suffix."""
        offset_rounded = np.round(offset, 2)
        if offset_rounded >= 0:
            return f"{offset_rounded} {pos_suffix}"
        else:
            return f"{abs(offset_rounded)} {neg_suffix}"
    
    def _adjust_label_position(self, 
                              PA: float, 
                              raoff: str, 
                              ralabel: float, 
                              declabel: float) -> Tuple[float, float]:
        """Adjust label position based on position angle."""
        if (30 < PA < 130) and raoff.endswith('E'):
            return LABEL_POSITIONS['bottom']
        elif (130 <= PA < 175) and raoff.endswith('E'):
            return LABEL_POSITIONS['right']
        elif (210 < PA < 270) and raoff.endswith('W'):
            return LABEL_POSITIONS['top']
        elif (PA > 350) and raoff.endswith('W'):
            return LABEL_POSITIONS['right']
        return ralabel, declabel
    
    def _add_star_visual_elements(self,
                                  fig: aplpy.FITSFigure,
                                  star_ra: float,
                                  star_dec: float,
                                  label: str,
                                  ralabel: float,
                                  declabel: float,
                                  ang_sep: float) -> None:
        """Add visual elements (rectangle and label) for a star."""
        # Add label
        fig.add_label(
            star_ra + ralabel * self.radius,
            star_dec + declabel * self.radius,
            label,
            color='blue',
            size='x-large',
            verticalalignment='center',
            family='serif'
        )
        
        # Determine rectangle size based on separation
        rectangle_size = 3 * self.radius if ang_sep >= 8.0 else 2 * self.radius
        
        # Show rectangle
        fig.show_rectangles(
            star_ra, star_dec,
            rectangle_size, rectangle_size,
            edgecolor='blue', lw=1
        )
    
    def _add_star_text_info(self,
                           fig: aplpy.FITSFigure,
                           label: str,
                           rastring: str,
                           decstring: str,
                           raoff: str,
                           decoff: str,
                           PA: float,
                           ang_sep: float,
                           ypos: float) -> None:
        """Add text information box for a star."""
        lines = [
            f"Star {label}",
            f"ra = {rastring}; dec = {decstring}",
            f"ra_off = {raoff} arcsec; dec_off = {decoff} arcsec",
            f"PA = {PA:.2f} deg; Ang. Sep.={ang_sep:.2f} arcsec"
        ]
        
        print("=" * 45)
        print("Pivot Star")
        print("\n".join(lines))
        
        # Add text box to plot
        ax = plt.gca()
        boxdict = dict(facecolor='white', alpha=0.5, edgecolor='none')
        ax.text(
            0.02, ypos, "\n".join(lines),
            transform=ax.transAxes,
            fontsize='small',
            bbox=boxdict
        )
    
    def plot_slit(self,
                  fig: aplpy.FITSFigure,
                  ra: float,
                  dec: float,
                  slit_params: Tuple[float, float, float]) -> None:
        """
        Plot a slit overlay on the finding chart.
        
        Parameters
        ----------
        fig : aplpy.FITSFigure
            The figure object
        ra : float
            Right ascension of the slit center (degrees)
        dec : float
            Declination of the slit center (degrees)
        slit_params : tuple
            (width, length, position_angle) in (arcsec, arcsec, degrees)
        """
        width, length, pa = slit_params
        slit_label = f'PA=${pa:.2f}$deg'
        print(f'Slit label: {slit_label}')
        
        # Convert to degrees
        width_deg = width / 3600.0
        length_deg = length / 3600.0
        
        # Show rectangle
        fig.show_rectangles(
            ra, dec, width_deg, length_deg,
            edgecolor='w', lw=1, angle=pa,
            coords_frame='world'
        )
        
        # Calculate label position
        ra_label, dec_label, rot_label = self._calculate_slit_label_position(
            ra, dec, pa
        )
        
        # Add label
        fig.add_label(
            ra_label, dec_label, slit_label,
            rotation=rot_label, size='large', color='w'
        )
    
    def _calculate_slit_label_position(self,
                                      ra: float,
                                      dec: float,
                                      pa: float) -> Tuple[float, float, float]:
        """Calculate optimal position for slit label."""
        rot_label = pa + 90.0
        dec_label = dec + 3 * self.radius
        ra_label = ra
        
        # Adjust based on position angle
        if (pa < 20) or (pa > 345):
            ra_label += self.radius * -4
        
        if 0 < pa < 90:
            rot_label += 180
            if pa > 20:
                ra_label += -2 * self.radius
        
        if 90 < pa < 180:
            ra_label += self.radius * -4
            rot_label += 180
            if pa > 160:
                ra_label += self.radius * -25
            elif 130 < pa < 160:
                ra_label += self.radius * -3
            elif 115 < pa < 130:
                dec_label += self.radius * 2
                ra_label += self.radius * 1
        
        if 180 < pa < 200:
            ra_label += self.radius * -4
        
        if pa > 350:
            ra_label += self.radius * -14
        
        return ra_label, dec_label, rot_label
    
    def make_chart(self,
                   image_path: str,
                   ra: List[float],
                   dec: List[float],
                   stars_file: Optional[str] = None,
                   output_name: str = 'fchart.png',
                   size_arcmin: Optional[float] = None,
                   extnum: Optional[int] = None,
                   show_slit: Optional[Tuple[float, float, float]] = None,
                   show_default_slit: bool = False) -> None:
        """
        Create a finding chart from a FITS image.
        
        Parameters
        ----------
        image_path : str
            Path to the FITS image
        ra : list of float
            Right ascension(s) of target object(s) in degrees
        dec : list of float
            Declination(s) of target object(s) in degrees
        stars_file : str, optional
            Path to file containing reference star information
        output_name : str
            Name for the output PNG file
        size_arcmin : float, optional
            Size of the chart in arcminutes (width=height)
        extnum : int, optional
            Extension number for multi-extension FITS files
        show_slit : tuple, optional
            Slit parameters (width, length, PA) in (arcsec, arcsec, degrees)
        show_default_slit : bool
            Show default slit at PA of star and target
        """
        # Load FITS data
        data, header = self._load_fits_data(image_path, extnum)
        
        # Create figure
        hdu = fits.PrimaryHDU(data, header)
        fig = aplpy.FITSFigure(hdu, north=True)
        
        # Set size if specified
        if size_arcmin is not None:
            fig.recenter(ra[0], dec[0], radius=size_arcmin / 60.0 * 0.5)
        
        # Display image with appropriate scaling
        self._apply_image_scaling(fig, data)
        
        # Add chart elements
        self._add_chart_elements(fig, ra, dec)
        
        # Add stars if specified
        if stars_file is not None:
            self.add_stars_info(stars_file, fig, ra[0], dec[0], show_default_slit)
        
        # Add slit if specified
        if show_slit is not None:
            self.plot_slit(fig, ra[0], dec[0], show_slit)
        
        # Finalize and save
        self._finalize_chart(fig, ra[0], dec[0], output_name)
    
    def _load_fits_data(self, 
                        image_path: str, 
                        extnum: Optional[int]) -> Tuple[np.ndarray, fits.Header]:
        """Load FITS data from file."""
        if extnum is not None:
            return fits.getdata(image_path, extnum, header=True)
        else:
            try:
                return fits.getdata(image_path, 1, header=True)
            except IndexError as e:
                print(f"Extension 1 not found: {e}")
                return fits.getdata(image_path, header=True)
    
    def _apply_image_scaling(self, fig: aplpy.FITSFigure, data: np.ndarray) -> None:
        """Apply appropriate image scaling (ZScale or percentile)."""
        try:
            zscale = ZScaleInterval()
            z1, z2 = zscale.get_limits(data)
            fig.show_grayscale(vmin=z1, vmax=z2)
            print(f'ZScale limits: z1={z1:.2f}, z2={z2:.2f}')
        except Exception as e:
            print(f"ZScale failed: {e}")
            print("Using percentile grayscale scaling")
            fig.show_grayscale(pmin=10, pmax=99)
    
    def _add_chart_elements(self, 
                           fig: aplpy.FITSFigure, 
                           ra: List[float], 
                           dec: List[float]) -> None:
        """Add standard chart elements (circles, scalebar, grid)."""
        # Add target circles
        if len(ra) > 1:
            fig.show_circles(xw=ra, yw=dec, radius=self.radius, 
                           edgecolor='red', alpha=1, lw=3)
        else:
            fig.show_circles(xw=ra[0], yw=dec[0], radius=self.radius, 
                           edgecolor='red', alpha=0.8, lw=2)
        
        # Add scalebar
        try:
            fig.add_scalebar(ARCMIN_IN_DEG, '1 arcmin', color='black', 
                           font='serif', linewidth=2)
        except:  # Compatibility with aplpy 1.0
            fig.add_scalebar(ARCMIN_IN_DEG, '1 arcmin', color='black')
        
        # Add grid
        try:
            fig.show_grid()
        except AttributeError:
            fig.add_grid()
            fig.grid.show()
    
    def _finalize_chart(self, 
                       fig: aplpy.FITSFigure, 
                       ra: float, 
                       dec: float, 
                       output_name: str) -> None:
        """Add final touches to the chart and save."""
        # Set axis labels
        fig.axis_labels.set_xtext('Right Ascension')
        fig.axis_labels.set_ytext('Declination')
        
        # Create title
        c = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree))
        title = f'RA= {c.ra.to_string(precision=3, sep=":", unit=u.hour)} ; ' \
                f'DEC = {c.dec.to_string(precision=3, sep=":", unit=u.degree, alwayssign=True)}'
        plt.title(title)
        
        # Set theme and save
        fig.set_theme('publication')
        fig.save(output_name, transparent=False)
        print(f'{output_name} created')


def get_uhs_info(hdu: fits.HDUList, 
                 ra: float, 
                 dec: float, 
                 size: float) -> Tuple[np.ndarray, fits.Header]:
    """
    Extract appropriate layer from UHS multi-layer images.
    
    Parameters
    ----------
    hdu : fits.HDUList
        HDU list from UHS image
    ra : float
        Right ascension in degrees
    dec : float
        Declination in degrees
    size : float
        Size of cutout in degrees
    
    Returns
    -------
    image : np.ndarray
        Image data
    header : fits.Header
        Image header
    """
    pixelscale = 0.400  # arcsec/pixel
    size_pixels = size / pixelscale
    
    # Try each layer until we find one with sufficient data
    for i, layer_idx in enumerate([1, 2, 3, 4]):
        image = hdu[layer_idx].data
        header = hdu[layer_idx].header
        wcs_img = wcs.WCS(header, relax=True)
        
        # Calculate cutout boundaries
        xc, yc = wcs_img.wcs_world2pix(ra, dec, 1)
        x1 = int(xc - 0.5 * size_pixels)
        x2 = int(xc + 0.5 * size_pixels)
        y1 = int(yc - 0.5 * size_pixels)
        y2 = int(yc + 0.5 * size_pixels)
        
        # Extract cutout
        imgcut = image[y1:y2, x1:x2]
        
        # Check if cutout has sufficient size
        if imgcut.size > 1000:
            print(f"Using Layer {layer_idx}")
            
            # Apply rotations as needed for layers 2 and 4
            if layer_idx == 2:
                image = np.rot90(image)
            elif layer_idx == 4:
                image = np.rot90(image, k=3)
            
            return image, header
    
    # If all layers fail, return the last one tried
    return image, header


def parse_slit_params(s: str) -> Tuple[float, float, float]:
    """
    Parse slit parameters from string.
    
    Parameters
    ----------
    s : str
        Comma-separated string "width,length,PA"
    
    Returns
    -------
    tuple
        (width, length, PA) as floats
    """
    try:
        width, length, pa = map(float, s.split(','))
        return width, length, pa
    except:
        raise argparse.ArgumentTypeError("Slit arguments must be width,length,PA")


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Create astronomical finding charts from FITS images. "
                    "The tool highlights objects of interest and reference stars, "
                    "with support for slit overlays and customizable appearance.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=EXAMPLES
    )
    
    parser.add_argument(
        'textfile',
        type=str,
        help='Text file with columns: ra, dec, image_path, [extnum]. '
             'The ra and dec should be in degrees.'
    )
    
    parser.add_argument(
        '--stars_info',
        type=str,
        default=None,
        help='Text file with reference stars to plot. '
             'Required columns: ra, dec, label'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='fchart.png',
        help='Output filename for the finding chart'
    )
    
    parser.add_argument(
        '--label_position',
        type=str,
        default='left',
        choices=["left", "right", "top", "bottom"],
        help='Position of star labels relative to the stars'
    )
    
    parser.add_argument(
        '--suffix',
        type=str,
        default="RA",
        help='Suffix for output image names (deprecated - use -o instead)'
    )
    
    parser.add_argument(
        '--size_stamp',
        type=float,
        default=None,
        help='Size of the stamp in arcminutes (width=height). '
             'If None, the whole image is used.'
    )
    
    parser.add_argument(
        '--show_slit',
        type=parse_slit_params,
        default=None,
        help='Slit parameters: width,length,PA (arcsec,arcsec,degrees). '
             'Example: 1.3,120,330.23'
    )
    
    parser.add_argument(
        '--show_default_slit',
        type=bool,
        default=False,
        help='Show default slit (3"x120") at PA between star and target'
    )
    
    return parser.parse_args()


def main():
    """Main execution function."""
    start_time = time.time()
    
    # Parse arguments
    args = parse_arguments()
    
    # Read input data
    data = ascii.read(args.textfile)
    print(f"Processing image: {data['img'][0]}")
    
    # Check for extension number
    extnum = None
    if 'extnum' in data.colnames:
        try:
            extnum = int(data['extnum'][0])
        except ValueError:
            print("Warning: Extension number must be an integer, ignoring.")
    
    # Create finding chart generator
    generator = FindingChartGenerator(label_position=args.label_position)
    
    # Generate the chart
    generator.make_chart(
        image_path=data['img'][0],
        ra=data['ra'],
        dec=data['dec'],
        stars_file=args.stars_info,
        output_name=args.output,
        size_arcmin=args.size_stamp,
        extnum=extnum,
        show_slit=args.show_slit,
        show_default_slit=args.show_default_slit
    )
    
    # Print elapsed time
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time:.2f} sec")


if __name__ == '__main__':
    main()