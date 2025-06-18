# Finding Chart Generator for Astronomical Observations


A Python tool for creating astronomical finding charts from FITS images, with support for Pan-STARRS (originally designed for high-z quasar panstarrs follow up; Banados+2014,2016,2016,2018,2023) and other survey data. This tool is designed to help astronomers prepare for observations by generating clear, annotated charts that highlight target objects and reference stars.

## Features

- **Multi-format Support**: Works with Pan-STARRS, UHS, and other FITS format images
- **Reference Star Annotation**: Automatically labels and highlights reference stars with position angles and separations
- **Slit Overlay**: Display spectrograph slit positions with customizable width, length, and position angle
- **Smart Labeling**: Intelligent label positioning to avoid overlaps and maintain clarity
- **Customizable Output**: Adjustable chart size, label positions, and visual styling
- **Multi-Extension FITS**: Support for multi-extension FITS files
- **Batch Processing**: Process multiple targets from a single input file

### Prerequisites

- Python 3.6 or higher
- pip package manager

### Required Dependencies

```bash
pip install numpy matplotlib astropy aplpy
```


## Usage

### Basic Usage

Create a simple finding chart:

```bash
./finding_chart.py targets.txt -o my_finding_chart.png
```

### Input File Format

The input text file should contain columns separated by spaces or tabs:

```
# ra(deg)    dec(deg)    image_path    [extnum]
245.89742    -17.92835   image1.fits   1
246.12456    -18.03421   image2.fits
```

- `ra`: Right ascension in degrees
- `dec`: Declination in degrees  
- `image_path`: Path to the FITS image
- `extnum`: (Optional) Extension number for multi-extension FITS files

### Reference Stars File Format

For highlighting reference stars, create a file with:

```
# ra(deg)    dec(deg)    label
245.89142    -17.92235   S1
245.90342    -17.93435   S2
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o, --output` | Output filename for the finding chart | `fchart.png` |
| `--stars_info` | Path to reference stars file | None |
| `--size_stamp` | Size of the chart in arcminutes (square) | Full image |
| `--show_slit` | Slit parameters: width,length,PA (arcsec,arcsec,deg) | None |
| `--show_default_slit` | Show default 3"×120" slit at star-target PA | False |
| `--label_position` | Position of star labels (left/right/top/bottom) | left |

## Examples

### Example 1: Basic Finding Chart

```bash
./finding_chart.py targets.txt -o NGC1234_finder.png
```

### Example 2: Chart with Reference Stars

```bash
./finding_chart.py qso_target.txt --stars_info reference_stars.txt -o qso_finder.png
```

### Example 3: Chart with Slit Overlay

Show a 1.5 arcsec wide, 120 arcsec long slit at PA=45°:

```bash
./finding_chart.py target.txt --show_slit 1.5,120,45 -o target_with_slit.png
```

### Example 4: Customized Chart Size

Create a 5 arcminute finding chart:

```bash
./finding_chart.py target.txt --size_stamp 5 -o small_finder.png
```

## Output

The tool generates PNG finding charts with:

- Target object(s) marked with red circles
- Reference stars marked with blue rectangles and labels
- 1 arcminute scale bar
- Coordinate grid
- North orientation indicator
- RA/Dec axis labels
- Title with precise coordinates

### Sample Output Features

Each reference star annotation includes:
- Star label
- RA/Dec coordinates
- Offset from target (E/W, N/S)
- Position angle (PA) in degrees
- Angular separation in arcseconds

## Advanced Features

### Multi-Layer Image Support

The tool automatically handles multi-layer images (e.g., UHS survey data) by selecting the appropriate layer containing the target.

### Intelligent Label Positioning

Labels are automatically positioned to avoid overlapping with chart features based on:
- Position angle between star and target
- Proximity to other objects
- Chart boundaries

### Custom Slit Positioning

When using `--show_default_slit`, the tool calculates the optimal slit position based on the reference star positions.

## Troubleshooting

### Common Issues

1. **WCS Headers**: If you encounter WCS-related errors, ensure your FITS files have proper WCS information in the headers.

2. **Missing Extensions**: For multi-extension FITS files, specify the correct extension number in the input file.

3. **Memory Issues**: For very large images, consider using the `--size_stamp` option to create smaller cutouts.

### Error Messages

- `Extension X not found`: The specified extension doesn't exist in the FITS file
- `ZScale failed`: The image scaling algorithm failed; the tool will automatically fall back to percentile scaling

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

### Development Setup

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



## Acknowledgments

- Built with [Astropy](https://www.astropy.org/) and [APLpy](https://aplpy.github.io/)

---

**Note**: This tool is under active development. Features and interfaces may change in future versions.
