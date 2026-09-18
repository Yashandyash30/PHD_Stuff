# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python (Astro Photometry)
#     language: python
#     name: astro_photometry
# ---

# %%
"""
daophot_pipeline_standalone.py
================================
Automated DAOPHOT II photometry pipeline using Stetson's original
standalone FORTRAN binaries (daophot + allstar) via subprocess.

Drop-in replacement for daophot_pipeline.py — all Python logic
(FWHM, WCS, Vizier, forced photometry, plots) is identical.
Only the D1–D7 reduction stages are changed.

Assumes you already have:
  • daophot  binary on PATH (or set DAOPHOT_BIN below)
  • allstar  binary on PATH (or set ALLSTAR_BIN below)
  • DAOPHOT.OPT in your working directory  (RE, GA, FW, HI, FI, PS set)
  • PHOTO.OPT in your working directory    (aperture radii set)
  • ALLSTAR.OPT in your working directory  (optional but recommended)

IMPORTANT — set WA = -2 in your DAOPHOT.OPT for non-interactive PSF
construction.  WA=-2 means DAOPHOT uses every PSF star without prompting
and writes minimal terminal output — essential for subprocess piping.

Pipeline stages
---------------
  D0  : FWHM, sky mean, sky std  (Gaussian fits to bright stars)
  D1  : FIND        → stem.coo
  D2  : PHOT        → stem.ap
  D3  : PICK        → stem.lst
  D4  : PSF         → stem.psf  (iterative neighbour-subtract loop)
  D5  : ALLSTAR     → stem.als  (full frame)
  D6  : Forced phot → stem_merged.ap → stem_merged.als  (when target missed)
  D7  : Source lookup  (RA/Dec → pixel → .als magnitude)
  D8  : Zero-point calibration  (Vizier / Pan-STARRS)
  D9  : Final calibrated magnitude
"""

# %%
# ============================================================
# CELL 1 : CONFIGURATION
# ============================================================

# --- Target ---
TARGET_NAME   = "NGC253"
TARGET_RA     =  None               # degrees (set to None for all-star photometry)
TARGET_DEC    =  None               # degrees (set to None for all-star photometry)
GRB_T0_JD     =  2458951.881027           # Burst trigger time in Julian Days (adjust as needed)

# --- Input image (WCS-solved, stacked FITS) ---
STACKED_IMAGE = "NGC253_R.fits"

# --- Astrometry / WCS Fallback ---
# If your primary STACKED_IMAGE lacks a WCS, the pipeline will borrow the WCS from this image.
# Set to None if your primary image already has WCS.
ASTROMETRY_REF_IMAGE = None

# How to determine FRAMES_AVG and FRAMES_SUM: 'auto_avg', 'auto_sum', or 'manual'
FRAMES_COUNT_MODE = 'auto_avg'
STACKED_IMAGE_WILDCARD = "al_*" # Wildcard pattern for auto counting frames (used if mode is auto)

# Pipeline mode :
PIPELINE_MODE    = 'semi-auto' # 'auto' = full auto, 'semi-auto' = run AT, FI, PH, PI automatically, 'manual' = full manual (skip all, assume manual run)

# --- Zero-point calibration ---
CALIBRATION_SYSTEM  = 'Pan-STARRS'    # Options: 'Pan-STARRS' or 'DIRECT_VEGA'
DIRECT_VEGA_CATALOG = 'APASS'          # Options: 'APASS', 'Landolt', or 'Gaia'
FILTER_BAND         = 'R'              # Target observing filter



# ============================================================

# Manual values (used if FRAMES_COUNT_MODE == 'manual')
FRAMES_AVG    =  3       # N frames that were averaged together
FRAMES_SUM    =  1       # N frames that were summed together
               # entry format sent to FIND: "FRAMES_AVG,FRAMES_SUM"
               # examples: single exp → "1,1"; avg of 5 → "5,1"; sum of 3 → "1,3"

# --- Detector ---
READ_NOISE    =  10.0     # e⁻ / ADU (one raw frame)
GAIN          =  5.0     # photons / ADU (one raw frame)
SATURATION    = 60000    # ADU — detector linearity limit

# --- Detection --- (Only required in `auto` mode)
DETECT_SIGMA  =  3.0     # σ threshold for FIND

# --- PSF --- (Only required in `auto` mode)
MAX_PSF_STARS    = 25     # PICK target number of candidates
PSF_FAINT_MAG    = 99     # faintest magnitude for PICK (99 = no limit)
PSF_ITERATIONS   =  3     # PSF rebuild iterations (≥2 recommended for crowded fields)
VARIABLE_PSF     =  0     # 0=constant, 1=linear, 2=quadratic, -1=pure analytic

# --- OPT File Parameters (Dynamic scaling based on FWHM) ---
DAOPHOT_FI_MULT = 1.0
DAOPHOT_PS_MULT = 3.0

ALLSTAR_FI_MULT = 1.0
ALLSTAR_IS_MULT = 6.0
ALLSTAR_OS_CONST = 10.0

PHOTO_A1_MULT = 1.0
PHOTO_APERTURE_MULT = 1.177
PHOTO_AC_MULT = 6.0
PHOTO_OS_CONST = 10.0

# --- Forced photometry ---
FORCE_COORDS   = None      # (ra_deg, dec_deg) or None
FORCE_CENTROID = True      # snap to local brightness peak before injecting

# --- Matching ---
MATCH_RADIUS_PIX    = 8.0
MATCH_RADIUS_ARCSEC = 2.0

# --- Stellarity / Quality Cuts (for all-star catalog) ---
FILTER_EXTENDED = True     # If True, filter out non-stellar objects (galaxies, cosmic rays)
MAX_CHI         = 3.0      # Reject objects with CHI > this value (bad PSF fit)
MAX_SHARP_ABS   = 1.0      # Reject objects with |SHARP| > this value (extended or cosmic ray)

# --- DYNAMIC QUALITY CONSTRAINTS ---
# Group 1: Dense, Deep, High-Precision Catalogs
if CALIBRATION_SYSTEM == 'Pan-STARRS' or DIRECT_VEGA_CATALOG == 'Gaia':
    QUERY_RADIUS_ARCMIN = 12               
    REF_MAG_BRIGHT      = 13.0
    REF_MAG_FAINT       = 19.5
    REF_MAG_ERR_MAX     = 0.08
    ZP_MATCH_RADIUS_AS  = 3.5             
# Group 2: Sparse, Shallow, or Classical Catalogs
elif CALIBRATION_SYSTEM == 'DIRECT_VEGA' and DIRECT_VEGA_CATALOG in ['APASS', 'Landolt']:
    QUERY_RADIUS_ARCMIN = 15                
    REF_MAG_BRIGHT      = 12.0          
    REF_MAG_FAINT       = 17.5
    REF_MAG_ERR_MAX     = 0.10
    ZP_MATCH_RADIUS_AS  = 5.0             

# Global settings for all catalogs
VIZIER_ROW_LIMIT    = 5000             
ZP_MIN_STARS        = 3

# --- ZP Calibration Star Quality (tighter than catalog-wide cuts) ---
ZP_MAX_CHI          = 2.0      # Only use stars with CHI < this for ZP (good PSF fit)
ZP_MAX_SHARP_ABS    = 0.5      # Only use stars with |SHARP| < this for ZP (point sources)
ZP_MAX_MERR         = 0.10     # Max instrumental mag error for ZP calibration stars
GALAXY_EXCL_RADIUS_ARCSEC = 0  # Set > 0 to exclude stars within this radius of image centre
EDGE_EXCL_MARGIN_PIXELS   = 50 # Set > 0 to exclude stars within this many pixels of image edges
ZP_MAX_SCATTER_MAG        = 0.20 # Reject stars whose ZP differs from the median ZP by more than this

# ============================================================



import glob
if FRAMES_COUNT_MODE == 'auto_avg':
    FRAMES_AVG = max(1, len(glob.glob(STACKED_IMAGE_WILDCARD)))
    FRAMES_SUM = 1
elif FRAMES_COUNT_MODE == 'auto_sum':
    FRAMES_AVG = 1
    FRAMES_SUM = max(1, len(glob.glob(STACKED_IMAGE_WILDCARD)))


# --- Directories ---
WORK_DIR  = f"daophot_work_{FILTER_BAND}_{CALIBRATION_SYSTEM}"
PLOTS_DIR = f"daophot_output_{FILTER_BAND}_{CALIBRATION_SYSTEM}"

# --- Subprocess timeouts (seconds) ---
TIMEOUT_DAOPHOT  = 600
TIMEOUT_ALLSTAR  = 1200

# --- Standalone binary paths ---
# If the binaries are on your PATH, just use "daophot" / "allstar".
# Otherwise supply full paths, e.g. "/usr/local/daophot/daophot"
DAOPHOT_BIN   = "distrobox enter astro-box -- bash -c '~/dao2/daophot'"
ALLSTAR_BIN   = "distrobox enter astro-box -- bash -c '~/dao2/allstar'"



# %%
# ============================================================
# CELL 2 : IMPORTS
# ============================================================

import os
import sys
import glob
import shutil
import warnings
import subprocess
import textwrap
import numpy as np
import matplotlib
try:
    ipy = get_ipython()
    try:
        ipy.run_line_magic('matplotlib', 'notebook')
    except Exception:
        ipy.run_line_magic('matplotlib', 'inline')
except NameError:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
from astropy.visualization import ZScaleInterval

from astropy.io    import fits
from astropy.wcs   import WCS
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u
from astroquery.vizier import Vizier
import warnings
from astropy.io.fits.verify import VerifyWarning
warnings.filterwarnings('ignore', category=VerifyWarning)
# %%
# CELL 3 : Functions
# ============================================================
# ============================================================

def get_catalog_columns(calib_system, vega_catalog, filter_band):
    if calib_system == 'Pan-STARRS':
        return (
            'II/349/ps1',
            ['RAJ2000', 'DEJ2000', 'gmag', 'e_gmag',
             'rmag', 'e_rmag', 'imag', 'e_imag'],
            None, None
        )
    elif calib_system == 'DIRECT_VEGA':
        catalog_map = {
            'APASS'  : ('II/336/apass9',  ['RAJ2000', 'DEJ2000'],
                        {'V': ('Vmag','e_Vmag'), 'B': ('Bmag','e_Bmag'),
                         'g': ("g'mag","e_g'mag", "g_mag","e_g_mag"), 
                         'r': ("r'mag","e_r'mag", "r_mag","e_r_mag"),
                         'i': ("i'mag","e_i'mag", "i_mag","e_i_mag")}),
            'Landolt': ('II/183A/table2', ['RAJ2000', 'DEJ2000'],
                        {'V': ('Vmag','e_Vmag')}),
            'Gaia'   : ('I/355/gaiadr3',  ['RA_ICRS', 'DE_ICRS'],
                        {'G':  ('Gmag','e_Gmag'),
                         'BP': ('BPmag','e_BPmag'),
                         'RP': ('RPmag','e_RPmag')}),
        }
        if vega_catalog not in catalog_map:
            raise ValueError(f"Unknown DIRECT_VEGA_CATALOG: {vega_catalog}")
        cat_id, coord_cols, valid_cols = catalog_map[vega_catalog]
        if filter_band not in valid_cols:
            raise ValueError(
                f"[CRITICAL] {vega_catalog} has no column for "
                f"filter '{filter_band}'. Available: {list(valid_cols)}"
            )
        mag_info = valid_cols[filter_band]
        if len(mag_info) == 2:
            q_mag, q_err = mag_info
            r_mag, r_err = mag_info
        else:
            q_mag, q_err, r_mag, r_err = mag_info
        return cat_id, coord_cols + [q_mag, q_err], r_mag, r_err
    raise ValueError(f"Unknown CALIBRATION_SYSTEM: {calib_system}")


VIZIER_CATALOG_ID, VIZIER_COLUMNS, VIZIER_MAG_COL, VIZIER_ERR_COL = (
    get_catalog_columns(CALIBRATION_SYSTEM, DIRECT_VEGA_CATALOG, FILTER_BAND)
)


# ============================================================
# CELL 4 : FILE-MANAGEMENT HELPERS
# ============================================================

def setup_dirs():
    for d in [WORK_DIR, PLOTS_DIR]:
        Path(d).mkdir(parents=True, exist_ok=True)
    for sub in ['zero_point', 'photometry', 'fwhm']:
        Path(os.path.join(PLOTS_DIR, sub)).mkdir(parents=True, exist_ok=True)
    print(f"[Setup] Work dir  : {os.path.abspath(WORK_DIR)}")
    print(f"[Setup] Plots dir : {os.path.abspath(PLOTS_DIR)}")


def workpath(filename):
    return os.path.join(os.path.abspath(WORK_DIR), filename)


def copy_image_to_workdir(image_path):
    """Copy input FITS image into WORK_DIR so all DAOPHOT outputs land there."""
    dest = workpath(os.path.basename(image_path))
    if not os.path.exists(dest):
        shutil.copy2(image_path, dest)
        print(f"[Setup] Copied image → {dest}")
    else:
        print(f"[Setup] Image already present: {dest}")
    return dest


def copy_opt_files_to_workdir(opt_source_dir='.'):
    """
    Copy .OPT files from opt_source_dir into WORK_DIR so DAOPHOT finds them.
    DAOPHOT looks for DAOPHOT.OPT, PHOTO.OPT in its current working directory.
    Call this if your .OPT files live somewhere other than WORK_DIR.
    """
    for opt in ['daophot.opt', 'photo.opt', 'allstar.opt',
                'DAOPHOT.OPT', 'PHOTO.OPT', 'ALLSTAR.OPT']:
        src = os.path.join(opt_source_dir, opt)
        if os.path.exists(src):
            dst = workpath(opt.lower())   # lowercase on Unix
            shutil.copy2(src, dst)
            print(f"[Setup] Copied {opt} → {dst}")


def get_stem(image_path):
    """Return bare filename without extension: stem for DAOPHOT default files."""
    return os.path.splitext(os.path.basename(image_path))[0]


def clean_previous_run(stem, pipeline_mode='auto'):
    """
    Remove output files from a previous run to prevent DAOPHOT's
    "file already exists — overwrite?" interactive prompt from firing.
    """
    patterns = [
        f"{stem}.als",
        f"{stem}s.fits", f"{stem}s",
        f"{stem}_nei*.als",  f"{stem}_nei*.fits",
        f"{stem}_clean*.fits",
        f"{stem}_merged.ap", f"{stem}_merged.als",
        f"{stem}_forced.coo", f"{stem}_forced.ap",
        f"{stem}_forced_off.ap",
    ]
    if pipeline_mode == 'auto':
        patterns.extend([
            f"{stem}.coo", f"{stem}.ap", f"{stem}.lst",
            f"{stem}.nei", f"{stem}.psf"
        ])
    for pat in patterns:
        for f in glob.glob(workpath(pat)):
            os.remove(f)
            print(f"[Clean] Removed: {os.path.basename(f)}")


# ============================================================
# CELL 5 : NATIVE FILE PARSERS
#           Standalone DAOPHOT writes plain ASCII — no IRAF overhead
# ============================================================

def parse_coo_header(coo_path):
    """
    Read the two-line header from a standalone DAOPHOT .coo file.
    Returns a dict: nx, ny, lowbad, highbad, thresh
    """
    defaults = dict(nx=0, ny=0, lowbad=0.0, highbad=50000.0, thresh=0.0)
    if not os.path.exists(coo_path):
        return defaults
    with open(coo_path) as fh:
        lines = [ln for ln in fh if ln.strip()]
    # Line 0: keyword names — skip
    # Line 1: numeric values in order: NL NX NY LOWBAD HIGHBAD THRESH ...
    if len(lines) < 2:
        return defaults
    try:
        parts = lines[1].split()
        return dict(
            nl      = int(float(parts[0])),
            nx      = int(float(parts[1])),
            ny      = int(float(parts[2])),
            lowbad  = float(parts[3]),
            highbad = float(parts[4]),
            thresh  = float(parts[5]) if len(parts) > 5 else 0.0,
        )
    except (ValueError, IndexError):
        return defaults


def parse_als_file(als_path):
    """
    Parse a standalone ALLSTAR .als file.

    Header (2 lines):
      NL NX NY LOWBAD HIGHBAD THRESH AP1 PH/ADU RNOISE FRAD
      1  284  492  400.0  24000.0  20.0  3.00  20.00  6.50  2.0

    Data (one line per star):
      ID  XCENTER  YCENTER  MAG  MERR  MSKY  NITER  CHI  SHARP
    """
    if not os.path.exists(als_path):
        print(f"[Parser] ❌ Not found: {als_path}")
        return None

    rows = []
    with open(als_path) as fh:
        header_lines = 0
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            # Header lines begin with non-digit characters (keyword or value line)
            if header_lines < 2 and not stripped[0].isdigit():
                header_lines += 1
                continue
            if header_lines < 2:
                header_lines += 1
                continue
            parts = stripped.split()
            if len(parts) < 9:
                continue
            try:
                rows.append({
                    'ID'     : int(parts[0]),
                    'XCENTER': float(parts[1]),
                    'YCENTER': float(parts[2]),
                    'MAG'    : float(parts[3]),
                    'MERR'   : float(parts[4]),
                    'MSKY'   : float(parts[5]),
                    'NITER'  : int(float(parts[6])),
                    'CHI'    : float(parts[7]),
                    'SHARP'  : float(parts[8]),
                })
            except (ValueError, IndexError):
                continue

    if not rows:
        print(f"[Parser] ⚠️  No valid rows in {os.path.basename(als_path)}")
        return None

    t = Table(rows=rows)
    print(f"[Parser] Read {len(t)} stars from {os.path.basename(als_path)}")
    return t


def create_forced_coo(forced_x, forced_y, coo_template_path, output_path,
                      id_offset=100000):
    """
    Write a minimal valid DAOPHOT .coo file containing one star
    at (forced_x, forced_y).
    """
    # Copy exact header from template to avoid DAOPHOT rejecting it
    with open(coo_template_path, 'r') as f:
        lines = [f.readline() for _ in range(3)]
        
    with open(output_path, 'w') as fh:
        fh.writelines(lines)
        # Use star_id=1 to fit in I5 format, OFFSET will add 100000 later
        star_id = 1
        fh.write(
            f" {star_id:5d} {forced_x:9.3f} {forced_y:9.3f} "
            f"  0.000   0.500   0.000   0.000\n"
        )
    print(f"[ForcedCoo] Written: {os.path.basename(output_path)}  "
          f"X={forced_x:.2f} Y={forced_y:.2f}")


# ============================================================
# CELL 6 : D0 — FWHM & SKY MEASUREMENT (Robust)
# ============================================================

def get_wcs(header, ref_image_path=None):
    """
    Get WCS from header. If missing/invalid, try to load from ref_image_path.
    """
    try:
        wcs = WCS(header, naxis=2)
        if wcs.wcs.has_cd() or wcs.wcs.has_pc():
            return wcs
    except Exception:
        pass
        
    if ref_image_path and os.path.exists(ref_image_path):
        print(f"  [WCS] Primary image lacks WCS. Falling back to {ref_image_path}")
        with fits.open(ref_image_path) as ref_hdul:
            return WCS(ref_hdul[0].header, naxis=2)
    return WCS(header, naxis=2)  # Will likely be invalid, but return it anyway


def measure_fwhm_and_sky(image_path):
    """
    Estimate FWHM (Gaussian fits) and sky statistics from the image.
    Robustly handles galaxy images by using shape-filtering and (if WCS available)
    Gaia cross-matching to ensure only true point sources are measured.
    """
    print("\n" + "="*60)
    print("  D0: FWHM & SKY MEASUREMENT")
    print("="*60)

    with fits.open(image_path) as hdul:
        data   = hdul[0].data.astype(float)
        header = hdul[0].header

    mean_sky, median_sky, std_sky = sigma_clipped_stats(
        data, sigma=3.0, maxiters=10
    )
    modal_sky = 3.0 * median_sky - 2.0 * mean_sky
    print(f"  Sky → mean={mean_sky:.2f}  median={median_sky:.2f}  "
          f"mode={modal_sky:.2f}  std={std_sky:.2f}")

    from astropy.modeling  import models, fitting
    from photutils.detection import find_peaks
    from astropy.nddata    import Cutout2D

    sat  = header.get('SATURATE', SATURATION)
    ny, nx = data.shape
    edge = int(min(nx, ny) * 0.05)

    # Find many peaks to build a robust statistical sample
    peaks = find_peaks(data - modal_sky, threshold=5.0 * std_sky, box_size=15)
    
    if peaks is None or len(peaks) == 0:
        raise RuntimeError("[D0] ❌ No sources found above 5 sigma.")

    valid = (
        (peaks['peak_value'] < sat - modal_sky) &
        (peaks['x_peak'] > edge) & (peaks['x_peak'] < nx - edge) &
        (peaks['y_peak'] > edge) & (peaks['y_peak'] < ny - edge)
    )
    peaks = peaks[valid]
    peaks.sort('peak_value')
    peaks.reverse()

    fitter = fitting.LevMarLSQFitter()
    box_sz = 21
    cy, cx = box_sz / 2.0, box_sz / 2.0
    yy, xx = np.mgrid[:box_sz, :box_sz]

    measured_sources = []
    
    # Process up to 500 bright peaks
    for star in peaks[:500]:
        xc, yc = int(star['x_peak']), int(star['y_peak'])
        try:
            cut = Cutout2D(data, (xc, yc), box_sz).data
        except Exception:
            continue
            
        g_init = (
            models.Const2D(amplitude=modal_sky) +
            models.Gaussian2D(
                amplitude=float(star['peak_value']),
                x_mean=cx, y_mean=cy,
                x_stddev=2.0, y_stddev=2.0
            )
        )
        g_init.x_stddev_1.bounds = (0.5, 7.0)
        g_init.y_stddev_1.bounds = (0.5, 7.0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            g_fit = fitter(g_init, xx, yy, cut)
            
        sigma_x = abs(g_fit.x_stddev_1.value)
        sigma_y = abs(g_fit.y_stddev_1.value)
        sigma = (sigma_x + sigma_y) / 2.0
        fwhm  = 2.355 * sigma
        
        # Shape filtering: Must be roughly round (reject galaxies) and right size
        ratio = sigma_x / sigma_y if sigma_y > 0 else 0
        if 1.5 < fwhm < 15.0 and 0.8 < ratio < 1.25:
            measured_sources.append({
                'x': xc, 'y': yc, 'fwhm': fwhm, 'flux': float(star['peak_value'])
            })

    if len(measured_sources) < 3:
        raise RuntimeError(
            "[D0] ❌ Not enough clean, round stars to measure FWHM. "
            "Check your image or SATURATION limit."
        )

    # Convert to table for cross-matching
    src_table = Table(rows=measured_sources)
    wcs = get_wcs(header, ASTROMETRY_REF_IMAGE)
    
    final_fwhms = []
    cross_match_successful = False

    # Attempt Gaia Cross-match if WCS is valid
    if wcs.wcs.has_cd() or wcs.wcs.has_pc():
        try:
            ra_cen, dec_cen = wcs.pixel_to_world_values(nx/2, ny/2)
            coord_cen = SkyCoord(ra_cen, dec_cen, unit='deg')
            radius = np.max([nx, ny]) / 2.0 * proj_plane_pixel_scales(wcs)[0] * 3600 * u.arcsec
            
            print(f"  [D0] WCS found. Querying Gaia DR3 for point sources within {radius.value:.1f}\"")
            v = Vizier(columns=['RA_ICRS', 'DE_ICRS', 'Gmag'], row_limit=5000)
            res = v.query_region(coord_cen, radius=radius, catalog='I/355/gaiadr3')
            
            if res and len(res) > 0:
                gaia_cat = res[0]
                gaia_coords = SkyCoord(gaia_cat['RA_ICRS'], gaia_cat['DE_ICRS'], unit='deg')
                src_coords = SkyCoord(*wcs.pixel_to_world_values(src_table['x'], src_table['y']), unit='deg')
                
                idx, d2d, _ = src_coords.match_to_catalog_sky(gaia_coords)
                match_limit = 2.0 * u.arcsec
                matched = d2d < match_limit
                
                matched_sources = src_table[matched]
                if len(matched_sources) >= 3:
                    final_fwhms = matched_sources['fwhm']
                    cross_match_successful = True
                    print(f"  [D0] Found {len(matched_sources)} Gaia point sources. Using these for FWHM.")
                else:
                    print("  [D0] Not enough Gaia matches. Falling back to stellar locus.")
        except Exception as e:
            print(f"  [D0] Gaia query failed: {e}. Falling back to stellar locus.")
    else:
        print("  [D0] No valid WCS found. Falling back to stellar locus.")

    if not cross_match_successful:
        # Fall back to the mode (densest region) of the FWHM distribution
        # Since we already shape-filtered, the median of the smallest objects is usually safe
        # but to be robust against a few small artifacts, we use sigma_clipped_stats
        _, fwhm_median, _ = sigma_clipped_stats(src_table['fwhm'], sigma=2.0)
        final_fwhms = src_table['fwhm']
        print(f"  [D0] Using shape-filtered locus ({len(src_table)} objects).")

    _, fwhm_median, _ = sigma_clipped_stats(final_fwhms, sigma=2.5)
    print(f"  FWHM → {fwhm_median:.2f} px")

    # --- SAVE FWHM DATA & PLOT ---
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches
        from astropy.visualization import ZScaleInterval, ImageNormalize
        fwhm_dir = os.path.join(PLOTS_DIR, 'fwhm')
        
        if cross_match_successful:
            final_table = matched_sources
        else:
            final_table = src_table
            
        txt_path = os.path.join(fwhm_dir, f"{TARGET_NAME}_fwhm_stars.txt")
        with open(txt_path, 'w') as f:
            f.write("X, Y, FWHM, FLUX\n")
            for row in final_table:
                f.write(f"{row['x']:.2f}, {row['y']:.2f}, {row['fwhm']:.3f}, {row['flux']:.2f}\n")
                
        img_path = os.path.join(fwhm_dir, f"{TARGET_NAME}_fwhm_stars.png")
        fig, ax = plt.subplots(figsize=(10, 10))
        norm = ImageNormalize(data, interval=ZScaleInterval())
        ax.imshow(data, cmap='gray', origin='lower', norm=norm)
        for row in final_table:
            ax.add_patch(patches.Circle((row['x'], row['y']), radius=15, edgecolor='red', facecolor='none', lw=2))
            ax.text(row['x']+18, row['y']+18, f"{row['fwhm']:.2f}", color='red', fontsize=8)
        ax.set_title(f"Stars used for FWHM Calculation (Median FWHM = {fwhm_median:.2f})")
        fig.savefig(img_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  [D0] Saved FWHM plot to {img_path}")
        print(f"  [D0] Saved FWHM data to {txt_path}")
    except Exception as e:
        print(f"  [D0] Warning: Could not save FWHM plot/data: {e}")

    return fwhm_median, modal_sky, std_sky, header, data


# ============================================================
# CELL 7 : SUBPROCESS WRAPPERS  ← replaces every iraf.* call
# ============================================================

def _run_binary(binary, commands, cwd, timeout, label):
    """
    Pipe a list of command strings to a DAOPHOT / ALLSTAR binary.

    Each item in `commands` is one input line.
    Empty string ''  →  bare carriage return  (accept default).

    Returns (stdout, stderr, returncode).
    Raises RuntimeError on timeout or if the binary is not found.
    """
    stdin_text = "\n".join(str(c) for c in commands) + "\n"

    try:
        result = subprocess.run(
            binary,                 # <-- Removed brackets [ ] so it passes as a single string
            shell          = True,  # <-- Added to evaluate the complex Distrobox command
            input          = stdin_text,
            capture_output = True,
            text           = True,
            cwd            = cwd,
            timeout        = timeout,
        )
    except FileNotFoundError:
        raise RuntimeError(
            f"[{label}] ❌ Binary '{binary}' not found. "
            "Install DAOPHOT II and add it to PATH, "
            f"or set {label}_BIN to the full path."
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"[{label}] ❌ Timed out after {timeout}s. "
            "Increase TIMEOUT_DAOPHOT / TIMEOUT_ALLSTAR if the frame is large."
        )

    # Surface any fatal DAOPHOT error messages
    if result.returncode != 0:
        print(f"[{label}] Non-zero return code: {result.returncode}")
        print(f"[{label}] STDERR:\n{result.stderr[-500:]}")

    import datetime
    log_path = os.path.join(cwd, "daophot_run.log")
    try:
        with open(log_path, "a") as f:
            f.write(f"\n{'='*60}\n[{label}] RUN at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n{'='*60}\n")
            f.write("--- INPUT ---\n")
            f.write(stdin_text)
            f.write("--- OUTPUT ---\n")
            f.write(result.stdout)
            if result.stderr:
                f.write("--- STDERR ---\n")
                f.write(result.stderr)
            f.write("\n")
    except Exception as e:
        print(f"Failed to write log: {e}")

    return result.stdout, result.stderr, result.returncode


def run_daophot_session(commands, label="DAOPHOT"):
    """
    Run one complete daophot session in WORK_DIR.
    commands : list of strings, one per input line.
    """
    stdout, stderr, rc = _run_binary(
        DAOPHOT_BIN, commands,
        cwd     = os.path.abspath(WORK_DIR),
        timeout = TIMEOUT_DAOPHOT,
        label   = label,
    )
    _print_session_summary(stdout, label)
    return stdout


def run_allstar_session(commands, label="ALLSTAR"):
    """Run one complete allstar session in WORK_DIR."""
    stdout, stderr, rc = _run_binary(
        ALLSTAR_BIN, commands,
        cwd     = os.path.abspath(WORK_DIR),
        timeout = TIMEOUT_ALLSTAR,
        label   = label,
    )
    _print_session_summary(stdout, label)
    return stdout


def _print_session_summary(stdout, label):
    """Print a condensed view of DAOPHOT/ALLSTAR terminal output."""
    if not stdout:
        return
    lines = stdout.splitlines()
    # Show first 5 and last 15 lines — catches startup params and final stats
    head = lines[:5]
    tail = lines[-15:] if len(lines) > 20 else lines[5:]
    print(f"\n  [{label} output — head]")
    for ln in head:
        print(f"    {ln}")
    if len(lines) > 20:
        print(f"    ... ({len(lines) - 20} lines) ...")
    print(f"  [{label} output — tail]")
    for ln in tail:
        print(f"    {ln}")


# ============================================================
# CELL 8 : OPTIONAL DAOPHOT.OPT UPDATE
#           Updates FW (FWHM) and RE/GA if measured dynamically.
#           Leaves every other parameter in your existing file intact.
# ============================================================

def update_daophot_opt(fwhm, std_sky=None):
    """
    Patch FW, FI, and PS in the copy of DAOPHOT.OPT that lives in WORK_DIR.
    """
    opt_path = workpath('daophot.opt')
    fi_val = fwhm * DAOPHOT_FI_MULT
    ps_val = fi_val * DAOPHOT_PS_MULT

    if not os.path.exists(opt_path):
        print("[OPT] ⚠️  daophot.opt not found in WORK_DIR — creating minimal version.")
        with open(opt_path, 'w') as fh:
            fh.write(f"RE = {READ_NOISE:.2f}\n")
            fh.write(f"GA = {GAIN:.2f}\n")
            fh.write(f"FW = {fwhm:.3f}\n")
            fh.write(f"HI = {float(SATURATION):.1f}\n")
            fh.write(f"WA = -2\n")
            fh.write(f"FI = {fi_val:.3f}\n")
            fh.write(f"PS = {ps_val:.3f}\n")
        print(f"[OPT] Created minimal daophot.opt with FW={fwhm:.3f}, FI={fi_val:.3f}, PS={ps_val:.3f}")
        return

    with open(opt_path) as fh:
        lines = fh.readlines()

    updated = {}
    new_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped or '=' not in stripped:
            new_lines.append(line)
            continue
        key = stripped.split('=')[0].strip().upper()[:2]
        if key == 'FW':
            new_lines.append(f"FW = {fwhm:.3f}\n")
            updated['FW'] = True
        elif key == 'FI':
            new_lines.append(f"FI = {fi_val:.3f}\n")
            updated['FI'] = True
        elif key == 'PS':
            new_lines.append(f"PS = {ps_val:.3f}\n")
            updated['PS'] = True
        else:
            new_lines.append(line)

    if 'FW' not in updated: new_lines.append(f"FW = {fwhm:.3f}\n")
    if 'FI' not in updated: new_lines.append(f"FI = {fi_val:.3f}\n")
    if 'PS' not in updated: new_lines.append(f"PS = {ps_val:.3f}\n")

    with open(opt_path, 'w') as fh:
        fh.writelines(new_lines)

    print(f"[OPT] Updated daophot.opt — FW={fwhm:.3f}, FI={fi_val:.3f}, PS={ps_val:.3f}")


def update_allstar_opt(fwhm):
    """Patch FI, IS, OS in ALLSTAR.OPT"""
    opt_path = workpath('allstar.opt')
    fi_val = fwhm * ALLSTAR_FI_MULT
    is_val = min(fi_val * ALLSTAR_IS_MULT, 35.0)
    os_val = is_val + ALLSTAR_OS_CONST

    if not os.path.exists(opt_path):
        print("[OPT] ⚠️  allstar.opt not found in WORK_DIR — creating minimal version.")
        with open(opt_path, 'w') as fh:
            fh.write(f"FI = {fi_val:.3f}\n")
            fh.write(f"IS = {is_val:.3f}\n")
            fh.write(f"OS = {os_val:.3f}\n")
            fh.write(f"RED = 1\n")
        print(f"[OPT] Created minimal allstar.opt with FI={fi_val:.3f}, IS={is_val:.3f}, OS={os_val:.3f}")
        return

    with open(opt_path) as fh:
        lines = fh.readlines()

    updated = {}
    new_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped or '=' not in stripped:
            new_lines.append(line)
            continue
        key = stripped.split('=')[0].strip().upper()[:2]
        if key == 'FI':
            new_lines.append(f"FI = {fi_val:.3f}\n")
            updated['FI'] = True
        elif key == 'IS':
            new_lines.append(f"IS = {is_val:.3f}\n")
            updated['IS'] = True
        elif key == 'OS':
            new_lines.append(f"OS = {os_val:.3f}\n")
            updated['OS'] = True
        else:
            new_lines.append(line)

    if 'FI' not in updated: new_lines.append(f"FI = {fi_val:.3f}\n")
    if 'IS' not in updated: new_lines.append(f"IS = {is_val:.3f}\n")
    if 'OS' not in updated: new_lines.append(f"OS = {os_val:.3f}\n")

    with open(opt_path, 'w') as fh:
        fh.writelines(new_lines)

    print(f"[OPT] Updated allstar.opt — FI={fi_val:.3f}, IS={is_val:.3f}, OS={os_val:.3f}")


def update_photo_opt(fwhm):
    """Patch Apertures and Sky in PHOTO.OPT"""
    opt_path = workpath('photo.opt')
    a1_val = fwhm * PHOTO_A1_MULT
    
    apertures = {'A1': a1_val}
    current_ap = a1_val
    for i in range(2, 10):
        current_ap *= PHOTO_APERTURE_MULT
        apertures[f'A{i}'] = current_ap
    for i in ['A', 'B']:
        current_ap *= PHOTO_APERTURE_MULT
        apertures[f'A{i}'] = current_ap

    ac_val = a1_val * PHOTO_AC_MULT
    is_val = ac_val
    os_val = is_val + PHOTO_OS_CONST

    if not os.path.exists(opt_path):
        print("[OPT] ⚠️  photo.opt not found in WORK_DIR — creating minimal version.")
        with open(opt_path, 'w') as fh:
            for k, v in apertures.items():
                fh.write(f"{k} = {v:.3f}\n")
            fh.write(f"AC = {ac_val:.3f}\n")
            fh.write(f"IS = {is_val:.3f}\n")
            fh.write(f"OS = {os_val:.3f}\n")
        print(f"[OPT] Created minimal photo.opt with A1={a1_val:.3f}, AC={ac_val:.3f}, IS={is_val:.3f}, OS={os_val:.3f}")
        return

    with open(opt_path) as fh:
        lines = fh.readlines()

    updated = {}
    new_lines = []
    for line in lines:
        stripped = line.strip()
        if not stripped or '=' not in stripped:
            new_lines.append(line)
            continue
        key = stripped.split('=')[0].strip().upper()[:2]
        if key in apertures:
            new_lines.append(f"{key} = {apertures[key]:.3f}\n")
            updated[key] = True
        elif key == 'AC':
            new_lines.append(f"AC = {ac_val:.3f}\n")
            updated['AC'] = True
        elif key == 'IS':
            new_lines.append(f"IS = {is_val:.3f}\n")
            updated['IS'] = True
        elif key == 'OS':
            new_lines.append(f"OS = {os_val:.3f}\n")
            updated['OS'] = True
        else:
            new_lines.append(line)

    for k, v in apertures.items():
        if k not in updated:
            new_lines.append(f"{k} = {v:.3f}\n")
    if 'AC' not in updated: new_lines.append(f"AC = {ac_val:.3f}\n")
    if 'IS' not in updated: new_lines.append(f"IS = {is_val:.3f}\n")
    if 'OS' not in updated: new_lines.append(f"OS = {os_val:.3f}\n")

    with open(opt_path, 'w') as fh:
        fh.writelines(new_lines)

    print(f"[OPT] Updated photo.opt — A1={a1_val:.3f}, AC={ac_val:.3f}, IS={is_val:.3f}, OS={os_val:.3f}")


# ============================================================
# CELL 9 : D1 — FIND  (standalone)
# ============================================================

def run_find(image_basename, fwhm):
    """
    Run DAOPHOT FIND.

    image_basename : just the filename (no path) — DAOPHOT runs in WORK_DIR
    fwhm           : used only for reporting; actual value comes from DAOPHOT.OPT

    DAOPHOT command sequence
    ------------------------
    AT  image              → attach image
    FI                     → FIND
    FRAMES_AVG,FRAMES_SUM  → frame averaging / summing info
    <CR>                   → accept default .coo output filename
    Y                      → accept detections
    EX                     → exit

    Returns path to the output .coo file.
    """
    print("\n" + "="*60)
    print("  D1: FIND — Source Detection")
    print("="*60)

    stem     = get_stem(image_basename)
    coo_path = workpath(f"{stem}.coo")
    if os.path.exists(coo_path): os.remove(coo_path)

    # Strip .fits extension for ATTACH — DAOPHOT finds the file itself
    attach_name = stem  # e.g. "stacked_science"

    stdout = run_daophot_session([
        f"AT {attach_name}",            # ATTACH
        "FI",                           # FIND
        f"{FRAMES_AVG},{FRAMES_SUM}",   # averaging/summing info
        "",                             # default .coo filename
        "Y",                            # happy with detections
        "EX",                           # exit DAOPHOT
    ], label="D1-FIND")

    # Count detections
    n = 0
    if os.path.exists(coo_path):
        with open(coo_path) as fh:
            n = sum(1 for ln in fh
                    if ln.strip() and ln.strip()[0].isdigit())
    print(f"  FIND → {n} sources  ({os.path.basename(coo_path)})")

    if n == 0:
        warnings.warn(
            "[D1] FIND detected 0 sources. Check DAOPHOT.OPT (FW, TH, LO, HI) "
            "and verify the image contains stars.",
            RuntimeWarning
        )
    return coo_path, n


# ============================================================
# CELL 10 : D2 — PHOT  (standalone)
# ============================================================

def run_phot(image_basename, coo_file=None, ap_out=None):
    """
    Run DAOPHOT PHOTOMETRY.

    Uses PHOTO.OPT from WORK_DIR for aperture / annulus settings.
    Pass coo_file / ap_out to override the defaults (used for forced phot).

    DAOPHOT command sequence
    ------------------------
    AT  stem
    PH
    <CR>          → accept PHOTO.OPT
    <CR>          → accept current table (no changes)
    [coo_name]    → coordinate file  ('' = default)
    [ap_name]     → output mag file  ('' = default)
    EX
    """
    print("\n" + "="*60)
    print("  D2: PHOT — Aperture Photometry")
    print("="*60)

    stem = get_stem(image_basename)

    # Decide input / output filenames — '' means DAOPHOT uses its default
    coo_arg = os.path.basename(coo_file) if coo_file else ""
    ap_arg  = os.path.basename(ap_out)   if ap_out  else ""

    default_ap = workpath(f"{stem}.ap")
    out_ap = workpath(ap_arg) if ap_arg else default_ap
    if os.path.exists(out_ap): os.remove(out_ap)

    stdout = run_daophot_session([
        f"AT {stem}",     # ATTACH
        "PH",             # PHOTOMETRY
        "",               # accept PHOTO.OPT
        "",               # accept PHO> table unchanged
        coo_arg,          # coordinate file ('' = default ?.coo)
        ap_arg,           # output file     ('' = default ?.ap)
        "EX",
    ], label="D2-PHOT")

    out_ap = workpath(ap_arg) if ap_arg else default_ap
    print(f"  PHOT → {os.path.basename(out_ap)}")
    return out_ap


# ============================================================
# CELL 11 : D3 — PICK  (standalone)
# ============================================================

def run_pick(image_basename, fwhm):
    """
    Run DAOPHOT PICK to select PSF-star candidates.

    DAOPHOT command sequence
    ------------------------
    AT  stem
    PICK
    <CR>                  → default .ap input
    MAX_PSF_STARS PSF_FAINT_MAG
    <CR>                  → default .lst output
    EX
    """
    print("\n" + "="*60)
    print("  D3: PICK — PSF Star Selection")
    print("="*60)

    stem     = get_stem(image_basename)
    lst_path = workpath(f"{stem}.lst")
    if os.path.exists(lst_path): os.remove(lst_path)

    stdout = run_daophot_session([
        f"AT {stem}",
        "PICK",
        "",                                       # default .ap
        f"{MAX_PSF_STARS} {PSF_FAINT_MAG}",       # n_stars  faintest_mag
        "",                                       # default .lst
        "EX",
    ], label="D3-PICK")

    n_psf = 0
    if os.path.exists(lst_path):
        with open(lst_path) as fh:
            n_psf = sum(1 for ln in fh
                        if ln.strip() and ln.strip()[0].isdigit())
    print(f"  PICK → {n_psf} PSF candidates  ({os.path.basename(lst_path)})")

    if n_psf < 3:
        warnings.warn(
            f"[D3] Only {n_psf} PSF stars selected. "
            "Consider increasing PSF_FAINT_MAG or reducing DETECT_SIGMA.",
            RuntimeWarning
        )
    return lst_path, n_psf


# ============================================================
# CELL 12 : D4 — PSF  (standalone, iterative)
# ============================================================

def run_psf_iterative(image_basename):
    """
    Build the empirical PSF with PSF_ITERATIONS of the
        PSF → ALLSTAR(neighbors) → SUBSTAR → PSF(cleaner) → …
    loop recommended in the DAOPHOT II manual (Appendix III).

    Iteration 0  : pure analytic (VA=-1) — avoids neighbor ghosts
                   contaminating the look-up table on first pass.
    Iteration 1+ : full analytic + look-up table (VA=VARIABLE_PSF).

    The DAOPHOT OPTIONS command is used within each session to set
    VA and WA without touching DAOPHOT.OPT.

    Returns path to the final .psf file.
    """
    print("\n" + "="*60)
    print(f"  D4: PSF BUILD  ({PSF_ITERATIONS} iterations)")
    print("="*60)

    stem       = get_stem(image_basename)
    psf_path   = workpath(f"{stem}.psf")
    nei_path   = workpath(f"{stem}.nei")

    # Current working image — starts as original, replaced after each SUBSTAR
    work_img   = stem            # DAOPHOT attach name (no extension, no path)
    work_img_f = image_basename  # full path for ALLSTAR / file-existence checks

    for iteration in range(PSF_ITERATIONS):
        va = -1 if iteration == 0 else VARIABLE_PSF
        print(f"\n  --- PSF iteration {iteration + 1}/{PSF_ITERATIONS} "
              f"  VA={va} ---")

        # ── Build PSF ──────────────────────────────────────────────────
        # OPTIONS: set WA=-2 (silent) and VA=va
        # Then PSF using default .ap, default .lst, default .psf
        
        if os.path.exists(psf_path): os.remove(psf_path)

        stdout = run_daophot_session([
            f"AT {work_img}",    # attach current work image
            "OP",                # OPTIONS
            "",                  # keyboard input
            "WA=-2",             # silent mode — essential for subprocess
            f"VA={va}",          # PSF variation degree
            "",                  # exit OPTIONS
            "PSF",               # PSF command
            "",                  # default .ap  (full aperture photometry)
            "",                  # default .lst (PSF candidate list)
            "",                  # default .psf (output)
            "EX",
        ], label=f"D4-PSF-iter{iteration+1}")

        if not os.path.exists(psf_path):
            raise RuntimeError(
                f"[D4] PSF file not created on iteration {iteration + 1}. "
                "Check DAOPHOT.OPT (WA must be -1 or -2 for non-interactive PSF). "
                "DAOPHOT stdout is printed above."
            )
        print(f"  PSF written: {psf_path}")

        # ── On all but the last iteration: reduce neighbours, then clean ──
        if iteration < PSF_ITERATIONS - 1 and os.path.exists(nei_path):
            print("  Reducing PSF-star neighbours (ALLSTAR on .nei) …")

            nei_als  = workpath(f"{stem}_nei_iter{iteration+1}.als")
            nei_sub  = workpath(f"{stem}_clean_iter{iteration+1}.fits")
            nei_sub_stem = f"{stem}_clean_iter{iteration+1}"

            # Remove previous outputs
            for f in [nei_als, nei_sub]:
                if os.path.exists(f): os.remove(f)

            # ALLSTAR on .nei file
            #   Input image : current work image
            #   PSF         : default ?.psf (matches current work_img stem — may differ)
            #   Input file  : stem.nei  (we specify explicitly)
            #   Results     : nei_als
            #   Sub image   : nei_sub_stem
            run_allstar_session([
                work_img_f,                       # full image path for ALLSTAR
                "",                               # default .psf
                os.path.basename(nei_path),       # .nei as input
                os.path.basename(nei_als),        # custom output .als
                nei_sub_stem,                     # custom subtracted image name
            ], label=f"D4-ALLSTAR-nei-iter{iteration+1}")

            if not os.path.exists(nei_als):
                print("  ⚠️  Neighbour ALLSTAR failed — proceeding with same image.")
                continue

            # SUBSTAR: subtract all neighbours but keep PSF stars
            #   File with photometry : nei_als
            #   Stars to leave in    : .lst  (our PSF candidates)
            #   Output image name    : nei_sub_stem  (explicit, avoid ?s default)
            if os.path.exists(nei_sub):
                os.remove(nei_sub)

            stdout = run_daophot_session([
                f"AT {work_img}",                 # original or previous clean image
                "SU",                             # SUBSTAR
                "",                               # default .psf
                os.path.basename(nei_als),        # neighbour photometry
                "Y",                              # yes, we have stars to leave in
                "",                               # default .lst (keep PSF stars)
                nei_sub_stem,                     # output image name
                "EX",
            ], label=f"D4-SUBSTAR-iter{iteration+1}")

            if os.path.exists(nei_sub):
                # Next PSF will be derived from the cleaned image
                work_img   = nei_sub_stem
                work_img_f = nei_sub
                print(f"  Next PSF iteration uses: {nei_sub_stem}")
            else:
                print("  ⚠️  SUBSTAR output not found — using same image.")

    print(f"\n  ✅ Final PSF: {psf_path}")
    return psf_path


# ============================================================
# CELL 13 : D5 — ALLSTAR  (standalone, full frame)
# ============================================================

def run_allstar_full(image_basename, ap_file=None):
    """
    Run standalone ALLSTAR on the complete star list.

    image_basename : filename of the original image (in WORK_DIR)
    ap_file        : .ap file to use; defaults to stem.ap

    ALLSTAR command sequence
    ------------------------
    image_name
    <CR>          → default .psf
    [ap_basename] → input photometry file
    <CR>          → default .als output
    <CR>          → default subtracted image

    Returns (als_path, sub_image_path, als_table).
    """
    print("\n" + "="*60)
    print("  D5: ALLSTAR — Full-Frame PSF Fitting")
    print("="*60)

    stem      = get_stem(image_basename)
    als_path  = workpath(f"{stem}.als")
    # ALLSTAR default subtracted image: ?s  →  stems (no fits extension noted)
    # Specify explicitly to be safe
    sub_stem  = f"{stem}s"
    sub_path  = workpath(f"{sub_stem}.fits")

    ap_arg = os.path.basename(ap_file) if ap_file else ""

    for f in [als_path, sub_path]:
        if os.path.exists(f): os.remove(f)

    run_allstar_session([
        "",               # Exit the OPT> menu that ALLSTAR opens automatically
        stem,             # bare filename without extension (ALLSTAR auto-adds .fits)
        "",               # default .psf
        ap_arg,           # .ap file  ('' = default ?.ap)
        "",               # default .als
        sub_stem,         # subtracted image name (explicit)
    ], label="D5-ALLSTAR")

    als_table = parse_als_file(als_path) if os.path.exists(als_path) else None
    n = len(als_table) if als_table is not None else 0
    print(f"  ALLSTAR fitted {n} stars  →  {os.path.basename(als_path)}")
    print(f"  Subtracted image         →  {sub_stem}.fits")
    return als_path, sub_path, als_table


# ============================================================
# CELL 14 : D6 — FORCED PHOTOMETRY  (standalone)
# ============================================================

def run_forced_photometry(image_basename, als_table, fwhm, wcs,
                           force_ra, force_dec):
    """
    Inject a manual RA/Dec position when DAOFIND misses the target.

    Workflow (all within standalone DAOPHOT)
    ----------------------------------------
    1. Convert RA/Dec → pixel via WCS; optionally re-centre
    2. Write a .coo file with the forced position
    3. DAOPHOT PHOT on that .coo  → forced.ap
    4. DAOPHOT OFFSET on forced.ap  (add +100 000 to IDs)  → forced_off.ap
    5. DAOPHOT APPEND  stem.ap + forced_off.ap  → stem_merged.ap
    6. ALLSTAR on merged .ap  → stem_merged.als

    Returns (merged_als_path, merged_table, forced_x, forced_y).
    """
    print("\n" + "="*60)
    print("  D6: FORCED PHOTOMETRY")
    print(f"  Target  RA={force_ra:.5f}  Dec={force_dec:.5f}")
    print("="*60)

    stem = get_stem(image_basename)

    # 1. RA/Dec → pixel
    forced_x, forced_y = wcs.wcs_world2pix(force_ra, force_dec, 0)
    forced_x, forced_y = float(forced_x), float(forced_y)
    print(f"  WCS pixel position: X={forced_x:.2f}, Y={forced_y:.2f}")

    # Optional re-centring on local peak
    if FORCE_CENTROID:
        forced_x, forced_y = _recenter(
            workpath(image_basename), forced_x, forced_y, fwhm, wcs
        )

    # 2. Write forced .coo
    coo_template = workpath(f"{stem}.coo")
    forced_coo   = workpath(f"{stem}_forced.coo")
    create_forced_coo(forced_x, forced_y, coo_template, forced_coo,
                      id_offset=100000)

    # 3. PHOT on forced position
    #    Use centroiding='none' by setting it in OPTIONS before PH.
    #    This prevents DAOPHOT from wandering off the forced position.
    forced_ap  = workpath(f"{stem}_forced.ap")
    if os.path.exists(forced_ap): os.remove(forced_ap)

    psf_file = workpath(f"{stem}.psf")
    psf_hidden = False
    if os.path.exists(psf_file):
        os.rename(psf_file, psf_file + ".bak")
        psf_hidden = True

    try:
        run_daophot_session([
            f"AT {stem}",
            "PH",
            "",          # default PHOTO.OPT
            "",          # accept table
            f"{stem}_forced.coo",       # forced coordinate file
            f"{stem}_forced.ap",        # output
            "EX",
        ], label="D6-FORCED-PHOT")
    finally:
        if psf_hidden:
            os.rename(psf_file + ".bak", psf_file)

    if not os.path.exists(forced_ap):
        print("  ❌ Forced PHOT failed — no .ap output.")
        return None, als_table, forced_x, forced_y

    # 4. OFFSET: add 100 000 to IDs so there's no clash
    forced_off = workpath(f"{stem}_forced_off.ap")
    if os.path.exists(forced_off): os.remove(forced_off)

    run_daophot_session([
        # OFFSET works without an attached image
        "OF",
        f"{stem}_forced.ap",        # input file
        "100000 0 0 0/",            # ID_offset DX DY DMAG  (/ terminates)
        f"{stem}_forced_off.ap",    # output file
        "EX",
    ], label="D6-OFFSET")

    # 5. APPEND: merge original .ap with forced_off .ap
    merged_ap = workpath(f"{stem}_merged.ap")
    if os.path.exists(merged_ap): os.remove(merged_ap)

    run_daophot_session([
        "AP",                        # APPEND (no image needed)
        f"{stem}.ap",                # first input
        f"{stem}_forced_off.ap",     # second input
        f"{stem}_merged.ap",         # output
        "EX",
    ], label="D6-APPEND")

    if not os.path.exists(merged_ap):
        print("  ❌ APPEND failed — using original .als.")
        return None, als_table, forced_x, forced_y

    # 6. ALLSTAR on merged list
    merged_als  = workpath(f"{stem}_merged.als")
    merged_sub  = workpath(f"{stem}_merged_sub.fits")
    merged_stem = f"{stem}_merged_sub"

    for f in [merged_als, merged_sub]:
        if os.path.exists(f): os.remove(f)

    run_allstar_session([
        "",                              # Exit OPT> menu
        stem,                            # bare filename without extension
        "",                              # default .psf
        f"{stem}_merged.ap",             # merged input
        f"{stem}_merged.als",            # output
        merged_stem,                     # subtracted image
    ], label="D6-ALLSTAR-FORCED")

    merged_table = parse_als_file(merged_als) if os.path.exists(merged_als) else None
    n = len(merged_table) if merged_table is not None else 0
    print(f"  Forced+merged ALLSTAR: {n} stars  →  "
          f"{os.path.basename(merged_als)}")
    return merged_als, merged_table, forced_x, forced_y


def _recenter(image_path, x, y, fwhm, wcs=None):
    """Snap (x, y) to the nearest Gaussian centroid within 2×FWHM."""
    try:
        from astropy.nddata import Cutout2D
        from photutils.centroids import centroid_2dg
        with fits.open(image_path) as hdul:
            data = hdul[0].data.astype(float)
        box  = max(int(2 * fwhm) * 2 + 1, 7)
        cut  = Cutout2D(data, (x, y), box)
        cx, cy = centroid_2dg(cut.data)
        nx, ny = cut.to_original_position((cx, cy))
        shift  = np.hypot(nx - x, ny - y)
        if not (np.isnan(nx) or np.isnan(ny)) and shift < 2.0 * fwhm:
            if wcs:
                ra, dec = wcs.pixel_to_world_values(nx, ny)
                print(f"  Re-centred {shift:.2f}px → "
                      f"X={nx:.2f} Y={ny:.2f}  "
                      f"(RA={float(ra):.5f} Dec={float(dec):.5f})")
            else:
                print(f"  Re-centred {shift:.2f}px → X={nx:.2f} Y={ny:.2f}")
            return float(nx), float(ny)
        else:
            print(f"  ⚠️  Re-centre shift too large ({shift:.2f}px) — "
                  "using WCS position.")
    except Exception as exc:
        print(f"  ⚠️  Re-centre failed ({exc}) — using WCS position.")
    return x, y


# ============================================================
# CELL 15 : D7 — SOURCE LOOKUP  (unchanged)
# ============================================================

def lookup_target_in_als(als_table, wcs, target_ra, target_dec,
                          match_radius_pix=MATCH_RADIUS_PIX):
    print("\n" + "="*60)
    print("  D7: TARGET LOOKUP IN .ALS")
    print("="*60)

    if als_table is None or len(als_table) == 0:
        print("  ❌ .als table is empty.")
        return _empty_result()

    target_x, target_y = wcs.wcs_world2pix(target_ra, target_dec, 0)
    target_x, target_y = float(target_x), float(target_y)
    print(f"  Target expected at X={target_x:.2f}, Y={target_y:.2f}")

    dx   = als_table['XCENTER'] - target_x
    dy   = als_table['YCENTER'] - target_y
    dist = np.sqrt(dx**2 + dy**2)
    idx  = int(np.argmin(dist))
    dmin = float(dist[idx])

    if dmin <= match_radius_pix:
        row = als_table[idx]
        print(f"  ✅ Matched!  ID={row['ID']}  "
              f"X={row['XCENTER']:.2f}  Y={row['YCENTER']:.2f}  "
              f"dist={dmin:.2f}px")
        print(f"  MAG={row['MAG']:.3f}±{row['MERR']:.3f}  "
              f"CHI={row['CHI']:.2f}  SHARP={row['SHARP']:.3f}")
        return {
            'detected': True,
            'mag'     : float(row['MAG']),
            'merr'    : float(row['MERR']),
            'chi'     : float(row['CHI']),
            'sharp'   : float(row['SHARP']),
            'pixel_x' : float(row['XCENTER']),
            'pixel_y' : float(row['YCENTER']),
            'dist_pix': dmin,
            'id'      : int(row['ID']),
        }
    else:
        print(f"  ❌ No match within {match_radius_pix:.1f}px  "
              f"(closest = {dmin:.1f}px)")
        return _empty_result()


def _empty_result():
    return dict(detected=False, mag=np.nan, merr=np.nan,
                chi=np.nan, sharp=np.nan, pixel_x=np.nan,
                pixel_y=np.nan, dist_pix=np.nan, id=-1)


# ============================================================
# CELL 16 : D8 — ZERO-POINT CALIBRATION  (unchanged)
# ============================================================

def calibrate_zero_point_daophot(als_table, wcs, image_shape, target_band=None, pipeline_mode='auto', image_path=None):
    if target_band is None:
        target_band = FILTER_BAND

    print("\n" + "="*60)
    print(f"  D8: ZERO-POINT CALIBRATION  "
          f"({VIZIER_CATALOG_ID}, band={target_band})")
    print("="*60)

    if als_table is None or len(als_table) == 0:
        print("  ❌ ALLSTAR table empty — cannot calibrate.")
        return np.nan, np.nan, None

    # Get instrumental magnitudes
    als_ra, als_dec = wcs.wcs_pix2world(
        als_table['XCENTER'], als_table['YCENTER'], 0
    )
    catalog_coords = SkyCoord(ra=als_ra * u.deg, dec=als_dec * u.deg)

    m_inst     = np.array(als_table['MAG'],  dtype=float)
    m_inst_err = np.array(als_table['MERR'], dtype=float)
    
    valid = (np.isfinite(m_inst) & (m_inst < 90) &
             np.isfinite(m_inst_err) & (m_inst_err < 1.0))
    m_inst[~valid] = np.nan
    m_inst_err[~valid] = np.nan

    ny, nx   = image_shape
    cra, cdec = wcs.wcs_pix2world(nx / 2, ny / 2, 0)
    center    = SkyCoord(cra, cdec, unit=(u.deg, u.deg))

    v = Vizier(columns=VIZIER_COLUMNS, row_limit=VIZIER_ROW_LIMIT)
    result = v.query_region(center, radius=QUERY_RADIUS_ARCMIN*u.arcmin, catalog=VIZIER_CATALOG_ID)

    if not result:
        print("  ❌ VizieR query returned nothing.")
        return np.nan, np.nan, None

    catalog_table = result[0]
    
    if CALIBRATION_SYSTEM == 'Pan-STARRS':
        has_valid_data = (
                        ~np.ma.getmaskarray(catalog_table['gmag']) & ~np.ma.getmaskarray(catalog_table['e_gmag']) &
                        ~np.ma.getmaskarray(catalog_table['rmag']) & ~np.ma.getmaskarray(catalog_table['e_rmag']) &
                        ~np.ma.getmaskarray(catalog_table['imag']) & ~np.ma.getmaskarray(catalog_table['e_imag'])
                        )
        good_stars = catalog_table[
            has_valid_data &
            (catalog_table['gmag'] > REF_MAG_BRIGHT) & (catalog_table['gmag'] < REF_MAG_FAINT) & (catalog_table['e_gmag'] < REF_MAG_ERR_MAX) &
            (catalog_table['rmag'] > REF_MAG_BRIGHT) & (catalog_table['rmag'] < REF_MAG_FAINT) & (catalog_table['e_rmag'] < REF_MAG_ERR_MAX) &
            (catalog_table['imag'] > REF_MAG_BRIGHT) & (catalog_table['imag'] < REF_MAG_FAINT) & (catalog_table['e_imag'] < REF_MAG_ERR_MAX)
        ]

        print(f"  VizieR: Downloaded {len(catalog_table)} stars, {len(good_stars)} pass quality cuts.")

        if len(good_stars) == 0:
            print("  ❌ No reference stars passed the quality cuts.")
            return np.nan, np.nan, None

        g = good_stars['gmag']
        r = good_stars['rmag']
        i = good_stars['imag']
        e_g = good_stars['e_gmag']
        e_r = good_stars['e_rmag']
        e_i = good_stars['e_imag']

        if target_band == 'V':
            true_mag = g - 0.5784 * (g - r) - 0.0038
            true_mag_err = np.sqrt((0.4216**2) * (e_g**2) + (0.5784**2) * (e_r**2))
        elif target_band == 'B':
            true_mag = g + 0.3130 * (g - r) + 0.2271
            true_mag_err = np.sqrt((1.3130**2) * (e_g**2) + (0.3130**2) * (e_r**2))
        elif target_band == 'R':
            true_mag = r - 0.2936 * (r - i) - 0.1439
            true_mag_err = np.sqrt((0.7064**2) * (e_r**2) + (0.2936**2) * (e_i**2))
        elif target_band == 'I':
            true_mag = r - 1.2444 * (r - i) - 0.3820
            true_mag_err = np.sqrt((0.2444**2) * (e_r**2) + (1.2444**2) * (e_i**2))
        elif target_band in ['g', 'r', 'i']:
            if target_band == 'g': true_mag, true_mag_err = g, e_g
            elif target_band == 'r': true_mag, true_mag_err = r, e_r
            elif target_band == 'i': true_mag, true_mag_err = i, e_i
        else:
            print(f"  ⚠️ Unknown target_band '{target_band}'. Transformation failed.")
            true_mag = np.full(len(good_stars), np.nan)
            true_mag_err = np.full(len(good_stars), np.nan)
            
        ra_col = 'RAJ2000'
        dec_col = 'DEJ2000'

    elif CALIBRATION_SYSTEM == 'DIRECT_VEGA':
        has_valid_data = (
            ~np.ma.getmaskarray(catalog_table[VIZIER_MAG_COL]) &
            ~np.ma.getmaskarray(catalog_table[VIZIER_ERR_COL])
        )
        good_stars = catalog_table[
            has_valid_data &
            (catalog_table[VIZIER_MAG_COL] > REF_MAG_BRIGHT) & 
            (catalog_table[VIZIER_MAG_COL] < REF_MAG_FAINT) & 
            (catalog_table[VIZIER_ERR_COL] < REF_MAG_ERR_MAX)
        ]
        print(f"  VizieR: Downloaded {len(catalog_table)} stars, {len(good_stars)} pass quality cuts in {VIZIER_MAG_COL}.")
        
        if len(good_stars) == 0:
            print("  ❌ No reference stars passed the quality cuts.")
            return np.nan, np.nan, None
            
        true_mag = good_stars[VIZIER_MAG_COL]
        true_mag_err = good_stars[VIZIER_ERR_COL]
        
        ra_col  = 'RA_ICRS' if 'RA_ICRS' in good_stars.colnames else 'RAJ2000'
        dec_col = 'DE_ICRS' if 'DE_ICRS' in good_stars.colnames else 'DEJ2000'

    else:
        print(f"  ❌ Unknown CALIBRATION_SYSTEM '{CALIBRATION_SYSTEM}'.")
        return np.nan, np.nan, None

    good_stars['true_mag'] = true_mag
    good_stars['true_mag_err'] = true_mag_err

    ref_coords = SkyCoord(ra=good_stars[ra_col], dec=good_stars[dec_col])
    idx, d2d, _ = ref_coords.match_to_catalog_sky(catalog_coords)

    # --- Galaxy exclusion zone ---
    ny, nx = image_shape
    img_center_ra, img_center_dec = wcs.wcs_pix2world(nx / 2, ny / 2, 0)
    img_center_coord = SkyCoord(img_center_ra, img_center_dec, unit=(u.deg, u.deg))

    # --- Extract ALLSTAR quality columns ---
    als_chi   = np.array(als_table['CHI'],   dtype=float)
    als_sharp = np.array(als_table['SHARP'], dtype=float)

    calculated_zps = []
    zp_indices = []
    reject_reasons = []  # track why each matched star was kept or auto-rejected
    
    n_matched = 0
    n_chi_rej = 0
    n_sharp_rej = 0
    n_merr_rej = 0
    n_gal_rej = 0
    n_edge_rej = 0
    
    for k in range(len(good_stars)):
        if d2d[k].arcsec < ZP_MATCH_RADIUS_AS:
            als_k = idx[k]
            inst_mag = m_inst[als_k]
            inst_err = m_inst_err[als_k]
            chi_val  = als_chi[als_k]
            sharp_val = als_sharp[als_k]
            
            if np.isnan(inst_mag):
                continue
            
            n_matched += 1
            reason = ''
            
            # Stellarity filters (tighter than catalog-wide)
            if chi_val > ZP_MAX_CHI:
                reason = f'CHI={chi_val:.2f}>{ZP_MAX_CHI}'
                n_chi_rej += 1
            elif abs(sharp_val) > ZP_MAX_SHARP_ABS:
                reason = f'|SHARP|={abs(sharp_val):.3f}>{ZP_MAX_SHARP_ABS}'
                n_sharp_rej += 1
            elif inst_err > ZP_MAX_MERR:
                reason = f'MERR={inst_err:.3f}>{ZP_MAX_MERR}'
                n_merr_rej += 1
            
            # Galaxy exclusion zone
            if not reason and GALAXY_EXCL_RADIUS_ARCSEC > 0:
                star_coord = SkyCoord(
                    ra=good_stars[ra_col][k], dec=good_stars[dec_col][k])
                sep = img_center_coord.separation(star_coord).arcsec
                if sep < GALAXY_EXCL_RADIUS_ARCSEC:
                    reason = f'galaxy_zone ({sep:.0f}"<{GALAXY_EXCL_RADIUS_ARCSEC}")'
                    n_gal_rej += 1
            
            # Edge exclusion filter
            if not reason and EDGE_EXCL_MARGIN_PIXELS > 0:
                x_c = als_table['XCENTER'][als_k]
                y_c = als_table['YCENTER'][als_k]
                if (x_c < EDGE_EXCL_MARGIN_PIXELS or x_c > nx - EDGE_EXCL_MARGIN_PIXELS or
                    y_c < EDGE_EXCL_MARGIN_PIXELS or y_c > ny - EDGE_EXCL_MARGIN_PIXELS):
                    reason = f'edge_zone (margin={EDGE_EXCL_MARGIN_PIXELS}px)'
                    n_edge_rej += 1
            
            zp_val = true_mag[k] - inst_mag
            calculated_zps.append(zp_val)
            zp_indices.append((k, als_k))
            reject_reasons.append(reason)  # '' = passed all filters

    # Report auto-filter results
    n_auto_pass = sum(1 for r in reject_reasons if r == '')
    print(f"  Matched {n_matched} catalog↔ALLSTAR pairs, "
          f"{n_auto_pass} pass stellarity filters")
    if n_chi_rej:   print(f"    ✂ CHI filter:   {n_chi_rej} rejected")
    if n_sharp_rej: print(f"    ✂ SHARP filter: {n_sharp_rej} rejected")
    if n_merr_rej:  print(f"    ✂ MERR filter:  {n_merr_rej} rejected")
    if n_gal_rej:   print(f"    ✂ Galaxy zone:  {n_gal_rej} rejected")
    if n_edge_rej:  print(f"    ✂ Edge zone:    {n_edge_rej} rejected")

    # Build enriched ZP table (including auto-rejected for visibility)
    als_ids, xcenters, ycenters = [], [], []
    matched_insts, matched_refs = [], []
    matched_chis, matched_sharps, matched_merrs = [], [], []
    auto_rejected = []

    for i, (ref_k, als_k) in enumerate(zp_indices):
        als_ids.append(als_table['ID'][als_k])
        xcenters.append(als_table['XCENTER'][als_k])
        ycenters.append(als_table['YCENTER'][als_k])
        matched_insts.append(m_inst[als_k])
        matched_refs.append(true_mag[ref_k])
        matched_chis.append(als_chi[als_k])
        matched_sharps.append(als_sharp[als_k])
        matched_merrs.append(m_inst_err[als_k])
        auto_rejected.append(reject_reasons[i])

    if len(calculated_zps) == 0:
        print("  ⚠️ No matches found for zero-point calibration.")
        return np.nan, np.nan, None

    # Sigma-clip ONLY the stars that passed auto filters
    calculated_zps = np.array(calculated_zps)
    auto_pass_mask = np.array([r == '' for r in reject_reasons])
    
    if np.sum(auto_pass_mask) < ZP_MIN_STARS:
        print(f"  ⚠️ Only {np.sum(auto_pass_mask)} stars passed stellarity filters "
              f"(need {ZP_MIN_STARS}). Relaxing to all matched stars.")
        auto_pass_mask[:] = True
        auto_rejected = ['' for _ in auto_rejected]

    # Sigma clip among passing stars
    pass_zps = calculated_zps[auto_pass_mask]
    sc = SigmaClip(sigma=2.5, maxiters=5) # Tightened from 3.0 to 2.5
    clipped_result = sc(pass_zps)
    clip_mask_pass = clipped_result.mask if hasattr(clipped_result, 'mask') else np.zeros(len(pass_zps), dtype=bool)

    # Hard absolute scatter cut
    if len(pass_zps[~clip_mask_pass]) > 0:
        temp_median = np.median(pass_zps[~clip_mask_pass])
        for i in range(len(pass_zps)):
            if not clip_mask_pass[i] and abs(pass_zps[i] - temp_median) > ZP_MAX_SCATTER_MAG:
                clip_mask_pass[i] = True

    # Build full clipped mask (auto-rejected stars are marked as clipped too)
    full_clipped = np.ones(len(calculated_zps), dtype=bool)
    j = 0
    for i in range(len(calculated_zps)):
        if auto_pass_mask[i]:
            full_clipped[i] = clip_mask_pass[j]
            j += 1

    good_zp = calculated_zps[~full_clipped]

    if len(good_zp) == 0:
        print("  ⚠️ Not enough valid matches for a reliable zero point after clipping.")
        return np.nan, np.nan, None

    master_zp = float(np.median(good_zp))
    master_zp_err = float(np.std(good_zp, ddof=1) / np.sqrt(len(good_zp))) if len(good_zp) > 1 else 0.05

    print(f"  ✅ Master Zero Point = {master_zp:.4f} ± {master_zp_err:.4f}  "
          f"(from {len(good_zp)}/{len(calculated_zps)} stars)")

    zp_table = Table({
        'als_id'       : als_ids,
        'xcenter'      : xcenters,
        'ycenter'      : ycenters,
        'm_inst'       : matched_insts,
        'm_ref'        : matched_refs,
        'zp_star'      : calculated_zps,
        'chi'          : matched_chis,
        'sharp'        : matched_sharps,
        'merr'         : matched_merrs,
        'auto_reject'  : auto_rejected,
        'clipped'      : full_clipped,
        'manual_reject': np.zeros(len(als_ids), dtype=bool),
    })

    # ── Interactive review (manual / semi-auto modes) ──────────────
    if pipeline_mode in ('manual', 'semi-auto'):

        # ── Show star map image BEFORE prompt ──────────────────
        if image_path is not None:
            try:
                with fits.open(image_path) as hdul:
                    _preview_data = hdul[0].data.astype(float)
                _zs = ZScaleInterval()
                _vmin, _vmax = _zs.get_limits(_preview_data)
                _fig, _ax = plt.subplots(1, 1, figsize=(12, 10))
                _ax.imshow(_preview_data, cmap='gray', origin='lower',
                           vmin=_vmin, vmax=_vmax)
                _ax.set_title(f"ZP Calibration Stars — review numbers below",
                              fontsize=14, color='white')
                _ax.set_facecolor('black')
                for i, row in enumerate(zp_table):
                    x, y = row['xcenter'], row['ycenter']
                    is_auto = row['auto_reject'] != ''
                    is_clip = bool(row['clipped'])
                    if is_auto:
                        clr = 'red'
                    elif is_clip:
                        clr = 'orange'
                    else:
                        clr = 'lime'
                    _ax.add_patch(patches.Circle((x, y), radius=15,
                                  edgecolor=clr, facecolor='none', lw=2))
                    _ax.text(x+18, y+18, f"{i+1}",
                             color=clr, fontsize=11, fontweight='bold')
                plt.tight_layout()
                plt.show()
            except Exception as e:
                print(f"  ⚠️ Could not generate preview: {e}")

        # ── Point source criteria guide ────────────────────────
        print("\n" + "─"*90)
        print("  📋 INTERACTIVE ZP STAR REVIEW")
        print("─"*90)
        print("  ┌─────────────────────────────────────────────────────────┐")
        print("  │  POINT SOURCE CRITERIA (ideal values for good stars):  │")
        print("  │    CHI   ≈ 1.0  (good PSF fit,  reject if > 2.0)      │")
        print("  │    SHARP ≈ 0.0  (point source,  reject if |val| > 0.5)│")
        print("  │    MERR  < 0.05 (well-measured,  caution if > 0.1)    │")
        print("  │    Resid ≈ 0.0  (consistent ZP,  suspect if > ±0.1)   │")
        print("  └─────────────────────────────────────────────────────────┘")
        print(f"\n  {'#':>3} {'ALS_ID':>7} {'m_inst':>8} {'m_ref':>8} "
              f"{'ZP':>8} {'Resid':>8} {'CHI':>6} {'SHARP':>7} {'MERR':>6} {'Status'}")
        print("  " + "-"*90)

        for i, row in enumerate(zp_table):
            resid = row['zp_star'] - master_zp
            # Flag suspicious values
            chi_flag   = '⚠' if row['chi'] > 2.0 else ' '
            sharp_flag = '⚠' if abs(row['sharp']) > 0.5 else ' '
            merr_flag  = '⚠' if row['merr'] > 0.05 else ' '
            resid_flag = '⚠' if abs(resid) > 0.1 else ' '

            status = '✓ USED'
            if row['auto_reject']:
                status = f'✂ AUTO ({row["auto_reject"]})'
            elif row['clipped']:
                status = '✂ SIGMA-CLIP'
            print(f"  {i+1:>3} {row['als_id']:>7} {row['m_inst']:>8.3f} "
                  f"{row['m_ref']:>8.3f} {row['zp_star']:>8.4f} "
                  f"{resid_flag}{resid:>+7.4f} {chi_flag}{row['chi']:>5.2f} {sharp_flag}{row['sharp']:>+6.3f} "
                  f"{merr_flag}{row['merr']:>5.3f}  {status}")

        print(f"\n  Current ZP = {master_zp:.4f} ± {master_zp_err:.4f}")
        print("\n  Enter star numbers to REJECT (comma-separated), or press Enter to accept all:")
        print("  Example: 3,5,7   |   Look for ⚠ flags above")
        
        try:
            user_input = input("  > ").strip()
        except EOFError:
            user_input = ''

        if user_input:
            try:
                reject_nums = [int(x.strip()) for x in user_input.split(',') if x.strip()]
                for num in reject_nums:
                    if 1 <= num <= len(zp_table):
                        zp_table['manual_reject'][num - 1] = True
                        print(f"    ❌ Rejected star #{num} (ALS_ID={zp_table['als_id'][num-1]})")
                    else:
                        print(f"    ⚠️ Star #{num} out of range, ignoring.")
            except ValueError:
                print("    ⚠️ Could not parse input, keeping all stars.")

            # Recompute ZP excluding manual rejections
            final_mask = (~zp_table['clipped']) & (~np.array(zp_table['manual_reject'])) & \
                         np.array([r == '' for r in zp_table['auto_reject']])
            final_zps = np.array(zp_table['zp_star'])[final_mask]

            if len(final_zps) >= 2:
                master_zp = float(np.median(final_zps))
                master_zp_err = float(np.std(final_zps, ddof=1) / np.sqrt(len(final_zps)))
            elif len(final_zps) == 1:
                master_zp = float(final_zps[0])
                master_zp_err = 0.05
            else:
                print("  ⚠️ No stars remaining after manual rejection!")
                return np.nan, np.nan, zp_table

            print(f"\n  ✅ Updated Zero Point = {master_zp:.4f} ± {master_zp_err:.4f}  "
                  f"(from {len(final_zps)} stars after manual review)")
        else:
            print("    ✓ All stars accepted.")

    return master_zp, master_zp_err, zp_table


# ============================================================
# CELL 17 : D9 — APPLY ZP  (unchanged)
# ============================================================

def compute_final_magnitude(target_result, zp, zp_err):
    r = dict(target_result)
    if not r.get('detected', False) or np.isnan(r.get('mag', np.nan)) or np.isnan(zp):
        r['cal_mag']     = np.nan
        r['cal_mag_err'] = np.nan
        return r
    r['cal_mag']     = r['mag'] + zp
    r['cal_mag_err'] = float(np.sqrt(r['merr']**2 + zp_err**2))
    return r

def compute_all_magnitudes(als_table, zp, zp_err, wcs=None):
    if als_table is None or len(als_table) == 0:
        return None, None
    
    all_cal_mags = []
    all_cal_mag_errs = []
    filt_cal_mags = []
    filt_cal_mag_errs = []
    
    for row in als_table:
        mag = float(row['MAG'])
        merr = float(row['MERR'])
        chi = float(row['CHI'])
        sharp = float(row['SHARP'])
        
        # Unfiltered magnitudes
        if np.isnan(mag) or mag > 89 or np.isnan(zp):
            all_cal_mags.append(np.nan)
            all_cal_mag_errs.append(np.nan)
        else:
            all_cal_mags.append(mag + zp)
            all_cal_mag_errs.append(np.sqrt(merr**2 + zp_err**2))
            
        # Filtered magnitudes
        if (chi > MAX_CHI or abs(sharp) > MAX_SHARP_ABS) or np.isnan(mag) or mag > 89 or np.isnan(zp):
            filt_cal_mags.append(np.nan)
            filt_cal_mag_errs.append(np.nan)
        else:
            filt_cal_mags.append(mag + zp)
            filt_cal_mag_errs.append(np.sqrt(merr**2 + zp_err**2))
            
    out_all = als_table.copy()
    out_all['CAL_MAG'] = all_cal_mags
    out_all['CAL_MAG_ERR'] = all_cal_mag_errs
    
    out_filt = als_table.copy()
    out_filt['CAL_MAG'] = filt_cal_mags
    out_filt['CAL_MAG_ERR'] = filt_cal_mag_errs
    
    if wcs is not None:
        ras, decs = wcs.wcs_pix2world(out_all['XCENTER'], out_all['YCENTER'], 0)
        out_all['RA'] = ras
        out_all['DEC'] = decs
        out_filt['RA'] = ras
        out_filt['DEC'] = decs
        
    return out_all, out_filt


# ============================================================
# CELL 18 : DIAGNOSTIC PLOTS  (unchanged)
# ============================================================

def plot_chi_sharp(als_table, target_name, save=True):
    if als_table is None: return
    mag   = np.array(als_table['MAG'],   dtype=float)
    chi   = np.array(als_table['CHI'],   dtype=float)
    sharp = np.array(als_table['SHARP'], dtype=float)
    valid = (mag < 90) & np.isfinite(chi) & np.isfinite(sharp)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    axes[0].scatter(mag[valid], chi[valid],   s=5, alpha=0.4, color='steelblue')
    axes[0].axhline(1.0, color='red', lw=1.5, ls='--', label='CHI=1 (ideal)')
    axes[0].set_xlabel("Instrumental Magnitude")
    axes[0].set_ylabel("CHI")
    axes[0].set_title(f"{target_name} — CHI Diagnostic")
    axes[0].legend()
    axes[0].set_ylim(0, min(chi[valid].max()*1.1, 5))

    axes[1].scatter(mag[valid], sharp[valid], s=5, alpha=0.4, color='darkorange')
    axes[1].axhline(0.0, color='red', lw=1.5, ls='--', label='SHARP=0 (star)')
    axes[1].set_xlabel("Instrumental Magnitude")
    axes[1].set_ylabel("SHARP")
    axes[1].set_title(f"{target_name} — SHARP Diagnostic")
    axes[1].legend()

    plt.tight_layout()
    if save:
        p = os.path.join(PLOTS_DIR, 'photometry', f"{target_name}_chi_sharp.png")
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"[Plot] CHI/SHARP → {p}")
    plt.close(fig)


def plot_calibration_diagnostics(image_path, zp_table, master_zp, master_zp_err, save=True):
    if zp_table is None or len(zp_table) == 0:
        return

    with fits.open(image_path) as hdul:
        clean_data = hdul[0].data.astype(float)

    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(clean_data)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 9))
    ax1.imshow(clean_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)

    # Categorise each star
    has_manual = 'manual_reject' in zp_table.colnames
    has_auto   = 'auto_reject'   in zp_table.colnames

    used_x, used_y, used_inst, used_ref, used_idx = [], [], [], [], []
    n_used = 0

    for i, row in enumerate(zp_table):
        x, y = row['xcenter'], row['ycenter']
        is_auto    = has_auto and row['auto_reject'] != ''
        is_clipped = bool(row['clipped'])
        is_manual  = has_manual and bool(row['manual_reject'])

        if is_auto:
            color, label_prefix = 'red', '✂'
        elif is_manual:
            color, label_prefix = 'magenta', '✗'
        elif is_clipped:
            color, label_prefix = 'orange', 'σ'
        else:
            color, label_prefix = 'lime', ''
            used_x.append(x); used_y.append(y)
            used_inst.append(row['m_inst']); used_ref.append(row['m_ref'])
            used_idx.append(i)
            n_used += 1

        ax1.add_patch(patches.Circle((x, y), radius=15,
                                      edgecolor=color, facecolor='none', lw=2))
        ax1.text(x+20, y+20, f"{label_prefix}{i+1}",
                 color=color, fontsize=11, fontweight='bold')

    # Legend patches for image panel
    from matplotlib.lines import Line2D
    legend_items = [
        Line2D([0], [0], marker='o', color='lime',   lw=0, markersize=10, label=f'Used ({n_used})'),
        Line2D([0], [0], marker='o', color='red',    lw=0, markersize=10, label='Auto-rejected'),
        Line2D([0], [0], marker='o', color='orange', lw=0, markersize=10, label='Sigma-clipped'),
    ]
    if has_manual:
        legend_items.append(
            Line2D([0], [0], marker='o', color='magenta', lw=0, markersize=10, label='Manual-rejected'))
    ax1.legend(handles=legend_items, loc='upper left', fontsize=10,
               facecolor='black', edgecolor='white', labelcolor='white')
    ax1.set_title(f"ZP Calibration Stars ({n_used} used / {len(zp_table)} matched)", fontsize=14)

    # Scatter plot: used stars
    if len(used_inst) > 0:
        x_vals = np.array([min(used_inst), max(used_inst)])
        ax2.scatter(used_inst, used_ref, color='blue', s=80, edgecolor='black',
                    zorder=5, label='Used')
        ax2.plot(x_vals, x_vals + master_zp, color='red', linestyle='--', lw=2,
                 label=f'ZP = {master_zp:.3f}')

    # Plot rejected stars as faded markers
    for i, row in enumerate(zp_table):
        is_auto    = has_auto and row['auto_reject'] != ''
        is_clipped = bool(row['clipped'])
        is_manual  = has_manual and bool(row['manual_reject'])
        if is_auto or is_clipped or is_manual:
            ax2.scatter(row['m_inst'], row['m_ref'], color='gray', s=50,
                        edgecolor='red', alpha=0.4, zorder=3, marker='x')

    ax2.set_xlabel(r"Instrumental Magnitude ($m_{inst}$)", fontsize=12)
    ax2.set_ylabel(f"{CALIBRATION_SYSTEM} True Magnitude", fontsize=12)
    ax2.set_title("Calibration Quality", fontsize=14)
    ax2.legend(fontsize=12)
    ax2.grid(True, linestyle=':', alpha=0.7)
    ax2.text(0.05, 0.95, f"ZP Uncertainty: ±{master_zp_err:.3f} mag",
             transform=ax2.transAxes, fontsize=12, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    if save:
        p = os.path.join(PLOTS_DIR, 'zero_point', f'{TARGET_NAME}_zp_calibration.png')
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"[Plot] Calibration Diagnostics → {p}")

        # --- Generate Second Plot (Used Only) ---
        fig_used, (ax1_u, ax2_u) = plt.subplots(1, 2, figsize=(20, 9))
        ax1_u.imshow(clean_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        
        for x, y, idx in zip(used_x, used_y, used_idx):
            ax1_u.add_patch(patches.Circle((x, y), radius=15, edgecolor='lime', facecolor='none', lw=2))
            ax1_u.text(x+20, y+20, f"{idx+1}", color='lime', fontsize=11, fontweight='bold')
            
        legend_items_u = [Line2D([0], [0], marker='o', color='lime', lw=0, markersize=10, label=f'Used ({n_used})')]
        ax1_u.legend(handles=legend_items_u, loc='upper left', fontsize=10, facecolor='black', edgecolor='white', labelcolor='white')
        ax1_u.set_title(f"ZP Calibration Stars ({n_used} used)", fontsize=14)
        
        if len(used_inst) > 0:
            x_vals = np.array([min(used_inst), max(used_inst)])
            ax2_u.scatter(used_inst, used_ref, color='blue', s=80, edgecolor='black', zorder=5, label='Used')
            ax2_u.plot(x_vals, x_vals + master_zp, color='red', linestyle='--', lw=2, label=f'ZP = {master_zp:.3f}')
            
        ax2_u.set_xlabel(r"Instrumental Magnitude ($m_{inst}$)", fontsize=12)
        ax2_u.set_ylabel(f"{CALIBRATION_SYSTEM} True Magnitude", fontsize=12)
        ax2_u.set_title("Calibration Quality (Used Only)", fontsize=14)
        ax2_u.legend(fontsize=12)
        ax2_u.grid(True, linestyle=':', alpha=0.7)
        ax2_u.text(0.05, 0.95, f"ZP Uncertainty: ±{master_zp_err:.3f} mag", transform=ax2_u.transAxes, fontsize=12, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        p_used = os.path.join(PLOTS_DIR, 'zero_point', f'{TARGET_NAME}_zp_calibration_used_only.png')
        fig_used.savefig(p_used, dpi=150, bbox_inches='tight')
        print(f"[Plot] Calibration Diagnostics (Used Only) → {p_used}")
        plt.close(fig_used)

    plt.close(fig)


def plot_target_confirmation(image_path, target_x, target_y, target_name, save=True):
    with fits.open(image_path) as hdul:
        clean_data = hdul[0].data.astype(float)
        
    zscale     = ZScaleInterval()
    vmin, vmax = zscale.get_limits(clean_data)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
    for ax in (ax1, ax2):
        ax.imshow(clean_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        ax.scatter(target_x, target_y, marker='+', color='red', s=2000, lw=2)
        
    ax1.set_title(f"Full Image: {target_name} Position", fontsize=14)
    ax1.add_patch(patches.Circle((target_x, target_y), radius=100,
                                  edgecolor='red', facecolor='none', lw=2, linestyle='--'))
    ax2.set_title(f"Zoom: {target_name} Site", fontsize=14)
    ax2.add_patch(patches.Circle((target_x, target_y), radius=6,
                                  edgecolor='red', facecolor='none', lw=1.5, linestyle='--'))
    zoom = 100
    ax2.set_xlim(target_x - zoom, target_x + zoom)
    ax2.set_ylim(target_y - zoom, target_y + zoom)
    ax2.grid(color='white', alpha=0.2, linestyle=':')
    plt.tight_layout()
    if save:
        p = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_confirmation.png')
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"[Plot] Target Confirmation → {p}")
    plt.close(fig)


def plot_source_detections(image_path, coo_file, target_name, save=True):
    if not os.path.exists(coo_file):
        return
        
    x_coords, y_coords = [], []
    with open(coo_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit(): continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    x_coords.append(float(parts[1]))
                    y_coords.append(float(parts[2]))
                except ValueError:
                    continue
                
    if len(x_coords) == 0:
        return
        
    with fits.open(image_path) as hdul:
        clean_data = hdul[0].data.astype(float)
        
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(clean_data)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(clean_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
    ax.scatter(x_coords, y_coords, marker='o', facecolors='none', edgecolors='cyan', s=50, lw=1)
    
    ax.set_title(f"Source Detections ({len(x_coords)} sources)", fontsize=14)
    ax.axis('off')
    plt.tight_layout()
    if save:
        p = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_source_detections.png')
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"[Plot] Source Detections → {p}")
    plt.close(fig)

def plot_all_stars_magnitudes(image_path, cal_table, target_name, save=True):
    if cal_table is None or len(cal_table) == 0:
        return
        
    with fits.open(image_path) as hdul:
        clean_data = hdul[0].data.astype(float)
        
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(clean_data)
    
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.imshow(clean_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
    
    valid = ~np.isnan(cal_table['CAL_MAG'])
    for row in cal_table[valid]:
        x = row['XCENTER']
        y = row['YCENTER']
        mag = row['CAL_MAG']
        
        ax.add_patch(patches.Circle((x, y), radius=10, edgecolor='yellow', facecolor='none', lw=1.5))
        ax.text(x + 12, y + 12, f"{mag:.2f}", color='yellow', fontsize=9, fontweight='bold')
        
    ax.set_title(f"All Stars Calibrated Magnitudes ({np.sum(valid)} stars)", fontsize=14)
    ax.axis('off')
    plt.tight_layout()
    if save:
        p = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_all_stars_mags.png')
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"[Plot] All Stars Magnitudes → {p}")
    plt.close(fig)


# ============================================================
# CELL 18.5 : COMPREHENSIVE PIPELINE LOG
# ============================================================

def write_comprehensive_log(final_results, pipeline_mode, image_path=STACKED_IMAGE, target_ra=TARGET_RA, target_dec=TARGET_DEC):
    import datetime

    log_path = os.path.join(PLOTS_DIR, f"{TARGET_NAME}_pipeline.log")
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    fwhm    = final_results.get('fwhm', np.nan)
    sky     = final_results.get('sky', np.nan)
    std_sky = final_results.get('std_sky', np.nan)
    als_table = final_results.get('als_table')
    num_sources = len(als_table) if als_table is not None else 0

    zp      = final_results.get('zp', np.nan)
    zp_err  = final_results.get('zp_err', np.nan)
    zp_table = final_results.get('zp_table')

    detected = final_results.get('detected', False)
    cal_mag  = final_results.get('cal_mag', np.nan)
    cal_mag_err = final_results.get('cal_mag_err', np.nan)
    inst_mag = final_results.get('inst_mag', np.nan)
    inst_merr = final_results.get('inst_merr', np.nan)
    chi = final_results.get('chi', np.nan)
    sharp = final_results.get('sharp', np.nan)

    detect_type = final_results.get('detect_type', 'None')
    sigma_th = final_results.get('sigma_th', np.nan)
    ap_raw = final_results.get('ap_raw')
    als_raw = final_results.get('als_raw')

    with open(log_path, 'a') as f:
        f.write("==================================================\n")
        f.write(f" PIPELINE RUN: {TARGET_NAME} @ {timestamp}\n")
        f.write("==================================================\n")
        f.write(f"Run Mode    : DAOPHOT II ({pipeline_mode.upper()})\n")
        f.write(f"Input file  : {image_path}\n")
        f.write(f"Calibration : {CALIBRATION_SYSTEM} ({VIZIER_CATALOG_ID})\n")
        
        with fits.open(image_path, ignore_missing_simple=True) as hdul:
            hdul.verify('fix')
            orig_header = hdul[0].header
        date_obs_str = orig_header.get('DATE-OBS', 'UNKNOWN')
        t_obs_jd = np.nan
        delta_t_days = np.nan
        if date_obs_str != 'UNKNOWN':
            date_obs_clean = date_obs_str.replace('"T"', 'T')
            try:
                t_obs = Time(date_obs_clean, format='isot', scale='utc')
                t_obs_jd = t_obs.jd
                delta_t_days = t_obs_jd - GRB_T0_JD
            except Exception as e:
                print(f"⚠️ Could not parse DATE-OBS: {e}")
                pass

        f.write(f"DATE-OBS    : {date_obs_str}\n")
        f.write(f"Obs JD      : {t_obs_jd:.4f}\n")
        f.write(f"T - T_0     : {delta_t_days:.4f} days\n")
        f.write(f"CCD Gain    : {GAIN:.2f} e⁻/ADU\n")
        f.write(f"CCD ReadN   : {READ_NOISE:.2f} e⁻\n\n")

        with fits.open(image_path) as hdul:
            wcs = get_wcs(hdul[0].header, ASTROMETRY_REF_IMAGE)
        pixel_scale_degrees = proj_plane_pixel_scales(wcs)[0]
        pixel_scale_arcsec = pixel_scale_degrees * 3600
        f.write(f"Pixel Scale : {pixel_scale_arcsec:.3f} arcsec/pix\n")
        f.write(f"GRB Target  : RA={target_ra}, Dec={target_dec}\n")
        
        pixel_x = final_results.get('pixel_x', np.nan)
        pixel_y = final_results.get('pixel_y', np.nan)
        if not np.isnan(pixel_x):
            recentered_ra, recentered_dec = wcs.pixel_to_world_values(pixel_x, pixel_y)
            f.write(f"Recentered  : RA={float(recentered_ra):.5f}, Dec={float(recentered_dec):.5f} (X={pixel_x:.1f}, Y={pixel_y:.1f})\n")
        f.write("\n")

        f.write("--- IMAGE QUALITY SUMMARY ---\n")
        f.write(f"FWHM        : {fwhm:.2f} px\n")
        f.write(f"Sky mean    : {sky:.2f}\n")
        f.write(f"Sky std     : {std_sky:.2f}\n")
        f.write(f"Num sources : {num_sources}\n\n")

        f.write("--- ZERO POINT CALIBRATION ---\n")
        f.write(f"Zero Point  : {zp:.3f} ± {zp_err:.3f}\n")
        
        if zp_table is not None:
            has_auto = 'auto_reject' in zp_table.colnames
            has_manual = 'manual_reject' in zp_table.colnames
            has_chi = 'chi' in zp_table.colnames

            f.write("\nCalibration Stars (All Matched):\n")
            f.write(f"{'#':<4} | {'m_inst':<10} | {'m_ref(cat)':<12} | "
                    f"{'ZP_star':<10} | {'ΔZP':<8} |")
            if has_chi:
                f.write(f" {'CHI':<6} | {'SHARP':<7} | {'MERR':<6} |")
            f.write(f" {'Status'}\n")
            f.write("-" * 105 + "\n")
            
            for i, row in enumerate(zp_table):
                zp_star = row['zp_star']
                delta_zp = zp_star - zp  # per-star ZP minus master ZP
                
                status = 'USED'
                if has_auto and row['auto_reject']:
                    status = f"AUTO-REJ ({row['auto_reject']})"
                elif has_manual and row['manual_reject']:
                    status = 'MANUAL-REJ'
                elif row['clipped']:
                    status = 'SIGMA-CLIP'
                
                f.write(f"{i+1:<4} | {row['m_inst']:<10.3f} | {row['m_ref']:<12.3f} | "
                        f"{zp_star:<10.4f} | {delta_zp:+.4f} |")
                if has_chi:
                    f.write(f" {row['chi']:<6.2f} | {row['sharp']:<+7.3f} | {row['merr']:<6.3f} |")
                f.write(f" {status}\n")

            f.write("\nCalibration Stars (ONLY USED):\n")
            f.write(f"{'#':<4} | {'m_inst':<10} | {'m_ref(cat)':<12} | "
                    f"{'ZP_star':<10} | {'ΔZP':<8} |")
            if has_chi:
                f.write(f" {'CHI':<6} | {'SHARP':<7} | {'MERR':<6} |")
            f.write(f" {'Status'}\n")
            f.write("-" * 105 + "\n")
            
            for i, row in enumerate(zp_table):
                zp_star = row['zp_star']
                delta_zp = zp_star - zp  # per-star ZP minus master ZP
                
                status = 'USED'
                if has_auto and row['auto_reject']: continue
                elif has_manual and row['manual_reject']: continue
                elif row['clipped']: continue
                
                f.write(f"{i+1:<4} | {row['m_inst']:<10.3f} | {row['m_ref']:<12.3f} | "
                        f"{zp_star:<10.4f} | {delta_zp:+.4f} |")
                if has_chi:
                    f.write(f" {row['chi']:<6.2f} | {row['sharp']:<+7.3f} | {row['merr']:<6.3f} |")
                f.write(f" {status}\n")
            
            f.write(f"\n  Master ZP = {zp:.4f} ± {zp_err:.4f}\n")
            f.write(f"  ΔZP = per-star ZP − Master ZP (scatter indicates calibration quality)\n")
        f.write("\n")

        f.write("--- PHOTOMETRY RESULTS ---\n")
        if target_ra is not None and target_dec is not None:
            if detected and not np.isnan(cal_mag):
                f.write("Result      : DETECTED\n")
                f.write(f"Detect Type : {detect_type}\n")
                f.write(f"Sigma Thresh: {sigma_th}\n")
                f.write(f"FINAL MAG   : {cal_mag:.3f} ± {cal_mag_err:.3f} mag ({FILTER_BAND}-band)\n")
                f.write(f"Instrumental: {inst_mag:.3f} ± {inst_merr:.3f}\n")
                f.write(f"CHI / SHARP : {chi:.2f} / {sharp:.3f}\n")
                if ap_raw:
                    f.write("\n  --- Raw .ap row data ---\n")
                    f.write(f"{ap_raw}\n")
                if als_raw:
                    f.write("\n  --- Raw .als row data ---\n")
                    f.write(f"{als_raw}\n")
            else:
                f.write("Result      : NON-DETECTION (Target not recovered by DAOPHOT)\n")
        else:
            f.write(f"Result      : ALL STARS PROCESSED\n")
            f.write(f"See {TARGET_NAME}_filtered_stars_mags.png for visual magnitudes.\n")
            cal_table_all = final_results.get('all_sources_calibrated')
            cal_table_filt = final_results.get('filtered_stars_calibrated')
            if cal_table_all is not None:
                valid_all = ~np.isnan(cal_table_all['CAL_MAG'])
                f.write(f"Total sources (unfiltered): {np.sum(valid_all)}\n")
            if cal_table_filt is not None:
                valid_filt = ~np.isnan(cal_table_filt['CAL_MAG'])
                f.write(f"Total clean stars (filtered): {np.sum(valid_filt)}\n")
        f.write("==================================================\n")

    print(f"\n[Save] Comprehensive Log → {os.path.abspath(log_path)}")


# ============================================================
# CELL 19 : MAIN ORCHESTRATOR
# ============================================================


def extract_raw_target_data(file_path, target_id):
    """
    Extracts header line, column names, and raw data lines for a target_id 
    from a DAOPHOT output file (.ap or .als).
    """
    import os
    if not file_path or not os.path.exists(file_path):
        return None
        
    with open(file_path, "r") as f:
        lines = f.readlines()
        
    if len(lines) < 3:
        return None
        
    header_str = lines[0].rstrip() + "\n" + lines[1].rstrip()
    
    try:
        nl = int(lines[1].split()[0])
    except:
        nl = 1
    
    data_str = ""
    for i in range(2, len(lines)):
        tokens = lines[i].split()
        if tokens:
            try:
                line_id = int(tokens[0])
                if line_id == target_id:
                    for j in range(nl):
                        if i + j < len(lines):
                            data_str += lines[i+j].rstrip() + "\n"
                    break
            except ValueError:
                pass
                
    if data_str:
        return f"{header_str}\n{data_str.rstrip()}"
    return None

def prepare_manual_workspace(image_path, fwhm_override=None, sky_override=None, std_sky_override=None, pipeline_mode='semi-auto'):
    """
    Prepare the workspace for manual DAOPHOT terminal sessions.
    This creates directories, copies the image, and updates OPT files based on FWHM.
    """
    print("\n" + "╔" + "═"*58 + "╗")
    print(f"║{'PREPARING MANUAL DAOPHOT WORKSPACE':^58}║")
    print("╚" + "═"*58 + "╝")
    
    setup_dirs()
    copy_opt_files_to_workdir('.')
    image_work = copy_image_to_workdir(image_path)
    stem       = get_stem(image_work)
    
    clean_previous_run(stem, pipeline_mode='manual')
    
    if fwhm_override is not None:
        fwhm    = fwhm_override
        sky     = sky_override     if sky_override     is not None else 0.0
        std_sky = std_sky_override if std_sky_override is not None else 1.0
        print(f"\n  D0: Using supplied FWHM={fwhm:.2f}  sky={sky:.2f}  std={std_sky:.2f}")
    else:
        fwhm, sky, std_sky, _, _ = measure_fwhm_and_sky(image_work)
        
    update_daophot_opt(fwhm)
    update_allstar_opt(fwhm)
    update_photo_opt(fwhm)
    
    if pipeline_mode in ['auto', 'semi-auto']:
        print("\n  [Manual Mode] Automating D1 (FIND), D2 (PHOT), and D3 (PICK)...")
        image_bn = os.path.basename(image_work)
        run_find(image_bn, fwhm)
        run_phot(image_bn)
        run_pick(image_bn, fwhm)
        
        print(f"\n  Workspace ready. You can now run DAOPHOT interactively in the terminal.")
        print(f"  Note: FIND, PHOT, and PICK are already done. You only need to run 'AT', 'PSF', and 'EX'.")
    else:
        print(f"\n  Workspace ready. You can now run DAOPHOT interactively in the terminal.")
        print(f"  Note: FULL MANUAL MODE. You need to run 'AT', 'FI', 'PH', 'PI', 'PSF', and 'EX'.")
    return fwhm

def run_daophot_pipeline(
    image_path       = STACKED_IMAGE,
    target_ra        = TARGET_RA,
    target_dec       = TARGET_DEC,
    target_name      = TARGET_NAME,
    filter_band      = FILTER_BAND,
    force_coords     = FORCE_COORDS,
    fwhm_override    = None,
    sky_override     = None,
    std_sky_override = None,
    pipeline_mode    = PIPELINE_MODE,
):
    """
    Run the full standalone DAOPHOT II pipeline.

    Parameters
    ----------
    image_path       : WCS-solved stacked FITS (absolute or relative path)
    target_ra/dec    : target celestial coordinates (degrees)
    target_name      : used in filenames and plot titles
    filter_band      : photometric band  ('r', 'V', 'B', …)
    force_coords     : (ra, dec) tuple or None — enable forced photometry
    fwhm_override    : pass your photutil FWHM to skip D0
    sky_override     : pass your photutil sky mode to skip D0
    std_sky_override : pass your photutil sky std to skip D0

    Returns
    -------
    dict with: detected, cal_mag, cal_mag_err, zp, zp_err,
               als_table, target_info, fwhm, sky, std_sky
    """
    print("\n" + "╔" + "═"*58 + "╗")
    print(f"║{'DAOPHOT II  (standalone)  PIPELINE':^58}║")
    print(f"║  Target  : {target_name:<46}║")
    print(f"║  Image   : {os.path.basename(image_path):<46}║")
    print(f"║  Filter  : {filter_band:<46}║")
    print("╚" + "═"*58 + "╝")

    # ── Setup ──────────────────────────────────────────────────────────
    setup_dirs()
    copy_opt_files_to_workdir('.')
    image_work = copy_image_to_workdir(image_path)
    stem       = get_stem(image_work)
    image_bn   = os.path.basename(image_work)   # basename used by DAOPHOT

    clean_previous_run(stem, pipeline_mode=pipeline_mode)

    with fits.open(image_path) as hdul:
        header      = hdul[0].header
        image_shape = hdul[0].data.shape
    wcs = get_wcs(header, ASTROMETRY_REF_IMAGE)

    # ── D0 : FWHM & sky ────────────────────────────────────────────────
    if fwhm_override is not None:
        fwhm    = fwhm_override
        sky     = sky_override     if sky_override     is not None else 0.0
        std_sky = std_sky_override if std_sky_override is not None else 1.0
        print(f"\n  D0: Using supplied FWHM={fwhm:.2f}  "
              f"sky={sky:.2f}  std={std_sky:.2f}")
    else:
        fwhm, sky, std_sky, _, _ = measure_fwhm_and_sky(image_work)

    # Update DAOPHOT.OPT, ALLSTAR.OPT, and PHOTO.OPT with the measured FWHM
    update_daophot_opt(fwhm)
    update_allstar_opt(fwhm)
    update_photo_opt(fwhm)

    if pipeline_mode == 'auto':
        # ── D1 : FIND ──────────────────────────────────────────────────────
        coo_file, n_sources = run_find(image_bn, fwhm)
    
        # ── D2 : PHOT ──────────────────────────────────────────────────────
        ap_file = run_phot(image_bn)
    
        # ── D3 : PICK ──────────────────────────────────────────────────────
        lst_file, n_psf = run_pick(image_bn, fwhm)
    
        # ── D4 : PSF (iterative) ───────────────────────────────────────────
        psf_file = run_psf_iterative(image_bn)
    else:
        print("\n  [Manual Mode] Skipping D1-D4. Assuming PSF is manually generated.")
        psf_file = workpath(f"{stem}.psf")
        if not os.path.exists(psf_file):
            raise RuntimeError(f"[Manual Mode] Expected PSF file at {psf_file} but none was found.")

    # ── D5 : ALLSTAR (full frame) ──────────────────────────────────────
    als_file, sub_image, als_table = run_allstar_full(image_bn)

    # ── D7 : Target lookup (Initial) ───────────────────────────────────
    target_info = None
    if target_ra is not None and target_dec is not None:
        target_info = lookup_target_in_als(
            als_table, wcs, target_ra, target_dec
        )

    # ── D6 : Forced photometry (Fallback) ─────────────────────────────
    final_als_table = als_table
    is_forced = False

    if target_ra is not None and target_dec is not None:
        if not target_info['detected'] or force_coords is not None:
            is_forced = True
            if not target_info["detected"] and force_coords is None:
                print("\n  ⚠️ Target not initially detected. Falling back to forced photometry...")
                frc_ra, frc_dec = target_ra, target_dec
            else:
                frc_ra, frc_dec = force_coords

            merged_als, merged_table, forced_x, forced_y = run_forced_photometry(
                image_bn, als_table, fwhm, wcs, frc_ra, frc_dec
            )
            if merged_table is not None:
                final_als_table = merged_table
                print(f"  Using merged table ({len(merged_table)} stars)")
                
                # ── D7 : Target lookup (Post-Forced) ───────────────────────
                target_info = lookup_target_in_als(
                    final_als_table, wcs, target_ra, target_dec
                )

    # ── D8 : Zero-point calibration ────────────────────────────────────
    zp, zp_err, zp_table = calibrate_zero_point_daophot(
        final_als_table, wcs, image_shape, filter_band,
        pipeline_mode=pipeline_mode, image_path=image_path
    )

    # ── D9 : Apply ZP ──────────────────────────────────────────────────
    final = {}
    if target_info is not None:
        final = compute_final_magnitude(target_info, zp, zp_err)
    final['detected'] = final.get('detected', False)

    all_sources_calibrated, filtered_stars_calibrated = compute_all_magnitudes(final_als_table, zp, zp_err, wcs)

    # ── Diagnostic plots ───────────────────────────────────────────────
    plot_chi_sharp(final_als_table, target_name)
    plot_calibration_diagnostics(image_path, zp_table, zp, zp_err)

    if zp_table is not None:
        reg_zp = os.path.join(PLOTS_DIR, 'zero_point', f'{target_name}_zp_calibration_stars.reg')
        try:
            with open(reg_zp, 'w') as f:
                f.write('global dashlist=8 3 width=2 font="helvetica 10 normal roman"\n')
                f.write('image\n')
                for i, row in enumerate(zp_table):
                    x, y = row['xcenter'], row['ycenter']
                    is_auto = row['auto_reject'] != ''
                    is_clip = bool(row['clipped'])
                    is_manual = 'manual_reject' in zp_table.colnames and bool(row['manual_reject'])
                    
                    if not (is_auto or is_clip or is_manual):
                        # Write used stars with a lime circle
                        f.write(f'circle({x:.2f},{y:.2f},15) # color=lime text="{{{i+1}}}"\n')
            print(f"  [Output] Calibration Stars Regions → {reg_zp}")
        except Exception as e:
            print(f"  ⚠️ Could not save ZP regions: {e}")

    if target_ra is not None and target_dec is not None:
        px = final.get('pixel_x', np.nan)
        py = final.get('pixel_y', np.nan)
        if np.isnan(px):
            px, py = wcs.wcs_world2pix(target_ra, target_dec, 0)
            px, py = float(px), float(py)
        plot_target_confirmation(image_path, px, py, target_name)
    
    coo_file = workpath(f"{get_stem(image_bn)}.coo")
    plot_source_detections(image_path, coo_file, target_name)
    plot_all_stars_magnitudes(image_path, filtered_stars_calibrated, target_name + "_filtered")

    # Extract raw rows if detected
    ap_raw = None
    als_raw = None
    sigma_th = 3.0  # default DAOPHOT sigma threshold
    opt_path = workpath("daophot.opt")
    if os.path.exists(opt_path):
        with open(opt_path, "r") as f:
            for ln in f:
                if "TH" in ln and "=" in ln:
                    try:
                        sigma_th = float(ln.split("=")[1].strip())
                    except:
                        pass

    if final.get('detected'):
        target_id = final.get("id", -1)
        # Determine the AP file used
        ap_f = workpath(f"{get_stem(image_bn)}_merged.ap") if is_forced else workpath(f"{get_stem(image_bn)}.ap")
        als_f = workpath(f"{get_stem(image_bn)}_merged.als") if is_forced else workpath(f"{get_stem(image_bn)}.als")
        
        ap_raw = extract_raw_target_data(ap_f, target_id)
        als_raw = extract_raw_target_data(als_f, target_id)

    # ── Summary ────────────────────────────────────────────────────────
    print("\n" + "="*60)
    print("  RESULT")
    print("="*60)
    
    if target_ra is not None and target_dec is not None:
        print(f"  Target      : {target_name}")
        print(f"  Detected    : {final['detected']}")
        
        if final['detected']:
            detection_type = "Forced Photometry" if is_forced else "Normal Intrinsic Find"
            print(f"  Detect Type : {detection_type}")
            print(f"  Sigma Thresh: {sigma_th}")

        if final['detected'] and not np.isnan(final.get('cal_mag', np.nan)):
            print(f"  Calibrated  : {final['cal_mag']:.3f} ± {final['cal_mag_err']:.3f}  ({filter_band}-band)")
            print(f"  Instrumental: {final['mag']:.3f} ± {final['merr']:.3f}")
            print(f"  Zero-point  : {zp:.4f} ± {zp_err:.4f}")
            print(f"  CHI / SHARP : {final['chi']:.2f} / {final['sharp']:.3f}")
            print(f"  Pixel pos   : X={final['pixel_x']:.2f}  Y={final['pixel_y']:.2f}  (dist={final['dist_pix']:.2f}px)")
                  
            if ap_raw:
                print("\n  --- Raw .ap row data ---")
                print(ap_raw)
            if als_raw:
                print("\n  --- Raw .als row data ---")
                print(als_raw)
        else:
            print("  ⚠️  Target not recovered by DAOPHOT.")
            print("     → Set  FORCE_COORDS = (ra, dec)  to inject position.")
            print(f"  Zero-point  : {zp:.4f} ± {zp_err:.4f}")
    else:
        print(f"  Processed all stars without specific target matching.")
        print(f"  Zero-point  : {zp:.4f} ± {zp_err:.4f}")
        
        if all_sources_calibrated is not None and filtered_stars_calibrated is not None:
            valid_all = ~np.isnan(all_sources_calibrated['CAL_MAG'])
            valid_filt = ~np.isnan(filtered_stars_calibrated['CAL_MAG'])
            print(f"  Total sources (unfiltered): {np.sum(valid_all)}")
            print(f"  Total clean stars (filtered): {np.sum(valid_filt)}")
            print(f"  See {target_name}_filtered_all_stars_mags.png for visual magnitudes.")
            
            # Save Unfiltered
            csv_all = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_all_sources_mags.csv')
            reg_all = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_all_sources_mags.reg')
            try:
                all_sources_calibrated.write(csv_all, format='csv', overwrite=True)
                with open(reg_all, 'w') as f:
                    f.write('global color=cyan dashlist=8 3 width=1 font="helvetica 10 normal roman"\n')
                    f.write('image\n')
                    for row in all_sources_calibrated[valid_all]:
                        f.write(f'circle({row["XCENTER"]:.2f},{row["YCENTER"]:.2f},10) # text="{{{row["ID"]} m:{row["CAL_MAG"]:.2f}±{row["CAL_MAG_ERR"]:.2f}}}"\n')
            except Exception as e: print(f"  ⚠️ Could not save unfiltered catalogs: {e}")

            # Save Filtered
            csv_filt = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_filtered_stars_mags.csv')
            reg_filt = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_filtered_stars_mags.reg')
            try:
                filtered_stars_calibrated.write(csv_filt, format='csv', overwrite=True)
                with open(reg_filt, 'w') as f:
                    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman"\n')
                    f.write('image\n')
                    for row in filtered_stars_calibrated[valid_filt]:
                        f.write(f'circle({row["XCENTER"]:.2f},{row["YCENTER"]:.2f},10) # text="{{{row["ID"]} m:{row["CAL_MAG"]:.2f}±{row["CAL_MAG_ERR"]:.2f}}}"\n')
                print(f"  Saved outputs to {PLOTS_DIR}/")
            except Exception as e: print(f"  ⚠️ Could not save filtered catalogs: {e}")
                
    print("="*60)

    final_output_dict = {
        'detected'   : final.get('detected', False),
        'cal_mag'    : final.get('cal_mag', np.nan),
        'cal_mag_err': final.get('cal_mag_err', np.nan),
        'inst_mag'   : final.get('mag', np.nan),
        'inst_merr'  : final.get('merr', np.nan),
        'zp'         : zp,
        'zp_err'     : zp_err,
        'chi'        : final.get('chi', np.nan),
        'sharp'      : final.get('sharp', np.nan),
        'pixel_x'    : final.get('pixel_x', np.nan),
        'pixel_y'    : final.get('pixel_y', np.nan),
        'id'         : final.get('id', -1),
        'dist_pix'   : final.get('dist_pix', np.nan),
        'fwhm'       : fwhm,
        'sky'        : sky,
        'std_sky'    : std_sky,
        'als_table'  : final_als_table if 'final_als_table' in locals() else None,
        'zp_table'   : zp_table if 'zp_table' in locals() else None,
        'all_sources_calibrated': all_sources_calibrated,
        'filtered_stars_calibrated': filtered_stars_calibrated,
        'detect_type': "Forced Photometry" if is_forced else "Normal Intrinsic Find",
        'sigma_th'   : sigma_th,
        'ap_raw'     : ap_raw,
        'als_raw'    : als_raw,
    }

    write_comprehensive_log(final_output_dict, pipeline_mode, image_path, target_ra, target_dec)

    return final_output_dict


# ============================================================
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# JUPYTER NOTEBOOK COPY-PASTE TEMPLATE
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The following blocks are formatted with standard `# %%` markers. 
# If you are using VS Code or Spyder, you can run them directly. 
# Otherwise, copy and paste each block into a separate Jupyter cell.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %% [markdown]
# # DAOPHOT II Pipeline (Automated / Manual)
# Follow these cells to seamlessly integrate photutils measurements,
# interactive terminal DAOPHOT sessions, and the automated PyRAF/subprocess pipeline.

# %%
# =====================================================================
# CELL 1: MANUAL WORKSPACE CONFIGURATION
# =====================================================================
manual_daophot_done = False   # Set to True after you finish the terminal steps

# Provide your photutils measurements here (or leave None for auto-measurement)
FWHM_VAL = None  
SKY_VAL  = None  
STD_VAL  = None  

# %%
# =====================================================================
# CELL 2: PREPARE WORKSPACE FOR MANUAL TERMINAL MODE
# =====================================================================
if __name__ == '__main__':
    if PIPELINE_MODE in ['manual', 'semi-auto']:
        if not manual_daophot_done:
            measured_fwhm = prepare_manual_workspace(
                image_path       = STACKED_IMAGE,
                fwhm_override    = FWHM_VAL,
                sky_override     = SKY_VAL,
                std_sky_override = STD_VAL,
                pipeline_mode    = PIPELINE_MODE
            )
        else:
            print("Manual workspace preparation skipped because manual_daophot_done is True.")
    else:
        print("Manual workspace preparation skipped because PIPELINE_MODE is 'auto'.")

# %% [markdown]
# ### 🛑 STOP: Open Your Terminal Now! (If PIPELINE_MODE in ['manual', 'semi-auto'] and manual_daophot_done is False)
# 1. `distrobox enter astro-box`
# 2. `cd <your_work_directory>/`
# 3. `~/dao2/daophot`
# 4. Interactive commands: `AT <image_name_without_extension>` -> `PSF` -> `EX`
#
# *Once you have created the `.psf` file and exited DAOPHOT, change `manual_daophot_done = True` in CELL 1 and re-run.*

# %%
if __name__ == '__main__':
    if PIPELINE_MODE in ['manual', 'semi-auto']:
        if manual_daophot_done:
            final_results = run_daophot_pipeline(
                image_path    = STACKED_IMAGE,
                target_ra     = TARGET_RA,
                target_dec    = TARGET_DEC,
                target_name   = TARGET_NAME,
                filter_band   = FILTER_BAND,
                pipeline_mode = PIPELINE_MODE
            )
            print("\nFinal Calibrated Magnitude:", final_results.get('cal_mag'))
        else:
            print("Waiting for manual DAOPHOT completion. Set manual_daophot_done = True when finished.")
    else:
        print("Skipping CELL 3 because PIPELINE_MODE is 'auto'.")

# %%
# =====================================================================
# FULLY AUTOMATED RUN
# =====================================================================
# If your frame is sparse and you trust DAOPHOT to build the PSF silently
# (ensure WA = -2 in daophot.opt), this will run everything.

if __name__ == '__main__':
    if PIPELINE_MODE == 'auto':
        final_results = run_daophot_pipeline(
            image_path    = STACKED_IMAGE,
            target_ra     = TARGET_RA,
            target_dec    = TARGET_DEC,
            target_name   = TARGET_NAME,
            filter_band   = FILTER_BAND,
            pipeline_mode = 'auto'
        )
        print("\nFinal Calibrated Magnitude:", final_results.get('cal_mag'))
    else:
        print("Skipping FULLY AUTOMATED RUN because PIPELINE_MODE is 'manual' or 'semi-auto'.")
