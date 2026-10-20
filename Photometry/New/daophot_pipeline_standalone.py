# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.6
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
TARGET_NAME   = "GRB261011A"
TARGET_RA     =  "01:27:35.39"               # Degrees float (e.g. 248.74775) OR sexagesimal string (e.g. "16:34:59.46" or "16 34 59.46")
TARGET_DEC    =  "-21:33:37.2"               # Degrees float (e.g. +26.58261) OR sexagesimal string (e.g. "+26:34:57.4" or "+26 34 57.4")
GRB_T0_JD     =  "2026-10-11 13:49:04"   # "YYYY-MM-DD HH:MM:SS", "19:17:03.16 UT on 10 October 2026", or JD float (2461324.303509)

# --- Input image (WCS-solved, stacked FITS) ---
STACKED_IMAGE = "R.fits"

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
FRAMES_AVG    =  12       # N frames that were averaged together
FRAMES_SUM    =  1       # N frames that were summed together
               # entry format sent to FIND: "FRAMES_AVG,FRAMES_SUM"
               # examples: single exp → "1,1"; avg of 5 → "5,1"; sum of 3 → "1,3"

# --- Detector ---
READ_NOISE    =  7.0     # e⁻ / ADU (one raw frame)
GAIN          =  2.0     # photons / ADU (one raw frame)
SATURATION    = 65535    # ADU — detector linearity limit

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
import re
import dateutil.parser
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
from scipy.ndimage import maximum_filter
import io
import urllib.request
from PIL import Image
# %%
# CELL 3 : Functions
# ============================================================
# ============================================================

def parse_target_coords(ra, dec):
    """
    Parse RA and Dec in either decimal degrees or sexagesimal formats:
    - RA: float degrees (e.g. 248.74775) or string (e.g. "16:34:59.46", "16 34 59.46", "16h34m59.46s")
    - Dec: float degrees (e.g. 26.58261) or string (e.g. "+26:34:57.4", "+26 34 57.4", "+26d34m57.4s")
    Returns (ra_deg, dec_deg) as floats.
    """
    if ra is None or dec is None:
        return None, None
    if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
        return float(ra), float(dec)

    ra_str = str(ra).strip().replace("h", ":").replace("m", ":").replace("s", "")
    dec_str = str(dec).strip().replace("d", ":").replace("m", ":").replace("s", "")
    try:
        if ":" in ra_str or " " in ra_str:
            c = SkyCoord(ra_str, dec_str, unit=(u.hourangle, u.deg))
        else:
            c = SkyCoord(float(ra_str), float(dec_str), unit=(u.deg, u.deg))
        return float(c.ra.deg), float(c.dec.deg)
    except Exception as e:
        print(f"  ⚠️ Could not parse coordinates with SkyCoord ({e}). Fallback to float conversion.")
        return float(ra), float(dec)


def format_hms_dms(ra, dec):
    """
    Format RA and Dec into standard sexagesimal strings:
    RA: HH:MM:SS.ss
    Dec: +DD:MM:SS.s (or -DD:MM:SS.s)
    Returns (ra_hms_str, dec_dms_str).
    """
    if ra is None or dec is None:
        return "N/A", "N/A"
    try:
        ra_deg, dec_deg = parse_target_coords(ra, dec)
        if ra_deg is None or dec_deg is None or np.isnan(ra_deg) or np.isnan(dec_deg):
            return "N/A", "N/A"
        c = SkyCoord(ra=float(ra_deg) * u.deg, dec=float(dec_deg) * u.deg)
        ra_hms = c.ra.to_string(unit=u.hour, sep=':', precision=2, pad=True)
        dec_dms = c.dec.to_string(unit=u.deg, sep=':', precision=1, pad=True, alwayssign=True)
        return ra_hms, dec_dms
    except Exception:
        return "N/A", "N/A"


# ====================================================================================================
# CELL 3B : ROBUST APERTURE PHOTOMETRY & PHOTOMETRIC UPPER LIMIT ENGINE
# ====================================================================================================
#
# ── THEORY & MATHEMATICAL FORMULATION OF PHOTOMETRIC DETECTION & UPPER LIMITS ──────────────────────
#
# 1. LIMITATION OF THE STANDARD CCD EQUATION IN REDUCED / STACKED DATA:
#    In idealized raw CCD frames, photon noise is Poissonian, governed by the classic CCD equation:
#        σ_theo² = (F_net / g) + N_ap * [ (B / g) + (σ_RN / g)² ]
#    where g is gain (e-/ADU), B is local sky background (ADU/px), and σ_RN is read noise (e-).
#
#    HOWEVER, in real astronomical science frames (co-added stacks, aligned images, flat-fielded frames):
#      a. Spatial resampling (sub-pixel shift, drizzling) introduces pixel-to-pixel covariance.
#      b. Residual flat-field non-uniformities, fringing, and sky gradients create spatial background
#         fluctuations that far exceed pure Poisson shot noise: σ_sky,emp >> sqrt(B / g).
#      c. Applying the raw Poisson formula underestimates background noise by factors of 3× to 10×,
#         which artificially inflates the signal-to-noise ratio (S/N) and produces spurious "detections"
#         out of random background ripples.
#
# 2. EMPIRICAL BACKGROUND VARIANCE & FLUX ERROR FORMULATION:
#    We measure the true empirical spatial variance σ_sky,emp directly from the local background annulus
#    using iterative sigma-clipping (rejecting stellar wings and cosmic rays).
#    The total empirical variance in an aperture of area N_ap with annulus area N_ann is:
#        σ_flux² = (max(0, F_net) / g) + σ_bkg,ap²
#    where the background-induced aperture flux uncertainty is:
#        σ_bkg,ap = σ_sky,emp * sqrt( N_ap * (1.0 + N_ap / N_ann) )
#    The term (1.0 + N_ap / N_ann) rigorously accounts for the statistical uncertainty in subtracting
#    the local sky mean derived from a finite annulus.
#
# 3. 4-GATE DETECTION VERIFICATION FRAMEWORK:
#    To prevent false detections on diffuse noise fluctuations or centroid jumps, a candidate must
#    pass all 4 gates to be classified as 'DETECTED':
#      Gate 1 (Peak Significance) : Central 5×5 peak pixel must exceed local background by >= 3.5σ:
#                                   Peak_sig = (max(stamp_5x5) - B_med) / σ_sky,emp >= 3.5
#      Gate 2 (Empirical S/N)     : Net integrated flux must satisfy S/N_emp = F_net / σ_flux >= 3.0
#      Gate 3 (Positive Flux)     : Net aperture flux F_net > 0
#      Gate 4 (Valid Calibration) : Photometric zero-point ZP must be finite.
#    If ANY gate fails, the measurement is declared a NON-DETECTION, and the pipeline calculates
#    rigorous photometric limiting magnitudes (Upper Limits).
#
# 4. MATHEMATICAL FORMULATION OF PHOTOMETRIC UPPER LIMITS:
#    In the absence of a statistically significant point source, we determine the maximum flux
#    that could remain undetected at the specified confidence level:
#        F_3σ = 3.0 * σ_bkg,ap   (3-sigma limiting flux, ~99.73% Gaussian confidence)
#        F_5σ = 5.0 * σ_bkg,ap   (5-sigma limiting flux, ~99.99994% Gaussian confidence)
#
#    Converting to calibrated instrumental and standard magnitudes:
#        m_inst,3σ = 25.0 - 2.5 * log10( F_3σ )
#        m_cal,3σ  = m_inst,3σ + ZP   ==>   Reported as:  Mag > m_cal,3σ  (3σ Upper Limit)
#
#        m_inst,5σ = 25.0 - 2.5 * log10( F_5σ )
#        m_cal,5σ  = m_inst,5σ + ZP   ==>   Reported as:  Mag > m_cal,5σ  (5σ Upper Limit)
# ====================================================================================================

def measure_aperture_photometry_at_pix(data, x, y, fwhm, zp, zp_err,
                                       gain=GAIN, readnoise=READ_NOISE,
                                       r_ap_mult=PHOTO_A1_MULT,
                                       min_peak_sigma=3.5,
                                       min_snr=3.0):
    """
    Perform rigorous circular aperture photometry directly at pixel (x, y) on the FITS array.
    Uses aperture radius r_ap = max(2.5, fwhm * r_ap_mult) and local background annulus.

    Features:
      • Directly samples local annulus pixels to derive TRUE empirical sky variance (σ_sky,emp),
        preventing inflated S/N from stacked/resampled image covariance.
      • Tests peak-pixel significance (must be >= min_peak_sigma) to reject diffuse background ripples.
      • Evaluates 4-gate detection criteria.
      • Computes both 3-sigma and 5-sigma calibrated limiting magnitudes (Upper Limits) in case of non-detection.

    Parameters:
      data           : 2D numpy float array of image counts (ADU)
      x, y           : 0-indexed pixel coordinates
      fwhm           : Seeing FWHM in pixels
      zp, zp_err     : Photometric zero-point and uncertainty
      gain           : Detector gain in e-/ADU
      readnoise      : Readout noise in e-
      r_ap_mult      : Aperture radius multiplier relative to FWHM
      min_peak_sigma : Minimum significance (above background std) required for peak pixel (default: 3.5)
      min_snr        : Minimum empirical S/N required for detection (default: 3.0)

    Returns:
      dict containing net_flux, flux_err, snr, peak_sig, m_cal, m_cal_err, ul_3sig_cal, ul_5sig_cal, status, r_ap
    """
    if data is None or np.isnan(x) or np.isnan(y):
        return {
            'net_flux': np.nan, 'flux_err': np.nan, 'bkg_flux_err': np.nan, 'snr': np.nan,
            'peak_val': np.nan, 'peak_sig': np.nan,
            'bkg_med': np.nan, 'bkg_std': np.nan,
            'm_inst': np.nan, 'merr': np.nan, 'm_cal': np.nan, 'm_cal_err': np.nan,
            'ul_3sig_cal': np.nan, 'ul_5sig_cal': np.nan, 'ul_cal': np.nan,
            'status': 'N/A', 'r_ap': np.nan
        }
    
    ny, nx = data.shape
    if x < 0 or x >= nx or y < 0 or y >= ny:
        return {
            'net_flux': np.nan, 'flux_err': np.nan, 'bkg_flux_err': np.nan, 'snr': np.nan,
            'peak_val': np.nan, 'peak_sig': np.nan,
            'bkg_med': np.nan, 'bkg_std': np.nan,
            'm_inst': np.nan, 'merr': np.nan, 'm_cal': np.nan, 'm_cal_err': np.nan,
            'ul_3sig_cal': np.nan, 'ul_5sig_cal': np.nan, 'ul_cal': np.nan,
            'status': 'OUT_OF_BOUNDS', 'r_ap': np.nan
        }

    try:
        from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry
        from astropy.stats import sigma_clipped_stats

        r_ap = max(2.5, float(fwhm * r_ap_mult))
        r_in = r_ap * PHOTO_AC_MULT
        r_out = r_in + PHOTO_OS_CONST

        ap = CircularAperture((x, y), r=r_ap)
        ann = CircularAnnulus((x, y), r_in=r_in, r_out=r_out)

        # ── 1. Local background extraction & empirical spatial variance ──────────
        ann_mask = ann.to_mask(method='center')
        ann_vals = ann_mask.get_values(data)
        if ann_vals is not None and len(ann_vals) > 5:
            _, bkg_med, bkg_std = sigma_clipped_stats(ann_vals, sigma=3.0, maxiters=5)
            bkg_med = float(bkg_med)
            bkg_std = float(bkg_std)
            n_ann   = float(len(ann_vals))
        else:
            tab_ann = aperture_photometry(data, ann)
            bkg_med = float(tab_ann['aperture_sum'][0] / ann.area)
            bkg_std = float(np.sqrt(max(0.0, bkg_med) / gain + (readnoise / gain)**2))
            n_ann   = float(ann.area)

        # ── 2. Net aperture flux ────────────────────────────────────────────────
        tab_ap = aperture_photometry(data, ap)
        ap_sum = float(tab_ap['aperture_sum'][0])
        n_ap   = float(ap.area)
        net_flux = float(ap_sum - bkg_med * n_ap)

        # ── 3. Empirical aperture flux noise budget ─────────────────────────────
        # Background uncertainty in aperture:
        # σ_bkg,ap = σ_sky,emp * sqrt( N_ap * (1 + N_ap / N_ann) )
        bkg_flux_err = float(bkg_std * np.sqrt(n_ap * (1.0 + n_ap / n_ann)))

        # Total flux error including source Poisson shot noise:
        source_shot_var = max(0.0, net_flux) / gain
        flux_err = float(np.sqrt(source_shot_var + bkg_flux_err**2))
        snr = float(net_flux / flux_err) if flux_err > 0 else 0.0

        # ── 4. Peak pixel significance in central 5x5 sub-stamp ─────────────────
        iy, ix = int(round(y)), int(round(x))
        sub_stamp = data[max(0, iy-2):min(ny, iy+3), max(0, ix-2):min(nx, ix+3)]
        peak_val = float(np.max(sub_stamp)) if sub_stamp.size > 0 else bkg_med
        peak_sig = float((peak_val - bkg_med) / bkg_std) if bkg_std > 0 else 0.0

        # ── 4b. Point-source PSF geometry and encircled energy fraction ──────
        sigma_psf = float(fwhm / 2.35482)
        area_psf  = float(2.0 * np.pi * (sigma_psf ** 2))   # ~ 1.1331 * fwhm^2
        f_enc = float(1.0 - np.exp(-0.5 * (r_ap / sigma_psf)**2))
        f_enc = max(0.5, min(1.0, f_enc))

        # Aperture correction in magnitudes (converts aperture flux to total PSF flux)
        ap_corr_mag = float(-2.5 * np.log10(f_enc))

        # Instrumental magnitude and error (corrected to total flux matching PSF Zero Point)
        m_inst = float(25.0 - 2.5 * np.log10(net_flux) - ap_corr_mag) if net_flux > 0 else np.nan
        merr   = float(1.0857 * (flux_err / net_flux)) if net_flux > 0 else np.nan

        # ── 5. Point-Source Limiting Fluxes & Realistic Upper Limits ───────────
        # 3σ and 5σ peak-detectable flux for a point source (Peak >= 3σ and 5σ above sky std)
        ul_3sig_flux = float(3.0 * bkg_std * area_psf)
        m_inst_ul3   = float(25.0 - 2.5 * np.log10(ul_3sig_flux)) if ul_3sig_flux > 0 else np.nan
        ul_3sig_cal  = float(m_inst_ul3 + zp) if np.isfinite(zp) and np.isfinite(m_inst_ul3) else np.nan

        ul_5sig_flux = float(5.0 * bkg_std * area_psf)
        m_inst_ul5   = float(25.0 - 2.5 * np.log10(ul_5sig_flux)) if ul_5sig_flux > 0 else np.nan
        ul_5sig_cal  = float(m_inst_ul5 + zp) if np.isfinite(zp) and np.isfinite(m_inst_ul5) else np.nan

        # Legacy alias for backward compatibility
        ul_cal = ul_3sig_cal

        # ── 6. 4-Gate Detection Decision ────────────────────────────────────────
        # A source is only classified as DETECTED if:
        #  1. Net flux > 0
        #  2. Empirical S/N >= min_snr (3.0σ)
        #  3. Peak pixel significance >= min_peak_sigma (3.5σ above local sky)
        #  4. Photometric zero-point is available
        is_detected = (net_flux > 0 and snr >= min_snr and peak_sig >= min_peak_sigma and np.isfinite(zp))

        if is_detected:
            m_cal     = float(m_inst + zp)
            m_cal_err = float(np.sqrt(merr**2 + zp_err**2))
            status    = 'DETECTED'
        else:
            m_cal     = np.nan
            m_cal_err = np.nan
            status    = f'LIMIT (<3σ, peak={peak_sig:.1f}σ)'

        return {
            'net_flux': net_flux, 'flux_err': flux_err, 'bkg_flux_err': bkg_flux_err,
            'snr': snr, 'peak_val': peak_val, 'peak_sig': peak_sig,
            'bkg_med': bkg_med, 'bkg_std': bkg_std,
            'm_inst': m_inst, 'merr': merr,
            'm_cal': m_cal, 'm_cal_err': m_cal_err,
            'ul_3sig_cal': ul_3sig_cal, 'ul_5sig_cal': ul_5sig_cal, 'ul_cal': ul_cal,
            'status': status, 'r_ap': r_ap
        }
    except Exception as e:
        print(f"  ⚠️ Aperture photometry calculation failed: {e}")
        return {
            'net_flux': np.nan, 'flux_err': np.nan, 'bkg_flux_err': np.nan, 'snr': np.nan,
            'peak_val': np.nan, 'peak_sig': np.nan,
            'bkg_med': np.nan, 'bkg_std': np.nan,
            'm_inst': np.nan, 'merr': np.nan, 'm_cal': np.nan, 'm_cal_err': np.nan,
            'ul_3sig_cal': np.nan, 'ul_5sig_cal': np.nan, 'ul_cal': np.nan,
            'status': 'ERROR', 'r_ap': np.nan
        }



def format_positions_summary_text(positions, target_name, filter_band, zp, zp_err, sigma_th=3.0):
    """
    Format a comprehensive multi-position photometry and coordinate comparison table.
    Covers:
      1. Exact Initial Location (GCN Input)
      2. Refined Centroid Location (FITS Sub-pixel Centroid)
      3. Matched ALLSTAR Detection (PSF-Fitted Source)
    Shows coordinates in both sexagesimal HH:MM:SS.ss, +DD:MM:SS.s and decimal degrees.
    """
    lines = []
    lines.append("=" * 125)
    lines.append(f"  TARGET COORDINATES & MULTI-POSITION PHOTOMETRY SUMMARY: {target_name}")
    lines.append("=" * 125)
    lines.append(f"  Photometric Filter : {filter_band}-band")
    lines.append(f"  Zero Point         : {zp:.4f} ± {zp_err:.4f}")
    lines.append("")

    # --- 1. Exact Initial (GCN) ---
    p_init = positions['initial']
    lines.append("  1. EXACT INITIAL LOCATION (GCN Input):")
    lines.append(f"     • Coordinates (Sexagesimal) : RA = {p_init['ra_hms']}, Dec = {p_init['dec_dms']}")
    lines.append(f"     • Coordinates (Decimal)     : RA = {p_init['ra_deg']:.6f}°, Dec = {p_init['dec_deg']:+.6f}°")
    lines.append(f"     • FITS Pixel (X, Y)         : X = {p_init['pixel_x']:.2f}, Y = {p_init['pixel_y']:.2f}")
    phot_i = p_init['phot']
    if phot_i.get('status') == 'DETECTED':
        lines.append(f"     • Detection Status          : ✅ DETECTED (S/N = {phot_i['snr']:.1f}σ, Peak = {phot_i.get('peak_sig', 0.0):.1f}σ)")
        lines.append(f"     • Calibrated Magnitude      : {phot_i['m_cal']:.3f} ± {phot_i['m_cal_err']:.3f} mag ({filter_band}-band, r={phot_i['r_ap']:.1f}px)")
        lines.append(f"     • Instrumental Photometry   : Inst Mag = {phot_i['m_inst']:.3f} ± {phot_i['merr']:.3f} (Net Flux = {phot_i['net_flux']:.1f} ± {phot_i['flux_err']:.1f} ADU)")
    elif np.isfinite(phot_i.get('ul_3sig_cal', np.nan)) or np.isfinite(phot_i.get('ul_cal', np.nan)):
        ul_3 = phot_i.get('ul_3sig_cal', phot_i.get('ul_cal', np.nan))
        ul_5 = phot_i.get('ul_5sig_cal', np.nan)
        lines.append(f"     • Detection Status          : ❌ NON-DETECTION (Peak={phot_i.get('peak_sig', 0.0):.1f}σ < 3.5σ threshold)")
        lines.append(f"     • 3σ Upper Limit            : Cal Mag > {ul_3:.2f} mag ({filter_band}-band, r={phot_i['r_ap']:.1f}px)")
        if np.isfinite(ul_5):
            lines.append(f"     • 5σ Upper Limit            : Cal Mag > {ul_5:.2f} mag ({filter_band}-band)")
        lines.append(f"     • Net Aperture Flux         : {phot_i['net_flux']:.1f} ± {phot_i['flux_err']:.1f} ADU (S/N = {phot_i['snr']:.1f}σ, Sky Std = {phot_i.get('bkg_std', 0.0):.1f} ADU)")
    else:
        lines.append("     • Aperture Photometry       : N/A")
    lines.append("")

    # --- 2. Refined Centroid (FITS) ---
    p_ref = positions['refined']
    lines.append("  2. REFINED CENTROID LOCATION (FITS Sub-Pixel Peak):")
    lines.append(f"     • Coordinates (Sexagesimal) : RA = {p_ref['ra_hms']}, Dec = {p_ref['dec_dms']}")
    lines.append(f"     • Coordinates (Decimal)     : RA = {p_ref['ra_deg']:.6f}°, Dec = {p_ref['dec_deg']:+.6f}°")
    lines.append(f"     • FITS Pixel (X, Y)         : X = {p_ref['pixel_x']:.2f}, Y = {p_ref['pixel_y']:.2f}")
    lines.append(f"     • Astrometric Offset (GCN)  : Δ = {p_ref['offset_arcsec']:.2f}\" ({p_ref['offset_pix']:.2f} px)")
    phot_r = p_ref['phot']
    if phot_r.get('status') == 'DETECTED':
        lines.append(f"     • Detection Status          : ✅ DETECTED (S/N = {phot_r['snr']:.1f}σ, Peak = {phot_r.get('peak_sig', 0.0):.1f}σ)")
        lines.append(f"     • Calibrated Magnitude      : {phot_r['m_cal']:.3f} ± {phot_r['m_cal_err']:.3f} mag ({filter_band}-band, r={phot_r['r_ap']:.1f}px)")
        lines.append(f"     • Instrumental Photometry   : Inst Mag = {phot_r['m_inst']:.3f} ± {phot_r['merr']:.3f} (Net Flux = {phot_r['net_flux']:.1f} ± {phot_r['flux_err']:.1f} ADU)")
    elif np.isfinite(phot_r.get('ul_3sig_cal', np.nan)) or np.isfinite(phot_r.get('ul_cal', np.nan)):
        ul_3 = phot_r.get('ul_3sig_cal', phot_r.get('ul_cal', np.nan))
        ul_5 = phot_r.get('ul_5sig_cal', np.nan)
        lines.append(f"     • Detection Status          : ❌ NON-DETECTION (Peak={phot_r.get('peak_sig', 0.0):.1f}σ < 3.5σ threshold)")
        lines.append(f"     • 3σ Upper Limit            : Cal Mag > {ul_3:.2f} mag ({filter_band}-band, r={phot_r['r_ap']:.1f}px)")
        if np.isfinite(ul_5):
            lines.append(f"     • 5σ Upper Limit            : Cal Mag > {ul_5:.2f} mag ({filter_band}-band)")
        lines.append(f"     • Net Aperture Flux         : {phot_r['net_flux']:.1f} ± {phot_r['flux_err']:.1f} ADU (S/N = {phot_r['snr']:.1f}σ, Sky Std = {phot_r.get('bkg_std', 0.0):.1f} ADU)")
    else:
        lines.append("     • Aperture Photometry       : N/A")
    lines.append("")

    # --- 3. Matched ALLSTAR PSF Detection ---
    p_mat = positions['matched']
    lines.append("  3. MATCHED ALLSTAR DETECTION (DAOPHOT PSF-Fitted Source):")
    if p_mat.get('detected'):
        lines.append(f"     • Coordinates (Sexagesimal) : RA = {p_mat['ra_hms']}, Dec = {p_mat['dec_dms']}")
        lines.append(f"     • Coordinates (Decimal)     : RA = {p_mat['ra_deg']:.6f}°, Dec = {p_mat['dec_deg']:+.6f}°")
        lines.append(f"     • FITS Pixel (X, Y)         : X = {p_mat['pixel_x']:.2f}, Y = {p_mat['pixel_y']:.2f} (Star ID #{p_mat['id']})")
        lines.append(f"     • Astrometric Offsets       : Δ(Initial) = {p_mat['offset_init_arcsec']:.2f}\" ({p_mat['offset_init_pix']:.2f} px) | Δ(Refined) = {p_mat['offset_ref_arcsec']:.2f}\" ({p_mat['offset_ref_pix']:.2f} px)")
        lines.append(f"     • Detection Type            : {p_mat['detect_type']} (Sigma Thresh: {sigma_th})")
        lines.append(f"     • PSF Calibrated Mag        : {p_mat['psf_cal_mag']:.3f} ± {p_mat['psf_cal_mag_err']:.3f} mag ({filter_band}-band)")
        lines.append(f"     • PSF Instrumental Mag      : {p_mat['psf_inst_mag']:.3f} ± {p_mat['psf_inst_merr']:.3f}")
        lines.append(f"     • PSF Profile Fit Quality   : CHI = {p_mat['chi']:.2f}, SHARP = {p_mat['sharp']:+.3f}")
        if p_mat['chi'] > 2.0 or abs(p_mat['sharp']) > 0.5:
            lines.append("     • ⚠️ Caution               : Fit parameters exceed stellar thresholds (CHI>2 or |SHARP|>0.5); likely noise or blend.")
    else:
        lines.append("     • Status                    : ❌ NON-DETECTION (Target not recovered in .als table)")
    lines.append("")

    # --- Side-by-Side Comparison Table ---
    lines.append("-" * 125)
    lines.append(f"{'POSITION':<28} | {'RA (HH:MM:SS.ss)  Dec (+DD:MM:SS.s)':<34} | {'RA (deg)':<11} {'Dec (deg)':<11} | {'Pixel (X, Y)':<15} | {'Offset (GCN)':<12} | {f'Photometry ({filter_band}-band)'}")
    lines.append("-" * 125)

    # Row 1: Initial
    c_i = f"{p_init['ra_hms']}   {p_init['dec_dms']}"
    pix_i = f"({p_init['pixel_x']:.1f}, {p_init['pixel_y']:.1f})"
    if phot_i.get('status') == 'DETECTED':
        phot_i_str = f"Aper {phot_i['m_cal']:.2f} ± {phot_i['m_cal_err']:.2f} mag"
    elif np.isfinite(phot_i.get('ul_3sig_cal', np.nan)) or np.isfinite(phot_i.get('ul_cal', np.nan)):
        ul_3 = phot_i.get('ul_3sig_cal', phot_i.get('ul_cal', np.nan))
        phot_i_str = f"Aper > {ul_3:.2f} mag (3σ)"
    else:
        phot_i_str = "N/A"
    lines.append(f"{'1. Exact Initial (GCN)':<28} | {c_i:<34} | {p_init['ra_deg']:<11.6f} {p_init['dec_deg']:<+11.6f} | {pix_i:<15} | {'-':<12} | {phot_i_str}")

    # Row 2: Refined
    c_r = f"{p_ref['ra_hms']}   {p_ref['dec_dms']}"
    pix_r = f"({p_ref['pixel_x']:.1f}, {p_ref['pixel_y']:.1f})"
    off_r = f"{p_ref['offset_arcsec']:.2f}\" ({p_ref['offset_pix']:.1f}p)"
    if phot_r.get('status') == 'DETECTED':
        phot_r_str = f"Aper {phot_r['m_cal']:.2f} ± {phot_r['m_cal_err']:.2f} mag"
    elif np.isfinite(phot_r.get('ul_3sig_cal', np.nan)) or np.isfinite(phot_r.get('ul_cal', np.nan)):
        ul_3 = phot_r.get('ul_3sig_cal', phot_r.get('ul_cal', np.nan))
        phot_r_str = f"Aper > {ul_3:.2f} mag (3σ)"
    else:
        phot_r_str = "N/A"
    lines.append(f"{'2. Refined Centroid (FITS)':<28} | {c_r:<34} | {p_ref['ra_deg']:<11.6f} {p_ref['dec_deg']:<+11.6f} | {pix_r:<15} | {off_r:<12} | {phot_r_str}")

    # Row 3: Matched
    if p_mat.get('detected'):
        c_m = f"{p_mat['ra_hms']}   {p_mat['dec_dms']}"
        pix_m = f"({p_mat['pixel_x']:.1f}, {p_mat['pixel_y']:.1f})"
        off_m = f"{p_mat['offset_init_arcsec']:.2f}\" ({p_mat['offset_init_pix']:.1f}p)"
        phot_m_str = f"PSF  {p_mat['psf_cal_mag']:.3f} ± {p_mat['psf_cal_mag_err']:.3f} mag (ID#{p_mat['id']})"
        lines.append(f"{'3. Matched ALLSTAR (PSF)':<28} | {c_m:<34} | {p_mat['ra_deg']:<11.6f} {p_mat['dec_deg']:<+11.6f} | {pix_m:<15} | {off_m:<12} | {phot_m_str}")
    else:
        lines.append(f"{'3. Matched ALLSTAR (PSF)':<28} | {'N/A':<34} | {'N/A':<11} {'N/A':<11} | {'N/A':<15} | {'N/A':<12} | NON-DETECTION")

    lines.append("-" * 125)

    # --- Ready-to-Publish GCN Circular Statement ---
    ul_report = phot_i.get('ul_3sig_cal', phot_i.get('ul_cal', np.nan))
    ul5_report = phot_i.get('ul_5sig_cal', np.nan)
    ul_str = f"> {ul_report:.2f}" if np.isfinite(ul_report) else "N/A"
    ul5_str = f"> {ul5_report:.2f}" if np.isfinite(ul5_report) else "N/A"

    lines.append("")
    lines.append("  📢 READY-TO-PUBLISH GCN CIRCULAR / ASTRONOMICAL REPORT STATEMENT:")
    lines.append("  " + "─" * 123)
    if phot_i.get('status') == 'DETECTED':
        lines.append(f"  \"We observed the optical field of {target_name} in {filter_band}-band. An optical counterpart is detected at")
        lines.append(f"   RA = {p_init['ra_hms']}, Dec = {p_init['dec_dms']} with a calibrated magnitude of {filter_band} = {phot_i['m_cal']:.2f} ± {phot_i['m_cal_err']:.2f} mag")
        lines.append(f"   (S/N = {phot_i['snr']:.1f}σ, calibrated against {VIZIER_CATALOG_ID}).\"")
    else:
        lines.append(f"  \"We observed the field of {target_name} in {filter_band}-band. No optical counterpart is detected at")
        lines.append(f"   the counterpart position (RA = {p_init['ra_hms']}, Dec = {p_init['dec_dms']}) above the 3-sigma background noise.")
        lines.append(f"   We derive a 3-sigma limiting magnitude of {filter_band} {ul_str} mag (5-sigma limit {filter_band} {ul5_str} mag,")
        lines.append(f"   calibrated against {VIZIER_CATALOG_ID}; Zero Point = {zp:.4f} ± {zp_err:.4f}).\"")
    lines.append("  " + "─" * 123)
    lines.append("")
    return "\n".join(lines)


def parse_trigger_time(t0, fits_header=None):
    """
    Parse trigger time T0 in either:
    - Full ISO timestamp (e.g. "2026-10-10 19:17:03.16" or "2026-10-10T19:17:03.16")
    - Text date & time (e.g. "19:17:03.16 UT on 10 October 2026" or "10 October 2026 19:17:03.16 UT")
    - Julian Date float (e.g. 2461324.303509)
    - Time-of-day only (e.g. "19:17:03.16 UT", uses DATE-OBS from FITS header for date)
    Returns Julian Date (float).
    """
    if t0 is None:
        return np.nan
    if isinstance(t0, (int, float)):
        return float(t0)

    s = str(t0).strip()
    try:
        val = float(s)
        if val > 2400000:
            return val
    except ValueError:
        pass

    clean_s = re.sub(r'\b(UT|UTC|on)\b', ' ', s, flags=re.IGNORECASE).strip()
    clean_s = re.sub(r'\s+', ' ', clean_s)

    # Check if time-of-day only (e.g. HH:MM:SS[.ss])
    time_only_match = re.fullmatch(r'(\d{1,2}):(\d{2})(?::(\d{2}(?:\.\d+)?))?', clean_s)
    if time_only_match:
        if fits_header is not None:
            date_obs = fits_header.get('DATE-OBS', fits_header.get('DATE', None))
            if date_obs:
                date_part = str(date_obs).replace('"T"', 'T').split('T')[0].split()[0].strip()
                try:
                    return float(Time(f"{date_part} {clean_s}", scale='utc').jd)
                except Exception:
                    pass
        return np.nan

    # Try standard astropy Time formats (iso, isot)
    for fmt in [None, 'iso', 'isot']:
        try:
            return float(Time(clean_s if fmt != 'isot' else clean_s.replace(' ', 'T'), format=fmt, scale='utc').jd)
        except Exception:
            pass

    # Try dateutil flexible date parser (handles '10 October 2026 19:17:03.16', '19:17:03.16 on 10 October 2026', etc.)
    try:
        dt = dateutil.parser.parse(clean_s, fuzzy=True)
        return float(Time(dt, scale='utc').jd)
    except Exception:
        pass

    try:
        return float(s)
    except Exception:
        return np.nan


def fetch_online_archival_cutout(ra, dec, size_arcsec, pixscale, survey="ps1"):
    """
    Fetch an archival image cutout from the Legacy Survey / Pan-STARRS viewer API.
    survey can be 'ps1' (Pan-STARRS1) or 'ls-dr10' (DESI Legacy Survey).
    Returns a PIL Image object or None.
    """
    size_pix = int(round(size_arcsec / pixscale))
    size_pix = max(80, min(512, size_pix))
    url = (
        f"https://www.legacysurvey.org/viewer/jpeg-cutout"
        f"?ra={ra:.6f}&dec={dec:.6f}&width={size_pix}&height={size_pix}"
        f"&layer={survey}&pixscale={pixscale:.4f}"
    )
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (Astro-Pipeline)"})
        with urllib.request.urlopen(req, timeout=8) as response:
            return Image.open(io.BytesIO(response.read()))
    except Exception:
        if survey == "ps1":
            try:
                url_fb = (
                    f"https://www.legacysurvey.org/viewer/jpeg-cutout"
                    f"?ra={ra:.6f}&dec={dec:.6f}&width={size_pix}&height={size_pix}"
                    f"&layer=ls-dr10&pixscale={pixscale:.4f}"
                )
                req_fb = urllib.request.Request(url_fb, headers={"User-Agent": "Mozilla/5.0 (Astro-Pipeline)"})
                with urllib.request.urlopen(req_fb, timeout=8) as response:
                    return Image.open(io.BytesIO(response.read()))
            except Exception:
                pass
        return None


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
        elif key == 'GA':
            new_lines.append(f"GA = {GAIN:.2f}\n")
            updated['GA'] = True
        elif key == 'RE':
            new_lines.append(f"RE = {READ_NOISE:.2f}\n")
            updated['RE'] = True
        else:
            new_lines.append(line)

    if 'FW' not in updated: new_lines.append(f"FW = {fwhm:.3f}\n")
    if 'FI' not in updated: new_lines.append(f"FI = {fi_val:.3f}\n")
    if 'PS' not in updated: new_lines.append(f"PS = {ps_val:.3f}\n")
    if 'GA' not in updated: new_lines.append(f"GA = {GAIN:.2f}\n")
    if 'RE' not in updated: new_lines.append(f"RE = {READ_NOISE:.2f}\n")

    with open(opt_path, 'w') as fh:
        fh.writelines(new_lines)

    print(f"[OPT] Updated daophot.opt — GA={GAIN:.2f}, RE={READ_NOISE:.2f}, FW={fwhm:.3f}, FI={fi_val:.3f}, PS={ps_val:.3f}")


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
    t_hms, t_dms = format_hms_dms(target_ra, target_dec)
    print(f"  Target expected at RA={target_ra:.6f}°, Dec={target_dec:+.6f}° ({t_hms}, {t_dms})")
    print(f"  Target expected at X={target_x:.2f}, Y={target_y:.2f}")

    dx   = als_table['XCENTER'] - target_x
    dy   = als_table['YCENTER'] - target_y
    dist = np.sqrt(dx**2 + dy**2)
    idx  = int(np.argmin(dist))
    dmin = float(dist[idx])

    if dmin <= match_radius_pix:
        row = als_table[idx]
        matched_ra, matched_dec = wcs.wcs_pix2world(row['XCENTER'], row['YCENTER'], 0)
        matched_ra, matched_dec = float(matched_ra), float(matched_dec)
        m_hms, m_dms = format_hms_dms(matched_ra, matched_dec)
        print(f"  ✅ Matched!  ID={row['ID']}  "
              f"X={row['XCENTER']:.2f}  Y={row['YCENTER']:.2f}  "
              f"dist={dmin:.2f}px")
        print(f"     Matched WCS : RA={matched_ra:.6f}°, Dec={matched_dec:+.6f}° ({m_hms}, {m_dms})")
        print(f"  MAG={row['MAG']:.3f}±{row['MERR']:.3f}  "
              f"CHI={row['CHI']:.2f}  SHARP={row['SHARP']:.3f}")
        return {
            'detected'   : True,
            'mag'        : float(row['MAG']),
            'merr'       : float(row['MERR']),
            'chi'        : float(row['CHI']),
            'sharp'      : float(row['SHARP']),
            'pixel_x'    : float(row['XCENTER']),
            'pixel_y'    : float(row['YCENTER']),
            'matched_ra' : matched_ra,
            'matched_dec': matched_dec,
            'matched_hms': m_hms,
            'matched_dms': m_dms,
            'dist_pix'   : dmin,
            'id'         : int(row['ID']),
        }
    else:
        print(f"  ❌ No match within {match_radius_pix:.1f}px  "
              f"(closest = {dmin:.1f}px)")
        return _empty_result()


def _empty_result():
    return dict(detected=False, mag=np.nan, merr=np.nan,
                chi=np.nan, sharp=np.nan, pixel_x=np.nan,
                pixel_y=np.nan, matched_ra=np.nan, matched_dec=np.nan,
                matched_hms='N/A', matched_dms='N/A',
                dist_pix=np.nan, id=-1)


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


def plot_target_confirmation(image_path, target_x, target_y, target_name, init_x=None, init_y=None, save=True):
    with fits.open(image_path) as hdul:
        clean_data = hdul[0].data.astype(float)
        
    zscale     = ZScaleInterval()
    vmin, vmax = zscale.get_limits(clean_data)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))
    for ax in (ax1, ax2):
        ax.imshow(clean_data, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
        if init_x is not None and init_y is not None and (abs(init_x - target_x) > 0.3 or abs(init_y - target_y) > 0.3):
            ax.scatter(init_x, init_y, marker='+', color='red', s=400, lw=2, label='Initial GCN', alpha=0.7)
        ax.scatter(target_x, target_y, marker='o', facecolors='none', edgecolors='lime', s=250, lw=2, label='Target / Centroid')
        ax.scatter(target_x, target_y, marker='+', color='lime', s=100, lw=1.5)
        
    ax1.set_title(f"Full Image: {target_name} Position", fontsize=14)
    ax1.add_patch(patches.Circle((target_x, target_y), radius=100,
                                  edgecolor='red', facecolor='none', lw=2, linestyle='--'))
    ax2.set_title(f"Zoom: {target_name} Site", fontsize=14)
    ax2.add_patch(patches.Circle((target_x, target_y), radius=6,
                                  edgecolor='lime', facecolor='none', lw=1.5, linestyle='--'))
    zoom = 60
    ax2.set_xlim(target_x - zoom, target_x + zoom)
    ax2.set_ylim(target_y - zoom, target_y + zoom)
    ax2.grid(color='white', alpha=0.2, linestyle=':')
    ax2.legend(loc='upper right', framealpha=0.85)
    plt.tight_layout()
    if save:
        p = os.path.join(PLOTS_DIR, 'photometry', f'{target_name}_confirmation.png')
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"[Plot] Target Confirmation → {p}")
    plt.close(fig)


def confirm_and_refine_target(
    image_path           = STACKED_IMAGE,
    target_ra            = TARGET_RA,
    target_dec           = TARGET_DEC,
    target_name          = TARGET_NAME,
    search_radius_arcsec = 5.0,
    cutout_size_arcsec   = 40.0,
    selected_candidate   = 1,
    ref_image            = ASTROMETRY_REF_IMAGE,
    plots_dir            = PLOTS_DIR,
    auto_update          = True,
    show_plot            = True,
    save_plot            = True,
):
    """
    Pre-confirmation and coordinate refinement for transient / GRB optical afterglow.

    1. Parses RA and Dec in either decimal degrees or sexagesimal HMS/DMS strings.
    2. Outputs direct clickable links to online archival sky surveys (Legacy Survey, Pan-STARRS, Aladin).
    3. Fetches pre-burst archival catalogue image from Pan-STARRS1 / Legacy Survey.
    4. Generates a 3-panel side-by-side comparison:
       - Panel 1: Your FITS File (Wide Field Context)
       - Panel 2: Online Catalogue (Pre-Burst Archival: Pan-STARRS1 / Legacy Survey)
       - Panel 3: Your FITS File (Zoomed Frame) with a VERY TRANSPARENT marker (alpha=0.3)
                 at the initial coordinate so underlying pixels are fully visible.
    5. Measures the sub-pixel centroid of candidate point sources and converts to refined WCS (RA, Dec).
    """
    parsed_ra, parsed_dec = parse_target_coords(target_ra, target_dec)
    if parsed_ra is None or parsed_dec is None:
        print("  [Pre-Check] ⚠️ No valid target coordinates provided (TARGET_RA / TARGET_DEC is None).")
        return None

    if not os.path.exists(image_path):
        print(f"  [Pre-Check] ❌ Image not found: {image_path}")
        return None

    with fits.open(image_path) as hdul:
        header = hdul[0].header
        data   = hdul[0].data.astype(float)

    wcs = get_wcs(header, ref_image)
    pix_scales = proj_plane_pixel_scales(wcs)
    pix_scale = float(pix_scales[0] * 3600.0)
    ny, nx = data.shape

    init_x, init_y = wcs.wcs_world2pix(parsed_ra, parsed_dec, 0)
    init_x, init_y = float(init_x), float(init_y)
    init_ra_hms, init_dec_dms = format_hms_dms(parsed_ra, parsed_dec)

    print("\n" + "╔" + "═"*66 + "╗")
    print(f"║{'TARGET PRE-CONFIRMATION & CENTROID REFINEMENT':^66}║")
    print(f"║  Target  : {target_name:<55}║")
    print(f"║  Initial : RA = {parsed_ra:.6f}°, Dec = {parsed_dec:+.6f}°{'':<21}║")
    print(f"║  Sexages : RA = {init_ra_hms}, Dec = {init_dec_dms}{'':<20}║")
    if str(target_ra) != f"{parsed_ra:.6f}" or str(target_dec) != f"{parsed_dec:+.6f}":
        print(f"║  Input   : RA = {str(target_ra):<20} Dec = {str(target_dec):<22}║")
    print(f"║  FITS Pix: X = {init_x:.2f}, Y = {init_y:.2f} (Scale: {pix_scale:.4f}\"/pix){'':<20}║")
    print("╚" + "═"*66 + "╝")

    # ── 1. Online Archival Sky Survey Links ────────────────────────────
    legacy_url = f"https://www.legacysurvey.org/viewer?ra={parsed_ra:.6f}&dec={parsed_dec:.6f}&layer=ls-dr10&zoom=16"
    ps1_url    = f"https://ps1images.stsci.edu/cgi-bin/ps1cutouts?pos={parsed_ra:.6f}+{parsed_dec:.6f}&filter=color"
    aladin_url = f"https://aladin.cds.unistra.fr/AladinLite/?target={parsed_ra:.6f}%20%2B{parsed_dec:.6f}&fov=0.03&survey=P%2FPanSTARRS%2FDR1%2Fcolor-r-i-z"

    print("\n  🌐 Direct Online Archival Check (Pre-explosion images):")
    print(f"     • Legacy Survey (DESI / PS1 / SDSS) : {legacy_url}")
    print(f"     • Pan-STARRS1 STScI Image Cutout    : {ps1_url}")
    print(f"     • CDS Aladin Lite Interactive View  : {aladin_url}")

    # ── 2. Local Candidate Detection in Search Area ───────────────────
    search_r_pix = search_radius_arcsec / pix_scale
    box_r = int(np.ceil(search_r_pix + 15))

    x0, x1 = max(0, int(round(init_x)) - box_r), min(nx, int(round(init_x)) + box_r + 1)
    y0, y1 = max(0, int(round(init_y)) - box_r), min(ny, int(round(init_y)) + box_r + 1)
    sub = data[y0:y1, x0:x1]

    _, med_bkg, std_bkg = sigma_clipped_stats(sub, sigma=3.0, maxiters=5)
    sub_s = sub - med_bkg

    # ── ROBUST CANDIDATE DETECTION & CENTROID REFINEMENT ────────────────
    # In pure 2D Gaussian background noise, extreme value statistics dictate that a 5" search box
    # contains dozens of random 2.0σ-2.5σ noise fluctuations.
    # We require a strict peak threshold of >= 3.5σ above the local background standard deviation
    # to prevent centroid refinement from locking onto noise ripples and creating artificial detections.
    CANDIDATE_SIGMA_THRESH = 3.5
    local_max = (sub_s == maximum_filter(sub_s, size=5)) & (sub_s >= CANDIDATE_SIGMA_THRESH * std_bkg)
    yy_loc, xx_loc = np.where(local_max)

    candidates = []
    for xi, yi in zip(xx_loc, yy_loc):
        xf = x0 + xi
        yf = y0 + yi
        d_px = np.hypot(xf - init_x, yf - init_y)
        if d_px <= search_r_pix:
            c_ra, c_dec = wcs.wcs_pix2world(xf, yf, 0)
            candidates.append({
                "x": float(xf), "y": float(yf),
                "peak_flux": float(sub_s[yi, xi]),
                "snr": float(sub_s[yi, xi] / std_bkg),
                "dist_pix": float(d_px),
                "dist_as": float(d_px * pix_scale),
                "ra": float(c_ra), "dec": float(c_dec)
            })

    # Prioritize candidates closest to the GCN position
    candidates.sort(key=lambda c: c["dist_pix"])

    # ── 3. Centroid Refinement ─────────────────────────────────────────
    detected = len(candidates) > 0
    if detected:
        cand_idx = max(1, min(selected_candidate, len(candidates))) - 1
        sel = candidates[cand_idx]
        px_c, py_c = sel["x"], sel["y"]

        # Sub-pixel centroid using 5x5 weighted box
        c_r = 2
        cx0, cx1 = max(0, int(round(px_c)) - c_r), min(nx, int(round(px_c)) + c_r + 1)
        cy0, cy1 = max(0, int(round(py_c)) - c_r), min(ny, int(round(py_c)) + c_r + 1)
        w_stamp = data[cy0:cy1, cx0:cx1] - med_bkg
        w_stamp[w_stamp < 0] = 0
        if w_stamp.sum() > 0:
            wy, wx = np.indices(w_stamp.shape)
            ref_x = float(cx0 + (wx * w_stamp).sum() / w_stamp.sum())
            ref_y = float(cy0 + (wy * w_stamp).sum() / w_stamp.sum())
        else:
            ref_x, ref_y = float(px_c), float(py_c)

        ref_ra, ref_dec = wcs.wcs_pix2world(ref_x, ref_y, 0)
        ref_ra, ref_dec = float(ref_ra), float(ref_dec)
        ref_ra_hms, ref_dec_dms = format_hms_dms(ref_ra, ref_dec)
        offset_as = float(SkyCoord(parsed_ra*u.deg, parsed_dec*u.deg).separation(SkyCoord(ref_ra*u.deg, ref_dec*u.deg)).arcsec)
        offset_px = float(np.hypot(ref_x - init_x, ref_y - init_y))

        print(f"\n  🔍 Detected {len(candidates)} candidate point source(s) >= {CANDIDATE_SIGMA_THRESH:.1f}σ within {search_radius_arcsec:.1f}\":")
        for i, c in enumerate(candidates):
            marker_str = "--> [SELECTED]" if i == cand_idx else "   "
            print(f"    {marker_str} #{i+1}: Pixel=({c['x']:.2f}, {c['y']:.2f})  "
                  f"Offset={c['dist_as']:.2f}\" ({c['dist_pix']:.2f}px)  "
                  f"S/N~{c['snr']:.1f}σ  RA={c['ra']:.6f}°, Dec={c['dec']:+.6f}°")

        print(f"\n  🎯 Refined Centroid (Candidate #{cand_idx+1}):")
        print(f"     Refined FITS Pixel : X = {ref_x:.2f}, Y = {ref_y:.2f}")
        print(f"     Refined WCS Coords : RA = {ref_ra:.6f}°, Dec = {ref_dec:+.6f}°")
        print(f"     Refined Sexagesimal: RA = {ref_ra_hms}, Dec = {ref_dec_dms}")
        print(f"     Astrometric Offset : {offset_as:.2f}\" ({offset_px:.2f} pixels)")
    else:
        print(f"\n  ⚠️ No point source found above {CANDIDATE_SIGMA_THRESH:.1f}σ within {search_radius_arcsec:.1f}\" of GCN position.")
        print("     Retaining initial coordinates for pipeline execution (preventing noise-hunting centroid drift).")
        ref_ra, ref_dec = parsed_ra, parsed_dec
        ref_ra_hms, ref_dec_dms = init_ra_hms, init_dec_dms
        ref_x, ref_y    = init_x, init_y
        offset_as       = 0.0
        offset_px       = 0.0

    # ── 4. 3-Panel Side-by-Side Comparison Plot ───────────────────────
    print("\n  📥 Fetching pre-burst archival catalogue image...")
    online_img = fetch_online_archival_cutout(parsed_ra, parsed_dec, cutout_size_arcsec, pix_scale, survey="ps1")

    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(data)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(21, 7))

    # Panel 1: Your FITS File (Wide Context)
    cut_pix = int(round(cutout_size_arcsec / pix_scale / 2.0))
    wx0, wx1 = max(0, int(round(init_x)) - cut_pix), min(nx, int(round(init_x)) + cut_pix)
    wy0, wy1 = max(0, int(round(init_y)) - cut_pix), min(ny, int(round(init_y)) + cut_pix)

    ax1.imshow(data, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)
    ax1.set_xlim(wx0, wx1)
    ax1.set_ylim(wy0, wy1)
    ax1.scatter([init_x], [init_y], marker="+", color="red", s=300, lw=2, label="Initial GCN")
    ax1.add_patch(patches.Circle((init_x, init_y), radius=search_r_pix,
                                 edgecolor="yellow", facecolor="none", lw=1.5, ls="--",
                                 label=f"Search ({search_radius_arcsec:.1f}\")"))
    ax1.set_title(f"1. Your FITS Image (Wide: {cutout_size_arcsec:.0f}\" × {cutout_size_arcsec:.0f}\")", fontsize=13)
    ax1.set_xlabel("FITS X (pixels)")
    ax1.set_ylabel("FITS Y (pixels)")
    ax1.legend(loc="upper right", framealpha=0.85)

    # Panel 2: Online Catalogue Archival Pre-burst
    if online_img is not None:
        ax2.imshow(online_img)
        w_on, h_on = online_img.size
        ax2.scatter([w_on / 2.0], [h_on / 2.0], marker="+", color="red", s=350, lw=2, alpha=0.85, label="Initial GCN Position")
        ax2.add_patch(patches.Circle((w_on / 2.0, h_on / 2.0), radius=(search_radius_arcsec / pix_scale),
                                     edgecolor="yellow", facecolor="none", lw=1.5, ls="--",
                                     label=f"Search ({search_radius_arcsec:.1f}\")"))
        ax2.set_title("2. Online Catalogue (Pre-Burst: Pan-STARRS1)", fontsize=13)
        ax2.set_xlabel("Sky Field (North Up, East Left)")
        ax2.legend(loc="upper right", framealpha=0.85)
        ax2.set_xticks([])
        ax2.set_yticks([])
    else:
        ax2.text(0.5, 0.5, "Archival Catalogue Image\n(Check links printed in log)",
                 ha="center", va="center", color="gray", fontsize=12, transform=ax2.transAxes)
        ax2.set_title("2. Online Catalogue (Offline)", fontsize=13)
        ax2.set_xticks([])
        ax2.set_yticks([])

    # Panel 3: Your FITS File (Zoomed Site with VERY TRANSPARENT marker at initial coord)
    zoom_arcsec = max(12.0, search_radius_arcsec * 2.5)
    zoom_pix = int(round(zoom_arcsec / pix_scale / 2.0))
    zx0, zx1 = max(0, int(round(init_x)) - zoom_pix), min(nx, int(round(init_x)) + zoom_pix)
    zy0, zy1 = max(0, int(round(init_y)) - zoom_pix), min(ny, int(round(init_y)) + zoom_pix)

    ax3.imshow(data, cmap="gray", origin="lower", vmin=vmin, vmax=vmax)
    ax3.set_xlim(zx0, zx1)
    ax3.set_ylim(zy0, zy1)

    # Very transparent markers for initial position (alpha=0.25 - 0.30)
    ax3.scatter([init_x], [init_y], marker="+", color="red", s=400, lw=1.5, alpha=0.30,
                label="Initial GCN (Transparent Reticle, α=0.3)")
    ax3.add_patch(patches.Circle((init_x, init_y), radius=search_r_pix,
                                 edgecolor="yellow", facecolor="none", lw=1.2, ls="--", alpha=0.35,
                                 label=f"Search Area ({search_radius_arcsec:.1f}\")"))
    ax3.add_patch(patches.Circle((init_x, init_y), radius=2.5,
                                 edgecolor="red", facecolor="none", lw=1.0, ls=":", alpha=0.25))

    for i, c in enumerate(candidates):
        is_sel = (detected and i == cand_idx)
        col = "cyan" if is_sel else "lime"
        s_sz = 260 if is_sel else 180
        ax3.scatter([c["x"]], [c["y"]], marker="o", facecolors="none", edgecolors=col, s=s_sz, lw=2, alpha=0.85)
        ax3.scatter([c["x"]], [c["y"]], marker="x", color=col, s=70, lw=1.5, alpha=0.85)
        lbl_text = f"#{i+1} ({c['dist_as']:.1f}\")"
        ax3.text(c["x"] + 1.2, c["y"] + 1.2, lbl_text, color=col, fontsize=11, fontweight="bold", alpha=0.9)

    if detected:
        ax3.scatter([ref_x], [ref_y], marker="*", color="yellow", s=180, alpha=0.9, label=f"Refined Centroid (#{cand_idx+1})")
        ax3.set_title(f"3. Your FITS Zoom (Refined Offset = {offset_as:.2f}\" / {offset_px:.2f}px)", fontsize=13)
    else:
        ax3.set_title(f"3. Your FITS Zoom Site", fontsize=13)

    ax3.set_xlabel("FITS X (pixels)")
    ax3.set_ylabel("FITS Y (pixels)")
    ax3.grid(color="white", alpha=0.20, ls=":")
    ax3.legend(loc="upper right", framealpha=0.85)
    plt.tight_layout()

    if save_plot:
        plot_dir = os.path.join(plots_dir, 'photometry')
        Path(plot_dir).mkdir(parents=True, exist_ok=True)
        p = os.path.join(plot_dir, f"{target_name}_pre_confirmation.png")
        fig.savefig(p, dpi=150, bbox_inches='tight')
        print(f"\n  [Plot] Target 3-Panel Pre-Confirmation Saved → {p}")

    if show_plot:
        try:
            plt.show()
        except Exception:
            pass
    plt.close(fig)

    print("═"*68 + "\n")
    return {
        "detected"       : detected,
        "init_ra"        : parsed_ra,
        "init_dec"       : parsed_dec,
        "init_ra_hms"    : init_ra_hms,
        "init_dec_dms"   : init_dec_dms,
        "init_x"         : init_x,
        "init_y"         : init_y,
        "refined_ra"     : ref_ra,
        "refined_dec"    : ref_dec,
        "refined_ra_hms" : ref_ra_hms,
        "refined_dec_dms": ref_dec_dms,
        "refined_x"      : ref_x,
        "refined_y"      : ref_y,
        "offset_arcsec"  : offset_as,
        "offset_pix"     : offset_px,
        "candidates"     : candidates,
        "urls"           : {"legacy_survey": legacy_url, "panstarrs": ps1_url, "aladin": aladin_url}
    }


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
    pixel_x = final_results.get('pixel_x', np.nan)
    pixel_y = final_results.get('pixel_y', np.nan)
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
        raw_date_obs = orig_header.get('DATE-OBS') or orig_header.get('DATE') or 'UNKNOWN'
        date_obs_str = str(raw_date_obs).strip()
        t_obs_jd = np.nan
        t0_jd = np.nan
        delta_t_days = np.nan

        try:
            t0_jd = parse_trigger_time(GRB_T0_JD, orig_header)
        except Exception as e:
            print(f"⚠️ Could not parse T0: {e}")

        if date_obs_str != 'UNKNOWN':
            date_obs_clean = date_obs_str.replace('"T"', 'T')
            for fmt in ['isot', 'iso', None]:
                try:
                    s = date_obs_clean if fmt != 'isot' else date_obs_clean.replace(' ', 'T')
                    t_obs = Time(s, format=fmt, scale='utc')
                    t_obs_jd = float(t_obs.jd)
                    break
                except Exception:
                    continue

        if np.isfinite(t_obs_jd) and np.isfinite(t0_jd):
            delta_t_days = t_obs_jd - t0_jd

        f.write(f"DATE-OBS    : {date_obs_str}\n")
        f.write(f"Obs JD      : {t_obs_jd:.4f}\n" if np.isfinite(t_obs_jd) else "Obs JD      : N/A\n")
        if np.isfinite(t0_jd):
            try:
                t0_iso = Time(t0_jd, format='jd', scale='utc').iso
                f.write(f"Trigger T_0 : {t0_jd:.4f} ({t0_iso} UT)\n")
            except Exception:
                f.write(f"Trigger T_0 : {t0_jd:.4f}\n")
            if np.isfinite(delta_t_days):
                f.write(f"T - T_0     : {delta_t_days:.4f} days ({delta_t_days * 24.0:.2f} hours)\n")
            else:
                f.write(f"T - T_0     : N/A\n")
        else:
            f.write(f"Trigger T_0 : N/A\n")
            f.write(f"T - T_0     : N/A\n")
        f.write(f"CCD Gain    : {GAIN:.2f} e⁻/ADU\n")
        f.write(f"CCD ReadN   : {READ_NOISE:.2f} e⁻\n\n")

        with fits.open(image_path) as hdul:
            wcs = get_wcs(hdul[0].header, ASTROMETRY_REF_IMAGE)
        pixel_scale_degrees = proj_plane_pixel_scales(wcs)[0]
        pixel_scale_arcsec = pixel_scale_degrees * 3600
        f.write(f"Pixel Scale : {pixel_scale_arcsec:.3f} arcsec/pix\n")
        parsed_log_ra, parsed_log_dec = parse_target_coords(target_ra, target_dec)
        hms_log_ra, dms_log_dec = format_hms_dms(parsed_log_ra, parsed_log_dec)
        f.write(f"GRB Target (GCN)  : RA={parsed_log_ra:.6f}°, Dec={parsed_log_dec:+.6f}° ({hms_log_ra}, {dms_log_dec})\n")
        if str(target_ra) != f"{parsed_log_ra:.6f}" or str(target_dec) != f"{parsed_log_dec:+.6f}":
            f.write(f"Input Target Coords: RA={target_ra}, Dec={target_dec}\n")
        
        pos_sum = final_results.get('positions_summary')
        if pos_sum:
            p_init = pos_sum.get('initial')
            p_ref  = pos_sum.get('refined')
            p_mat  = pos_sum.get('matched')
            if p_init:
                f.write(f"Exact Initial Pos : RA={p_init['ra_deg']:.6f}°, Dec={p_init['dec_deg']:+.6f}° ({p_init['ra_hms']}, {p_init['dec_dms']}) [X={p_init['pixel_x']:.2f}, Y={p_init['pixel_y']:.2f}]\n")
            if p_ref and np.isfinite(p_ref.get('ra_deg', np.nan)):
                f.write(f"Refined Centroid  : RA={p_ref['ra_deg']:.6f}°, Dec={p_ref['dec_deg']:+.6f}° ({p_ref['ra_hms']}, {p_ref['dec_dms']}) [X={p_ref['pixel_x']:.2f}, Y={p_ref['pixel_y']:.2f}, Δ={p_ref['offset_arcsec']:.2f}\"]\n")
            if p_mat and p_mat.get('detected'):
                f.write(f"Matched ALLSTAR   : RA={p_mat['ra_deg']:.6f}°, Dec={p_mat['dec_deg']:+.6f}° ({p_mat['ra_hms']}, {p_mat['dec_dms']}) [X={p_mat['pixel_x']:.2f}, Y={p_mat['pixel_y']:.2f}, ID#{p_mat['id']}, Δ={p_mat['offset_init_arcsec']:.2f}\"]\n")
        elif not np.isnan(pixel_x):
            recentered_ra, recentered_dec = wcs.pixel_to_world_values(pixel_x, pixel_y)
            recentered_hms, recentered_dms = format_hms_dms(recentered_ra, recentered_dec)
            f.write(f"Recentered        : RA={float(recentered_ra):.6f}°, Dec={float(recentered_dec):+.6f}° ({recentered_hms}, {recentered_dms}) [X={pixel_x:.2f}, Y={pixel_y:.2f}]\n")
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
        pos_sum = final_results.get('positions_summary')
        if target_ra is not None and target_dec is not None:
            if pos_sum:
                f.write(format_positions_summary_text(pos_sum, TARGET_NAME, FILTER_BAND, zp, zp_err, sigma_th) + "\n\n")
            elif detected and not np.isnan(cal_mag):
                f.write("Result      : DETECTED\n")
                f.write(f"Detect Type : {detect_type}\n")
                f.write(f"Sigma Thresh: {sigma_th}\n")
                f.write(f"FINAL MAG   : {cal_mag:.3f} ± {cal_mag_err:.3f} mag ({FILTER_BAND}-band)\n")
                f.write(f"Instrumental: {inst_mag:.3f} ± {inst_merr:.3f}\n")
                f.write(f"CHI / SHARP : {chi:.2f} / {sharp:.3f}\n")
            else:
                f.write("Result      : NON-DETECTION (Target not recovered by DAOPHOT)\n")

            if ap_raw:
                f.write("  --- Raw .ap row data ---\n")
                f.write(f"{ap_raw}\n\n")
            if als_raw:
                f.write("  --- Raw .als row data ---\n")
                f.write(f"{als_raw}\n\n")
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
    image_path            = STACKED_IMAGE,
    target_ra             = TARGET_RA,
    target_dec            = TARGET_DEC,
    target_name           = TARGET_NAME,
    filter_band           = FILTER_BAND,
    force_coords          = FORCE_COORDS,
    fwhm_override         = None,
    sky_override          = None,
    std_sky_override      = None,
    pipeline_mode         = PIPELINE_MODE,
    refine_target_coords  = False,
    initial_target_coords = None,
    refined_target_coords = None,
):
    """
    Run the full standalone DAOPHOT II pipeline.

    Parameters
    ----------
    image_path            : WCS-solved stacked FITS (absolute or relative path)
    target_ra/dec         : target celestial coordinates (degrees or sexagesimal HMS/DMS)
    target_name           : used in filenames and plot titles
    filter_band           : photometric band  ('r', 'V', 'B', …)
    force_coords          : (ra, dec) tuple or None — enable forced photometry
    fwhm_override         : pass your photutil FWHM to skip D0
    sky_override          : pass your photutil sky mode to skip D0
    std_sky_override      : pass your photutil sky std to skip D0
    refine_target_coords  : boolean — refine target position via centroid search before DAOPHOT
    initial_target_coords : (ra, dec) tuple of original GCN / user coordinates
    refined_target_coords : (ra, dec) tuple of pre-refined centroid coordinates

    Returns
    -------
    dict with: detected, cal_mag, cal_mag_err, zp, zp_err,
               als_table, target_info, fwhm, sky, std_sky, positions_summary
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
        header          = hdul[0].header
        image_shape     = hdul[0].data.shape
        clean_fits_data = hdul[0].data.astype(float)
    wcs = get_wcs(header, ASTROMETRY_REF_IMAGE)

    # ── Parse Initial & Refined Target Coordinates ───────────────────────
    init_ra, init_dec = parse_target_coords(target_ra, target_dec)
    if initial_target_coords is not None:
        init_ra, init_dec = parse_target_coords(*initial_target_coords)

    init_x, init_y = np.nan, np.nan
    init_ra_hms, init_dec_dms = "N/A", "N/A"
    if init_ra is not None and init_dec is not None:
        init_x, init_y = wcs.wcs_world2pix(init_ra, init_dec, 0)
        init_x, init_y = float(init_x), float(init_y)
        init_ra_hms, init_dec_dms = format_hms_dms(init_ra, init_dec)

    ref_ra, ref_dec = init_ra, init_dec
    ref_x, ref_y = init_x, init_y
    ref_ra_hms, ref_dec_dms = init_ra_hms, init_dec_dms
    ref_offset_as = 0.0
    ref_offset_px = 0.0

    if refined_target_coords is not None:
        ref_ra, ref_dec = parse_target_coords(*refined_target_coords)
        if ref_ra is not None and ref_dec is not None:
            ref_x, ref_y = wcs.wcs_world2pix(ref_ra, ref_dec, 0)
            ref_x, ref_y = float(ref_x), float(ref_y)
            ref_ra_hms, ref_dec_dms = format_hms_dms(ref_ra, ref_dec)
            if np.isfinite(init_x) and np.isfinite(ref_x):
                ref_offset_px = float(np.hypot(ref_x - init_x, ref_y - init_y))
                ref_offset_as = float(SkyCoord(init_ra*u.deg, init_dec*u.deg).separation(SkyCoord(ref_ra*u.deg, ref_dec*u.deg)).arcsec)
    elif init_ra is not None and init_dec is not None:
        ref_check = confirm_and_refine_target(
            image_path           = image_path,
            target_ra            = init_ra,
            target_dec           = init_dec,
            target_name          = target_name,
            ref_image            = ASTROMETRY_REF_IMAGE,
            plots_dir            = PLOTS_DIR,
            show_plot            = False,
            save_plot            = True,
        )
        if ref_check and ref_check.get('detected'):
            ref_ra        = ref_check['refined_ra']
            ref_dec       = ref_check['refined_dec']
            ref_x         = ref_check['refined_x']
            ref_y         = ref_check['refined_y']
            ref_offset_as = ref_check['offset_arcsec']
            ref_offset_px = ref_check['offset_pix']
            ref_ra_hms, ref_dec_dms = format_hms_dms(ref_ra, ref_dec)

    search_ra, search_dec = ref_ra, ref_dec

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
    if search_ra is not None and search_dec is not None:
        target_info = lookup_target_in_als(
            als_table, wcs, search_ra, search_dec
        )
        if not target_info['detected'] and (search_ra != init_ra or search_dec != init_dec):
            info_init = lookup_target_in_als(als_table, wcs, init_ra, init_dec)
            if info_init['detected']:
                target_info = info_init

    # ── D6 : Forced photometry (Fallback) ─────────────────────────────
    final_als_table = als_table
    is_forced = False

    if search_ra is not None and search_dec is not None:
        if not target_info['detected'] or force_coords is not None:
            is_forced = True
            if not target_info["detected"] and force_coords is None:
                print("\n  ⚠️ Target not initially detected. Falling back to forced photometry...")
                frc_ra, frc_dec = search_ra, search_dec
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
                    final_als_table, wcs, search_ra, search_dec
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

    # ── Multi-Position Photometry (Exact Initial / Refined / Matched) ───
    positions_summary = None
    if init_ra is not None and init_dec is not None:
        eff_gain = float(header.get('GAIN', GAIN * FRAMES_AVG))
        eff_rn   = float(header.get('RDNOISE', header.get('READNOI', READ_NOISE / np.sqrt(max(1, FRAMES_AVG)))))
        phot_exact = measure_aperture_photometry_at_pix(clean_fits_data, init_x, init_y, fwhm, zp, zp_err, gain=eff_gain, readnoise=eff_rn)
        phot_ref   = measure_aperture_photometry_at_pix(clean_fits_data, ref_x, ref_y, fwhm, zp, zp_err, gain=eff_gain, readnoise=eff_rn)
        
        match_detected = target_info.get('detected', False) if target_info else False
        if match_detected:
            m_x = float(target_info['pixel_x'])
            m_y = float(target_info['pixel_y'])
            m_ra = float(target_info.get('matched_ra', np.nan))
            m_dec = float(target_info.get('matched_dec', np.nan))
            m_hms = target_info.get('matched_hms', 'N/A')
            m_dms = target_info.get('matched_dms', 'N/A')
            if np.isnan(m_ra):
                m_ra, m_dec = wcs.wcs_pix2world(m_x, m_y, 0)
                m_ra, m_dec = float(m_ra), float(m_dec)
                m_hms, m_dms = format_hms_dms(m_ra, m_dec)
            m_off_init_as = float(SkyCoord(init_ra*u.deg, init_dec*u.deg).separation(SkyCoord(m_ra*u.deg, m_dec*u.deg)).arcsec)
            m_off_init_px = float(np.hypot(m_x - init_x, m_y - init_y))
            m_off_ref_as  = float(SkyCoord(ref_ra*u.deg, ref_dec*u.deg).separation(SkyCoord(m_ra*u.deg, m_dec*u.deg)).arcsec)
            m_off_ref_px  = float(np.hypot(m_x - ref_x, m_y - ref_y))
            phot_match = measure_aperture_photometry_at_pix(clean_fits_data, m_x, m_y, fwhm, zp, zp_err, gain=eff_gain, readnoise=eff_rn)
        else:
            m_x, m_y = np.nan, np.nan
            m_ra, m_dec = np.nan, np.nan
            m_hms, m_dms = "N/A", "N/A"
            m_off_init_as, m_off_init_px = np.nan, np.nan
            m_off_ref_as, m_off_ref_px = np.nan, np.nan
            phot_match = None

        positions_summary = {
            'initial': {
                'desc': 'Exact Initial Location (GCN)',
                'ra_deg': init_ra, 'dec_deg': init_dec,
                'ra_hms': init_ra_hms, 'dec_dms': init_dec_dms,
                'pixel_x': init_x, 'pixel_y': init_y,
                'phot': phot_exact,
            },
            'refined': {
                'desc': 'Refined Centroid Location (FITS)',
                'ra_deg': ref_ra, 'dec_deg': ref_dec,
                'ra_hms': ref_ra_hms, 'dec_dms': ref_dec_dms,
                'pixel_x': ref_x, 'pixel_y': ref_y,
                'offset_arcsec': ref_offset_as, 'offset_pix': ref_offset_px,
                'phot': phot_ref,
            },
            'matched': {
                'desc': 'Matched ALLSTAR Detection (PSF Fit)',
                'detected': match_detected,
                'id': target_info.get('id', -1) if target_info else -1,
                'ra_deg': m_ra, 'dec_deg': m_dec,
                'ra_hms': m_hms, 'dec_dms': m_dms,
                'pixel_x': m_x, 'pixel_y': m_y,
                'offset_init_arcsec': m_off_init_as, 'offset_init_pix': m_off_init_px,
                'offset_ref_arcsec': m_off_ref_as, 'offset_ref_pix': m_off_ref_px,
                'psf_inst_mag': target_info.get('mag', np.nan) if target_info else np.nan,
                'psf_inst_merr': target_info.get('merr', np.nan) if target_info else np.nan,
                'psf_cal_mag': final.get('cal_mag', np.nan),
                'psf_cal_mag_err': final.get('cal_mag_err', np.nan),
                'chi': target_info.get('chi', np.nan) if target_info else np.nan,
                'sharp': target_info.get('sharp', np.nan) if target_info else np.nan,
                'detect_type': "Forced Photometry" if is_forced else "Normal Intrinsic Find",
                'phot': phot_match,
            }
        }

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

    if init_ra is not None and init_dec is not None:
        px = final.get('pixel_x', np.nan)
        py = final.get('pixel_y', np.nan)
        if np.isnan(px):
            px, py = ref_x, ref_y
        plot_target_confirmation(image_path, px, py, target_name, init_x=init_x, init_y=init_y)
    
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
        ap_f = workpath(f"{get_stem(image_bn)}_merged.ap") if is_forced else workpath(f"{get_stem(image_bn)}.ap")
        als_f = workpath(f"{get_stem(image_bn)}_merged.als") if is_forced else workpath(f"{get_stem(image_bn)}.als")
        
        ap_raw = extract_raw_target_data(ap_f, target_id)
        als_raw = extract_raw_target_data(als_f, target_id)

    # ── Summary Output ──────────────────────────────────────────────────
    print("\n" + "="*60)
    print("  RESULT")
    print("="*60)
    
    if init_ra is not None and init_dec is not None:
        if positions_summary is not None:
            print("\n" + format_positions_summary_text(positions_summary, target_name, filter_band, zp, zp_err, sigma_th))
            if ap_raw:
                print("\n  --- Raw .ap row data ---")
                print(ap_raw)
            if als_raw:
                print("\n  --- Raw .als row data ---")
                print(als_raw)
        elif final['detected'] and not np.isnan(final.get('cal_mag', np.nan)):
            print(f"  Target      : {target_name}")
            print(f"  Calibrated  : {final['cal_mag']:.3f} ± {final['cal_mag_err']:.3f}  ({filter_band}-band)")
            print(f"  Instrumental: {final['mag']:.3f} ± {final['merr']:.3f}")
            print(f"  Zero-point  : {zp:.4f} ± {zp_err:.4f}")
            print(f"  CHI / SHARP : {final['chi']:.2f} / {final['sharp']:.3f}")
            print(f"  Pixel pos   : X={final['pixel_x']:.2f}  Y={final['pixel_y']:.2f}  (dist={final['dist_pix']:.2f}px)")
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
        'detected'         : final.get('detected', False),
        'cal_mag'          : final.get('cal_mag', np.nan),
        'cal_mag_err'      : final.get('cal_mag_err', np.nan),
        'inst_mag'         : final.get('mag', np.nan),
        'inst_merr'        : final.get('merr', np.nan),
        'zp'               : zp,
        'zp_err'           : zp_err,
        'chi'              : final.get('chi', np.nan),
        'sharp'            : final.get('sharp', np.nan),
        'pixel_x'          : final.get('pixel_x', np.nan),
        'pixel_y'          : final.get('pixel_y', np.nan),
        'id'               : final.get('id', -1),
        'dist_pix'         : final.get('dist_pix', np.nan),
        'fwhm'             : fwhm,
        'sky'              : sky,
        'std_sky'          : std_sky,
        'als_table'        : final_als_table if 'final_als_table' in locals() else None,
        'zp_table'         : zp_table if 'zp_table' in locals() else None,
        'all_sources_calibrated': all_sources_calibrated,
        'filtered_stars_calibrated': filtered_stars_calibrated,
        'detect_type'      : "Forced Photometry" if is_forced else "Normal Intrinsic Find",
        'sigma_th'         : sigma_th,
        'ap_raw'           : ap_raw,
        'als_raw'          : als_raw,
        'positions_summary': positions_summary,
    }

    write_comprehensive_log(final_output_dict, pipeline_mode, image_path, init_ra, init_dec)

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
# CELL 0 : TARGET PRE-CONFIRMATION & WCS REFINEMENT
# =====================================================================
# 1. Click the online archival links (Legacy Survey, Pan-STARRS, Aladin).
# 2. Inspect the 3-panel side-by-side comparison:
#      Panel 1: Your FITS File (Wide Context)
#      Panel 2: Online Catalogue (Pre-Burst Archival: Pan-STARRS1)
#      Panel 3: Your FITS File (Zoomed Frame) with faint transparent marker (alpha=0.3)
# 3. Set ADOPT_REFINED_COORDS = False to keep your original GCN coordinates.

# Save initial coordinates before any potential refinement update
INITIAL_TARGET_RA  = TARGET_RA
INITIAL_TARGET_DEC = TARGET_DEC
REFINED_TARGET_RA  = None
REFINED_TARGET_DEC = None
REFINE_RESULT      = None

ADOPT_REFINED_COORDS = True   # Set to False if you want to keep initial GCN coordinates!

if __name__ == '__main__':
    if TARGET_RA is not None and TARGET_DEC is not None:
        refine_result = confirm_and_refine_target(
            image_path           = STACKED_IMAGE,
            target_ra            = TARGET_RA,
            target_dec           = TARGET_DEC,
            target_name          = TARGET_NAME,
            search_radius_arcsec = 5.0,     # Search radius around GCN position (arcsec)
            cutout_size_arcsec   = 40.0,    # Field of view for context (arcsec)
            selected_candidate   = 1,       # 1 = closest peak, or choose #2, #3 from candidate list
            show_plot            = True,
            save_plot            = True,
        )
        REFINE_RESULT = refine_result
        if refine_result and refine_result.get('refined_ra') is not None:
            REFINED_TARGET_RA  = refine_result['refined_ra']
            REFINED_TARGET_DEC = refine_result['refined_dec']

        if ADOPT_REFINED_COORDS and REFINED_TARGET_RA is not None:
            TARGET_RA  = REFINED_TARGET_RA
            TARGET_DEC = REFINED_TARGET_DEC
            print(f"  [Updated Target] TARGET_RA  = {TARGET_RA:.6f}°")
            print(f"  [Updated Target] TARGET_DEC = {TARGET_DEC:+.6f}°")
        else:
            print(f"  [Preserved Initial Target] TARGET_RA  = {TARGET_RA}")
            print(f"  [Preserved Initial Target] TARGET_DEC = {TARGET_DEC}")

# %%
# =====================================================================
# CELL 1: MANUAL WORKSPACE CONFIGURATION
# =====================================================================
manual_daophot_done = True   # Set to True after you finish the terminal steps

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
                image_path            = STACKED_IMAGE,
                target_ra             = TARGET_RA,
                target_dec            = TARGET_DEC,
                target_name           = TARGET_NAME,
                filter_band           = FILTER_BAND,
                pipeline_mode         = PIPELINE_MODE,
                initial_target_coords = (INITIAL_TARGET_RA, INITIAL_TARGET_DEC) if 'INITIAL_TARGET_RA' in locals() else None,
                refined_target_coords = (REFINED_TARGET_RA, REFINED_TARGET_DEC) if 'REFINED_TARGET_RA' in locals() and REFINED_TARGET_RA is not None else None,
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
            image_path            = STACKED_IMAGE,
            target_ra             = TARGET_RA,
            target_dec            = TARGET_DEC,
            target_name           = TARGET_NAME,
            filter_band           = FILTER_BAND,
            pipeline_mode         = 'auto',
            initial_target_coords = (INITIAL_TARGET_RA, INITIAL_TARGET_DEC) if 'INITIAL_TARGET_RA' in locals() else None,
            refined_target_coords = (REFINED_TARGET_RA, REFINED_TARGET_DEC) if 'REFINED_TARGET_RA' in locals() and REFINED_TARGET_RA is not None else None,
        )
        print("\nFinal Calibrated Magnitude:", final_results.get('cal_mag'))
    else:
        print("Skipping FULLY AUTOMATED RUN because PIPELINE_MODE is 'manual' or 'semi-auto'.")
