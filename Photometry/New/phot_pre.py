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
#     display_name: PyRAF
#     language: python
#     name: pyraf
# ---

# %% [markdown]
# # 🔭 CCD Photometric Reduction Pipeline (PyRAF / IRAF)
# ### Unified notebook — handles both 3-D FITS cubes *and* pre-unpacked science frames
#
# **Pipeline steps:**
# 1. ⚙️  Configuration  
# 2. 📦  Package setup & output-directory creation  
# 3. 🔪  (Optional) Unpack 3-D FITS cubes with `imslice`  
# 4. 🧮  Master Bias creation  
# 5. 🔧  Bias-correction → Science frames  
# 6. 🔧  Bias-correction → Flat frames  
# 7. 🎛️  Master Flat creation  
# 8. 📐  Master Flat normalisation  
# 9. ✨  Flat-field correction → Science frames  
#
# ---

# %% [markdown]
# ##  Cell 1 — Configuration
# > **Edit only this cell.** All other cells read from these variables.

# %%
# ╔══════════════════════════════════════════════════════════════╗
# ║              USER CONFIGURATION — EDIT HERE ONLY            ║
# ╚══════════════════════════════════════════════════════════════╝

import os

# ── 1. DATA LOCATIONS ──────────────────────────────────────────
# Directory that contains ALL your raw .fits files.
# Use '.' for the current working directory.
RAW_DATA_DIR = '.'

# ── 2. OUTPUT DIRECTORY ROOT ───────────────────────────────────
# Sub-folders  bias/ | flats/ | science/  will be created here.
OUTPUT_DIR = './reduced'

# ── 3. INPUT MODE ──────────────────────────────────────────────
# 'cube'      → raw files are 3-D FITS cubes; imslice will unpack them.
# 'unstacked' → files are already individual 2-D science frames.
SCIENCE_INPUT_MODE = 'cube'    # <-- change to 'unstacked' if already sliced

# ── 4. FILE GLOB PATTERNS ──────────────────────────────────────
# Patterns are relative to RAW_DATA_DIR (.fit or .fits supported).
BIAS_PATTERN       = 'Bias_*.fits*'      # raw bias frames 
FLAT_PATTERN       = 'FLAT_R*.fits*'     # raw flat frames
# For cubes: pattern that selects the 3-D cube files to slice
CUBE_PATTERN       = 'GRB_R_*.fits*'
# For unstacked: pattern that selects the already-sliced science frames
UNSTACKED_PATTERN  = 'GRB261010A_P10_r_*.fit*' 

# ── 5. INSTRUMENT / DETECTOR ───────────────────────────────────
RDNOISE  = 7.0    # e⁻  — read noise of your detector
GAIN     = 2.0   # e⁻/ADU — gain of your detector

# ── 6. ZEROCOMBINE PARAMETERS ──────────────────────────────────
BIAS_COMBINE  = 'average'   # 'average' | 'median' | 'sum'
BIAS_REJECT   = 'minmax'    # 'minmax'  | 'sigclip' | 'avsigclip'
BIAS_NLOW     = 0
BIAS_NHIGH    = 1

# ── 7. FLATCOMBINE PARAMETERS ──────────────────────────────────
FLAT_COMBINE  = 'median'    # 'median'  | 'average' | 'sum'
FLAT_REJECT   = 'avsigclip' # 'avsigclip' | 'sigclip' | 'minmax'
FLAT_LSIGMA   = 3.0
FLAT_HSIGMA   = 3.0
FLAT_SCALE    = 'mode'      # 'mode' | 'mean' | 'median' | 'none'

# ── 8. FILE PREFIX CONVENTIONS ─────────────────────────────────
PREFIX_BIAS_CORR  = 'b_'    # bias-subtracted files
PREFIX_FLAT_CORR  = 'f'     # flat-fielded files  (stacks with PREFIX_BIAS_CORR → fb_)

# ── 9. OUTPUT FILE NAMES ───────────────────────────────────────
MASTER_BIAS_NAME      = 'MasterBias.fits'
MASTER_FLAT_NAME      = 'MasterFlat.fits'
NORM_MASTER_FLAT_NAME = 'NormMasterFlat.fits'
BPM_NAME              = 'BPM.fits'          # Bad pixel mask output name

# ── 10. IMSLICE SETTINGS (cube mode only) ──────────────────────
SLICE_DIMENSION = 3   # 3 = slice along Z-axis (standard for data cubes)

# ── 11. BAD PIXEL & 2D ILLUMINATION CORRECTION ─────────────────
# WHEN TO USE:
# • Keep APPLY_ILLUM_FLATTEN = True (Recommended for all optical datasets).
#   Removes asymmetric 2D sky tilts, airmass gradients, and residual vignetting
#   (e.g., +5% corner brightnings) while strictly preserving star counts and true flux.
# • ILLUM_BOX_SIZE: 64 px is standard for 2048x2048 detectors. Use 32-48 px for
#   steep gradients; use 128 px for very smooth fields.
BPM_FLAT_THRESH       = 0.35   # Normalized flat threshold; pixels below this are dead/masked
APPLY_ILLUM_FLATTEN   = True   # Flatten residual 2D sky gradients & vignetting tilts
ILLUM_BOX_SIZE        = 64     # Mesh block size (pixels) for 2D background surface

# ── 12. DONUT & ARTIFACT MASKING IN BPM ────────────────────────
# WHEN TO USE:
# • Use when dust shifted between flats and science (due to filter wheel micro-
#   rotation or gravity sag), leaving uncorrected bright/dark donut rings.
# • Daophot and PyRAF will ignore these masked pixels during PSF and aperture fitting.
# • Format: (x_center, y_center, radius_pixels). In FITS, (1,1) is bottom-left.
DONUT_MASK_CIRCLES    = [
    (659, 1072, 28),   # Prominent dust donut in field
    (1534, 647, 22),   # Secondary donut
]
AUTO_DETECT_DONUTS    = True   # Automatically detect strong ring anomalies in flat
DONUT_SIGMA_THRESH    = 4.5    # Threshold (sigma) for auto donut detection

# ── 13. NIGHT-SKY SUPERFLAT (ADVANCED) ─────────────────────────
# WHEN TO USE:
# • ONLY if you have >= 15-30 WIDELY DITHERED science frames in a SPARSE field.
# • In near-IR (J/H/K) where night-sky background dominates over twilight.
# WHEN NOT TO USE (WARNING):
# • DO NOT USE if you have few frames (< 10) or UNDITHERED frames (dither < 20 px).
# • Why? Stacking few science frames injects high Poisson noise into your images
#   and imprints unrejected faint "ghost stars" as artificial dark potholes.
# • For small datasets (e.g. 2-6 frames), keep ENABLE_SUPERFLAT = False and use
#   standard twilight flats + 2D illumination flattening (Section 11).
ENABLE_SUPERFLAT          = False  # Set True ONLY if >= 15 dithered frames available
SUPERFLAT_MIN_FRAMES      = 5      # Minimum science frames required (>=5 recommended)
USE_SUPERFLAT_FOR_SCIENCE = False  # Set True to flat-field with Superflat instead of twilight flat
SUPERFLAT_NAME            = 'MasterSuperFlat.fits'

# ── 14. EDGE TRIMMING PARAMETERS (Pre-Alignment) ───────────────
# WHEN TO USE:
# • Set > 0 if detector has overscan columns/rows, dead borders, or severe edge
#   rolloff that should be permanently cut before image registration/stacking.
TRIM_INPUT_WILDCARD = "fb_*.fit*"
TRIM_BOTTOM         = 0  # Pixels to cut from bottom (white overscan strip: rows 0..99)
TRIM_LEFT           = 0    # Pixels to cut from left edge (set > 0 if you want to trim left)
TRIM_TOP            = 0    # Pixels to cut from top edge (0 = keep all)
TRIM_RIGHT          = 0  # Pixels to cut from right edge (overscan strip: cols 2048..2147; set 0 if not needed)
TRIM_OUTPUT_PREFIX  = "t"

# ╔══════════════════════════════════════════════════════════════╗
# ║                 END OF USER CONFIGURATION                    ║
# ╚══════════════════════════════════════════════════════════════╝

# ── Derived absolute paths (do NOT edit below this line) ───────
RAW_DATA_DIR = os.path.abspath(RAW_DATA_DIR)

DIR_BIAS    = os.path.abspath(os.path.join(OUTPUT_DIR, 'bias'))
DIR_FLATS   = os.path.abspath(os.path.join(OUTPUT_DIR, 'flats'))
DIR_SCIENCE = os.path.abspath(os.path.join(OUTPUT_DIR, 'science'))

MASTER_BIAS_PATH      = os.path.join(DIR_BIAS,  MASTER_BIAS_NAME)
MASTER_FLAT_PATH      = os.path.join(DIR_FLATS, MASTER_FLAT_NAME)
NORM_MASTER_FLAT_PATH = os.path.join(DIR_FLATS, NORM_MASTER_FLAT_NAME)
SUPERFLAT_PATH        = os.path.join(DIR_FLATS, SUPERFLAT_NAME)
BPM_PATH              = os.path.join(DIR_FLATS, BPM_NAME)
BPM_SCIENCE_PATH      = os.path.join(DIR_SCIENCE, BPM_NAME)

print("✅ Configuration loaded.")
print(f"   RAW_DATA_DIR  : {RAW_DATA_DIR}")
print(f"   OUTPUT_DIR    : {os.path.abspath(OUTPUT_DIR)}")
print(f"   Mode          : {SCIENCE_INPUT_MODE}")
print(f"   Bias dir      : {DIR_BIAS}")
print(f"   Flats dir     : {DIR_FLATS}")
print(f"   Science dir   : {DIR_SCIENCE}")
print(f"   BPM Path      : {BPM_PATH}")
print(f"   Superflat     : {SUPERFLAT_PATH if ENABLE_SUPERFLAT else '(disabled)'}")


# %% [markdown]
# ##  Cell 2 — Package Setup & Output Directories

# %%
import glob
import shutil
from pyraf import iraf

# ── Load IRAF packages ─────────────────────────────────────────
iraf.noao()
iraf.imred()
iraf.ccdred()
iraf.images()
iraf.imutil()
print("✅ IRAF packages loaded.")

# ── Create output directories ──────────────────────────────────
for label, path in [('bias', DIR_BIAS), ('flats', DIR_FLATS), ('science', DIR_SCIENCE)]:
    os.makedirs(path, exist_ok=True)
    print(f"   📁 {label:8s} → {path}")

print("\n✅ Output directories ready.")

# ── Helper: resolve glob relative to RAW_DATA_DIR ──────────────
def find_files(pattern, base=RAW_DATA_DIR):
    """Glob *pattern* inside *base* and return a sorted list of absolute paths.
    Supports both .fit and .fits extensions."""
    results = sorted(glob.glob(os.path.join(base, pattern)))
    if not results and (pattern.endswith('.fit') or pattern.endswith('.fits')):
        alt_pattern = os.path.splitext(pattern)[0] + '.fit*'
        results = sorted(glob.glob(os.path.join(base, alt_pattern)))
    return results

# ── Helper: write an IRAF @-list file ──────────────────────────
def write_list(paths, list_path):
    """Write *paths* (one per line) to *list_path* for IRAF @-file syntax."""
    with open(list_path, 'w') as fh:
        for p in paths:
            fh.write(p + '\n')
    print(f"   📝 List written: {list_path}  ({len(paths)} entries)")



# %% [markdown]
# ##  Cell 3 — (Optional) Unpack 3-D FITS Cubes with `imslice`
# > **Skipped automatically** when `SCIENCE_INPUT_MODE = 'unstacked'`.

# %%
# This cell only runs when SCIENCE_INPUT_MODE == 'cube'.
# Sliced frames are written directly into DIR_SCIENCE so that
# subsequent cells always look in the same place regardless of mode.

if SCIENCE_INPUT_MODE == 'cube':
    cube_files = find_files(CUBE_PATTERN)
    if not cube_files:
        print(f"⚠️  No cube files matched pattern '{CUBE_PATTERN}' in {RAW_DATA_DIR}")
    else:
        print(f"Found {len(cube_files)} cube file(s) to slice:\n")
        for cube in cube_files:
            # e.g.  /data/01_12.fits  →  <DIR_SCIENCE>/01_12_slice
            base     = os.path.splitext(os.path.basename(cube))[0]
            out_root = os.path.join(DIR_SCIENCE, base + '_slice')
            print(f"  Slicing  {os.path.basename(cube)}  →  {out_root}NNN.fits")
            try:
                iraf.imslice(
                    input          = cube,
                    output         = out_root,
                    slice_dimension= SLICE_DIMENSION,
                    verbose        = 'yes',
                    mode           = 'h'
                )
            except iraf.IrafError as exc:
                print(f"  ⚠️  imslice failed for {cube}:\n     {exc}")
        print("\n✅ All cubes sliced into:", DIR_SCIENCE)

elif SCIENCE_INPUT_MODE == 'unstacked':
    print("ℹ️  Mode = 'unstacked' — imslice step skipped.")
    print("   Pre-sliced science frames will be read from RAW_DATA_DIR in the next steps.")
else:
    raise ValueError(f"SCIENCE_INPUT_MODE must be 'cube' or 'unstacked', got '{SCIENCE_INPUT_MODE}'")


# %% [markdown]
# ##  Cell 4 — Master Bias Creation

# %%
bias_files = find_files(BIAS_PATTERN)

if not bias_files:
    raise FileNotFoundError(
        f"No bias frames found matching '{BIAS_PATTERN}' in {RAW_DATA_DIR}\n"
        "  ↳ Check BIAS_PATTERN and RAW_DATA_DIR in the Configuration cell."
    )

print(f"Found {len(bias_files)} bias frame(s).")

# Write @-list into the bias output directory
bias_list_path = os.path.join(DIR_BIAS, 'bias.list')
write_list(bias_files, bias_list_path)

# Remove stale output
if os.path.exists(MASTER_BIAS_PATH):
    os.remove(MASTER_BIAS_PATH)
    print(f"   🗑️  Removed existing {MASTER_BIAS_PATH}")

iraf.unlearn('zerocombine')
try:
    iraf.zerocombine(
        input   = f'@{bias_list_path}',
        output  = MASTER_BIAS_PATH,
        combine = BIAS_COMBINE,
        reject  = BIAS_REJECT,
        ccdtype = '',
        process = 'no',
        delete  = 'no',
        clobber = 'no',
        scale   = 'none',
        statsec = '',
        nlow    = BIAS_NLOW,
        nhigh   = BIAS_NHIGH,
        nkeep   = 1,
        mclip   = 'yes',
        lsigma  = 3.0,
        hsigma  = 3.0,
        rdnoise = RDNOISE,
        gain    = GAIN
    )
    print(f"\n✅ Master Bias saved → {MASTER_BIAS_PATH}")
except iraf.IrafError as exc:
    print(f"❌ zerocombine failed:\n   {exc}")


# %% [markdown]
# ##  Cell 5 — Bias Correction: Science Frames

# %%
# Locate science frames depending on input mode
if SCIENCE_INPUT_MODE == 'cube':
    # Slices were written to DIR_SCIENCE in Cell 3
    science_files = sorted(glob.glob(os.path.join(DIR_SCIENCE, '*_slice*.fit*')))
    science_files = [f for f in science_files if not os.path.basename(f).startswith((PREFIX_BIAS_CORR, PREFIX_FLAT_CORR))]
else:  # unstacked
    science_files = find_files(UNSTACKED_PATTERN)

if not science_files:
    raise FileNotFoundError(
        "No science frames found.\n"
        f"  Mode = '{SCIENCE_INPUT_MODE}'\n"
        f"  Searched: {'DIR_SCIENCE' if SCIENCE_INPUT_MODE == 'cube' else 'RAW_DATA_DIR'}\n"
        "  ↳ Check UNSTACKED_PATTERN / CUBE_PATTERN in the Configuration cell."
    )

print(f"Found {len(science_files)} science frame(s) to bias-correct.\n")
iraf.unlearn('ccdproc')

for src in science_files:
    bname  = PREFIX_BIAS_CORR + os.path.basename(src)
    dest   = os.path.join(DIR_SCIENCE, bname)
    if os.path.exists(dest):
        os.remove(dest)
    print(f"  {os.path.basename(src):40s}  →  {bname}")
    try:
        iraf.ccdproc(
            images   = src,
            output   = dest,
            ccdtype  = '',
            interactive = 'no',
            overscan = 'no',
            trim     = 'no',
            fixpix   = 'no',
            zerocor  = 'yes',        # ← Bias subtraction ON
            darkcor  = 'no',
            flatcor  = 'no',
            illumcor = 'no',
            fringecor= 'no',
            readcor  = 'no',
            scancor  = 'no',
            zero     = MASTER_BIAS_PATH
        )
    except iraf.IrafError as exc:
        print(f"  ⚠️  ccdproc failed for {src}:\n     {exc}")

print(f"\n✅ Bias-corrected science frames saved → {DIR_SCIENCE}")


# %% [markdown]
# ##  Cell 6 — Bias Correction: Flat Frames

# %%
flat_files = find_files(FLAT_PATTERN)

if not flat_files:
    raise FileNotFoundError(
        f"No flat frames found matching '{FLAT_PATTERN}' in {RAW_DATA_DIR}\n"
        "  ↳ Check FLAT_PATTERN and RAW_DATA_DIR in the Configuration cell."
    )

print(f"Found {len(flat_files)} flat frame(s) to bias-correct.\n")
iraf.unlearn('ccdproc')

for src in flat_files:
    bname = PREFIX_BIAS_CORR + os.path.basename(src)
    dest  = os.path.join(DIR_FLATS, bname)
    if os.path.exists(dest):
        os.remove(dest)
    print(f"  {os.path.basename(src):40s}  →  {bname}")
    try:
        iraf.ccdproc(
            images   = src,
            output   = dest,
            ccdtype  = '',
            interactive = 'no',
            overscan = 'no',
            trim     = 'no',
            fixpix   = 'no',
            zerocor  = 'yes',        # ← Bias subtraction ON
            darkcor  = 'no',
            flatcor  = 'no',
            illumcor = 'no',
            fringecor= 'no',
            readcor  = 'no',
            scancor  = 'no',
            zero     = MASTER_BIAS_PATH
        )
    except iraf.IrafError as exc:
        print(f"  ⚠️  ccdproc failed for {src}:\n     {exc}")

print(f"\n✅ Bias-corrected flat frames saved → {DIR_FLATS}")


# %% [markdown]
# ##  Cell 7 — Master Flat Creation

# %%
flat_stem = os.path.splitext(FLAT_PATTERN)[0].lstrip('*')
ext = '.fit*' if os.path.splitext(FLAT_PATTERN)[1].lower().startswith('.fit') else '*.fit*'
b_flat_pattern = PREFIX_BIAS_CORR + (flat_stem if (not flat_stem or flat_stem.endswith('*')) else flat_stem + '*') + ext
b_flat_files   = sorted(glob.glob(os.path.join(DIR_FLATS, b_flat_pattern)))

# Filter to only bias-corrected flats (exclude any master files)
b_flat_files = [f for f in b_flat_files
                if MASTER_FLAT_NAME not in f and NORM_MASTER_FLAT_NAME not in f]

if not b_flat_files:
    raise FileNotFoundError(
        f"No bias-corrected flat frames found in {DIR_FLATS}\n"
        "  ↳ Make sure Cell 6 completed without errors."
    )

print(f"Found {len(b_flat_files)} bias-corrected flat frame(s) to combine.\n")

flat_list_path = os.path.join(DIR_FLATS, 'flat.list')
write_list(b_flat_files, flat_list_path)

if os.path.exists(MASTER_FLAT_PATH):
    os.remove(MASTER_FLAT_PATH)
    print(f"   🗑️  Removed existing {MASTER_FLAT_PATH}")

# ── Quality report on raw flats before combine ─────────────────
from astropy.io import fits
import numpy as np

flat_meds = []
for f in b_flat_files:
    with fits.open(f) as h:
        m = float(np.nanmedian(h[0].data))
        flat_meds.append(m)
        print(f"   Flat: {os.path.basename(f):35s}  median = {m:8.1f} ADU")

if flat_meds:
    group_med = float(np.median(flat_meds))
    for f, m in zip(b_flat_files, flat_meds):
        dev_pct = (m / group_med - 1.0) * 100.0
        if abs(dev_pct) > 35.0:
            print(f"   ⚠️  WARNING: {os.path.basename(f)} deviates {dev_pct:+.1f}% from median level ({group_med:.1f})")

iraf.unlearn('flatcombine')
try:
    iraf.flatcombine(
        input   = f'@{flat_list_path}',
        output  = MASTER_FLAT_PATH,
        combine = FLAT_COMBINE,
        reject  = FLAT_REJECT,
        ccdtype = '',
        process = 'no',
        subsets = 'no',
        delete  = 'no',
        clobber = 'no',
        scale   = FLAT_SCALE,
        statsec = '',
        nlow    = 0,
        nhigh   = 0,
        nkeep   = 1,
        mclip   = 'yes',
        lsigma  = FLAT_LSIGMA,
        hsigma  = FLAT_HSIGMA,
        rdnoise = RDNOISE,
        gain    = GAIN,
        snoise  = '0.0',
        pclip   = '-0.5',
        blank   = '1.0',
        mode    = 'ql'
    )
    print(f"\n✅ Master Flat saved → {MASTER_FLAT_PATH}")
except iraf.IrafError as exc:
    print(f"❌ flatcombine failed:\n   {exc}")


# %% [markdown]
# ##  Cell 8 — Normalise Master Flat & Generate Bad Pixel Mask (BPM)
# > Divides the Master Flat by the central illuminated region median,
# > flags dead/cold pixels into BPM.fits, and repairs bad pixels via local median
# > to prevent division-by-near-zero spikes in science frames.

# %%
if os.path.exists(NORM_MASTER_FLAT_PATH):
    os.remove(NORM_MASTER_FLAT_PATH)
    print(f"   🗑️  Removed existing {NORM_MASTER_FLAT_PATH}")
if os.path.exists(BPM_PATH):
    os.remove(BPM_PATH)
    print(f"   🗑️  Removed existing {BPM_PATH}")

print(f"Normalising {MASTER_FLAT_PATH} & generating Bad Pixel Mask ...")

try:
    from astropy.io import fits
    import numpy as np
    from scipy import ndimage

    with fits.open(MASTER_FLAT_PATH) as hdul:
        flat_data = np.squeeze(hdul[0].data).astype(np.float64)
        ny, nx = flat_data.shape
        cy, cx = ny // 2, nx // 2

        # ── Normalise by CENTRE-BOX median (not global) ──────────
        CBOX = 200  # half-size of centre box
        centre_box = flat_data[cy - CBOX:cy + CBOX, cx - CBOX:cx + CBOX]
        norm_val = float(np.nanmedian(centre_box))
        print(f"   Master Flat centre-box median = {norm_val:.4f} ADU")
        print(f"   (Calculated on {2*CBOX}×{2*CBOX} box centred at ({cx},{cy}))")

        norm_data = flat_data / norm_val

        # ── Bad Pixel Mask (BPM) detection ───────────────────────
        # Flag pixels with dead/negative response, non-finite, or < BPM_FLAT_THRESH
        bpm = (norm_data < BPM_FLAT_THRESH) | (~np.isfinite(norm_data)) | (norm_data <= 0.0)
        n_bad = int(np.sum(bpm))
        print(f"   Bad pixels flagged (<{BPM_FLAT_THRESH} or non-positive): {n_bad} ({100*n_bad/norm_data.size:.3f}%)")

        # ── Heal bad pixels in NormMasterFlat via local median ───
        if n_bad > 0:
            local_med = ndimage.median_filter(norm_data, size=5)
            norm_data[bpm] = local_med[bpm]
            norm_data = np.clip(norm_data, BPM_FLAT_THRESH, None)
            print(f"   ✅ Healed {n_bad} bad pixel(s) in NormMasterFlat using 5×5 local median.")

        # ── Donut & Custom Defect Masking ─────────────────────────
        Y, X = np.ogrid[:ny, :nx]
        n_custom_donut_px = 0
        if DONUT_MASK_CIRCLES:
            for cx_d, cy_d, rad_d in DONUT_MASK_CIRCLES:
                circle_mask = ((X - cx_d)**2 + (Y - cy_d)**2) <= (rad_d**2)
                bpm[circle_mask] = 1
                n_custom_donut_px += int(np.sum(circle_mask))
            print(f"   ⭕ Masked {len(DONUT_MASK_CIRCLES)} custom donut region(s) in BPM: {n_custom_donut_px} px")

        # Auto-detect residual ring anomalies if enabled
        if AUTO_DETECT_DONUTS:
            smooth_flat = ndimage.uniform_filter(norm_data, size=51)
            flat_res = norm_data - smooth_flat
            res_sigma = float(np.nanstd(flat_res))
            auto_donut_mask = np.abs(flat_res) > (DONUT_SIGMA_THRESH * res_sigma)
            labeled_donuts, n_feat = ndimage.label(auto_donut_mask)
            if n_feat > 0:
                sizes = ndimage.sum(auto_donut_mask, labeled_donuts, range(1, n_feat + 1))
                auto_donut_px = 0
                for comp_id, sz in enumerate(sizes, start=1):
                    if sz >= 20:
                        bpm[labeled_donuts == comp_id] = 1
                        auto_donut_px += int(sz)
                if auto_donut_px > 0:
                    print(f"   🍩 Auto-detected & masked {auto_donut_px} donut/artifact px in BPM (>{DONUT_SIGMA_THRESH}σ)")

        # Save BPM.fits (1 = bad, 0 = good) to flats dir and sync to science dir
        bpm_hdr = hdul[0].header.copy()
        if 'NAXIS3' in bpm_hdr:
            del bpm_hdr['NAXIS3']
        bpm_hdr['NAXIS']   = 2
        bpm_hdr['BPM_SRC'] = 'MasterFlat threshold + Donut masks'
        bpm_hdr['THRESH']  = BPM_FLAT_THRESH
        fits.writeto(BPM_PATH, bpm.astype(np.uint8), bpm_hdr, overwrite=True, output_verify='silentfix')
        print(f"   📝 Bad Pixel Mask saved → {BPM_PATH}")
        fits.writeto(BPM_SCIENCE_PATH, bpm.astype(np.uint8), bpm_hdr, overwrite=True, output_verify='silentfix')
        print(f"   📝 Bad Pixel Mask synced → {BPM_SCIENCE_PATH}")

        # Save NormMasterFlat.fits
        hdr = hdul[0].header.copy()
        if 'NAXIS3' in hdr:
            del hdr['NAXIS3']
        hdr['NAXIS']   = 2
        hdr['CCDMEAN'] = 1.0  # Reset to 1.0 so IRAF ccdproc does not multiply science by old flat mean
        hdr['HISTORY'] = f'Normalised by centre-box median = {norm_val:.4f}'
        hdr['HISTORY'] = f'Bad pixels healed and floor clamped at {BPM_FLAT_THRESH}'
        fits.writeto(NORM_MASTER_FLAT_PATH, norm_data.astype(np.float32),
                     hdr, overwrite=True, output_verify='silentfix')

    print(f"\n✅ Normalised Master Flat saved → {NORM_MASTER_FLAT_PATH}")

except Exception as exc:
    print(f"❌ Error during flat normalisation: {exc}")


# %% [markdown]
# ##  Cell 8b — (Optional) Night-Sky Superflat Generation
# > **WHEN TO USE:**
# > • Only when you have >= 15-30 WIDELY DITHERED science frames in a sparse star field.
# > • In infrared bands (J/H/K) where sky brightness dominates twilight.
# >
# > **WHEN NOT TO USE (WARNING):**
# > • DO NOT USE if you have few frames (< 10) or undithered exposures (dither < 20 px).
# > • Why? Stacking few science frames injects high Poisson noise and creates fake dark
# >   "potholes" around unrejected faint ghost stars.
# > • For small or undithered datasets, use standard twilight flats + 2D illumination flattening (Cell 9).

# %%
FLAT_CORR_PATH = NORM_MASTER_FLAT_PATH  # Default: normalised twilight master flat

if ENABLE_SUPERFLAT:
    # Look for bias-corrected science frames (ready in DIR_SCIENCE)
    b_sci_for_superflat = sorted(glob.glob(os.path.join(DIR_SCIENCE, PREFIX_BIAS_CORR + '*.fit*')))
    b_sci_for_superflat = [f for f in b_sci_for_superflat if not os.path.basename(f).startswith(PREFIX_FLAT_CORR)]

    if len(b_sci_for_superflat) < SUPERFLAT_MIN_FRAMES:
        print(f"ℹ️ Found {len(b_sci_for_superflat)} science frame(s).")
        print(f"   Night-Sky Superflat requires >= {SUPERFLAT_MIN_FRAMES} frames.")
        print(f"   Falling back to standard Master Flat: {NORM_MASTER_FLAT_PATH}")
    else:
        print(f"🌌 Building Night-Sky Superflat from {len(b_sci_for_superflat)} science frame(s) ...")
        from astropy.io import fits
        import numpy as np
        from scipy import ndimage

        sf_frames = []
        for sf_src in b_sci_for_superflat:
            with fits.open(sf_src) as h:
                d = h[0].data.astype(np.float64)
                med = float(np.median(d))
                p16, p84 = np.percentile(d, [16, 84])
                std = float((p84 - p16) / 2.0)
                # Mask stars (sigma-clipping + dilation)
                star_mask = ndimage.binary_dilation(d > (med + 2.5 * std), iterations=2)
                d_norm = d / med
                d_norm[star_mask] = np.nan
                sf_frames.append(d_norm)

        # Median combine across unmasked science pixels
        superflat = np.nanmedian(np.array(sf_frames), axis=0)

        # Fill remaining all-masked pixels using 9x9 local median filter
        nan_mask = np.isnan(superflat)
        if np.sum(nan_mask) > 0:
            fill_map = ndimage.median_filter(np.nan_to_num(superflat, nan=1.0), size=9)
            superflat[nan_mask] = fill_map[nan_mask]

        # Normalise by centre box
        ny, nx = superflat.shape
        cy, cx = ny // 2, nx // 2
        CBOX = 200
        sf_norm = float(np.median(superflat[cy-CBOX:cy+CBOX, cx-CBOX:cx+CBOX]))
        superflat /= sf_norm

        # Save Superflat
        with fits.open(b_sci_for_superflat[0]) as h:
            sf_hdr = h[0].header.copy()
        sf_hdr['CCDMEAN'] = 1.0  # Reset to 1.0 so IRAF ccdproc does not multiply science by old flat mean
        sf_hdr['HISTORY'] = f'Night-Sky Superflat created from {len(b_sci_for_superflat)} frames'
        fits.writeto(SUPERFLAT_PATH, superflat.astype(np.float32), sf_hdr, overwrite=True, output_verify='silentfix')
        print(f"✅ Night-Sky Superflat saved → {SUPERFLAT_PATH}")

        if USE_SUPERFLAT_FOR_SCIENCE:
            FLAT_CORR_PATH = SUPERFLAT_PATH
            print(f"   👉 Using Superflat for science frame calibration in Cell 9.")
else:
    print("ℹ️ Night-Sky Superflat is disabled (ENABLE_SUPERFLAT = False). Using twilight Master Flat.")


# %% [markdown]
# ##  Cell 9 — Flat-Field Correction & 2D Illumination Flattening
# > Applies flat-field division, then fits and removes residual 2D sky/vignetting gradients
# > while strictly preserving stellar flux and positive ADU levels.

# %%
# Bias-corrected science frames are always in DIR_SCIENCE
bc_science_files = sorted(glob.glob(os.path.join(DIR_SCIENCE, PREFIX_BIAS_CORR + '*.fit*')))
bc_science_files = [f for f in bc_science_files if not os.path.basename(f).startswith(PREFIX_FLAT_CORR)]

if not bc_science_files:
    raise FileNotFoundError(
        f"No bias-corrected science frames found in {DIR_SCIENCE}\n"
        "  ↳ Make sure Cell 5 completed without errors."
    )

print(f"Found {len(bc_science_files)} bias-corrected science frame(s) to flat-field.\n")
iraf.unlearn('ccdproc')

from astropy.io import fits
import numpy as np
from scipy import ndimage

for src in bc_science_files:
    fname = PREFIX_FLAT_CORR + os.path.basename(src)   # e.g.  fb_b_01_12_slice001.fits
    dest  = os.path.join(DIR_SCIENCE, fname)
    if os.path.exists(dest):
        os.remove(dest)
    print(f"  {os.path.basename(src):45s}  →  {fname}")
    try:
        iraf.ccdproc(
            images   = src,
            output   = dest,
            ccdtype  = '',
            interactive = 'no',
            overscan = 'no',
            trim     = 'no',
            fixpix   = 'no',
            zerocor  = 'no',         # ← Already done in Cell 5
            darkcor  = 'no',
            flatcor  = 'yes',        # ← Flat-field correction ON
            flat     = FLAT_CORR_PATH,
            illumcor = 'no',
            fringecor= 'no',
            readcor  = 'no',
            scancor  = 'no'
        )

        # ── 2D Illumination / Residual Gradient Flattening ───
        if APPLY_ILLUM_FLATTEN and os.path.exists(dest):
            with fits.open(dest) as h:
                sci_data = np.squeeze(h[0].data).astype(np.float64)
                sci_hdr  = h[0].header.copy()

            med_val = float(np.median(sci_data))
            p16, p84 = np.percentile(sci_data, [16, 84])
            std_est = float((p84 - p16) / 2.0)

            # Mask point sources (stars) so they do not bias the background mesh
            star_mask = sci_data > (med_val + 2.5 * std_est)
            star_mask = ndimage.binary_dilation(star_mask, iterations=2)

            box = ILLUM_BOX_SIZE
            ny, nx = sci_data.shape
            grid_y, grid_x = ny // box, nx // box
            mesh = np.zeros((grid_y, grid_x), dtype=np.float64)

            for i in range(grid_y):
                for j in range(grid_x):
                    block = sci_data[i*box:(i+1)*box, j*box:(j+1)*box]
                    m = star_mask[i*box:(i+1)*box, j*box:(j+1)*box]
                    valid = block[~m]
                    if len(valid) > 10:
                        mesh[i, j] = np.median(valid)
                    else:
                        mesh[i, j] = med_val

            # Interpolate mesh to full image size via cubic spline
            smooth_bg = ndimage.zoom(mesh, (ny / grid_y, nx / grid_x), order=3)

            # Flatten 2D residual gradient while strictly conserving flux & positive ADU
            flat_sci = sci_data - smooth_bg + med_val

            sci_hdr['HISTORY'] = f'2D illumination/sky gradient flattened (mesh={box}px)'
            fits.writeto(dest, flat_sci.astype(np.float32), sci_hdr, overwrite=True, output_verify='silentfix')
            print(f"    ✨ Flattened residual 2D gradient on {fname} (mesh={box}px, sky={med_val:.1f} ADU)")

    except iraf.IrafError as exc:
        print(f"  ⚠️  ccdproc failed for {src}:\n     {exc}")
    except Exception as exc:
        print(f"  ⚠️  Illumination flattening error on {dest}:\n     {exc}")

print(f"\n✅ Fully-reduced science frames saved → {DIR_SCIENCE}")


# %% [markdown]
# # Cell 10 : Permanent Edge Trimming (Pre-Alignment)

# %%
# =====================================================================
# 10. Permanent Edge Trimming (Pre-Alignment)
# =====================================================================
import os
import glob
import numpy as np
from astropy.io import fits

def trim_science_frames(file_pattern, trim_bottom=0, trim_left=0, trim_top=0, trim_right=0, out_prefix="t"):
    files = sorted(glob.glob(file_pattern))
    if not files:
        print(f"⚠️ No files found matching '{file_pattern}'. Check your prefix.")
        return

    print(f"✂️ Trimming edges for {len(files)} files: [Bottom: {trim_bottom}px, Left: {trim_left}px, Top: {trim_top}px, Right: {trim_right}px]...")

    for file_path in files:
        dir_name = os.path.dirname(file_path)
        base_name = os.path.basename(file_path)
        out_path = os.path.join(dir_name, f"{out_prefix}{base_name}")

        with fits.open(file_path) as hdul:
            # Force 64-bit float to prevent byte-order errors later
            data = np.squeeze(hdul[0].data).astype(np.float64) 
            header = hdul[0].header.copy()

            # PERMANENT ASYMMETRIC TRIM:
            # Slicing: [Y_start:Y_end, X_start:X_end]
            # In FITS: Y=0 is bottom, X=0 is left
            y_start = trim_bottom
            y_end = data.shape[0] - trim_top if trim_top > 0 else None
            x_start = trim_left
            x_end = data.shape[1] - trim_right if trim_right > 0 else None

            trimmed_data = data[y_start:y_end, x_start:x_end]

            # Update header to reflect the new dimensions
            header['NAXIS1'] = trimmed_data.shape[1]
            header['NAXIS2'] = trimmed_data.shape[0]

            # Update WCS reference pixel if present
            if 'CRPIX1' in header and x_start > 0:
                header['CRPIX1'] = float(header['CRPIX1']) - x_start
            if 'CRPIX2' in header and y_start > 0:
                header['CRPIX2'] = float(header['CRPIX2']) - y_start

            header['HISTORY'] = f"Trimmed edges: bottom={trim_bottom}, left={trim_left}, top={trim_top}, right={trim_right}"

            # Save the cleanly trimmed file
            fits.writeto(out_path, trimmed_data, header, overwrite=True, output_verify='silentfix')

        print(f"  [Trimmed] {base_name}  ->  {os.path.basename(out_path)}  shape={trimmed_data.shape}")

    print(f"\n✅ Trimming complete! Your images are now safe to align.")

# Run the trimming process
trim_science_frames(
    os.path.join(DIR_SCIENCE, TRIM_INPUT_WILDCARD),
    trim_bottom=TRIM_BOTTOM,
    trim_left=TRIM_LEFT,
    trim_top=TRIM_TOP,
    trim_right=TRIM_RIGHT,
    out_prefix=TRIM_OUTPUT_PREFIX
)

# %% [markdown]
# ##  Cell 11 — Reduction Summary

# %%
print("=" * 60)
print("         REDUCTION PIPELINE COMPLETE")
print("=" * 60)

for label, directory, pattern in [
    ("BIAS",    DIR_BIAS,    '*.fit*'),
    ("FLATS",   DIR_FLATS,   '*.fit*'),
    ("SCIENCE", DIR_SCIENCE, 'fb_*.fit*'),
    ("TRIMMED", DIR_SCIENCE, f"{TRIM_OUTPUT_PREFIX}*.fit*"),
]:
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    print(f"\n📁 {label}  ({directory})")
    for f in files:
        size_kb = os.path.getsize(f) / 1024
        print(f"   {os.path.basename(f):50s}  {size_kb:8.1f} KB")
    if not files:
        print("   (empty)")

print("\n" + "=" * 60)
print("Key calibration files:")
print(f"  Master Bias         : {MASTER_BIAS_PATH}")
print(f"  Master Flat         : {MASTER_FLAT_PATH}")
print(f"  Norm. Master Flat   : {NORM_MASTER_FLAT_PATH}")
print(f"  Bad Pixel Mask      : {BPM_PATH}")
if ENABLE_SUPERFLAT:
    print(f"  Night-Sky Superflat : {SUPERFLAT_PATH}")
print("=" * 60)


