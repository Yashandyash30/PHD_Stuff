# GRB Photometry Pipeline — In-Depth Technical Documentation

**Target:** GRB 241002A (configuration uses GRB 200412 as example)  
**Instrument:** Optical telescope + CCD camera  
**Language:** Python (Jupyter Notebook)

---

## Table of Contents

1. [Pipeline Overview](#1-pipeline-overview)
2. [Cell 1 — Imports & Dependencies](#2-cell-1--imports--dependencies)
3. [Cell 2 — Configuration & Parameters](#3-cell-2--configuration--parameters)
4. [Cell 3 — Output Folder Setup](#4-cell-3--output-folder-setup)
5. [Cell 4 — Core Function Definitions](#5-cell-4--core-function-definitions)
   - [F0: Image Stacking](#f0-stack_images--deep-image-co-addition)
   - [F1: Local Plate Solving](#f1-run_local_plate_solve)
   - [F1.5: Online Plate Solving](#f15-run_online_plate_solve)
   - [F2: FWHM & Background Estimation](#f2-measure_image_fwhm)
   - [F4: DAOStarFinder + Morphology Filtering](#f4-run_daofinder)
   - [F5: Multi-Aperture Photometry](#f5-run_multi_aperture_photometry)
   - [F6: WCS Astrometry & Target Matching](#f6-match_target_to_catalog)
   - [F6.5: Aperture Recentering](#f65-recenter_target_coordinates)
   - [F7: Zero-Point Calibration](#f7-calibrate_zero_point)
   - [F8: Forced Photometry & Upper Limit](#f8-forced_photometry_upper_limit)
   - [F13: PSF Photometry](#f13-run_psf_photometry)
   - [F14: Extract PSF Magnitude](#f14-extract_psf_magnitude)
   - [Plotting Functions (F9–F12, F15)](#plotting-functions-f9f12-f15)
6. [Execution Cells (5–16): Pipeline Flow](#6-execution-cells-pipeline-flow)
7. [Mathematical Reference](#7-mathematical-reference)
8. [Catalog Selection Guide](#8-catalog-selection-guide)
9. [Tuning & Troubleshooting Guide](#9-tuning--troubleshooting-guide)

---

## 1. Pipeline Overview

This pipeline performs **differential aperture photometry** on optical FITS images to detect and measure the brightness (magnitude) of a Gamma-Ray Burst (GRB) afterglow. The pipeline proceeds through these major stages:

```
Raw FITS frames
      │
      ▼
[F0]  Image Stacking (optional, mode='stack')
      │  Align multiple frames → sigma-clipped mean
      ▼
[F1]  Plate Solving
      │  Associates pixel coordinates with RA/Dec (WCS)
      ▼
[F2]  Background Estimation + FWHM Measurement
      │  2D sky background model; Gaussian star fits
      ▼
[F4]  Source Detection (DAOStarFinder + morphology filter)
      │  Finds all stellar-like objects in the image
      ▼
[F5]  Multi-Aperture Photometry (Curve of Growth)
      │  Measures flux in rings of increasing radius
      ▼
[F6]  Astrometry: Match catalog to target RA/Dec
      │  Identifies which (if any) detected source is the GRB
      ▼
[F7]  Zero-Point Calibration via Pan-STARRS / APASS
      │  Converts raw ADU flux → calibrated magnitudes
      ▼
[F8]  Forced Photometry & Upper Limit
      │  Always measures at the GRB position, even if faint
      ▼
[F13] PSF Photometry (Advanced)
      │  Fits analytical PSF models; handles blended sources
      ▼
      Final magnitude or 3σ upper limit
```

The pipeline supports three run modes:
- **`single`** — processes one FITS frame
- **`stack`** — co-adds multiple frames then processes the deep stack
- **`batch`** — loops through all frames matching a glob pattern

---

## 2. Cell 1 — Imports & Dependencies

This cell loads all required Python libraries. Key packages and their roles:

| Package | Role |
|---|---|
| `astropy` | FITS I/O, WCS, coordinates, time, statistics, modeling |
| `photutils` | Background estimation, source detection, aperture photometry, PSF photometry |
| `astroquery.vizier` | Remote catalog queries (Pan-STARRS, APASS, Gaia, Landolt) |
| `astroquery.astrometry_net` | Online plate solving via astrometry.net API |
| `astroalign` | Geometric image registration for stacking |
| `astroscrappy` | L.A.Cosmic cosmic-ray rejection (imported, stub for future use) |
| `matplotlib` | All diagnostic and output plots |
| `numpy` | Array math, masking, sigma clipping |
| `scipy` (via astropy) | Gaussian fitting via `LevMarLSQFitter` |

---

## 3. Cell 2 — Configuration & Parameters

All pipeline tuning lives here. Nothing is hardcoded deeper in the code. Key parameter groups:

### 2A. Pipeline Mode & Paths

```python
RUN_MODE            = 'stack'     # 'single', 'stack', or 'batch'
PLATE_SOLVE_MODE    = 'online'    # 'local', 'online', 'copy', or 'external'
SINGLE_FRAME_FILE   = "fbGRB200412B_R_2148x2148_18022.fits"
STACK_PATTERN       = "fbGRB200412B_R*.fits"
STACKED_FILENAME    = "MASTER_STACKED_fbGRB200412B_R.fits"
```

### 2B. Target Coordinates

```python
TARGET_RA    = 278.313375   # degrees (J2000)
TARGET_DEC   = +62.532528   # degrees (J2000)
GRB_T0_JD   = 2458951.881027  # Julian Date of GRB trigger
```

These are the known coordinates of the GRB, obtained from GCN circulars or Swift/Fermi catalogs. The pipeline projects these onto pixels using the WCS solution to search for a counterpart.

### 2C. Detection & Photometry Parameters

| Parameter | Value | Meaning |
|---|---|---|
| `DETECT_SIGMA` | 4.0 | Source must be 4× the sky noise to be detected |
| `EDGE_BUFFER_PIXELS` | 100 | Sources within 100px of edge are rejected |
| `BG_BOX_SIZE` | (25, 25) | Size of tiles for 2D background model |
| `SATURATION_LIMIT` | 55000 | ADU value above which a pixel is saturated |
| `EFFECTIVE_GAIN` | 5.0 | e⁻/ADU — converts ADU to electrons for Poisson noise |
| `READ_NOISE` | 10.0 | e⁻ — CCD read noise floor |
| `APERTURE_MULTIPLIERS` | [1.25, 1.5, 1.75, 2.0, 2.5, 3.0] | Aperture radii as multiples of FWHM |
| `OPTIMAL_APERTURE` | 1.5 | Which multiplier to use for the final magnitude |
| `ANNULUS_INNER_MULT` | 2.0 | Inner sky annulus = 2 × FWHM |
| `ANNULUS_OUTER_MULT` | 4.0 | Outer sky annulus = 4 × FWHM |

### 2D. Morphology Filters (The "Bouncer" Rules)

DAOStarFinder assigns each detected source three shape metrics. Only sources satisfying all three pass:

| Parameter | Range | Rejects |
|---|---|---|
| `ROUNDNESS1_RANGE` | (-0.6, 0.6) | Vertically/horizontally elongated objects (e.g., satellites) |
| `ROUNDNESS2_RANGE` | (-0.6, 0.6) | Diagonally elongated objects |
| `SHARPNESS_RANGE` | (0.2, 1.0) | Fuzzy galaxies (low) and hot pixels/cosmic rays (high) |

### 2E. Catalog Calibration

```python
CALIBRATION_SYSTEM  = 'Pan-STARRS'   # or 'DIRECT_VEGA'
DIRECT_VEGA_CATALOG = 'APASS'        # only used if DIRECT_VEGA
FILTER_BAND         = 'R'
```

When `Pan-STARRS` is selected, the pipeline queries the PS1 catalog and applies Lupton (2005) color transformations to synthesize Bessell R or I magnitudes from the native g/r/i PS1 passbands.

---

## 4. Cell 3 — Output Folder Setup

```python
def setup_output_dirs(input_name, mode):
```

Creates a hierarchical output directory tree rooted at `./pipeline_output/`. The structure is:

```
pipeline_output/
├── fits_solved/              ← WCS-solved FITS files (shared across runs)
└── <name>_<mode>_<calib>/
    ├── plots/
    │   ├── detection/        ← DAOFind source maps
    │   ├── calibration/      ← ZP scatter plots
    │   └── target/           ← GRB confirmation stamps
    ├── fits/
    │   ├── reduced/          ← Background-subtracted FITS
    │   └── stacked/          ← Co-added master FITS
    ├── catalogs/             ← Photometry tables (.fit, .csv)
    └── logs/                 ← Detailed text run logs
```

The directory structure is stored in a global dict `DIRS` so all functions can find the right path without hardcoded strings.

---

## 5. Cell 4 — Core Function Definitions

This is the heart of the pipeline. All reusable logic is encapsulated here.

---

### F0: `stack_images` — Deep Image Co-addition

**Purpose:** Combines multiple short exposures into one deep image, suppressing cosmic rays and noise.

**Mathematical steps:**

1. **Load & normalize:** All frames are cast to `float64`. NaN pixels (from bad columns, overscan strips) are replaced with the frame median so they don't corrupt alignment math.

2. **Alignment (astroalign):** Each frame is registered to a reference (frame 0) using a geometric transform. The transform is found on a central crop (`[200:-200, 200:-200]`) to avoid edge artifacts:

   ```
   T, _ = astroalign.find_transform(source_crop, reference_crop)
   registered = astroalign.apply_transform(T, source, reference)
   ```

   Internally, astroalign detects star-like sources (via a LoG-based finder), builds a k-d tree of triangle descriptors, and finds a similarity transform (rotation + translation + scale) via RANSAC.

3. **Sigma-clipped mean stacking:** All aligned frames are stacked into a 3D cube `(N_frames, height, width)`. A sigma clip (σ=3.0, 3 iterations) is applied along the frame axis (axis=0). This masks outliers — cosmic rays, satellite trails, bad pixels — that appear in only a subset of frames. The unmasked mean is then computed:

   $$\text{Stack}[y, x] = \frac{1}{N_{\text{valid}}} \sum_{k=1}^{N_{\text{valid}}} I_k[y, x]$$

   The signal-to-noise ratio improves as √N relative to a single frame.

4. **Output:** Saved as a FITS file with the reference frame's header.

---

### F1: `run_local_plate_solve`

**Purpose:** Assigns a World Coordinate System (WCS) to a FITS image by matching stars against an astrometric reference catalog using the `astrometry.net` `solve-field` command-line tool.

**Key subprocess call:**

```bash
solve-field image.fits \
  --config /path/to/local.cfg \
  --ra TARGET_RA --dec TARGET_DEC --radius 2.0 \
  --scale-units arcsecperpix --scale-low 0.50 --scale-high 1.00 \
  --downsample 2 --cpulimit 300 --no-plots --overwrite
```

The result is a `.new` FITS file containing a full WCS (TAN projection with SIP distortion coefficients). If found, this is moved to `fits_solved/` and a `WCS` object is extracted.

**Caching:** If a solved file already exists for that input, it skips the computationally expensive solve and loads the cached result directly.

---

### F1.5: `run_online_plate_solve`

**Purpose:** Same as F1 but uploads the image to `nova.astrometry.net` via its REST API. Useful when local index files are not available.

**Steps:**
1. Initialises `astroquery.astrometry_net.AstrometryNet` with the API key.
2. Sends the image with a RA/Dec hint and scale bounds.
3. Waits for the server to return a WCS FITS header.
4. Writes a new FITS file combining the original image data with the new WCS header.

---

### F2: `measure_image_fwhm`

**Purpose:** Estimates the image **Point Spread Function (PSF) size** as a Full Width at Half Maximum (FWHM) in pixels, and computes a spatially-varying 2D sky background model.

This function does four things:

#### Step 1: Iterative Source Masking

Before fitting the background, the code runs three iterations of source masking to prevent the stars themselves from biasing the background estimate:

```
Iteration 1–3:
  threshold = detect_threshold(data, nsigma=2.0, sigma_clip=...)
  seg_map   = detect_sources(data, threshold, npixels=10)
  source_mask = seg_map.make_source_mask(footprint=circular_footprint(radius=10))
```

Each iteration masks all previously detected sources. This progressively reveals fainter sources that were previously hidden in the "noise" next to bright stars.

#### Step 2: 2D Background Model (`Background2D`)

The unmasked pixels are divided into tiles of size `BG_BOX_SIZE = (25, 25)`. Within each tile, the background is estimated using `SExtractorBackground`, which computes:

$$\text{BKG} = 2.5 \times \text{median} - 1.5 \times \text{mean}$$

(This formula is robust to source contamination inside a tile.) A smoothing filter of size `BG_FILTER_SIZE = (3, 3)` is applied across tiles to produce a smooth spatial background surface `bkg.background`. The background-subtracted image is:

$$\text{clean\_data}[y, x] = \text{data}[y, x] - \text{bkg.background}[y, x]$$

The background RMS map `bkg.background_rms` provides a spatially-varying noise estimate.

#### Step 3: Peak Finding

Bright, unsaturated, non-edge stellar peaks are found using `photutils.detection.find_peaks` with a threshold of `FWHM_PEAK_SIGMA × std_bg`. The top 40 brightest valid peaks are selected for FWHM fitting.

#### Step 4: Gaussian Fitting for FWHM

For each candidate star, a 25×25 pixel cutout is extracted. A 2D model:

$$f(x, y) = C + A \cdot \exp\!\left(-\frac{(x - x_0)^2}{2\sigma_x^2} - \frac{(y - y_0)^2}{2\sigma_y^2}\right)$$

is fitted using the Levenberg-Marquardt least-squares fitter (`LevMarLSQFitter`). The σ bounds are constrained to `(0.5, 6.0)` pixels to reject non-stellar fits.

The FWHM for each star is then:

$$\text{FWHM} = 2.355 \times \frac{\sigma_x + \sigma_y}{2}$$

The constant 2.355 comes from the Gaussian relationship: $\text{FWHM} = 2\sqrt{2\ln 2}\,\sigma \approx 2.355\,\sigma$.

Only FWHM values in the range (1.5, 15.0) pixels are accepted. A sigma-clipped median over all accepted values gives the **final FWHM**, which represents the typical seeing disk size.

---

### F4: `run_daofinder`

**Purpose:** Detects all point-source objects in the background-subtracted image.

**Algorithm — DAOStarFinder (DAOPHOT method):**

1. Convolves the image with a Gaussian kernel whose sigma is derived from the FWHM.
2. Detects local maxima above the threshold: `DETECT_SIGMA × std_bg`.
3. For each candidate, computes three morphological statistics:

**Roundness 1** (horizontal-vertical asymmetry):
$$\text{roundness1} = \frac{2 \cdot (I_h - I_v)}{I_h + I_v}$$
where $I_h$ and $I_v$ are the sums along the horizontal and vertical directions of the source profile. A perfect circular star gives 0.0.

**Roundness 2** (diagonal asymmetry, based on a Gaussian fit):
$$\text{roundness2} = \frac{2(\sigma_y^2 - \sigma_x^2)}{\sigma_x^2 + \sigma_y^2}$$

**Sharpness** (central peak relative to the fitted Gaussian amplitude):
$$\text{sharpness} = \frac{I_{\text{peak}} - \text{background}}{G_{\text{peak}}}$$
where $G_\text{peak}$ is the amplitude of the fitted Gaussian. Real stars have sharpness ≈ 0.5–0.8. Cosmic rays have sharpness >> 1.

**Edge filter:** Sources within `EDGE_BUFFER_PIXELS = 100` of the image boundary are removed (often contain diffraction spikes, vignetting artifacts).

**Result:** A cleaned catalog of point sources — predominantly stars — with their pixel centroids.

---

### F5: `run_multi_aperture_photometry`

**Purpose:** For every detected source, measures the flux inside multiple circular apertures to build a **Curve of Growth** and choose the optimal aperture.

#### Step 1: Error Array

A pixel-by-pixel noise map is computed via:

$$\sigma[y,x] = \sqrt{\sigma_\text{sky}^2 + \frac{I[y,x]}{\text{gain}}}$$

where:
- $\sigma_\text{sky}$ = sky background noise (= `std_bg`, per pixel)
- $I[y,x]$ = source flux in ADU (Poisson shot noise term)
- `gain` = `EFFECTIVE_GAIN` in e⁻/ADU

This is implemented via `photutils.utils.calc_total_error`.

#### Step 2: Local Background Subtraction (Annulus)

A circular annulus is placed around each source:
- Inner radius: $r_\text{in} = 2.0 \times \text{FWHM}$
- Outer radius: $r_\text{out} = 4.0 \times \text{FWHM}$

The median pixel value in the annulus (after 3σ clipping) is used as the **local background** estimate. Using the median (rather than mean) makes it robust to nearby faint sources that partially overlap the annulus.

#### Step 3: Aperture Photometry Loop

For each multiplier $m$ in `[1.25, 1.5, 1.75, 2.0, 2.5, 3.0]`, a circular aperture of radius $r = m \times \text{FWHM}$ is constructed. The background-corrected flux is:

$$F_m = \sum_{i \in \text{aperture}} \left(I[i] - B_\text{local}\right)$$

The flux error is propagated via:

$$\sigma_{F_m} = \sqrt{\sum_{i \in \text{aperture}} \sigma[i]^2}$$

All aperture results are stored in columns `flux_1.25x`, `flux_1.5x`, ..., `flux_3.0x` in the master photometry table.

**Why multiple apertures?** The optimal aperture balances two effects:
- Too small → misses flux in the outer PSF wings (under-counts)
- Too large → includes more sky noise and neighbour contamination (over-noisy)

The Curve of Growth plot (F12) visualises this trade-off. The pipeline uses `OPTIMAL_APERTURE = 1.5` by default, which captures ~95% of the PSF flux for a Gaussian PSF while minimising noise.

---

### F6: `match_target_to_catalog`

**Purpose:** Uses the WCS to project the photometry catalog onto the sky and check whether any source lands at the GRB's known position.

#### Step 1: Pixel → Sky Coordinates

```python
ra, dec = wcs.all_pix2world(xcentroid, ycentroid, 0)
```

`all_pix2world` applies the full WCS transformation including any SIP polynomial distortion coefficients.

#### Step 2: Sky Match

The GRB's known coordinates are expressed as an `astropy.coordinates.SkyCoord`. The function calls:

```python
idx, d2d, _ = target_coord.match_to_catalog_sky(catalog_coords)
```

This finds the single nearest source in the catalog using a k-d tree on the unit sphere. The angular separation $d_{2d}$ is compared against `MATCH_RADIUS_ARCSEC = 3.0"`.

$$\text{Detected} \iff d_{2d} \leq 3.0''$$

#### Step 3: Inverse Projection

```python
grb_pixel_x, grb_pixel_y = wcs.wcs_world2pix(TARGET_RA, TARGET_DEC, 0)
```

This gives the expected pixel position of the GRB, which is used even when no catalog source is detected (for forced photometry).

---

### F6.5: `recenter_target_coordinates`

**Purpose:** Refines the GRB pixel position by fitting a 2D Gaussian centroid to the local image data. This corrects for small WCS residuals (typically < 1–2 pixels).

**Method:**
1. A square cutout of size `2 × RECENTER_RADIUS_PIX + 1` is extracted around the WCS-predicted position.
2. `photutils.centroids.centroid_2dg` fits a 2D Gaussian to the cutout and returns the sub-pixel centroid.
3. If the shift is ≤ `RECENTER_RADIUS_PIX`, the new coordinates are accepted. Otherwise, the original WCS coordinates are retained.

This prevents the recentering from jumping to a nearby unrelated source.

---

### F7: `calibrate_zero_point`

**Purpose:** Converts instrumental magnitudes (arbitrary log-flux units) to a physically meaningful standard photometric system by computing the **photometric zero point** ZP.

The fundamental relation is:

$$m_\text{calib} = m_\text{inst} + \text{ZP}$$

where the instrumental magnitude of a source with flux $F$ ADU is:

$$m_\text{inst} = -2.5 \log_{10}(F)$$

#### Step 1: Instrumental Magnitudes

For all detected sources with positive flux:

$$m_\text{inst} = -2.5 \log_{10}(F_\text{opt})$$
$$\sigma_{m_\text{inst}} = \frac{1.0857 \cdot \sigma_F}{F}$$

(The factor 1.0857 = $\frac{1}{\ln 10}$ converts fractional flux uncertainty to magnitude uncertainty.)

#### Step 2: Reference Star Query

The pipeline queries VizieR (an online astronomical catalog service) for stars within `QUERY_RADIUS_ARCMIN` of the image center. Depending on the calibration system:

**Pan-STARRS (PS1):** Downloads native PS1 g, r, i magnitudes. Stars are quality-filtered on magnitude range and photometric error.

**DIRECT_VEGA (APASS/Landolt/Gaia):** Downloads the native filter magnitude directly (e.g. Vmag from APASS for V-band observations).

#### Step 3: Passband Transformations (Pan-STARRS only)

Because PS1 passbands differ from the Bessel R/I passbands of the CCD, Lupton (2005) color transformations are applied:

**Bessel R from PS1:**
$$R = r_\text{PS1} - 0.2936(r_\text{PS1} - i_\text{PS1}) - 0.1439$$

**Bessel I from PS1:**
$$I = r_\text{PS1} - 1.2444(r_\text{PS1} - i_\text{PS1}) - 0.3820$$

**Bessel V from PS1:**
$$V = g_\text{PS1} - 0.5784(g_\text{PS1} - r_\text{PS1}) - 0.0038$$

**Bessel B from PS1:**
$$B = g_\text{PS1} + 0.3130(g_\text{PS1} - r_\text{PS1}) + 0.2271$$

These transformations synthesize standard Vega-system magnitudes from the AB-system PS1 magnitudes. The associated errors are propagated with standard error propagation (sum in quadrature of partial derivatives × individual catalog errors).

#### Step 4: Cross-Matching

For each reference star, its RA/Dec is matched to the nearest detected source in the pipeline catalog within `ZP_MATCH_RADIUS_AS`. The individual zero-point contribution of each match is:

$$\text{ZP}_k = m_{\text{true},k} - m_{\text{inst},k}$$

#### Step 5: Master Zero Point

After collecting all individual ZP estimates, a **3σ sigma-clipped mean** is computed:

$$\text{ZP}_\text{master} = \text{clipped\_mean}(\text{ZP}_k), \quad \sigma_\text{ZP} = \text{clipped\_std}(\text{ZP}_k)$$

Sigma-clipping removes the effect of mismatched pairs (blended sources, stars near saturation). A minimum of `ZP_MIN_STARS = 5` matches is required for a valid calibration.

---

### F8: `forced_photometry_upper_limit`

**Purpose:** Regardless of whether the GRB was catalogued, always perform photometry at its exact expected position. If no source is detected, reports a statistical upper limit.

#### Detection Branch

Aperture and annulus photometry are run at `(grb_pixel_x, grb_pixel_y)` using the same parameters as F5.

If the measured flux $F$ exceeds the $N\sigma$ threshold:

$$F > N_\sigma \cdot \sigma_F \quad (N_\sigma = 3 \text{ by default})$$

then the calibrated magnitude is:

$$m = \text{ZP} - 2.5\log_{10}(F)$$
$$\sigma_m = \sqrt{\left(\frac{1.0857\,\sigma_F}{F}\right)^2 + \sigma_\text{ZP}^2}$$

#### Non-Detection Branch (Upper Limit)

If the flux does not exceed the threshold, a 3σ upper limit is reported:

$$m_\text{UL} = \text{ZP} - 2.5\log_{10}(3\sigma_F)$$

This is the standard convention in GRB afterglow papers: the GRB was fainter than this magnitude at the time of observation.

---

### F13: `run_psf_photometry`

**Purpose:** Fits analytical PSF models to all sources simultaneously. Superior to aperture photometry for:
- **Blended sources:** Can separate two overlapping stars ("deblending")
- **Low SNR sources:** Extracts more signal by using the known PSF shape as a matched filter

**Model:** `CircularGaussianPRF` — a Pixel-Response Function (PRF), which is the PSF convolved with the pixel response (a square top-hat function). This accurately represents how a star's light is actually distributed across CCD pixels.

**Grouper:** `SourceGrouper(min_separation = 2.0 × FWHM)` — groups sources closer than 2 FWHM together so they are fitted simultaneously (the "peanut" rule). Sources in the same group compete for flux, enabling deblending.

**Fitting box:** `PSF_FIT_SHAPE = (5, 5)` — the fitting box must be smaller than the distance between stars but large enough to cover the high-SNR PSF core.

**Process:**
1. `DAOStarFinder` detects seed positions.
2. `LocalBackground` estimates local background within an annulus.
3. `PSFPhotometry` performs the simultaneous grouped least-squares fit.
4. A residual image (data minus all PSF models) is produced. Residuals near stars should be featureless noise; structure indicates poor fits.

---

### F14: `extract_psf_magnitude`

**Purpose:** Locates the PSF model that was placed at (or near) the GRB coordinates and reads out its fitted flux.

The PSF table entry closest to `(grb_pixel_x, grb_pixel_y)` within `PSF_MATCH_RADIUS_PIXELS = 5` pixels is selected:

$$\text{dist}_i = \sqrt{(x_\text{fit,i} - x_\text{GRB})^2 + (y_\text{fit,i} - y_\text{GRB})^2}$$

The PSF flux is then calibrated using the same zero point:

$$m_\text{PSF} = \text{ZP} - 2.5\log_{10}(F_\text{PSF})$$
$$\sigma_{m_\text{PSF}} = \sqrt{\left(\frac{1.0857\,\sigma_{F_\text{PSF}}}{F_\text{PSF}}\right)^2 + \sigma_\text{ZP}^2}$$

If no PSF model is found within the match radius (source too faint for PSF to detect), the function falls back gracefully and signals non-detection.

---

### Plotting Functions (F9–F12, F15)

| Function | Output | Purpose |
|---|---|---|
| `plot_detection_map` (F9) | `daofind_detection.png` | Side-by-side: raw detections (red) vs. filtered star-like sources (green) |
| `plot_target_confirmation` (F10) | `<TARGET>_confirmation.png` | Full image + 200px zoom with red crosshair on GRB position |
| `plot_calibration_diagnostics` (F11) | `zp_calibration.png` | Reference star positions on image + instrumental vs. true magnitude scatter with ZP line |
| `plot_aperture_diagnostic` (F12) | `<TARGET>_aperture_diagnostic.png` | Visual aperture rings + Curve of Growth (flux vs. aperture multiplier) |
| `plot_grb_psf` (F15) | `<TARGET>_PSF_profile.png` | 3-panel postage stamp: Original | PSF Model | Residual |

All plots are saved at 150 DPI. The Curve of Growth plot shades `DIAG_CONTAM_ZONE = (2.0, 3.0)` in red to flag aperture sizes where neighbour contamination typically begins.

---

## 6. Execution Cells: Pipeline Flow

### Cell 4.5 — Pre-Processing Switchboard

Decides what to process based on `RUN_MODE`:

```python
if RUN_MODE == 'stack':
    FINAL_IMAGE_TO_PROCESS = stack_images(STACK_PATTERN, stacked_path)
elif RUN_MODE == 'single':
    FINAL_IMAGE_TO_PROCESS = SINGLE_FRAME_FILE
```

### Cell 5 — Plate Solving (Step 0.5)

Dispatches to the correct plate-solving strategy:

| `PLATE_SOLVE_MODE` | Action |
|---|---|
| `'online'` | `run_online_plate_solve()` — uploads to nova.astrometry.net |
| `'local'` | `run_local_plate_solve()` — runs `solve-field` locally |
| `'copy'` | `run_copy_wcs_plate_solve()` — copies WCS from a pre-solved reference |
| `'external'` | Loads a pre-solved FITS file directly |

Output: `solved_path`, `data`, `header`, `wcs` — the WCS-equipped image.

### Cell 6 — FWHM & Background (Steps 1–2)

```python
clean_data, final_fwhm, background_rms, mean_bg, median_bg, std_bg, header = \
    measure_image_fwhm(solved_path)
```

`clean_data` is used for all subsequent detection and photometry steps.

### Cell 8 — Source Detection (Step 5)

```python
sources       = DAOStarFinder(fwhm=final_fwhm, threshold=DETECT_SIGMA * std_bg)(clean_data)
clean_sources = run_daofinder(clean_data, std_bg, final_fwhm)
plot_detection_map(data, sources, clean_sources)
```

Two runs are shown: raw DAOFind detections and filtered star-like sources.

### Cell 9 — Multi-Aperture Photometry (Step 4)

```python
master_phot_table, error_array = run_multi_aperture_photometry(clean_data, clean_sources, final_fwhm, std_bg)
```

Produces the master table with flux columns for all aperture sizes.

### Cell 10 — Astrometry & GRB Verdict (Step 7)

Calls `match_target_to_catalog()` then `recenter_target_coordinates()` and `plot_target_confirmation()`.

After this cell, `grb_pixel_x` / `grb_pixel_y` are the best-estimate pixel coordinates of the GRB, and `detected` is a boolean verdict.

### Cell 11 — Zero-Point Calibration (Step 8)

```python
master_phot_table, master_zp, master_zp_err, good_stars, zp_idx, zp_d2d = \
    calibrate_zero_point(master_phot_table, wcs, clean_data, catalog_coords)
```

Populates `mag_actual` and `mag_actual_err` columns in `master_phot_table`.

### Cell 12 — Calibration Diagnostics (Step 8.5)

Calls `plot_calibration_diagnostics()`. Useful for visually inspecting the quality of the ZP fit and spotting outlier reference stars.

### Cell 13 — Forced Photometry (Step 9)

```python
result = forced_photometry_upper_limit(clean_data, error_array, grb_pixel_x, grb_pixel_y,
                                        final_fwhm, master_zp, master_zp_err)
```

Prints the **final result**: either `magnitude ± mag_err` or `upper limit > X mag`.

### Cell 10.5 — Curve of Growth Diagnostic

Calls `plot_aperture_diagnostic()` to visualise how flux grows with aperture size at the GRB position. This confirms the chosen `OPTIMAL_APERTURE` is in the flat part of the growth curve.

### Cell 10.8 — PSF Photometry

```python
psf_phot_table, resid_image = run_psf_photometry(clean_data, error_array, final_fwhm, std_bg, TARGET_NAME)
```

### Cell 15 — PSF Magnitude Extraction

```python
psf_result = extract_psf_magnitude(psf_phot_table, grb_pixel_x, grb_pixel_y, master_zp, master_zp_err, TARGET_NAME)
```

### Final Cell — Save All Outputs

Saves:
1. **Clean FITS** (background-subtracted) with provenance history keyword
2. **Aperture catalog** as `.fit` (FITS binary table) and `.csv`
3. **PSF catalog** as `.csv` (if PSF photometry was run)
4. **Comprehensive log** (`.log` file) with all parameters, ZP cross-check table, and final magnitude

---

## 7. Mathematical Reference

### Magnitude System

All magnitudes in this pipeline follow the **Pogson** logarithmic scale:

$$m = -2.5 \log_{10}(F) + \text{ZP}$$

where $F$ is the source flux in ADU. The zero point ZP anchors the instrumental scale to the standard Vega or AB magnitude system via reference stars.

### Error Propagation Summary

| Quantity | Formula |
|---|---|
| Instrumental magnitude error | $\sigma_m = \frac{1.0857\,\sigma_F}{F}$ |
| Calibrated magnitude error | $\sigma_{m_\text{cal}} = \sqrt{\sigma_{m_\text{inst}}^2 + \sigma_\text{ZP}^2}$ |
| Total flux error (sky + Poisson) | $\sigma_F = \sqrt{\sigma_\text{sky}^2 + F/G}$  where $G$ = gain |
| FWHM from Gaussian σ | $\text{FWHM} = 2.355\,\sigma$ |
| 3σ Upper Limit magnitude | $m_\text{UL} = \text{ZP} - 2.5\log_{10}(3\sigma_F)$ |

### Signal-to-Noise Ratio

The detection significance of a source with flux $F$ in an aperture of $N_\text{pix}$ pixels with local background $B$ per pixel:

$$\text{SNR} = \frac{F}{\sqrt{F/G + N_\text{pix}(\sigma_\text{sky}^2 + B/G)}}$$

The pipeline detects sources at `DETECT_SIGMA = 4.0`, i.e., SNR ≥ 4.

---

## 8. Catalog Selection Guide

| Target Filter | Calibration System | Catalog | Notes |
|:---:|---|---|---|
| **Bessel U** | SDSS (manual) | `V/147/sdss12` | Neither APASS nor PS1 covers UV; use SDSS if in footprint |
| **Bessel B** | `DIRECT_VEGA` / `APASS` | `II/336/apass9` | Native Johnson B, no transformation needed |
| **Bessel V** | `DIRECT_VEGA` / `APASS` | `II/336/apass9` | Native Johnson V, no transformation needed |
| **Bessel R** | `Pan-STARRS` | `II/349/ps1` | Lupton R = r − 0.2936(r−i) − 0.1439 |
| **Bessel I** | `Pan-STARRS` | `II/349/ps1` | Lupton I = r − 1.2444(r−i) − 0.3820 |

**Decision rules:**
- B and V band: prefer APASS (native Vega system, no colour terms needed)
- R and I band: always use Pan-STARRS (deeper, more stars, well-calibrated transformations)
- U band: requires SDSS query and custom transformations (not automated in this pipeline)
- If APASS stars are too faint for a very deep V-band image, fall back to Pan-STARRS V

---

## 9. Tuning & Troubleshooting Guide

### Too few / no sources detected

| Symptom | Cause | Fix |
|---|---|---|
| 0 sources after DAOFind | Detection threshold too high | Lower `DETECT_SIGMA` (e.g., 3.5 → 3.0) |
| Good stars rejected by filters | Stars slightly trailed | Widen `ROUNDNESS1_RANGE` to (-0.8, 0.8) |
| Tight, bright stars rejected | Sharpness upper bound too low | Raise `SHARPNESS_RANGE` upper bound to 1.5 |
| FWHM fails | No peaks above `FWHM_PEAK_SIGMA × std` | Lower `FWHM_PEAK_SIGMA` |

### Zero point fails

| Symptom | Cause | Fix |
|---|---|---|
| < 5 ZP stars | Image is off-catalog footprint | Increase `QUERY_RADIUS_ARCMIN` |
| < 5 ZP stars | Stars are too faint/bright | Widen `REF_MAG_BRIGHT` / `REF_MAG_FAINT` range |
| Large ZP scatter | Blended reference stars | Decrease `ZP_MATCH_RADIUS_AS` |
| Wrong catalog for filter | Transformation mismatch | Check Table in Section 8 and update `CALIBRATION_SYSTEM` |

### Plate solving fails

| Symptom | Fix |
|---|---|
| Local solve times out | Increase `--cpulimit`; check index files cover the field |
| Online solve fails | Check API key; ensure image has stars (not empty field) |
| WCS looks wrong | Use `'copy'` mode from a previously solved frame of same field |

### GRB not detected (expected detection)

1. Check `plot_target_confirmation.png` — is the crosshair in the right place?
2. Check `MATCH_RADIUS_ARCSEC` — the WCS residual may exceed 3"
3. Run the Curve of Growth diagnostic — is there a flux excess at the GRB position even if it didn't make the catalog?
4. Check `forced_photometry_upper_limit` — the source may be below the detection threshold; a meaningful upper limit is still scientifically useful

### Image stacking issues

| Problem | Fix |
|---|---|
| `astroalign` alignment failure | Check that images overlap sufficiently; lower `detection_sigma` in `find_transform` |
| Stars trailed in stack | Images may not be all the same field; verify all frames are pointed at the same target |
| NaN regions in stack | Normal at image edges where not all frames overlap; add `EDGE_BUFFER_PIXELS` |
