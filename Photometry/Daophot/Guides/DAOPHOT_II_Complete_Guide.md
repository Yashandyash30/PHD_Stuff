# DAOPHOT II — Complete Practical Guide

> Based on Peter B. Stetson's official DAOPHOT II manual (Dominion Astrophysical Observatory).
> This guide covers every command, every parameter, what each value means, what to watch for, pitfalls to avoid, and what to do after the core reductions.

---

## Table of Contents

1. [Pre-flight Checklist — Before You Run DAOPHOT II](#1-pre-flight-checklist)
2. [Launching DAOPHOT II](#2-launching-daophot-ii)
3. [OPTIONS — Setting Reduction Parameters](#3-options)
4. [ATTACH — Loading Your Image](#4-attach)
5. [SKY — Estimating Sky Background](#5-sky)
6. [FIND — Detecting Stars](#6-find)
7. [PHOTOMETRY — Aperture Photometry](#7-photometry)
8. [PICK — Choosing PSF Candidate Stars](#8-pick)
9. [PSF — Building the Point-Spread Function](#9-psf)
10. [GROUP — Grouping Stars for Multi-Star Fitting](#10-group)
11. [NSTAR — Multiple Profile Fitting](#11-nstar)
12. [SUBSTAR — Subtracting Stars from the Image](#12-substar)
13. [ALLSTAR — All-in-One Alternative to GROUP+NSTAR+SUBSTAR](#13-allstar)
14. [Additional / Utility Commands](#14-additional-utility-commands)
15. [Full Reduction Workflow — Step by Step](#15-full-reduction-workflow)
16. [After DAOPHOT — What Comes Next](#16-after-daophot)
17. [Parameter Quick-Reference Table](#17-parameter-quick-reference-table)
18. [Common Pitfalls and How to Avoid Them](#18-common-pitfalls)

---

## 1. Pre-flight Checklist

Before launching DAOPHOT II, ensure the following are done **outside** of DAOPHOT:

### Image Preparation (do in IRAF / Midas / your pipeline)
- [ ] **Bias subtraction** — subtract the electronic zero level
- [ ] **Flat-fielding** — correct for pixel-to-pixel sensitivity variations
- [ ] **Trim** — remove bad rows/columns at frame edges
- [ ] Verify that pixel values are **linearly proportional to true intensity** (no non-linear regions in your data range)
- [ ] If combining multiple exposures, know **how** they were combined (average? sum? median?)

### Know These Numbers (write them down)
| Quantity | Where to find it | Why DAOPHOT needs it |
|---|---|---|
| **FWHM** (pixels) | Measure on a bright star | Optimizes the star-finding convolution kernel |
| **Read noise** (e⁻ or ADU per single frame) | Detector spec sheet / header | Noise model for detection and fitting |
| **Gain** (photons or e⁻ per ADU, per single frame) | Detector spec sheet / header | Noise model — Poisson statistics |
| **Saturation / linearity limit** (ADU) | Detector spec sheet | Sets `HIGH GOOD DATUM` |
| **Sky background level** (ADU) | Run `SKY` command | Sets detection threshold |

In daophot.opt we have to update following parmeters according to the FWHM
1. FW (FWHM)
2. FI ((equal to FW)Fitting radius (pixel))
3. PS ((3 times FI)PSF radius (in pixel))

In allstar.opt we have to update 
1. FI (fitting radius (equal to FWHN))
2. IS (Inner Sky Radius(6×FI))
3. OS (Outer Sky Radius( IS + 10))

In photo.opt we have to update

1. A1 (Aperture 1 radius (equal to FWHM))
2. A2-AB (A2,A3,A4,A5,A6,A7,A8,A9,AA,AB) (Intermediate apertures (each increasing geometrically by a multiple of ~1.177 from the previous))
3. AC (Maximum Aperture Radius (6×A1))
4. IS (Inner Sky Radius (equal to AC, or 6×A1))
5. OS (Outer Sky Radius (IS + 10))

Instead of hardcoding the constant and multiple factors can we declare these in configuration cell for each .opt files (in configuration take the default values of constant and multipiiers to that i have given here)

### Create a DAOPHOT.OPT File (Highly Recommended)
Place this file in **every working directory**. It is read automatically at startup, like a `.bashrc` for DAOPHOT. Minimum recommended contents:

```
RE = 6.50        # Read noise in ADU (for a single frame)
GA = 20.00       # Gain in photons/ADU (for a single frame)
FW = 3.5         # FWHM of stars in pixels
HI = 32000.0     # High good datum (saturation level in ADU)
FI = 3.0         # Fitting radius (≈ 1× FWHM is a good start)
PS = 11.0        # PSF radius (should exceed visible star halo)
```

> ⚠️ **Critical:** `RE` and `GA` default to **0**, which is deliberately invalid. If you forget to set them, DAOPHOT will nag you repeatedly and your noise model will be wrong.

---

## 2. Launching DAOPHOT II

```bash
daophot
```

On startup, DAOPHOT:
1. Reads `DAOPHOT.OPT` from your current directory (if it exists)
2. May display a "news" message — answer `Y` to read it if it has changed since your last run
3. Prints the current values of all optional parameters — **check them now**
4. Presents the `Command:` prompt

> Every command can be abbreviated to its **first two characters**:
> `FIND` → `FI`, `PHOTOMETRY` → `PH`, `OPTIONS` → `OP`, etc.
> Commands are **case-insensitive**.

---

## 3. OPTIONS

**Purpose:** Change any of the 19 user-definable reduction parameters.

```
Command: OP
```

DAOPHOT prints the current parameter table, then asks:

```
Parameter file (default KEYBOARD INPUT): <CR>
```

- Press `<CR>` to enter values from the keyboard
- Or provide a filename (same format as `DAOPHOT.OPT`) to read from a file

**Syntax for setting a parameter:**
```
OPT> FW=3.5
OPT> HS=0.9
OPT> <CR>     ← exit OPTIONS
```

Only the **first two characters** of the parameter name matter. Anything between the two-char ID and the `=` sign is ignored. So all of these are equivalent:
```
FW=3.0
fwhm = 3.0
Fwhatever = 3.0
```

Exit by pressing `<CR>` at an empty `OPT>` prompt, or `CTRL-Z` / `CTRL-D` (EOF on your system).

---

### Parameter-by-Parameter Explanation

#### (1) `RE` — Read Noise
```
RE = 6.50
```
- **Units:** ADU (data numbers), for **one single exposure**
- **What it is:** The electronic noise floor of your detector, expressed in ADU
- **If your frame is an average or sum of N frames:** DAOPHOT corrects for this automatically when FIND asks "Number of frames averaged, summed"
- **What to watch:** If you leave this at 0 (the default), your noise model is wrong and star detection thresholds will be meaningless
- **Typical values:** 2–15 ADU depending on detector

#### (2) `GA` — Gain
```
GA = 20.00
```
- **Units:** Photons (or electrons) per ADU, for **one single exposure**
- **What it is:** The conversion factor between detector counts (ADU) and actual photons detected; used to compute Poisson noise
- **What to watch:** Same as RE — defaults to 0, which breaks everything. Use the value for a **single raw frame**, not for your combined image
- **Typical values:** 1–30 e⁻/ADU

#### (3) `LO` — Low Good Datum
```
LO = 7.00
```
- **Units:** Standard deviations below the mean sky
- **What it is:** Any pixel more than this many σ below the sky level is flagged as **bad** (dead pixels, charge traps, cosmic ray shadows). Bad pixels are ignored by FIND and all subsequent routines
- **Flat background (Galactic field):** Use 5–7σ — a 5σ outlier in a 300×500 image has only ~1-in-20 chance of occurring naturally
- **Non-flat background (globular cluster, galaxy, nebula):** Use 10σ or more — otherwise real sky pixels near the bright center will be wrongly flagged as bad
- **What to watch:** If set too tight on a crowded field, PHOTOMETRY will refuse to proceed because fewer than 20 sky pixels survive rejection. This causes `σsky = −1`

#### (4) `HI` — High Good Datum
```
HI = 32000.0
```
- **Units:** Absolute ADU value (not relative to sky)
- **What it is:** Any pixel above this level is treated as saturated/non-linear and ignored. This is your detector's **linearity limit** in ADU
- **Important:** This is NOT 32767 (the integer maximum). After bias subtraction and flat-fielding, the effective saturation level will be somewhat lower. Check your detector's specification
- **What to watch:** Saturated stars will get `magnitude = 99.999` from PHOTOMETRY. They appear at the bottom of sorted lists, which PICK will correctly skip. Stars with **any bad pixel inside the fitting radius** get rejected by PSF

#### (5) `FW` — FWHM of Objects
```
FW = 3.50
```
- **Units:** Pixels
- **What it is:** The full-width at half-maximum of unresolved stellar profiles in your image. FIND convolves your image with a Gaussian of this width to detect stars
- **How to measure it:** Use your image display system; measure 5–10 moderately bright, isolated stars and take the average
- **What to watch:** FIND is most sensitive to stars whose FWHM matches this value. If your PSF varies significantly across the frame (e.g., due to guiding errors or field curvature), use a typical or median value
- **Range:** 0.2–15.0 pixels

#### (6) `TH` — Threshold
```
TH = 4.00
```
- **Units:** Standard deviations above the local sky background
- **What it is:** The minimum significance level for a brightness enhancement to be accepted as a probable star. FIND will reject anything fainter
- **Typical starting value:** 4σ — in a 300×500 frame this gives roughly 40 false positives from pure noise
- **Too low:** Floods your star list with noise spikes; massively increases CPU time for NSTAR (scales as ~4th power of star density); these ghost stars can migrate and corrupt real photometry
- **Too high:** You miss faint real stars
- **Best practice:** Run FIND, answer "N" to "Are you happy with this?", try a new threshold. Build a curve of N(detected) vs. threshold on paper — find the "elbow" where counts start exploding. Set threshold just above that

#### (7) `LS` — Low Sharpness Cutoff
```
LS = 0.20
```
- **What it is:** Minimum acceptable value of the SHARP index (see FIND section). Eliminates cold pixels (isolated dim pixels) that create spurious detections near them
- **Default is usually fine.** Only change if your PSF is highly unusual (extreme aberrations, severe undersampling)

#### (8) `HS` — High Sharpness Cutoff
```
HS = 1.00
```
- **What it is:** Maximum acceptable SHARP. Eliminates hot pixels (single bright pixels = cosmic rays)
- **SHARP > 1:** The object is sharper than a star = probably a cosmic ray or bad pixel
- **SHARP ∼ 0:** Probably a cold pixel artifact
- **Default is usually fine**

#### (9) `LR` — Low Roundness Cutoff
```
LR = -1.00
```
- **What it is:** Minimum ROUND. Rejects objects elongated along the **x-axis** (charge-bleed columns, bad rows)
- **ROUND < 0:** Elongated in x
- **ROUND > 0:** Elongated in y
- **Note:** ROUND only rejects objects elongated along cardinal directions. Objects elongated at oblique angles are not preferentially rejected

#### (10) `HR` — High Roundness Cutoff
```
HR = 1.00
```
- **What it is:** Maximum ROUND. Rejects objects elongated along the **y-axis**
- **When to change:** If your stars have genuine elongation (bad guiding in one direction), you may need to relax this. But be careful — you may let in real artifacts

#### (11) `WA` — Watch Progress
```
WA = 1.00
```
- **Values:**
  - `-2`: Silent — no terminal output at all. **Use this in batch jobs** (avoids filling logfiles, especially on Unix)
  - `-1`: Minimal output; PSF uses stars without prompting, just reports bad pixels
  - `0`: No real-time results on screen
  - `1`: Results printed to terminal as they are computed (default)
  - `2`: Results + gray-scale ASCII images of each star during PSF construction. Very slow; use only for debugging

#### (12) `FI` — Fitting Radius
```
FI = 3.00
```
- **Units:** Pixels
- **What it is:** The radius of the circular area around each star actually used in PSF fitting (PEAK, NSTAR) and in building the PSF itself
- **Rule of thumb:** ≈ 1× FWHM under normal conditions. For severe crowding, use slightly smaller (loses some accuracy but reduces neighbor contamination). For variable PSF across the frame, use slightly larger
- **Hard minimum:** 1.6 pixels (fewer than 4 pixels in the fit is forbidden by the program)
- **What to watch:** This parameter is written into file headers and flows forward into all subsequent routines. Setting it wrong at the beginning means re-running everything

#### (13) `PS` — PSF Radius
```
PS = 11.00
```
- **Units:** Pixels
- **What it is:** The radius within which the PSF look-up table is defined. This determines how far out the PSF profile extends
- **How to set it:** On your image display, look at a bright unsaturated star. The PSF radius should slightly **exceed** the radius at which the stellar profile merges into the background
- **Too small:** After subtraction, your bright stars will be surrounded by luminous rings with sharp inner edges (PSF halo not subtracted)
- **Too large:** CPU time for PSF fitting scales steeply with this value — don't make it unnecessarily large
- **Typical:** 2–4× FWHM is a reasonable starting point

#### (14) `VA` — Variable PSF
```
VA = 0.00
```
- **What it is:** Degree of spatial variation allowed in the PSF model
  - `-1`: Purely analytic PSF (no look-up table). Useful for the first pass in crowded fields since neighbors don't contaminate the model
  - `0`: Spatially constant PSF (analytic + look-up table, constant across frame). **Start here**
  - `1`: PSF varies linearly with position. Requires ≥3 PSF stars per degree of freedom (so ≥9 PSF stars total for a linear model with 3 DOF)
  - `2`: PSF varies quadratically. Requires ≥6 PSF stars per DOF
- **What to watch:** Don't use VA=1 or VA=2 until you have enough well-distributed PSF stars to actually constrain the variation. An under-constrained variable PSF is worse than a constant one
- **Workflow tip:** Start with VA=-1 (pure analytic) in the first PSF pass on a crowded field. Increase to 0, then 1, then 2 as you iterate and the frame cleans up

#### (15) `FR` — Fractional Pixel Expansion
- **Not implemented.** Leave it at 0.

#### (16) `AN` — Analytic Model PSF
```
AN = 1.00
```
- **What it is:** Which analytic function to use as the first-order approximation to the PSF before the empirical corrections are added
  - `1`: **Gaussian** — 2 free parameters (half-width in x and y). Fastest. Adequate for well-sampled ground-based data
  - `2`: **Moffat (β=1.5)** — 3 free parameters (half-widths + position angle of ellipse). Better for under-sampled data or broader wings
  - `3`: **Moffat (β=2.5)** — 3 free parameters. Even broader wings
  - `4`: **Lorentz** — 3 free parameters. Very broad wings
  - `5`: **Penny 1** — Gaussian + Lorentz, 4 free parameters. Good for HST-like PSFs
  - `6`: **Penny 2** — Gaussian + tilted Lorentz, 5 free parameters. Most flexible
- **When to change from Gaussian:** If your Chi value from the PSF output is large (>0.05 or so) and your PSF has notably non-Gaussian wings (check on the image display), try Moffat or Penny functions

#### (17) `EX` — Extra PSF Cleaning Passes
```
EX = 5.00
```
- **What it is:** Number of times PSF re-examines each pixel and downweights obvious outliers (cosmic rays in PSF stars, imperfectly subtracted neighbors)
- `0`: Every pixel taken at full face value — not recommended unless field is extremely sparse
- `5`: A reasonable default; PSF has usually converged by 5 passes
- **What to watch:** Each pass adds time. Values above 5 rarely help but do cost CPU

#### (18) `PE` — Percent Error
```
PE = 0.75
```
- **Units:** Percent
- **What it is:** Uncertainty contribution from flat-field fine structure — the "graininess" of the mismatch between your flat-field and your science frames, on sub-seeing-disk scales. This error grows linearly with intensity (star + sky)
- **Leave at default** (0.75%) until you've characterized your flat-field quality. Only adjust if you have strong reason to believe your flat-field is worse or better than typical

#### (19) `PR` — Profile Error
```
PR = 5.00
```
- **Units:** Percent
- **What it is:** Uncertainty from the PSF not being perfectly known — interpolation errors, residual PSF variation not captured by the model. Grows linearly with star brightness and inversely as FWHM⁴ (becomes more important with better seeing)
- **Leave at default** until late in your analysis when you are trying to characterize your photometric error budget

---

## 4. ATTACH

**Purpose:** Tell DAOPHOT which image file to work on.

```
Command: AT filename
```

Or just:
```
Command: AT
Enter file name: filename
```

- DAOPHOT prints the image header comment (if any) and the picture dimensions (`nnn × nnn`)
- File extension (`.DST`, `.imh`, `.bdf`) may be **omitted** if it is the standard extension for your system
- On **Unix:** filenames are case-sensitive. `OBS137` ≠ `obs137`
- On **Unix:** use the VMS-style directory alias syntax: if `setenv ccd /path/to/data` is in your `.cshrc`, you can write `ccd:obs137`
- **Filename limit:** 30 characters
- If you try to write a file that already exists on Unix, DAOPHOT warns you and gives you a chance to rename the output or overwrite. In automated scripts, use `rm` beforehand to prevent unexpected prompts

> Commands that operate only on data files (`SORT`, `OFFSET`, `APPEND`, `OPTIONS`, `MONITOR`) can be run **without** a prior ATTACH. Commands that need the image (`FIND`, `PHOTOMETRY`, `PSF`, `NSTAR`, `SUBSTAR`) will fail without it.

---

## 5. SKY

**Purpose:** Estimate the average sky brightness of the currently attached frame.

```
Command: SK
```

**Output:**
```
Approximate sky value for this frame = 156.8
Standard deviation of sky brightness = 4.16

Clipped mean and median = 157.9  157.5
Number of pixels used (after clip) = 8673
```

- Samples ~10,000 uniformly scattered pixels, sigma-clips the tails, then estimates the **mode** = 3×median − 2×mean
- The standard deviation is the width of the sky histogram after clipping — NOT the standard deviation of the mode
- **Note the mode value** — you'll need it for your notes and to sanity-check the FIND threshold
- **This step is optional** — FIND runs SKY automatically if you skip it
- **When to use SKY explicitly:** If you want to understand your frame's noise properties before committing to detection parameters

---

## 6. FIND

**Purpose:** Detect star-like objects above a significance threshold using matched-filter convolution.

```
Command: FI

Approximate sky value for this frame = 156.8
Standard deviation of sky brightness = 4.16

Relative error = 1.14

Number of frames averaged, summed: 5,1
```

### The "frames averaged, summed" question

This is one of the most important and confusing inputs. DAOPHOT needs this to correctly adjust the noise model for your image.

| Scenario | Entry | Meaning |
|---|---|---|
| Single exposure | `1,1` | 1 average of 1 frame |
| Average of N frames | `N,1` | e.g., average of 5 → `5,1` |
| Sum of N frames | `1,N` | e.g., sum of 5 → `1,5` |
| Average of K groups each summed N times | `K,N` | |
| Median of N frames | `0.667N,1` | Median ≈ average of 2/3 as many. Median of 3 → `2,1`. Median of 5 → `3.3,1` |

The program uses your stated `RE` and `GA` plus the averaging/summing scheme to compute the correct effective noise per pixel. Getting this wrong means your thresholds are wrong.

### Output prompt

```
File for the positions (default ?.COO): <CR>
```

Accept the default `.COO` filename — keeping default filenames reduces bookkeeping errors. Press `<CR>`.

FIND then counts through the image rows, performing the convolution. It reports candidate stars and computes:

### SHARP Statistic

```
SHARP = height of best-fitting delta function / height of best-fitting Gaussian
```

| SHARP value | Interpretation |
|---|---|
| ∼0 | Cold pixel artifact (isolated dark pixel creates spurious peaks nearby) |
| 0.2–1.0 | Normal star (default acceptance range) |
| >1.0 | Hot pixel / cosmic ray (sharper than any real star) |

### ROUND Statistic

```
ROUND = (height_x_Gaussian − height_y_Gaussian) / average(height_x, height_y)
```

| ROUND value | Interpretation |
|---|---|
| < −1 | Object elongated strongly in x-direction |
| −1.0 to +1.0 | Acceptable round-ish object (default range) |
| > +1 | Object elongated strongly in y-direction (e.g., charge bleed) |

> ⚠️ ROUND only catches objects elongated **along rows or columns** — diagonally elongated objects pass through undetected.

### Final question

```
Are you happy with this? Y
```

- Answer `Y` to accept and write the `.COO` file
- Answer `N` to try a new threshold **without re-running the convolution** (fast — use this to experiment)
- Build a curve of N(detected) vs threshold; look for the "elbow" and set threshold just above it

### Output file: `.COO`
```
ID   X        Y        Mag      Sharp    Round1   Round2
1    6.953    2.612    -0.053   0.333   -0.546   -0.532
```
- Magnitudes are **relative to the threshold** — a value of 0.000 means exactly at the detection limit. Stars brighter than threshold have **negative** magnitudes (FIND convention only)

### Things to be careful of with FIND

1. **Don't set threshold too low** — ghost detections cost CPU time that scales roughly as the **4th power** of surface density, and they corrupt real photometry by wandering around the frame
2. Stars at the image edge have positions returned to the **nearest pixel only**, and their ROUND index is not computed
3. Badly saturated stars may be detected multiple times. The program "continues to have some trouble with grossly saturated objects"
4. Closely blended pairs (separation < 1 FWHM) will often be detected as one object — handle with the iterative subtract-and-re-find approach (see Section 15)

---

## 7. PHOTOMETRY

**Purpose:** Compute sky brightness and concentric aperture photometry for all detected stars.

### Step 1 — Create or modify PHOTO.OPT

A template `PHOTO.OPT` at the DAO looks like:
```
A1 = 3
A2 = 4
A3 = 5
A4 = 6
A5 = 7
A6 = 8
A7 = 10
IS = 20
OS = 35
```

### Step 2 — Run PHOTOMETRY

```
Command: PH

Enter table name (default PHOTO.OPT): <CR>
```

The parameter table is displayed:
```
A1 RADIUS OF APERTURE 1  =  3.00   A2 RADIUS OF APERTURE 2  =  4.00
...
IS INNER SKY RADIUS      = 20.00   OS OUTER SKY RADIUS       = 35.00
PHO>
```

At the `PHO>` prompt, modify any values:
```
PHO> A1=2.5
PHO> IS=10
PHO> OS=20
PHO> <CR>
```

Then:
```
File with the positions (default ?.COO): <CR>
File for the magnitudes (default ?.AP): <CR>
```

### Aperture and Sky Parameters

| Parameter | Meaning | How to choose |
|---|---|---|
| `A1`–`A12` | Radii of up to 12 concentric apertures in pixels | Span from ∼FWHM to several×FWHM for growth-curve analysis |
| `IS` | Inner sky annulus radius (pixels) | Should clear the faint outer wings of the PSF; typically 3–5× FWHM |
| `OS` | Outer sky annulus radius (pixels) | Large enough for statistical sky estimate; typically IS + 10–20 pixels |

- The **first aperture** (`A1`) magnitude is used as the starting guess for PSF fitting — choose it to be a meaningful measurement (∼1× FWHM radius is common)
- A zero or negative radius terminates the aperture list
- If using DAOGROW for growth-curve analysis, aperture radii **must increase monotonically**
- Up to **12 apertures** are allowed

### Sky estimation method

DAOPHOT uses the **mode** (not mean or median) of pixel values in the sky annulus. The mode = 3×median − 2×mean after iterative sigma-clipping. This is a maximum-likelihood estimate of the "most probable sky value" — not just the diffuse background but the full contribution of all unresolved sources in that region.

### Output: `.AP` file

Two lines per star:
```
Line 1: ID   X       Y       Mag_A1  Mag_A2  Mag_A3  Mag_A4 ...
Line 2:      Sky_mode  Sky_σ  Sky_skew  σ_A1  σ_A2  σ_A3  σ_A4 ...
```

### When PHOTOMETRY returns 99.999 ± 9.999

A star gets this flag if:
- The aperture extends **outside the frame boundary**
- A **bad pixel** (above HI or below LO) falls inside the star aperture
- Star + sky is fainter than sky alone (negative star brightness)
- Sky mode could not be determined (flagged additionally by `σsky = −1`)

### Why PHOTOMETRY might abort: "fewer than 20 sky pixels"

If after sigma-clipping the sky annulus fewer than 20 pixels remain, PHOTOMETRY refuses to proceed. Causes:
- Sky annulus too small
- LO set too tight — real sky pixels flagged as bad
- Short exposures where sky ≈ read noise → many pixels falsely flagged as bad. Solution: set LO to a moderate-sized **negative** value (e.g., −5σ) for short-exposure frames

---

## 8. PICK

**Purpose:** Automatically select good PSF star candidates from the aperture photometry file.

```
Command: PICK

Input file name (default ?.AP): <CR>
Desired number of stars, faintest magnitude: 25 18.0
Output file name (default ?.LST): <CR>
```

### How PICK works

1. Sorts stars by apparent magnitude (brightest first)
2. Rejects stars within **one FITTING RADIUS** of the frame edge (profile would be truncated)
3. Rejects stars within **one PSF RADIUS + one FITTING RADIUS** of any brighter or saturated star (neighbor contamination)
4. Writes surviving candidates in order of increasing apparent magnitude up to your requested number

### Inputs to specify

- **Desired number of stars:** Depends on your `VA` setting:
  - VA=0 (constant PSF): minimum 3, but more is better; ≥10 recommended
  - VA=1 (linear): minimum 9; ≥20–25 recommended, distributed across the frame
  - VA=2 (quadratic): minimum 18; ≥30–50 recommended
- **Faintest magnitude:** Enter `99` if you don't know yet. Stars fainter than your limit are included only if you'd otherwise have too few to constrain the PSF variation degree you selected

### What to watch

- Even after PICK, **visually inspect** your PSF star candidates on the image display before committing. PICK uses simple geometric criteria — it cannot check for cosmic rays, asymmetric halos, or diffraction spikes
- Stars that are too bright (near saturation) will have `mag = 99.999` from PHOTOMETRY and will appear at the end of the sorted list; PICK skips them
- Stars flagged with `*` or `?` in later PSF output (scatter > 3× or 2× average respectively) should be removed from the `.LST` file and PSF re-run

---

## 9. PSF

**Purpose:** Derive the point-spread function for the frame by fitting an analytic model plus empirical corrections to a set of PSF stars.

```
Command: PSF

File with aperture results (default ?.AP): <CR>
File with PSF stars (default ?.LST): <CR>
File for the PSF (default ?.PSF): <CR>
```

### What PSF does

1. If `WA ≥ 1`: Shows each PSF star as a gray-scale ASCII image and asks `Y/N` whether to include it — **examine each carefully**
2. Fits the chosen analytic function (set by `AN`) simultaneously to all PSF stars' inner profiles (within one fitting radius)
3. Reports the **Chi** value = RMS residual / peak height (fraction). This tells you how well your analytic model approximates the true PSF
4. Reports per-star profile scatter; flags outliers with `*` (>3× mean) and `?` (>2× mean)
5. Builds an empirical look-up table of corrections from the analytic model to the observed profile
6. Writes a `.NEI` file listing PSF stars and their recognized neighbors

### Interpreting the output

```
Chi Parameters...
0.0282  0.71847  0.83218  0.30568

Profile errors:
180 0.023   555 0.026   122 0.026  ...
456 0.096 *     (flagged as outlier)
```

| What you see | What it means | What to do |
|---|---|---|
| Chi ∼ 0.01–0.03 | Excellent fit | Proceed |
| Chi ∼ 0.03–0.08 | Acceptable | Check analytic model; try Moffat if Gaussian |
| Chi > 0.10 | Poor fit | Try different AN value; check for blends; consider more PSF stars |
| Stars flagged `*` or `?` | High individual scatter | Remove from `.LST` and re-run PSF |

### After PSF

The `.NEI` file (PSF stars + their neighbors) should be run through **ALLSTAR** or **GROUP + NSTAR** to fit out the neighbors, then **SUBSTAR** used to remove those neighbors from the original image, leaving the PSF stars relatively isolated for the next iteration (see Section 15).

### Important: PSF is iterative

For crowded fields, **one PSF pass is almost never sufficient**. Plan on 3–5 iterations of:
1. Make PSF → 2. Reduce neighbors → 3. Subtract neighbors → 4. New PSF on cleaner image → repeat

Use `VA=-1` (pure analytic) on the first pass to avoid neighbor ghosts contaminating the look-up table. Increase VA progressively as the field cleans up.

---

## 10. GROUP

**Purpose:** Divide the star list into natural groups for simultaneous multi-star profile fitting with NSTAR. Two stars belong to the same group if one can significantly influence the other's fit.

```
Command: GR

File with the photometry (default ?.AP): <CR>
File with the PSF (default ?.PSF): <CR>
Critical overlap: 0.5
File for stellar groups (default ?.GRP): <CR>
```

### The "Critical Overlap" parameter

For two stars i (fainter) and j (brighter):
- GROUP calculates the brightness of star j (using the PSF model) at a distance of `FITTING RADIUS + 1 pixel` from the centroid of star i
- If this brightness > `critical_overlap × random_noise_per_pixel` → they are grouped together

| Critical overlap | Regime | Consequence |
|---|---|---|
| 0.1 | Crowding negligible relative to random noise | All groups < 60 stars; proceed with NSTAR |
| 1.0 | Crowding ≈ random noise | Some large groups possible |
| >> 1.0 | Crowding dominates noise | Many stars cannot be well-measured |

### Interpreting the group-size table

After GROUP runs, a table appears:
```
Group size   Number of groups
1            150
2            43
3            20
...
60           2
87           1    ← problem!
```

NSTAR can handle **maximum 60 stars per group**. If any group exceeds this:
- Increase critical overlap so that weakly interacting pairs are separated → smaller groups
- Or use **SELECT** to carve out the oversized groups, re-run GROUP with larger critical overlap on just those stars
- Or — best choice — just use **ALLSTAR**, which has no group size limit

---

## 11. NSTAR

**Purpose:** Multiple simultaneous PSF profile fitting for stars within the same group.

```
Command: NS

File with the PSF (default ?.PSF): <CR>
File with stellar groups (default ?.GRP): <CR>
File for NSTAR results (default ?.NST): <CR>
```

### How NSTAR works

- Simultaneously fits all stars in a group by iterative nonlinear least squares
- Each iteration updates all positions and magnitudes together, properly accounting for neighbor contamination
- Converges when: magnitude changes < 0.0001 mag or 0.05σ (whichever is larger) AND centroid changes < 0.001 pixel, for **all** stars in the group
- Gives up after **50 iterations** (but results may still be usable)
- Maximum **60 stars per group**

### Star rejection logic in NSTAR

NSTAR automatically eliminates stars during iteration that are:

1. **Too close together** (< 0.37 × FWHM separation) → positions merged; fainter star eliminated
2. **Too faint** (> 12.5 mag fainter than PSF reference star; < 0.2 ADU/pixel at center)
3. **Insignificant** (brightness < 1σ above zero after iteration 5; < 1.5σ after iter 10; < 2.0σ after iter 15)
4. **Faint and close** (within 1.0× FWHM and fainter than the significance threshold)

When a star is eliminated, the iteration counter is backed up by 1 so the remaining stars can settle.

> ⚠️ The rejection thresholds depend on your `RE` and `GA` values. Wrong noise parameters mean wrong rejection decisions.

### Output: `.NST` file

One line per star:
```
ID   X(improved)  Y(improved)  Mag   σ_Mag  Sky  N_iter  CHI   SHARP
```

### Interpreting CHI and SHARP in output

**CHI** — ratio of observed pixel scatter to expected pixel scatter from noise model:
- Most stars should cluster around **CHI ∼ 1.0**
- Systematic trend (CHI rising with brightness) → `PE` and/or `PR` too low; or RE/GA incorrect
- Individual stars with very high CHI → unrecognized blends, cosmic rays, or bad pixels surviving through FIND
- Stars with low CHI at bright end → those stars may be saturated

**SHARP** — diagnostic of intrinsic angular size:
- Should be **near 0** for point sources
- Significantly > 0 → probably a galaxy or unresolved double star
- Significantly < 0 → probably a cosmic ray or image artifact that survived FIND

---

## 12. SUBSTAR

**Purpose:** Subtract fitted star profiles from the image to create a residual image.

### Subtract ALL stars:

```
Command: SU

File with the PSF (default ?.PSF): <CR>
File with photometry (default ?.NST): <CR>
Do you have stars to leave in? N
Name for subtracted image (default ?s.DST): <CR>
```

### Subtract ALL EXCEPT selected stars (leave PSF stars in, for re-PSF):

```
Command: SU

File with the PSF (default ?.PSF): <CR>
File with photometry (default ?.NST): <CR>
Do you have stars to leave in? Y
File with star list (default ?.LST): <CR>    ← typically your PSF stars list
Name for subtracted image (default ?s.DST): <CR>
```

- Stars are matched by **ID number** — if you have renumbered or APPENDed files where IDs are duplicated, wrong stars may be left in. Be careful.
- The original image is **not modified** — a new image file is created
- Uses the positions and magnitudes from the photometry file, not from any fits performed live

### When to use SUBSTAR

1. **Quality check:** Examine the subtracted image on your display. Well-fit stars should leave flat residuals. Systematic patterns → PSF problems
2. **Second-pass FIND:** Subtract all known stars, then run FIND on the residual to find previously hidden stars
3. **PSF improvement:** Subtract all neighbors (not PSF stars) from the image; re-derive PSF from the cleaner image
4. **Aperture photometry zero-point:** Subtract all stars except a few bright isolated ones; do large-aperture photometry on those isolated stars

---

## 13. ALLSTAR

**Important:** ALLSTAR is a **separate stand-alone program**, not a command within DAOPHOT. Run it from the operating system prompt.

```bash
allstar
```

ALLSTAR displays and lets you modify its parameter table, then asks:

```
Input image name: filename
File with the PSF (default ?.PSF): <CR>
Input file (default ?.AP): <CR>         ← can feed .AP directly, no GROUP needed
File for results (default ?.ALS): <CR>
Name for subtracted image (default ?s.DST): <CR>   ← CTRL-Z to skip
```

### ALLSTAR vs NSTAR: Key differences

| Feature | NSTAR | ALLSTAR |
|---|---|---|
| Input requirement | Must be pre-grouped (.GRP) | Any .AP file directly |
| Star limit | 60 per group | 15,000 total |
| Group size limit | Hard limit of 60 | No hard limit; auto-manages internally |
| Subtracted image | Need separate SUBSTAR | Produced automatically |
| Sky re-determination | No | Yes (if OS > IS > 0) |
| Speed on sparse fields | Similar | Similar |
| Speed on crowded fields | Slow (many GROUP+SELECT iterations) | Faster overall |
| Recommendation | Historical; use for specific NSTAR features | **Preferred for most use cases** |

### ALLSTAR-specific parameters

```
FITTING RADIUS = 2.50
REDETERMINE CENTROIDS = 1.00    ← 1=yes, 0=no (fix positions from input)
WATCH PROGRESS = 1.00
MAXIMUM GROUP SIZE = 50.00      ← internal auto-grouping maximum
PERCENT ERROR = 0.75
PROFILE ERROR = 5.00
IS (INNER SKY RADIUS) = 0.00   ← set > 0 to enable sky re-determination
OS (OUTER SKY RADIUS) = 0.00   ← set > IS to enable sky re-determination
CLIPPING EXPONENT = 6.00        ← b in bad-pixel rejection formula
CLIPPING RANGE = 2.50           ← a in bad-pixel rejection formula
```

**REDETERMINE CENTROIDS:**
- `1`: Full astrometric + photometric solution (normal mode)
- `0`: Fix positions from input file; solve only for magnitudes. Use when positions from a high-quality frame are being transferred to poorer frames. **Caution:** Different bandpasses or airmasses introduce positional shifts, rotations, and atmospheric dispersion that a simple zero-point offset doesn't capture

**Sky re-determination (IS and OS):**
- Set `OS > IS > 0` to have ALLSTAR iteratively re-determine sky brightness for each star every 3rd iteration
- ALLSTAR subtracts all known stars from a working copy, then takes the **median** of the annulus around each position
- This is "the best way currently available to get at the sky under the star"
- Recommended: `IS` ∼ fitting radius or smaller; `OS` ∼ PSF radius or somewhat larger

**Clipping parameters:**
- `CLIPPING RANGE (a)`: A residual of a×σ gets half weight in the fit
- `CLIPPING EXPONENT (b)`: Must be a non-negative integer; higher = sharper clipping cutoff
- Set CLIPPING EXPONENT = 0 to disable clipping entirely (conservative mode)
- Default (RANGE=2.5, EXP=6) is a reasonable starting point

**Create ALLSTAR.OPT** in your working directory (same format as DAOPHOT.OPT) to set your preferred ALLSTAR defaults automatically.

---

## 14. Additional / Utility Commands

### MONITOR / NOMONITOR

Toggle real-time terminal output without going through OPTIONS:

```
Command: NO    ← suppress output (sets WA=0)
Command: MO    ← restore output (sets WA=1)
```

For WA values of -2, -1, or 2, use OPTIONS explicitly.

### SORT

Reorder any DAOPHOT output file:

```
Command: SO

Which do you want?
  +/- 1   By increasing/decreasing star ID
  +/- 2   By increasing/decreasing X coordinate
  +/- 3   By increasing/decreasing Y coordinate
  +/- 4   By increasing/decreasing magnitude
  +/- n   By increasing/decreasing column n (n ≤ 30)

Which do you want? 4
Input file name: file.AP
Output file name (default ?): <CR>
Do you want the stars renumbered? Y
```

- `+4` = sort by increasing magnitude (faintest first)
- `-4` = sort by decreasing magnitude (brightest first)
- `+/- n` allows sorting by any column — see Appendix IV for column numbers
- Renumbering: first star in output gets ID=1, second gets ID=2, etc.
- ⚠️ Sorting a `.GRP` file destroys the group structure

### SELECT

Extract groups within a specific size range from a `.GRP` file:

```
Command: SE

Input group file: file.GRP
Minimum, maximum group size: 1 60
Output group file (default ?.GRP): <CR>
```

- Use twice: once for groups 1–60 (feed to NSTAR), once for groups 61–9999 (re-run GROUP with larger critical overlap on these)
- Repeat until all stars are in manageable groups, or just use ALLSTAR

### OFFSET

Add constants to ID numbers, X, Y coordinates, or magnitudes in any data file:

```
Command: OF

Input file name: file.COO
Additive offsets ID, DX, DY, DMAG: 50000  0  0  0
Output file name (default ?.OFF): <CR>
```

- Numbers can be separated by spaces or commas; terminate with `/` to set remaining offsets to zero
- **Why add to ID?** After running FIND on a subtracted image for a second pass, add e.g. 50000 to new-star IDs before APPENDing, so no ID collisions with first-pass stars
- **Why add to magnitudes?** To shift artificial star magnitudes when inserting the same stars into multiple frames taken at different exposures

### APPEND

Concatenate two DAOPHOT data files, handling headers correctly:

```
Command: AP

First input file: file1.GRP
Second input file: file2.GRP
Output file (default ?.CMB): <CR>
```

- ⚠️ DAOPHOT does NOT check that the two files are the same type — appending mismatched types (e.g., .COO onto .AP) produces a corrupt file that will crash subsequent routines
- GROUP leaves a blank line at the end of every group, so groups remain distinct after APPENDing. If you edit group files and delete that last blank line, the last group of file 1 and first group of file 2 will merge

### DUMP

Display raw pixel values in a rectangular subregion:

```
Command: DU

Box size: 9
Coordinates of central pixel: 201,378
```

Output is a labeled grid of pixel values. Enter `0,0` or CTRL-Z to exit.

- Useful for: checking if pixel values make sense; verifying saturation levels; diagnosing why DAOPHOT isn't working on your data (e.g., all pixels = -32768)
- Maximum box that fits an 80-column terminal: **12×12**. Maximum for 132-column terminal: **21×21**

### FUDGE

Replace pixels in a rectangular region with interpolated or constant values:

```
Command: FU

Name for output picture (default ?f.DST): <CR>
Border (pixels): 3
Polynomial order (0=constant, 1=plane, etc.): 1
First, last column number: 267,268
First, last row number: 381,381
First, last column number: CTRL-Z    ← exit FUDGE
```

- **Constant value mode:** specify 0 for border, then enter a brightness value directly
- **Interpolated mode:** specify a border width and polynomial order. For each pixel inside the rectangle, a different polynomial is fit to the border pixels using 1/r⁴ weights → smooth surface that flows naturally between the borders
- **Use case:** Cosmic ray hitting a PSF star's outer wings. FUDGE the bad pixels, derive a clean PSF, then **delete the fudged image** before publication

> ⚠️ "Although it is morally wrong..." — Use this sparingly and only when there is truly no alternative. Document it in your reduction notes.

### ADDSTAR

Add synthetic stars to your image for completeness testing:

**Random placement:**
```
Command: AD

File with the PSF (default ?.PSF): <CR>
Seed (any integer): 12345
Photons per ADU: 20.0
Input data file (default RANDOM STARS): <CR>
Magnitude of PSF star is 12.033
Minimum, maximum magnitudes desired: 15.0 18.0
Number of stars to add to each frame: 100
Number of new frames: 5
File-name stem: FAKE
```

Creates: `FAKE01.DST` through `FAKE05.DST` (images) and `FAKE01.ADD` through `FAKE05.ADD` (input star positions/magnitudes)

**Deliberate placement:**
```
Input data file (default RANDOM STARS): ARTSTAR.MAG
```

Where `ARTSTAR.MAG` is a file you created with custom positions and magnitudes.

- **Seed:** Use the same seed on any machine to reproduce identical results
- **Photons per ADU:** Used to add correct Poisson noise to synthetic stars. Use the correct gain. For noiseless stars (testing only): use 99999
- **Best practice:** Add only a few stars per frame (∼ 5–10% of detected star count) to avoid artificially crowding the frame. Run many frames with different random seeds

### LIST (VMS only)

Examine image header fields (not implemented on Unix):
```
Command: LI
LIST> OBS.RA
LIST> OBS.DEC
```

### HELP

Print a list of all recognized commands:
```
Command: HE
```

### EXIT

Exit cleanly from DAOPHOT:
```
Command: EX
```

---

## 15. Full Reduction Workflow — Step by Step

### Pass 1: First-pass detection and photometry

```bash
# 1. Start DAOPHOT
daophot

# 2. Set parameters (or rely on DAOPHOT.OPT)
Command: OP
OPT> RE=6.5
OPT> GA=20.0
OPT> FW=3.5
OPT> HI=32000
OPT> FI=3.0
OPT> PS=11.0
OPT> <CR>

# 3. Load image
Command: AT myimage

# 4. Optional: check sky
Command: SK

# 5. Find stars
Command: FI
Number of frames averaged, summed: 1,1     [adjust as needed]
File for the positions (default ?.COO): <CR>
Are you happy with this? Y

# 6. Aperture photometry
Command: PH
PHO> <CR>
File with the positions (default ?.COO): <CR>
File for the magnitudes (default ?.AP): <CR>

# 7. Sort by magnitude (optional but convenient)
Command: SO
Which do you want? 4
Input file: myimage.AP
Do you want the stars renumbered? Y

# 8. Select PSF candidates
Command: PICK
Desired number of stars, faintest magnitude: 30 18
Output file name (default ?.LST): <CR>
```

### Pass 1: Build the PSF (iterative)

```bash
# First PSF — pure analytic to avoid neighbor contamination
Command: OP
OPT> VA=-1
OPT> <CR>

Command: PSF
File with aperture results (default ?.AP): <CR>
File with PSF stars (default ?.LST): <CR>
File for the PSF (default ?.PSF): <CR>
[Answer Y/N for each star shown]
[Note: .NEI file is created]

# Reduce PSF stars and neighbors
Command: EX
```

Now run ALLSTAR on the `.NEI` file:
```bash
allstar
Input image name: myimage
File with the PSF: myimage.PSF
Input file: myimage.NEI
File for results: myimage_nei.ALS
Name for subtracted image: myimage_neis.DST
```

Back in DAOPHOT:
```bash
daophot

# Subtract neighbors, leave PSF stars in
Command: AT myimage
Command: SU
File with the PSF: myimage.PSF
File with photometry: myimage_nei.ALS
Do you have stars to leave in? Y
File with star list: myimage.LST
Name for subtracted image: myimage_clean.DST

# Now build improved PSF from cleaner image
Command: OP
OPT> VA=0     [now use full look-up table]
OPT> <CR>

Command: AT myimage_clean.DST
Command: PSF
[Repeat until PSF residuals are flat and Chi is minimized]
```

### Pass 1: Full reduction

```bash
# Run ALLSTAR on complete star list
allstar
Input image name: myimage
File with the PSF: myimage.PSF
Input file: myimage.AP
File for results: myimage.ALS
Name for subtracted image: myimage_s.DST
```

### Pass 2: Find hidden stars in subtracted image

```bash
daophot
Command: AT myimage_s.DST    [the subtracted image]
Command: FI
Number of frames averaged, summed: 1,1
File for the positions: myimage2.COO
Are you happy with this? Y

# Get photometry for new stars from ORIGINAL image
Command: AT myimage             [back to original!]
Command: PH
File with the positions: myimage2.COO
File for the magnitudes: myimage2.AP

# Add 50000 to new-star IDs to avoid collisions
Command: OF
Input file: myimage2.AP
Additive offsets ID, DX, DY, DMAG: 50000 0 0 0/
Output file: myimage2_offset.AP

# Combine old and new star lists
Command: AP
First input file: myimage.AP
Second input file: myimage2_offset.AP
Output file: myimage_combined.AP

# Run ALLSTAR again on combined list
Command: EX
```

```bash
allstar
Input image name: myimage
File with the PSF: myimage.PSF
Input file: myimage_combined.AP
File for results: myimage_combined.ALS
Name for subtracted image: myimage_combined_s.DST
```

### Pass 3 (optional): Visual inspection pass

```bash
# Display myimage_combined_s.DST on your image display system
# Add any hand-identified stars to a new .COO file: missed_stars.COO

daophot
Command: AT myimage
Command: PH
File with the positions: missed_stars.COO
File for the magnitudes: missed_stars.AP

Command: OF
Input file: missed_stars.AP
Additive offsets ID, DX, DY, DMAG: 100000 0 0 0/
Output file: missed_stars_offset.AP

Command: AP
First: myimage_combined.AP
Second: missed_stars_offset.AP
Output: myimage_final.AP

Command: EX
```

```bash
allstar
Input file: myimage_final.AP
File for results: myimage_final.ALS
Name for subtracted image: myimage_final_s.DST
```

### Establishing the photometric zero-point

```bash
# In the final subtracted image, identify a few bright, isolated, unsaturated stars
# Use an editor to delete those bright stars from the .ALS file → myimage_nozp.ALS
# Or: create a file containing ONLY those bright stars → myimage_zpstars.LST

daophot
Command: AT myimage
Command: SU
File with photometry: myimage_final.ALS
Do you have stars to leave in? Y
File with star list: myimage_zpstars.LST   [stars to KEEP — subtract everything else]
Name for subtracted image: myimage_zp.DST

Command: AT myimage_zp.DST
Command: PH
[Use many apertures, A1=2 through A12=30 or whatever the frame allows]
File with the positions: myimage_zpstars.COO
File for the magnitudes: myimage_zpstars.AP
Command: EX
```

The growth curve from these aperture magnitudes establishes your aperture correction and absolute zero-point.

---

## 16. After DAOPHOT — What Comes Next

### 1. Aperture Correction (Growth Curve Analysis)

Profile-fitting magnitudes are measured within a finite fitting radius. To place them on an absolute scale:
- Use the isolated bright stars measured with large apertures (step above)
- Compute the difference: `aperture_mag(large) − aperture_mag(small)` as a function of aperture radius — the **growth curve**
- The aperture correction to apply to all PSF-fit magnitudes = the asymptotic value of the growth curve
- Tools: **DAOGROW** (Stetson 1990), or your own growth-curve analysis script

### 2. Calibration to Standard System

DAOPHOT produces **instrumental magnitudes** (relative to your PSF star's brightness). You need:
- Standard star frames reduced in the same way
- Transformation equations accounting for: zero-point offset, color terms, extinction
- Tools: **DAOMASTER** and **DAOMATCH** for frame-to-frame matching; **CCDSTD** or similar for atmospheric calibration

### 3. Frame Matching with DAOMASTER / DAOMATCH

If you have multiple frames of the same field (multiple epochs, filters, exposures):
- **DAOMATCH**: finds the coordinate transformation (shift + rotation + scale) between frames
- **DAOMASTER**: combines matched star lists across frames, computes mean magnitudes and coordinates
- Input: your `.ALS` or `.NST` files from all frames

### 4. Artificial Star Tests (Completeness)

Use ADDSTAR to determine:
- **Detection completeness** as a function of magnitude: what fraction of added stars does FIND recover?
- **Photometric bias**: is the recovered magnitude systematically offset from the input?
- **Photometric scatter**: what is the actual scatter in recovered magnitudes vs. input magnitude?

Best practice:
- Add only ∼5–10% of the detected star count per frame
- Use many frames with different seeds (∼20–100 frames) for statistical robustness
- Reduce each artificial frame identically to the science frame

### 5. Plotting Diagnostics (Do These Before Publishing)

Plot for every star from your final `.ALS` output:

| Plot | What to look for |
|---|---|
| CHI vs. magnitude | Should cluster around 1.0 with no trend. Rising CHI at bright end = saturation. Rising CHI at all magnitudes = wrong RE or GA |
| SHARP vs. magnitude | Point stars cluster near 0. High SHARP = galaxies or blends. Low SHARP = cosmic rays or artifacts |
| Magnitude residual vs. position (X, Y) | Check for spatially coherent PSF errors = may need higher VA |
| σ_mag vs. magnitude | Should follow Poisson + read noise curve. If scatter >> curve: crowding errors dominate |

### 6. Quality Flags to Apply

Stars to exclude or flag in final catalog:
- `CHI > some threshold` (e.g., > 2.0 or > 3.0 depending on field)
- `|SHARP| > some threshold` (e.g., > 0.5)
- `σ_mag > maximum acceptable error`
- Stars with `mag = 99.999` from any stage
- Stars with `N_iter = 50` (non-converged) should be inspected individually

---

## 17. Parameter Quick-Reference Table

| ID | Full name | Affected routines | Typical value | Range |
|---|---|---|---|---|
| `RE` | Read noise (ADU, 1 frame) | FIND (+ all via header) | Detector-specific | >0 |
| `GA` | Gain (photons/ADU, 1 frame) | FIND (+ all via header) | Detector-specific | >0 |
| `LO` | Low good datum (σ below sky) | FIND (+ all via header) | 7 | ≥0 |
| `HI` | High good datum (ADU absolute) | FIND (+ all via header) | Near saturation level | ≥0 |
| `FW` | FWHM of objects (pixels) | FIND | 2–5 typical | 0.2–15.0 |
| `TH` | Detection threshold (σ) | FIND | 4.0 | ≥0 |
| `LS` | Low sharpness cutoff | FIND | 0.20 | 0.0–1.0 |
| `HS` | High sharpness cutoff | FIND | 1.00 | 0.6–2.0 |
| `LR` | Low roundness cutoff | FIND | -1.00 | -2.0–0.0 |
| `HR` | High roundness cutoff | FIND | 1.00 | 0.0–2.0 |
| `WA` | Watch progress | Most routines | 1 (interactive), -2 (batch) | -2 to 2 |
| `FI` | Fitting radius (pixels) | PSF, PEAK, GROUP, NSTAR | ∼1× FWHM | 1.0–10.0 |
| `PS` | PSF radius (pixels) | PSF | ∼2–4× FWHM | 1.0–35.0 |
| `VA` | Variable PSF degree | PSF | 0 (start here) | -1 to 2 |
| `FR` | Fractional pixel expansion | — | 0 (not implemented) | — |
| `AN` | Analytic PSF model | PSF | 1 (Gaussian) | 1–6 |
| `EX` | Extra PSF cleaning passes | PSF | 5 | 0–9 |
| `PE` | Percent error (flat-field) | PEAK, NSTAR | 0.75 | 0–100 |
| `PR` | Profile error | PEAK, NSTAR | 5.0 | 0–100 |

---

## 18. Common Pitfalls and How to Avoid Them

### Read noise and gain left at default (0)

**Symptom:** Detection thresholds are meaningless; NSTAR rejection behaves erratically; CHI values are all wrong.

**Fix:** Always have `RE` and `GA` set correctly in `DAOPHOT.OPT` before starting. These are for **a single raw exposure**, not your combined frame.

---

### Threshold set too low

**Symptom:** Thousands of detections; NSTAR/ALLSTAR runs for hours; photometry scatter is catastrophic; spurious stars migrate and corrupt neighbors.

**Fix:** Use FIND's threshold-search feature (answer "N" repeatedly). Build N(detected) vs threshold and find the elbow. 4–5σ is usually the right region.

---

### PSF radius too small

**Symptom:** Subtracted bright stars have luminous rings with sharp inner edges on the residual image.

**Fix:** Increase `PS` to cover the entire visible halo of bright stars. Check on your image display before running PSF.

---

### Fitting radius too large on a crowded field

**Symptom:** Profile fits are contaminated by neighbors; CHI values are elevated; positions and magnitudes are systematically wrong.

**Fix:** Reduce `FI` slightly. The fitting radius controls how far out the profile is used — in a crowded field, including the outer wings drags in contaminated pixels.

---

### One PSF pass is not enough in a crowded field

**Symptom:** Subtracted PSF stars have systematic patterns in residuals; outer wings imperfectly subtracted; PSF look-up table has "ghost" features.

**Fix:** Iterate. Use VA=-1 first, subtract neighbors, make new PSF, increase VA progressively. Plan on 3–5 iterations.

---

### PHOTOMETRY aborts: "fewer than 20 sky pixels"

**Symptom:** PHOTOMETRY returns to Command: mode immediately.

**Fix:** (a) Increase outer sky radius; (b) relax LO (lower the sigma threshold); (c) for short-exposure frames with low sky, set LO to a negative number (e.g., LO=-5) so that legitimate sky pixels below the mean are not flagged as bad.

---

### APPEND of mismatched files

**Symptom:** Subsequent routines crash or produce garbage.

**Fix:** Only APPEND files of identical type. Double-check: `.AP` onto `.AP`, `.GRP` onto `.GRP`, never `.COO` onto `.AP`.

---

### Group size exceeds 60 for NSTAR

**Symptom:** NSTAR refuses to reduce some groups or crashes.

**Fix:** Use SELECT to extract the offending groups, run GROUP with larger critical overlap on them, feed subsets back to NSTAR. Or just switch to ALLSTAR, which has no group size limit.

---

### SUBSTAR leaves wrong stars in the image

**Symptom:** After SUBSTAR with "stars to leave in", the residual image has unexpected holes or leftover stars.

**Fix:** Star matching in SUBSTAR is done by ID number only. Ensure that when you APPENDed or sorted files, you used unique ID offsets (via OFFSET) and that the `.LST` file contains the correct IDs corresponding to the stars you actually want to keep.

---

### CHI systematically different from 1.0

**If CHI >> 1.0 for all stars:**
- RE and/or GA are underestimated → expected noise too low → observed scatter seems too high
- Or your PSF is a poor fit (wrong analytic model, too few PSF stars, unconverged PSF)

**If CHI << 1.0 for all stars:**
- RE and/or GA are overestimated → expected noise too high

**If CHI rises steeply with brightness:**
- PE and PR are too low; fine-scale flat-field errors or PSF mismatch dominate at high counts

---

### Artificial star tests showing large bias or scatter

**Symptom:** Recovered magnitudes differ systematically from input, or scatter much more than expected from random noise alone.

**Fix:** (a) Sky under the star is being overestimated in crowded regions — consider using ALLSTAR's sky re-determination option; (b) PSF is not well-constrained; (c) completeness limit is higher than thought — adjust your science goals accordingly.

---

*Guide compiled from: Stetson, P.B. — "User's Manual for DAOPHOT II: The Next Generation", Dominion Astrophysical Observatory, Herzberg Institute of Astrophysics, 2000 June 5.*
