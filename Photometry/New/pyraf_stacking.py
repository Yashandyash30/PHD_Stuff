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
# ## 1. Configuration

# %%
import glob
import os
import numpy as np
import astroalign as aa
from astropy.io import fits
from pyraf import iraf

# =====================================================================
# 1. USER CONFIGURATION
# =====================================================================

# ----------------- Align Configuration -----------------
input_wildcard_align = "ctfb_GRB261010A_P10_i_*.fit*"
output_prefix_align = "al_"

# ----------------- Imcombine Configuration -------------
input_wildcard_imcombine = f"{output_prefix_align}*.fit*"
output_file_imcombine = "GRB261010A_i.fit"

# Detector characteristics (REQUIRED for ccdclip to work!)
detector_gain = "1.0"      # e-/ADU
detector_rdnoise = "6.0"   # e-
exptime_keyword = "EXPTIME" # FITS header keyword for exposure time

# Combination settings
combine_option = "average"  # 'average', 'sum', or 'median'
reject_option = "ccdclip"   # 'ccdclip' (best), 'avsigclip', or 'none'


# %% [markdown]
# ## 2. Functions

# %%
def align_for_imcombine(file_pattern, output_prefix="aligned_", detector_gain=None, detector_rdnoise=None):
    """
    Aligns a set of FITS images to the first image in the sequence using astroalign.
    Saves the registered images to disk so they can be stacked by PyRAF's imcombine.
    Supports both .fit and .fits extensions and ensures GAIN/RDNOISE are in headers.
    """
    image_list = sorted(glob.glob(file_pattern))
    if not image_list and (file_pattern.endswith('.fit') or file_pattern.endswith('.fits')):
        alt_pattern = os.path.splitext(file_pattern)[0] + '.fit*'
        image_list = sorted(glob.glob(alt_pattern))
    
    if len(image_list) < 2:
        print("[Alignment] WARNING: Need at least 2 images to align. Exiting.")
        return

    ref_image_path = image_list[0]
    print(f"\n--- Starting Alignment Pipeline ---")
    print(f"Reference image selected: {ref_image_path}")

    # =====================================================================
    # 1. LOAD & PREP REFERENCE IMAGE
    # =====================================================================
    with fits.open(ref_image_path) as hdul:
        # Force 64-bit float to prevent byte-swapping errors
        ref_data = hdul[0].data.astype(np.float64) 
        ref_header = hdul[0].header.copy()
        
    # Clean any NaN values 
    ref_data[np.isnan(ref_data)] = np.nanmedian(ref_data)
    
    # Ensure detector characteristics are stored in header
    if detector_gain is not None:
        ref_header['GAIN'] = (float(detector_gain), 'Detector gain (e-/ADU)')
    if detector_rdnoise is not None:
        ref_header['RDNOISE'] = (float(detector_rdnoise), 'Detector read noise (e-)')

    # Crop the reference image for the calculation phase
    c_ref = ref_data[200:-200, 200:-200]

    # Save the reference image with the new prefix so it gets included in the stack
    ref_out_path = f"{output_prefix}{os.path.basename(ref_image_path)}"
    fits.writeto(ref_out_path, ref_data, ref_header, overwrite=True, output_verify='silentfix')
    print(f"[KEEP] {ref_image_path} -> {ref_out_path} (Master Reference)")

    successful_alignments = 1

    # =====================================================================
    # 2. ALIGN REMAINING IMAGES
    # =====================================================================
    for img_file in image_list[1:]:
        with fits.open(img_file) as hdul:
            source_data = hdul[0].data.astype(np.float64)
            source_header = hdul[0].header.copy()
            
        source_data[np.isnan(source_data)] = np.nanmedian(source_data)
        
        # Ensure detector characteristics are stored in header
        if detector_gain is not None:
            source_header['GAIN'] = (float(detector_gain), 'Detector gain (e-/ADU)')
        if detector_rdnoise is not None:
            source_header['RDNOISE'] = (float(detector_rdnoise), 'Detector read noise (e-)')

        try:
            # Trick 1: Central Cropping to bypass edge artifacts (amplifier glow, vignetting)
            c_source = source_data[200:-200, 200:-200]
            
            # Trick 2: Bulletproof Astroalign parameters for faint fields
            t, _ = aa.find_transform(c_source, c_ref, min_area=9, detection_sigma=2.0)
            
            # Apply the transform to the FULL source image (not the cropped one!)
            # This ensures your final aligned image is the original size.
            registered_data, _ = aa.apply_transform(t, source_data, ref_data)
            
            # Save the aligned image with its header
            out_path = f"{output_prefix}{os.path.basename(img_file)}"
            fits.writeto(out_path, registered_data, source_header, overwrite=True, output_verify='silentfix')
            
            successful_alignments += 1
            print(f"[KEEP] {img_file} -> {out_path} (Aligned successfully)")
            
        except Exception as e:
            # If clouds ruined the frame or it's too blurry to match stars, we drop it safely
            print(f"[REJECT] ⚠️ Could not align {img_file}: {e}")

    # =====================================================================
    # 3. SUMMARY
    # =====================================================================
    print(f"\n--- Alignment Complete ---")
    print(f"Successfully aligned {successful_alignments} out of {len(image_list)} images.")
    print(f"Ready for PyRAF: Run imcombine on '{output_prefix}*.fit*'")

def run_imcombine(input_wildcard, output_file, combine_option, reject_option, detector_gain, detector_rdnoise, exptime_keyword):
    file_list = sorted(glob.glob(input_wildcard))
    if not file_list and (input_wildcard.endswith('.fit') or input_wildcard.endswith('.fits')):
        alt_wildcard = os.path.splitext(input_wildcard)[0] + '.fit*'
        file_list = sorted(glob.glob(alt_wildcard))
        if file_list:
            input_wildcard = alt_wildcard

    N = len(file_list)

    print(f"Found {N} files matching '{input_wildcard}'.")
    if N == 0:
        print("Error: No files found. Check your directory or wildcard pattern.")
        return

    print(f"Stacking using '{combine_option}' with '{reject_option}' rejection...")

    # Write input list file for clean, robust IRAF ingestion
    list_file = "imcombine_input.list"
    with open(list_file, "w") as fh:
        for f in file_list:
            fh.write(f"{f}\n")

    # Remove existing output file if present
    if os.path.exists(output_file):
        os.remove(output_file)
        print(f"Removed existing output file: {output_file}")

    # Load the immatch package which contains imcombine
    iraf.immatch()

    # RUN IMCOMBINE
    try:
        iraf.imcombine(
            input=f"@{list_file}",
            output=output_file,
            headers="", 
            bpmasks="", 
            rejmasks="", 
            nrejmasks="", 
            expmasks="", 
            sigmas="",
            imcmb="$I", 
            logfile="STDOUT",
            
            combine=combine_option, 
            reject=reject_option, 
            project="no", 
            outtype="real",
            outlimits="", 
            offsets="none", 
            masktype="none", 
            maskvalue=0, 
            blank=0.0,
            
            scale="exposure", 
            zero="mode", 
            weight="none", 
            statsec="", 
            expname=exptime_keyword,
            
            lthreshold="INDEF", 
            hthreshold="INDEF", 
            nlow=1, 
            nhigh=1, 
            nkeep=1,
            
            mclip="yes", 
            lsigma=3.0, 
            hsigma=3.0,
            
            rdnoise=detector_rdnoise, 
            gain=detector_gain,
            snoise=0.0,
            
            sigscale=0.1, 
            pclip=-0.5, 
            grow=0.0, 
            mode="al"
        )
    finally:
        if os.path.exists(list_file):
            os.remove(list_file)
    
    print("\nStacking complete. Calculating new noise model and updating headers...")

    # CALCULATE NEW EFFECTIVE GAIN & READ NOISE NUMERICALLY
    gain_val = float(detector_gain)
    rn_val = float(detector_rdnoise)

    if combine_option == "average":
        eff_gain = N * gain_val
        eff_rn = rn_val / np.sqrt(N)
    elif combine_option == "sum":
        eff_gain = gain_val
        eff_rn = rn_val * np.sqrt(N)
    elif combine_option == "median":
        effective_N = (2.0 / 3.0) * N
        eff_gain = effective_N * gain_val
        eff_rn = rn_val / np.sqrt(effective_N)
    else:
        print(f"Unrecognized combine option '{combine_option}'! Using raw detector values.")
        eff_gain = gain_val
        eff_rn = rn_val

    eff_gain = round(float(eff_gain), 4)
    eff_rn = round(float(eff_rn), 4)

    # UPDATE THE OUTPUT FILE'S HEADER USING ASTROPY FITS
    with fits.open(output_file, mode='update') as hdul:
        hdul[0].header['GAIN'] = (eff_gain, 'Effective gain (e-/ADU)')
        hdul[0].header['RDNOISE'] = (eff_rn, 'Effective read noise (e-)')
        hdul.flush()

    # Also update via IRAF hedit with computed numeric literals
    try:
        iraf.hedit(images=output_file, fields="gain", value=eff_gain, add="yes", addonly="no", delete="no", verify="no", show="no", update="yes")
        iraf.hedit(images=output_file, fields="rdnoise", value=eff_rn, add="yes", addonly="no", delete="no", verify="no", show="no", update="yes")
    except Exception:
        pass

    print(f"Header updated in {output_file}: GAIN = {eff_gain} e-/ADU, RDNOISE = {eff_rn} e-")
    print("\nPipeline finished successfully! Your image is ready for DAOPHOT.")



# %% [markdown]
# ## 3. Execute align

# %%
align_for_imcombine(
    file_pattern=input_wildcard_align,
    output_prefix=output_prefix_align,
    detector_gain=detector_gain,
    detector_rdnoise=detector_rdnoise
)

# %% [markdown]
# ## 4. Execute imcombine

# %%
run_imcombine(
    input_wildcard=input_wildcard_imcombine,
    output_file=output_file_imcombine,
    combine_option=combine_option,
    reject_option=reject_option,
    detector_gain=detector_gain,
    detector_rdnoise=detector_rdnoise,
    exptime_keyword=exptime_keyword
)
