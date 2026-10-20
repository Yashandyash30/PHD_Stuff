# ---
# jupyter:
#   jupytext:
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
#================================================================
# Cosmic Ray Cleaning Script for GRB Frames
# This script uses the Astro-SCRAPPY package to identify and remove cosmic rays from your GRB frames. 
# Make sure to adjust the parameters (sigclip, readnoise, gain) according to your specific data and camera characteristics for optimal results.
#================================================================

import glob
from astropy.io import fits
import astroscrappy
import numpy as np

flat_sci_files = glob.glob('tfb_*.fit')

for file in flat_sci_files:
    out_name = f"c{file}" # 'c' for cosmic-ray cleaned (e.g., ctfb_image.fits)
    print(f"Cleaning cosmic rays: {file} -> {out_name}")
    
    # Open the FITS file
    with fits.open(file) as hdul:
        data = hdul[0].data.astype(np.float64) # Force data into a clean, standard 64-bit float array
        header = hdul[0].header
        
        # Run Astro-SCRAPPY (L.A.Cosmic)
        # Note: readnoise and gain should be set to your camera's actual values!
        crmask, clean_data = astroscrappy.detect_cosmics(
            data, 
            sigclip=4.5,     # Adjust depending on how aggressive you want the clipping
            readnoise=5.0,   # Replace with your CCD's Read Noise (e.g., from header 'RDNOISE')
            gain=10.0         # Replace with your CCD's Gain (e.g., from header 'GAIN')
        )
        
        # Save the cleaned data to a new file
        fits.writeto(out_name, clean_data, header, overwrite=True)

        # Save the mask as a separate FITS file
        mask_out_name = f"mask_{file}"
        fits.writeto(mask_out_name, crmask.astype(np.uint8), header, overwrite=True)

print("Cosmic ray cleaning complete! Your frames are fully reduced.")
