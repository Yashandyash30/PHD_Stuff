Fire up your Arch terminal. DAOPHOT II is an old-school beast, but for isolating a fading GRB optical counterpart or a supernova from a messy host galaxy, its iterative profile-fitting is exactly what you need.

Based on your terminal screenshot, your parameters are successfully loaded (using that Penny function with 5 free parameters for the PSF is a solid choice for complex profiles). Here is the exact, step-by-step CLI workflow to get you through the actual photometry and point-spread function (PSF) generation.

Whenever you see `<CR>`, it means simply hit **Enter** to accept the default file name the program offers. This is highly recommended to keep your file extensions cleanly synced.

### Phase 1: Detection & Aperture Photometry

This executes the initial finding and basic photometry based on the `daophot.opt` and `photo.opt` parameters you've already configured.

1. Attach your FITS file 

```text
Command: AT
Enter file name: your_image.fits

```

2. Find the sources 

```text
Command: FI
Number of frames averaged, summed: 1,1  
File for the positions (default ?.COO): <CR>

```

(Note: Enter `1,1` if this is a single exposure. The terminal will spit out the detection threshold. If it asks "Are you happy with this?", type `Y`.)

3. Run Aperture Photometry 

```text
Command: PH
Enter table name (default PHOTO.OPT): <CR>
File with the positions (default ?.COO): <CR>
File for the magnitudes (default ?.AP): <CR>

```

*This calculates rough sky backgrounds and magnitudes for everything it just found, saving it to a `.ap` file.*

---

### Phase 2: Building the Point-Spread Function (PSF)

This is the critical step. You need isolated, bright (but unsaturated) stars to build a mathematical model of how light spreads in your specific image.

4. Pick Candidate PSF Stars 

```text
Command: PI
Input file name (default ?.AP): <CR>
Desired number of stars, faintest magnitude: 25, 99
Output file name (default ?.LST): <CR>

```

(DAOPHOT will find 25 good candidates and reject ones too close to the edge or too close to brighter stars. Entering `99` just tells it to use whatever magnitude is necessary to get your 25 stars.)

5. Generate the PSF Model 

```text
Command: PS
File with aperture results (default ?.AP): <CR>
File with PSF stars (default ?.LST): <CR>
File for the PSF (default ?.PSF): <CR>

```

Because your `WATCH PROGRESS` is set to `1.00`, the terminal will now draw ASCII-art representations of each candidate star.

* For each star, it will ask if you want to include it.
* Look at the ASCII plot: if it looks like a clean, single star, type `Y`. If it looks like a cosmic ray, a doublet, or a smudge, type `N`.



---

### Phase 3: Profile-Fitting Photometry (ALLSTAR)

While DAOPHOT has internal commands (`NSTAR`, `PEAK`), the manual strongly recommends using the standalone program **ALLSTAR** for the final rigorous simultaneous least-squares fit.

6. Exit DAOPHOT 

```text
Command: EX

```

**7. Run ALLSTAR**
Back in your standard Linux shell, launch ALLSTAR (it will read your `.opt` files just like DAOPHOT did):

```bash
allstar

```

```text
Input image name: your_image.fits
File with the PSF (default ?.PSF): <CR>
Input file (default ?.AP): <CR>
File for results (default ?.ALS): <CR>
Name for subtracted image (default ?s.fits): <CR>

```

### The Output

ALLSTAR will crunch the matrix and give you two highly useful files:

1. 
**The `.als` file:** This contains your final, precise profile-fitted magnitudes and coordinates.


2. 
**The `?s.fits` (Subtracted Image):** This is your original image with all the known stars digitally erased using the PSF model.



**Pro-Tip for Crowded Fields:** Open that subtracted FITS file in a viewer (like ds9). If you see faint stars left behind that the initial `FIND` missed because they were hiding in the glare of brighter objects, you can run DAOPHOT again, `ATTACH` the subtracted image, run `FIND` to grab their coordinates, and feed them back through the loop.
