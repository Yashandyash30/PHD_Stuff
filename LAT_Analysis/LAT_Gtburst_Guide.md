Viewed Fermi_LAT_objectives.txt:1-97
Viewed LAT_Analyisi_guide_raw.txt:1-590

Here is a comprehensive and detailed guide for performing Fermi-LAT analysis on GRB 250919A, based on the specific objectives and GUI guidelines provided in your documents.

### 1. Likelihood Analysis GUI Parameters

When working within the `gtburst` GUI, selecting the correct parameters is crucial for both detecting the burst and generating accurate spectral files.

**Geometry Parameters**

* **rad (Radius of Interest):** Set to **$12^\circ$** for the Prompt phase (using Transient class) to capture spreading photons. Set to **$10^\circ$** for Extended/Afterglow phase (using Source class).
* **zmax (Zenith Angle Cut):** Set to **$100^\circ$**. This is the safety limit to avoid Earth-limb atmospheric gamma-ray contamination.

**Physics Parameters**

* **irf (Event Class):**
  * **`P8R3_TRANSIENT020_V3`**: Use for the **Prompt Phase** (short timescales $< 100$s). It accepts more photons to maximize statistics during the high-flux interval. *(Note: Avoid the "E" extended version unless you have manually symlinked the missing isotropic templates in your Conda environment).*
  * **`P8R3_SOURCE`**: Use for **Extended Emission/Afterglow** (e.g., $100$s to $1000$s). It is stricter and filters out more background during the lower-flux, late-time emission.
* **Strategy:** Change from `time` to **`events`**. The `time` strategy discards the *entire* time interval if even a tiny fraction touches the Earth limit. The `events` strategy only removes specific bad photons (acting like a Swiss cheese filter), which prevents the fatal "Zero Exposure" error.

**Likelihood Model Parameters (for Localization / Final Fits)**

* **particle_model:** `isotr template` (Models uniform background).
* **galactic_model:** `template` (Best for bright bursts, allows the software to calculate background scaling rather than forcing a fixed value).
* **source_model:** `powerlaw2` (Fits a simple Power Law $N(E) \propto E^{-\gamma}$, standard for GRBs).
* **fgl_mode:** `complete` (Includes known steady sources from the Fermi catalog in the background fit).
* **tsmin:** `20` (Minimum Test Statistic required for detection).
* **optimizeposition:** `yes` (Fine-tunes RA/Dec using the physical model).
* **spectralfiles:** `yes` (Critical: Generates the `.pha` and `.rsp` files you need for XSPEC).
* **liketype:** `unbinned` (Superior for low-count GRB data compared to binned analysis).

---

### 2. High-Energy Localization in Gtburst GUI

To refine the burst coordinates using the LAT's superior angular resolution compared to the GBM:

1. **Extract the Prompt Phase:** In the `gtburst` GUI, select the primary emission episode (e.g., $T_0$ to $T_0 + 50$s).
2. **Make Counts Map:** Run the tool with `rad=12`, `zmax=100`, `irf=p8_transient020`, and `strategy=events`. Validate the map to ensure you see a cluster of photons around the center cross.
3. **Run Localization Tool:** Go to `Tasks -> Localize` (this runs `gtfindsrc` in the background).
4. **Set Model Parameters:** Use `particle_model=isotr template`, `galactic_model=template`, `source_model=powerlaw2`, and `fgl_mode=complete`. Click **Run**.
5. **Accept and Update:** Once the likelihood maximization finishes, it will give you a new precise RA/Dec (e.g., RA 298.767, Dec -48.970). **Accept this new position** and manually update the RA and Dec fields in the main GUI window so your Region of Interest (ROI) is perfectly re-centered for all subsequent analyses.

---

### 3. Prompt Emission Analysis (Time-Resolved via Bayesian Blocks)

You want to run Time-Resolved Spectroscopy using a Bayesian Blocks algorithm to segment the data objectively.

**How to make the filtered data in Gtburst GUI:**
The analysis strictly calls for using data *before* `gtmktime` deletes entire Earth-limb time intervals.

* **In the GUI:** By setting your **`Strategy`** parameter to **`events`** instead of `time`, you force `gtburst` to behave this way. It skips the `gtmktime` step that deletes whole time blocks, and instead, directly cuts out bad individual photons using your Zenith $\le 100^\circ$ cut. The underlying data structure this uses is saved as the `_filt.fit` file in your working directory.
* **Applying Bayesian Blocks:** You cannot run the Bayesian block algorithm directly inside the standard `gtburst` GUI interface for LAT data. Instead, you run the provided Python script (from your guide) on the `gll_ft1_tr_bn..._filt.fit` file. The script parses the relative times of photons within $T_0 - 100$s to $+500$s, applying `zenith_angles <= 100`, and calculates the exact $T_{start}$ and $T_{stop}$ bins (using $p_0=0.05$).
* **Returning to GUI:** Once the script prints out the Bayesian Block time intervals (e.g., Bin 1: 0.0s to 4.5s), you manually drag the time-selection bars in the `gtburst` GUI to match those exact intervals, and generate `.pha` files for each slice with `spectralfiles=yes`.

**Python Code for Time Resolved Analysis Using Bayesian Block Binning:**

```python
import numpy as np
from astropy.io import fits
from astropy.stats import bayesian_blocks
import sys

# 1. Define your file and trigger time
# CRITICAL FIX: We are using the _filt.fit file. This file contains the data
# BEFORE gtmktime deleted the Earth-limb time intervals.
filename = 'gll_ft1_tr_bn250919020_v00_filt.fit'

# Trigger time from your gtburst GUI
trigger_time = 779934537.284

print(f"Loading data from {filename}...")

# 2. Open the FITS file and extract the data
try:
    with fits.open(filename) as hdul:
        # LAT event data is stored in the 1st extension ('EVENTS')
        events = hdul[1].data

        # Extract the absolute arrival times
        absolute_times = events['TIME']

        # NEW: Extract the Zenith angle of every individual photon
        zenith_angles = events['ZENITH_ANGLE']
except FileNotFoundError:
    print(f"\n⚠️ ERROR: '{filename}' not found.")
    print("Make sure this script is saved inside your 'Filtered_event_file' folder!")
    sys.exit(1)

# 3. Convert to time relative to the trigger (T - T0)
relative_times = absolute_times - trigger_time

# 4. Filter the data (The Python equivalent of the 'events' strategy)
# - Time window: -100 to 500 seconds (Captures the whole prompt and early tail)
# - Zenith cut: <= 100 degrees (Throws out individual Earth-limb photons, NOT the time)
prompt_mask = (relative_times > -100) & (relative_times < 500) & (zenith_angles <= 100)
prompt_times = relative_times[prompt_mask]

print(f"Found {len(prompt_times)} photons in the clean -100s to 500s window.")

# Safety check
if len(prompt_times) == 0:
    print("\n⚠️ ERROR: No photons found! Check your trigger time or filename.")
    sys.exit(1)

print("Calculating Bayesian Blocks...\n")

# 5. Run the Bayesian Blocks algorithm
# p0 is the false-alarm probability. 0.05 is standard.
edges = bayesian_blocks(prompt_times, fitness='events', p0=0.05)

# 6. Print the intervals formatted for your report and gtburst
print("--- BAYESIAN BLOCK TIME INTERVALS ---")
for i in range(len(edges) - 1):
    t_start = edges[i]
    t_stop = edges[i+1]

    # Filter the printout to only show the relevant prompt/early tail bins
    if t_stop > -10 and t_start < 200:
        print(f"Bin {i+1}: {t_start:>7.2f}  to  {t_stop:>7.2f} seconds")
print("-------------------------------------")
```

---

### 4. Temporally Extended > 100 MeV Emission

After the prompt phase, you must analyze the long-lasting late-time emission.

**Which Intervals to Select:**

* You should select **logarithmically spaced time bins** (e.g., $50$-$100$s, $100$-$300$s, $300$-$1000$s) or isolate secondary pulses as seen in the LAT light curve.
* *Integrated extended interval:* You can also select the entire extended range (e.g., $T_0 + 50$s to $T_0 + 1000$s) to calculate total late-time fluence.

**How to Proceed in Gtburst GUI for Likelihood Analysis:**

1. **Update Coordinates:** Ensure your RA/Dec are set to the optimized LAT localization coordinates found in step 2.
2. **Adjust Parameters:** Because the flux is lower, you must mitigate background contamination. Change **`irf`** to `Source` (e.g., `P8R3_SOURCE`) and shrink **`rad`** to $10^\circ$. Keep **`strategy=events`**.
3. **Run Likelihood Analysis:**
   * Select your extended time interval on the light curve.
   * Set `spectralfiles=yes`, `liketype=unbinned`.
   * Click Run. Note the Test Statistic (TS), Photon Index, and Flux for the interval.
4. **Save/Rename Files:** Immediately rename the generated `.pha` and `.rsp` files (e.g., `LAT_100_300s.pha`) so they are not overwritten when you analyze the next logarithmic bin.

---

### 5. High-Energy Photon Association and Jet Kinematics

To understand jet kinematics, you need to extract the highest-energy photons temporally and spatially coincident with the GRB and plot them.

**How to extract high-energy photons & Plot Energy vs. Time:**
You can achieve both by extending the Python script provided in your guide. Save the following code as a python file (e.g. `plot_photons.py`) in the same directory as your `_filt.fit` file and run it. It filters for $100^\circ$ Zenith and $0.5^\circ$ spatial separation, prints the highest energy photons, and then generates the Energy vs. Time scatter plot.

```python
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
import matplotlib.pyplot as plt

# 1. The File and Trigger Time
filename = 'gll_ft1_tr_bn250919020_v00_filt.fit'
trigger_time = 779934537.284 

# Your optimized GRB coordinates
grb_ra = 298.719
grb_dec = -49.133
grb_coord = SkyCoord(ra=grb_ra*u.degree, dec=grb_dec*u.degree, frame='icrs')

print(f"Scanning {filename} for extreme photons...\n")

try:
    with fits.open(filename) as hdul:
        events = hdul[1].data
        times = events['TIME'] - trigger_time
        energies_mev = events['ENERGY']
        ras = events['RA']
        decs = events['DEC']
        zeniths = events['ZENITH_ANGLE']
except FileNotFoundError:
    print(f"⚠️ ERROR: '{filename}' not found.")
    sys.exit(1)

# 2. Filter the data
# We want photons within the time window, avoiding the Earth limb, 
# and ONLY those within a 0.5 degree radius of the exact GRB center.
photon_coords = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')
separations = grb_coord.separation(photon_coords).degree

# Mask: Valid time, valid zenith, and tightly clustered to the GRB
mask = (times > -10) & (times < 1000) & (zeniths <= 100) & (separations <= 0.5)

filtered_times = times[mask]
filtered_energies_gev = energies_mev[mask] / 1000.0  # Convert MeV to GeV
filtered_separations = separations[mask]

if len(filtered_times) == 0:
    print("No photons found matching the criteria.")
    sys.exit(1)

# 3. Sort by Energy (Highest first)
sort_indices = np.argsort(filtered_energies_gev)[::-1]

top_times = filtered_times[sort_indices]
top_energies = filtered_energies_gev[sort_indices]
top_separations = filtered_separations[sort_indices]

# 4. Print the Top 5 Record Breakers
print("--- TOP 5 HIGHEST ENERGY PHOTONS ---")
print(f"{'Rank':<6} | {'Time (s)':<10} | {'Energy (GeV)':<14} | {'Distance (deg)':<14}")
print("-" * 52)

for i in range(min(5, len(top_times))):
    print(f"#{i+1:<5} | {top_times[i]:<10.2f} | {top_energies[i]:<14.2f} | {top_separations[i]:<14.3f}")
print("-" * 52)

# 5. Plotting Energy vs Time
plt.figure(figsize=(10, 6))

# Use a scatter plot. High energy photons are larger and colored differently.
plt.scatter(filtered_times, filtered_energies_gev, 
            c=filtered_energies_gev, cmap='viridis', 
            s=50, edgecolors='black', alpha=0.8)

# Formatting the plot
plt.yscale('log') # Log scale is standard for energy ranges covering MeV to GeV
plt.xscale('symlog', linthresh=10) # Symmetric log handles negative pre-trigger times well
plt.colorbar(label='Energy (GeV)')

plt.axvline(x=0, color='red', linestyle='--', alpha=0.6, label='Trigger Time ($T_0$)')

plt.title('LAT High-Energy Photons vs Time (GRB 250919A)')
plt.xlabel('Time since trigger (s)')
plt.ylabel('Photon Energy (GeV)')
plt.grid(True, which="both", ls="--", alpha=0.3)
plt.legend()

# Save and Show
plt.savefig('High_Energy_Photons_vs_Time.png', dpi=300, bbox_inches='tight')
print("\nPlot saved as 'High_Energy_Photons_vs_Time.png'")
plt.show()
```
