import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import bayesian_blocks

# 1. Load the GBM TTE Event Data
# Replace with your actual file path for GRB250919A
file_path = '/home/void/GRB_Data/GRB250919A/current/glg_tte_n6_bn250919020_v00.fit' 
with fits.open(file_path) as hdul:
    data = hdul['EVENTS'].data
    arrival_times = data['TIME']
    trigtime = hdul['PRIMARY'].header.get('TRIGTIME', 0)

# Zero the times relative to the trigger
time = arrival_times - trigtime

# 2. Filter the time range (e.g., -10 to +100 seconds)
mask = (time > -10) & (time < 100)
filtered_time = time[mask]

# ---------------------------------------------------------
# THE FIX: PRE-BIN THE DATA
# ---------------------------------------------------------
dt = 0.05  # Set fine bin width (e.g., 50 milliseconds)
fine_bins = np.arange(-10, 100 + dt, dt)
counts, _ = np.histogram(filtered_time, bins=fine_bins)

# 3. Run Bayesian Blocks on the PRE-BINNED data
# Note: For fitness='events', we pass bin CENTERS as t (from edges),
# and 'counts' as x (integer counts per bin)
bin_centers_fine = (fine_bins[:-1] + fine_bins[1:]) / 2
edges = bayesian_blocks(t=bin_centers_fine, x=counts, fitness='events', p0=0.01)

# 4. Create the rate arrays for plotting
bin_widths = edges[1:] - edges[:-1]
block_counts, _ = np.histogram(filtered_time, bins=edges)
rates = block_counts / bin_widths
bin_centers = (edges[1:] + edges[:-1]) / 2

# 5. Plot the result
plt.figure(figsize=(10, 5))

# Plot the background fine-binned lightcurve for reference (convert counts to rate)
plt.step(fine_bins[:-1], counts/dt, where='post', color='grey', alpha=0.5, label=f'{dt}s Bins')

# Plot the Bayesian Blocks
plt.step(bin_centers, rates, where='mid', color='red', linewidth=2, label='Bayesian Blocks')

plt.xlabel('Time since trigger (s)')
plt.ylabel('Rate (counts/s)')
plt.xlim(-10, 100)
plt.legend()
plt.tight_layout()
plt.show()

print(f"Success: Reduced from {len(filtered_time)} individual events down to {len(edges)-1} Bayesian Blocks.")

# Print the intervals for XSPEC/gtburst extraction
for i in range(len(edges)-1):
    print(f"Bin {i}: T_start = {edges[i]:.3f}, T_stop = {edges[i+1]:.3f}")