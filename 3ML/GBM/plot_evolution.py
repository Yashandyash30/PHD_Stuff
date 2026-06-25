# %% [markdown]
# # Evolution Plot
# This script generates a 7-panel plot showing the evolution of spectral parameters.
# It reads from `fit_params.csv` and the raw data to extract the background NaI lightcurve.

# %%
import os
import astropy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.ticker import AutoMinorLocator
from astropy.table import Table
import pandas as pd
from astropy.io import fits

# Set styling
font = {'family' : 'serif', 'weight' : 'normal', 'size' : 16}
plt.style.use('classic')
matplotlib.rc('font', family='serif', size='22', weight='normal')
matplotlib.rc('text', usetex=False) # Keep False if LaTeX is not installed, set to True if it is.

# Configure models paths
base_dir = 'gbm_pipeline_output'

def load_csv(model_name):
    path = os.path.join(base_dir, model_name, 'fit_params.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    else:
        print(f"Warning: {path} not found.")
        return None

# Load all model data
band = load_csv('band')
bandbb = load_csv('bandbb')
cpl = load_csv('cutoffpl')
cplbb = load_csv('cutoffbb')
syn = load_csv('synchrotron')

# We'll use band's time axis as the standard (they should all be identical)
if band is not None:
    t = (band['t_start'] + band['t_stop']) / 2.0
else:
    # fallback
    t = np.arange(31)

# Extract raw LC from NaI data
# Assuming 'glg_tte_n6_bn200412381_v00.fit.gz' is in current directory
tte_file = 'glg_tte_n6_bn200412381_v00.fit.gz'
if os.path.exists(tte_file):
    try:
        events = fits.getdata(tte_file, 'EVENTS')
        event_times = events['TIME']
        t_trigger = fits.getheader(tte_file, 'PRIMARY').get('TRIGTIME', event_times[0])
        event_times_rel = event_times - t_trigger
        
        # Bin the lightcurve (e.g., from -10 to +40 s with 0.5s resolution)
        lc_bins = np.arange(-10, 40, 0.5)
        lcr, lct_edges = np.histogram(event_times_rel, bins=lc_bins)
        lct = lct_edges[:-1]
    except Exception as e:
        print("Error loading raw lightcurve:", e)
        lct, lcr = np.array([]), np.array([])
else:
    lct, lcr = np.array([]), np.array([])

# %%
# Setup Figure
fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, sharex=True, figsize=(9, 15))
fig.subplots_adjust(wspace=0.3, hspace=0)

# PANEL 1: alpha (Band and BandBB)
if band is not None and 'B_alpha' in band.columns:
    ax1.errorbar(t, band['B_alpha'], yerr=(band['B_alpha_n_err'].abs(), band['B_alpha_p_err']), 
                 color="#1f77b4", marker="s", mfc='#1f77b4', mec='#1f77b4', ms=12, ls="--", label='Band')
if bandbb is not None and 'B_alpha' in bandbb.columns:
    ax1.errorbar(t, bandbb['B_alpha'], yerr=(bandbb['B_alpha_n_err'].abs(), bandbb['B_alpha_p_err']), 
                 color="#6baed6", marker="o", mfc='#6baed6', mec='k', ms=11, ls="--", label='Band+BB')

# PANEL 2: E_p (Band and BandBB)
if band is not None and 'E_peak' in band.columns:
    ax2.errorbar(t, band['E_peak'], yerr=(band['E_peak_n_err'].abs(), band['E_peak_p_err']), 
                 color="#1f77b4", marker="s", mfc='#1f77b4', mec='#1f77b4', ms=12, ls="--")
if bandbb is not None and 'E_peak' in bandbb.columns:
    ax2.errorbar(t, bandbb['E_peak'], yerr=(bandbb['E_peak_n_err'].abs(), bandbb['E_peak_p_err']), 
                 color="#6baed6", marker="o", mfc='#6baed6', mec='k', ms=11, ls="--")

# PANEL 3: beta (Band and BandBB)
if band is not None and 'beta' in band.columns:
    ax3.errorbar(t, band['beta'], yerr=(band['beta_n_err'].abs(), band['beta_p_err']), 
                 color="#1f77b4", marker="s", mfc='#1f77b4', mec='#1f77b4', ms=12, ls="--")
if bandbb is not None and 'beta' in bandbb.columns:
    ax3.errorbar(t, bandbb['beta'], yerr=(bandbb['beta_n_err'].abs(), bandbb['beta_p_err']), 
                 color="#6baed6", marker="o", mfc='#6baed6', mec='k', ms=11, ls="--")

# PANEL 4: kT (BandBB)
if bandbb is not None and 'kT' in bandbb.columns:
    ax4.errorbar(t, bandbb['kT'], yerr=(bandbb['kT_n_err'].abs(), bandbb['kT_p_err']), 
                 color="#6baed6", marker="o", mfc='#6baed6', mec='k', ms=11, ls="--")

# PANEL 5: B (Synchrotron)
if syn is not None and 'B' in syn.columns:
    ax5.errorbar(t, syn['B'], yerr=(syn['B_n_err'].abs(), syn['B_p_err']), 
                 color="#9ecae1", marker="s", mfc='#9ecae1', mec='k', ms=12, ls="--", label='Synchrotron')

# PANEL 6: p (Synchrotron index)
if syn is not None and 'index' in syn.columns:
    ax6.errorbar(t, syn['index'], yerr=(syn['index_n_err'].abs(), syn['index_p_err']), 
                 color="#9ecae1", marker="s", mfc='#9ecae1', mec='k', ms=12, ls="--")

# Configuration placeholders for other models:
# You can uncomment and modify these if you want to plot CutoffPL or CutoffBB parameters:
# if cpl is not None:
#     ax1.errorbar(t, cpl['B_alpha'], yerr=(cpl['B_alpha_n_err'].abs(), cpl['B_alpha_p_err']), color="green", marker="^", ls="--", label='CPL')
#     ax2.errorbar(t, cpl['E_c'], yerr=(cpl['E_c_n_err'].abs(), cpl['E_c_p_err']), color="green", marker="^", ls="--")
# if cplbb is not None:
#     ax1.errorbar(t, cplbb['B_alpha'], yerr=(cplbb['B_alpha_n_err'].abs(), cplbb['B_alpha_p_err']), color="lightgreen", marker="^", ls="--", label='CPL+BB')
#     ax2.errorbar(t, cplbb['E_c'], yerr=(cplbb['E_c_n_err'].abs(), cplbb['E_c_p_err']), color="lightgreen", marker="^", ls="--")
#     ax4.errorbar(t, cplbb['kT'], yerr=(cplbb['kT_n_err'].abs(), cplbb['kT_p_err']), color="lightgreen", marker="^", ls="--")

# Set Limits (Adjust based on your data)
ax1.set_xlim(-6, 33)
ax1.set_ylim(-1.7, -0.001)
ax2.set_ylim(11, 699)
ax3.set_ylim(-3.9, -1.6)
ax4.set_ylim(2, 29)
ax5.set_ylim(-9, 129)
ax6.set_ylim(2.9, 5.9)

# Formatting axes
for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7]:
    ax.tick_params(axis="both", direction='in', which='major', length=7, width=2, top=True, right=False, colors='k', grid_alpha=0.1, labelsize=20)
    ax.tick_params(axis="both", direction='in', which='minor', length=4, width=2, top=True, right=False, colors='k', grid_alpha=0.1, labelsize=20)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

ax1.legend(loc='upper right', fontsize=15, ncol=2, frameon=True, numpoints=1)
ax5.legend(loc='upper right', fontsize=15, ncol=2, frameon=True, numpoints=1)

ax1.set_ylabel(r"$\alpha$", size=20, color='black')
ax2.set_ylabel(r"$E_{p}$ (keV)", size=20, color='black')
ax3.set_ylabel(r"$\beta$", size=20, color='black')
ax4.set_ylabel(r"$kT$ (keV)", size=20, color='black')
ax5.set_ylabel(r"$B$ (G)", size=20, color='black')
ax6.set_ylabel(r"$p$", size=20, color='black')
ax7.set_xlabel(r"Time$-T_{0}$ (s)", size=20, color='black')

# PANEL 7: DIC differences
if band is not None and bandbb is not None and 'DIC' in band.columns and 'DIC' in bandbb.columns:
    ax7.errorbar(t, bandbb['DIC'] - band['DIC'], yerr=None, marker="s", color="tab:blue", mfc='tab:blue', mec='k', ms=5, ls="-", lw=2, label=r'DIC$_{Band+BB}-$DIC$_{Band}$')

if cpl is not None and cplbb is not None and 'DIC' in cpl.columns and 'DIC' in cplbb.columns:
    ax7.errorbar(t, cplbb['DIC'] - cpl['DIC'], yerr=None, marker="s", color="tab:cyan", mfc='tab:cyan', mec='k', ms=5, ls="-", lw=2, label=r'DIC$_{CPL+BB}-$DIC$_{CPL}$')

ax7.tick_params(axis="both", direction='in', which='major', length=7, width=2, top=True, right=True, colors='k', grid_alpha=0.1, labelsize=20)
ax7.tick_params(axis="both", direction='in', which='minor', length=4, width=2, top=True, right=True, colors='k', grid_alpha=0.1, labelsize=20)
ax7.legend(loc='lower left', fontsize=15, ncol=1, frameon=True, numpoints=1)
ax7.set_ylabel(r"$\Delta$ DIC", size=20, color='black')
ax7.axhline(-10, color='k', ls='--', lw=2)
ax7.set_xlim(-6, 33)
ax7.set_ylim(-159, 19)

# Background Twin Axes for lightcurve
if len(lct) > 0:
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7]:
        ax_twin = ax.twinx()
        ax_twin.step(lct, lcr, color="0.8", ls="-", zorder=1, alpha=0.3, where='post')
        ax_twin.fill_between(lct, 0, lcr, step='post', color='0.8', zorder=1, alpha=0.3)
        ax_twin.tick_params(direction='in', which='both', length=0, width=0, left=False, right=False, labelright=False)
        ax_twin.set_ylim(1, 1e4)

plt.savefig("lightcurve_evolution.pdf", format="pdf", bbox_inches='tight')
print("Saved plot to lightcurve_evolution.pdf")
