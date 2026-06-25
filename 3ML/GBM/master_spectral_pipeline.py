# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python (threeML)
#     language: python
#     name: threeml
# ---

# %% [markdown]
# # Fermi-GBM Time-Resolved Spectral Analysis Pipeline
# Each cell corresponds to one pipeline stage — run them top to bottom, or re-run any single stage independently.

# %% [markdown]
# ## Cell 1 — Imports & Setup

# %%
import os
import logging
import numpy as np
import matplotlib
matplotlib.use('Agg')   # swap to 'inline' if running interactively
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u

from threeML import (
    silence_logs, download_GBM_trigger_data, TimeSeriesBuilder,
    OGIPLike, DataList, Model, PointSource, BayesianAnalysis, plot_spectra,
)
from threeML.catalogs import FermiGBMBurstCatalog as FermiGBMBurst
from astromodels import Cutoff_powerlaw, Band, Log_uniform_prior, Truncated_gaussian, Blackbody           # add model name here after adding cell for new model
import pynchrotron
silence_logs()
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%H:%M:%S',
)
log = logging.getLogger('GBM_pipeline')
print('Imports OK')


# %%
from astropy.coordinates import SkyCoord
import astropy.units as u

# Your input string (note: replacing the special degree symbol with 'd' or 'deg' is sometimes 
# safer, but astropy can usually handle standard unicode degree symbols)
ra_str = "18h33m15.21s"
dec_str = "+62d31m57.1s" # Swapped symbols to standard d, m, s for safety

# Create the SkyCoord object
c = SkyCoord(ra=ra_str, dec=dec_str, frame='icrs')

# Extract the decimal degrees
ra_deg = c.ra.degree
dec_deg = c.dec.degree

print("--- Astropy Conversion ---")
print(f"RA in Degrees:  {ra_deg:.6f}")
print(f"DEC in Degrees: {dec_deg:.6f}")

# %% [markdown]
# ## Cell 2 — Configuration

# %%
# ── Edit this block only ────────────────────────────────────────────
CONFIG = dict(
    grb_name       = 'GRB200412381',  # catalog key      (with 'GRB' prefix)
    # Note: trigger_name, ra, and dec will be automatically fetched from the catalog
    # You can optionally define them here to override the catalog defaults!
    #ra=             # in degrees
    #dec=            # in degrees

    # detector selection
    max_rov        = 60.0,            # keep detectors with ROV < this [deg]

    # Analysis Type Selection
    analysis_type  = 'time_integrated',          # 'time_resolved', 'time_integrated', or 'both'

    # Time-Integrated settings
    ti_t_start     = 3.1,             # start of time-integrated window [s relative to trigger]
    ti_t_stop      = 21.5,            # end of time-integrated window

    # Bayesian-blocks settings (Time-resolved)
    bb_t_start     = -10.0,           # start of BB window [s relative to trigger]
    bb_t_stop      =  50.0,           # end   of BB window
    bb_p0          =  0.01,           # false-alarm probability for BB
    min_bin_width  =  0.05,           # seconds; narrower bins are merged away

    # spectral energy ranges
    nai_range      = '8-33, 37-900',  # keV  (gap skips K-edge ~34-36 keV)
    bgo_range      = '200-40000',     # keV

    # ── Add new model priors here ────────────────────────────────────────────────────────────
    
    # CutoffPL-function priors
    priors_cutoffpl = dict(
        K     = dict(lower=1e-5, upper=1e1),
        index = dict(lower=-2.0, upper=1.0,  mu=-1.3,  sigma=0.5),
        xc    = dict(lower=10.0, upper=1000, mu=300.0, sigma=200.0)
    ),

    # Band-function priors
    priors_band = dict(
        K     = dict(lower=1e-4, upper=1e1),
        alpha = dict(lower=-1.7, upper=1.0,  mu=-1.0,  sigma=0.5),
        xp    = dict(lower=10.0, upper=1000, mu=300.0, sigma=200.0),
        beta  = dict(lower=-4.0, upper=-1.8, mu=-2.5,  sigma=1.0),
    ),

    # Band+BB-function priors
    priors_bandbb = dict(
        alpha_1 = dict(lower=-1.7, upper=1.0,  mu=-1.0,  sigma=0.5),
        K_1     = dict(lower=1e-5, upper=1e1),
        xp_1    = dict(lower=10.0, upper=1000, mu=300.0, sigma=200.0),
        beta_1  = dict(lower=-4.0, upper=-1.8, mu=-2.5,  sigma=1.0),
        K_2     = dict(lower=1e-10, upper=1e-1),
        kT_2    = dict(lower=0.01,  upper=100.0, mu=10.0,  sigma=5.0),
    ),

    # Cutoff+BB-function priors
    priors_cutoffbb = dict(
        K_1     = dict(lower=1e-6, upper=1e1),
        index_1 = dict(lower=-2.0, upper=1.0,  mu=-1.0,  sigma=0.5),
        xc_1    = dict(lower=10.0, upper=1000, mu=300.0, sigma=200.0),
        K_2     = dict(lower=1e-9, upper=1e-2),
        kT_2    = dict(lower=0.1,  upper=50.0, mu=10.0,  sigma=7.0),
    ),

    # Synchrotron-function priors
    priors_synchrotron = dict(
        index_bounds = (1, 5),
        gamma_min = 1e5,
        gamma_max = 1e8,
        bulk_gamma = 300,
        K = dict(lower=1e-7, upper=1e-1),
        B = dict(lower=1.0, upper=100.0),
        index = dict(lower=2.0, upper=5.0, mu=3.0, sigma=1.0),
        gamma_cool = dict(lower=1e5, upper=1e8),
    ),

    # ──────────────────────────────────────────────────────────────────────────────────────────
    
    # MultiNest settings
    n_live_points  = 50,             # increase to >=1000 for publication quality

    # output
    outdir         = 'gbm_pipeline_output',
)

# ── Create common output directories ──────────────────────────────────
COMMON_PATHS = {
    'root': CONFIG['outdir'],                 # For common .h5 detector files
    'lc':   os.path.join(CONFIG['outdir'], 'lightcurves'), # For shared lightcurves
}

for d in COMMON_PATHS.values():
    os.makedirs(d, exist_ok=True)

# ── Helper to dynamically create model subdirectories ─────────────────
def get_model_paths(model_name):
    """Creates and returns specific paths for a given spectral model."""
    m_dir = os.path.join(COMMON_PATHS['root'], model_name)
    m_paths = {
        'root':    m_dir,
        'spectra': os.path.join(m_dir, 'spectra'),
        'corners': os.path.join(m_dir, 'corners'),
    }
    for d in m_paths.values():
        os.makedirs(d, exist_ok=True)
    return m_paths

print('Config loaded. Common output directories created at:', CONFIG['outdir'])


# %% [markdown]
# ## Cell 3 — Helper Function Definitions

# %%
def pick_reference_detector(detectors):
    """Select the first NaI detector as the BB reference, fallback to first overall."""
    nai = [d for d in detectors if d.startswith('n')]
    return nai[0] if nai else detectors[0]


def build_cutoffpl_model(ra, dec, priors_cutoffpl):
    """Return a fresh Band-function Model with the configured priors."""
    band  = Cutoff_powerlaw()
    model = Model(PointSource('GBM', ra, dec, spectral_shape=band))
    p = priors_cutoffpl
    model.GBM.spectrum.main.Cutoff_powerlaw.K.prior = Log_uniform_prior(
        lower_bound=p['K']['lower'], upper_bound=p['K']['upper'])
    model.GBM.spectrum.main.Cutoff_powerlaw.index.prior = Truncated_gaussian(
        lower_bound=p['index']['lower'], upper_bound=p['index']['upper'],
        mu=p['index']['mu'], sigma=p['index']['sigma'])
    model.GBM.spectrum.main.Cutoff_powerlaw.xc.prior = Truncated_gaussian(
        lower_bound=p['xc']['lower'], upper_bound=p['xc']['upper'],
        mu=p['xc']['mu'], sigma=p['xc']['sigma'])
    return model

def build_band_model(ra, dec, priors_band):
    """Return a fresh Band-function Model with the configured priors."""
    band  = Band()
    model = Model(PointSource('GBM', ra, dec, spectral_shape=band))
    p = priors_band
    model.GBM.spectrum.main.Band.K.prior = Log_uniform_prior(
        lower_bound=p['K']['lower'], upper_bound=p['K']['upper'])
    model.GBM.spectrum.main.Band.alpha.bounds = (
        p['alpha']['lower'], p['alpha']['upper'])
    model.GBM.spectrum.main.Band.alpha.prior = Truncated_gaussian(
        lower_bound=p['alpha']['lower'], upper_bound=p['alpha']['upper'],
        mu=p['alpha']['mu'], sigma=p['alpha']['sigma'])
    model.GBM.spectrum.main.Band.xp.prior = Truncated_gaussian(
        lower_bound=p['xp']['lower'], upper_bound=p['xp']['upper'],
        mu=p['xp']['mu'], sigma=p['xp']['sigma'])
    model.GBM.spectrum.main.Band.beta.prior = Truncated_gaussian(
        lower_bound=p['beta']['lower'], upper_bound=p['beta']['upper'],
        mu=p['beta']['mu'], sigma=p['beta']['sigma'])
    return model

def build_bandbb_model(ra, dec, p):
    mod = Band() + Blackbody()
    model = Model(PointSource('GBM', ra, dec, spectral_shape=mod))
    model.GBM.spectrum.main.composite.alpha_1.bounds = (-2, 1)
    model.GBM.spectrum.main.composite.K_1.prior = Log_uniform_prior(lower_bound=p['K_1']['lower'], upper_bound=p['K_1']['upper'])
    model.GBM.spectrum.main.composite.alpha_1.prior = Truncated_gaussian(lower_bound=p['alpha_1']['lower'], upper_bound=p['alpha_1']['upper'], mu=p['alpha_1']['mu'], sigma=p['alpha_1']['sigma'])
    model.GBM.spectrum.main.composite.xp_1.prior = Truncated_gaussian(lower_bound=p['xp_1']['lower'], upper_bound=p['xp_1']['upper'], mu=p['xp_1']['mu'], sigma=p['xp_1']['sigma'])
    model.GBM.spectrum.main.composite.beta_1.prior = Truncated_gaussian(lower_bound=p['beta_1']['lower'], upper_bound=p['beta_1']['upper'], mu=p['beta_1']['mu'], sigma=p['beta_1']['sigma'])
    model.GBM.spectrum.main.composite.K_2.prior = Log_uniform_prior(lower_bound=p['K_2']['lower'], upper_bound=p['K_2']['upper'])
    model.GBM.spectrum.main.composite.kT_2.prior = Truncated_gaussian(lower_bound=p['kT_2']['lower'], upper_bound=p['kT_2']['upper'], mu=p['kT_2']['mu'], sigma=p['kT_2']['sigma'])
    return model

def build_cutoffbb_model(ra, dec, p):
    band = Cutoff_powerlaw() + Blackbody()
    model = Model(PointSource('GBM', ra, dec, spectral_shape=band))
    model.GBM.spectrum.main.composite.K_1.prior = Log_uniform_prior(lower_bound=p['K_1']['lower'], upper_bound=p['K_1']['upper'])
    model.GBM.spectrum.main.composite.index_1.prior = Truncated_gaussian(lower_bound=p['index_1']['lower'], upper_bound=p['index_1']['upper'], mu=p['index_1']['mu'], sigma=p['index_1']['sigma'])
    model.GBM.spectrum.main.composite.xc_1.prior = Truncated_gaussian(lower_bound=p['xc_1']['lower'], upper_bound=p['xc_1']['upper'], mu=p['xc_1']['mu'], sigma=p['xc_1']['sigma'])
    model.GBM.spectrum.main.composite.K_2.prior = Log_uniform_prior(lower_bound=p['K_2']['lower'], upper_bound=p['K_2']['upper'])
    model.GBM.spectrum.main.composite.kT_2.prior = Truncated_gaussian(lower_bound=p['kT_2']['lower'], upper_bound=p['kT_2']['upper'], mu=p['kT_2']['mu'], sigma=p['kT_2']['sigma'])
    return model

def build_synchrotron_model(ra, dec, p):
    mod = pynchrotron.SynchrotronNumerical()
    model = Model(PointSource('Synchrotron', ra, dec, spectral_shape=mod))
    model.Synchrotron.spectrum.main.SynchrotronNumerical.index.bounds = p['index_bounds']
    model.Synchrotron.spectrum.main.SynchrotronNumerical.gamma_min.fix = p['gamma_min']
    model.Synchrotron.spectrum.main.SynchrotronNumerical.gamma_max.fix = p['gamma_max']
    model.Synchrotron.spectrum.main.SynchrotronNumerical.bulk_gamma.fix = p['bulk_gamma']
    model.Synchrotron.spectrum.main.SynchrotronNumerical.K.prior = Log_uniform_prior(lower_bound=p['K']['lower'], upper_bound=p['K']['upper'])
    model.Synchrotron.spectrum.main.SynchrotronNumerical.B.prior = Log_uniform_prior(lower_bound=p['B']['lower'], upper_bound=p['B']['upper'])
    model.Synchrotron.spectrum.main.SynchrotronNumerical.index.prior = Truncated_gaussian(lower_bound=p['index']['lower'], upper_bound=p['index']['upper'], mu=p['index']['mu'], sigma=p['index']['sigma'])
    model.Synchrotron.spectrum.main.SynchrotronNumerical.gamma_cool.prior = Log_uniform_prior(lower_bound=p['gamma_cool']['lower'], upper_bound=p['gamma_cool']['upper'])
    return model

print('Helper functions defined.')


# %% [markdown]
# ## Cell 4 — Detector Selection (ROV < 60°)

# %%
log.info('Querying Fermi-GBM catalog for %s ...', CONFIG['grb_name'])

cat = FermiGBMBurst()
cat.query_sources(CONFIG['grb_name'])
grb_info = cat.get_detector_information()[CONFIG['grb_name']]

# ── Auto-derive missing configs ───────────────────────────────────────
if 'trigger_name' not in CONFIG:
    CONFIG['trigger_name'] = CONFIG['grb_name'].replace('GRB', 'bn')
if 'ra' not in CONFIG:
    CONFIG['ra'] = float(grb_info['ra'])
if 'dec' not in CONFIG:
    CONFIG['dec'] = float(grb_info['dec'])

log.info('Auto-derived CONFIG -> trigger: %s, ra: %.3f, dec: %.3f', 
         CONFIG['trigger_name'], CONFIG['ra'], CONFIG['dec'])

# ── Fetch catalog source / background intervals ───────────────────────
source_interval     = grb_info['source']['fluence']
background_interval = grb_info['background']['full']
log.info('Source interval    : %s', source_interval)
log.info('Background interval: %s', background_interval)

# ── Manually Specify Detectors ────────────────────────────────────────
# Valid NaI detectors: 'n0' through 'n9', 'na', 'nb'
# Valid BGO detectors: 'b0', 'b1'
detectors = ['n6', 'n7', 'n8', 'b1']  # <-- Change these to your preferred detectors

log.info('Manually selected detectors done: %s', detectors)


# ── Quick sanity check ────────────────────────────────────────────────
assert len(detectors) > 0, 'No detectors selected. Please add them to the list.'
assert any(d.startswith('n') for d in detectors), 'Need at least one NaI detector for the fit.'


# %% [markdown]
# ## Cell 5 — Download GBM Data

# %%
log.info('Downloading TTE / CSPEC / RSP for: %s', detectors)
dload = download_GBM_trigger_data(CONFIG['trigger_name'], detectors=detectors)

# ── Verify all expected file types are present ────────────────────────
for det in detectors:
    for ftype in ('tte', 'cspec', 'rsp'):
        assert ftype in dload[det], f'Missing {ftype} file for detector {det}'

print('Download complete. Files available for:', list(dload.keys()))


# %% [markdown]
# ## Cell 6 — Build Time Series & Fit Backgrounds

# %%
time_series = {}

for det in detectors:
    log.info('Processing detector %s ...', det)

    # ── Step 6a: fit background polynomial on CSPEC ───────────────────
    ts_cspec = TimeSeriesBuilder.from_gbm_cspec_or_ctime(
        det,
        cspec_or_ctime_file = dload[det]['cspec'],
        rsp_file            = dload[det]['rsp'],
    )
    bkg_segs = [s.strip() for s in background_interval.split(',')]
    ts_cspec.set_background_interval(*bkg_segs)

    bkg_file = os.path.join(COMMON_PATHS['root'], f'{det}_bkg.h5')
    ts_cspec.save_background(bkg_file, overwrite=True)
    log.info('  Background saved to %s', bkg_file)

    # ── Step 6b: load TTE with restored background ────────────────────
    ts_tte = TimeSeriesBuilder.from_gbm_tte(
        det,
        tte_file           = dload[det]['tte'],
        rsp_file           = dload[det]['rsp'],
        restore_background = bkg_file,
    )
    ts_tte.set_active_time_interval(source_interval)
    time_series[det] = ts_tte
    log.info('  TTE loaded and active interval set.')

print('Time series built for:', list(time_series.keys()))

# ── Quick sanity: view lightcurve of first NaI ────────────────────────
first_nai = next(d for d in detectors if d.startswith('n'))
fig = time_series[first_nai].view_lightcurve(-40, 100)

fig.savefig(os.path.join(COMMON_PATHS['lc'], f'{first_nai}_raw_lc.png'), dpi=100, bbox_inches='tight')
plt.close('all')
print(f'Raw lightcurve saved for {first_nai}')


# %% [markdown]
# ## Cell 7 — Bayesian-Block Time Binning

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    ref_det = pick_reference_detector(detectors)
    ts_ref  = time_series[ref_det]
    log.info('Reference detector for BB: %s', ref_det)

    # ── Run Bayesian blocks ────────────────────────────────────────────────
    ts_ref.create_time_bins(
        CONFIG['bb_t_start'], CONFIG['bb_t_stop'],
        method='bayesblocks', use_background=True, p0=CONFIG['bb_p0'],
    )

    print('Raw BB bins:')
    print('  starts:', ts_ref.bins.starts)
    print('  stops :', ts_ref.bins.stops)
    print('  widths:', ts_ref.bins.widths)

    fig_raw = ts_ref.view_lightcurve(use_binner=True)
    fig_raw.savefig(os.path.join(COMMON_PATHS['lc'], f'{ref_det}_bb_raw.png'), dpi=120, bbox_inches='tight')
    plt.close('all')


# %% [markdown]
# ## Cell 8 — Merge Sub-Threshold Bins

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    min_w = CONFIG['min_bin_width']

    bad_bins = {i for i, w in enumerate(ts_ref.bins.widths) if w < min_w}
    print(f'Bins narrower than {min_w} s (to be merged): {sorted(bad_bins)}')

    edges = [ts_ref.bins.starts[0]]
    for i, b in enumerate(ts_ref.bins):
        if i not in bad_bins:
            edges.append(b.stop)

    bb_starts = np.array(edges[:-1])
    bb_stops  = np.array(edges[1:])
    n_bins    = len(bb_starts)

    ts_ref.create_time_bins(bb_starts, bb_stops, method='custom')

    print(f'Clean bins ({n_bins} total):')
    print('  starts:', bb_starts)
    print('  stops :', bb_stops)
    print('  widths:', bb_stops - bb_starts)

    fig_clean = ts_ref.view_lightcurve(use_binner=True)
    fig_clean.savefig(os.path.join(COMMON_PATHS['lc'], f'{ref_det}_bb_clean.png'), dpi=120, bbox_inches='tight')
    plt.close('all')

    assert n_bins > 0, 'Zero bins after merging — check bb_t_start / bb_t_stop or p0'
    print(f'\nReady to fit {n_bins} time bin(s).')


# %% [markdown]
# ## Cell 9 — Sync Detectors & Build Time-Resolved Plugins

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    time_resolved_plugins = {}   # {det: [plugin_bin_0, plugin_bin_1, ...]}

    for det in detectors:
        ts = time_series[det]
        ts.read_bins(ts_ref)                         # sync to reference bins
        plugin_list = ts.to_spectrumlike(from_bins=True)
        time_resolved_plugins[det] = plugin_list
        log.info('%s  ->  %d time-bin plugins', det, len(plugin_list))

    # ── Sanity: all detectors should have the same number of bins ─────────
    bin_counts = {d: len(v) for d, v in time_resolved_plugins.items()}
    print('Plugin counts per detector:', bin_counts)
    assert len(set(bin_counts.values())) == 1, 'Bin count mismatch across detectors!'
    print('All detectors in sync. Ready for spectral fitting.')


# %% [markdown]
# ## Cell 10 — Time-Resolved Spectral Fitting (CutoffPLModel)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    MODEL_PATHS = get_model_paths('cutoffpl')                             # Change this for each new model

    PARAM_NAMES = [
        'B_K',     'B_K_n_err',     'B_K_p_err',     'B_K_avg_err',
        'B_alpha', 'B_alpha_n_err', 'B_alpha_p_err',  'B_alpha_avg_err',
        'E_c',  'E_c_n_err',  'E_c_p_err',   'E_c_avg_err',
        'flux',    'flux_low',      'flux_hi', 'DIC', 'PDIC', 'BIC'
    ]

    fit_rows = []

    for i in range(n_bins):
        log.info('==== Time bin %d / %d ====', i + 1, n_bins)

        # ── Assemble per-detector plugins for this bin ────────────────────
        det_plugins = []
        for det in detectors:
            p = time_resolved_plugins[det][i]
            if det.startswith('b'):
                p.set_active_measurements(CONFIG['bgo_range'])
            else:
                ranges = [r.strip() for r in CONFIG['nai_range'].split(',')]
                p.set_active_measurements(*ranges)
            det_plugins.append(p)

        data_list = DataList(*det_plugins)

        # ── Fresh Band model per bin ──────────────────────────────────────
        model = build_cutoffpl_model(CONFIG['ra'], CONFIG['dec'], CONFIG['priors_cutoffpl'])

        # ── Bayesian fit with MultiNest ───────────────────────────────────
        bayes = BayesianAnalysis(model, data_list)
        bayes.set_sampler('multinest')
        bayes.sampler.setup(n_live_points=CONFIG['n_live_points'])
        bayes.sample()

        # ── νFν spectrum plot ─────────────────────────────────────────────
        fig_spec = plot_spectra(
            bayes.results,
            ene_min=1.0, ene_max=1e5, num_ene=200,
            contour_colors='white',
            contour_style_kwargs=dict(alpha=0.2),
            use_components=True, components_to_use=['total'],
            flux_unit='keV2/(cm2 s keV)',
        )
        ax = plt.gca()
        ax.set_xlabel('Energy (keV)', fontsize=16)
        ax.set_ylabel(r'$\nu\,F_\nu\;[\mathrm{keV\,cm^{-2}\,s^{-1}}]$', fontsize=16)
        ax.set_title(f"{CONFIG['trigger_name']}  bin {i+1}/{n_bins}  "
                     f"[{bb_starts[i]:.3f}, {bb_stops[i]:.3f}] s", fontsize=13)
        fig_spec.savefig(
            os.path.join(MODEL_PATHS['spectra'], f'cutoffpl_spec_bin{i+1:03d}.png'),
            dpi=120, bbox_inches='tight')
        plt.close('all')

        # ── Corner plot ───────────────────────────────────────────────────
        fig_corner = bayes.results.corner_plot()
        fig_corner.savefig(
            os.path.join(MODEL_PATHS['corners'], f'cutoffpl_corner_bin{i+1:03d}.png'),
            dpi=100, bbox_inches='tight')
        plt.close('all')

        # ── Integrated flux ───────────────────────────────────────────────
        flux_df = bayes.results.get_point_source_flux(10 * u.keV, 10 * u.MeV)
        log.info('  10 keV-10 MeV flux = %s', flux_df)

        # Extract the raw floats by stripping the astropy units (.value)
        f_val = flux_df['flux'].iloc[0].value
        f_low = flux_df['low bound'].iloc[0].value
        f_hi  = flux_df['hi bound'].iloc[0].value

        # ── Collect parameters ────────────────────────────────────────────
        df  = bayes.results.get_data_frame('equal tail')
        df  = df.drop('unit', axis=1, errors='ignore')
        param_array = df.to_numpy().flatten()

        stat_df = bayes.results.get_statistic_measure_frame()
        dic_val = stat_df.loc['DIC'].iloc[0] if 'DIC' in stat_df.index else np.nan
        pdic_val = stat_df.loc['PDIC'].iloc[0] if 'PDIC' in stat_df.index else np.nan
        bic_val = stat_df.loc['BIC'].iloc[0] if 'BIC' in stat_df.index else np.nan

        full_row = np.append(param_array, [f_val, f_low, f_hi, dic_val, pdic_val, bic_val])
        fit_rows.append(full_row)
        log.info('  Bin %d / %d complete.', i + 1, n_bins)

    print(f'\nAll {n_bins} bins fitted successfully.')


# %% [markdown]
# ## Cell 11 — Save Results to CSV

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    fit_array = np.array(fit_rows)
    cutoffpl_fit_param = Table(fit_array, names=PARAM_NAMES)

    # ── Attach bin edges as extra columns ────────────────────────────────
    cutoffpl_fit_param['t_start'] = bb_starts
    cutoffpl_fit_param['t_stop']  = bb_stops
    cutoffpl_fit_param['dt']      = bb_stops - bb_starts

    csv_path = os.path.join(MODEL_PATHS['root'], 'fit_params.csv')
    cutoffpl_fit_param.write(csv_path, overwrite=True)

    print(f'Results saved to: {csv_path}')
    print(f'Shape: {fit_array.shape}  ({n_bins} bins x 16 parameters)')
    cutoffpl_fit_param


# %% [markdown]
# ## Cell 10 — Time-Resolved Spectral Fitting (Band Model)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    MODEL_PATHS = get_model_paths('band')

    PARAM_NAMES = [
        'B_K',     'B_K_n_err',     'B_K_p_err',     'B_K_avg_err',
        'B_alpha', 'B_alpha_n_err', 'B_alpha_p_err',  'B_alpha_avg_err',
        'E_peak',  'E_peak_n_err',  'E_peak_p_err',   'E_peak_avg_err',
        'beta',    'beta_n_err',    'beta_p_err',      'beta_avg_err',
        'flux',    'flux_low',      'flux_hi', 'DIC', 'PDIC', 'BIC'
    ]

    fit_rows = []

    for i in range(n_bins):
        log.info('==== Time bin %d / %d ====', i + 1, n_bins)

        # ── Assemble per-detector plugins for this bin ────────────────────
        det_plugins = []
        for det in detectors:
            p = time_resolved_plugins[det][i]
            if det.startswith('b'):
                p.set_active_measurements(CONFIG['bgo_range'])
            else:
                ranges = [r.strip() for r in CONFIG['nai_range'].split(',')]
                p.set_active_measurements(*ranges)
            det_plugins.append(p)

        data_list = DataList(*det_plugins)

        # ── Fresh Band model per bin ──────────────────────────────────────
        model = build_band_model(CONFIG['ra'], CONFIG['dec'], CONFIG['priors_band'])

        # ── Bayesian fit with MultiNest ───────────────────────────────────
        bayes = BayesianAnalysis(model, data_list)
        bayes.set_sampler('multinest')
        bayes.sampler.setup(n_live_points=CONFIG['n_live_points'])
        bayes.sample()

        # ── νFν spectrum plot ─────────────────────────────────────────────
        fig_spec = plot_spectra(
            bayes.results,
            ene_min=1.0, ene_max=1e5, num_ene=200,
            contour_colors='white',
            contour_style_kwargs=dict(alpha=0.2),
            use_components=True, components_to_use=['total'],
            flux_unit='keV2/(cm2 s keV)',
        )
        ax = plt.gca()
        ax.set_xlabel('Energy (keV)', fontsize=16)
        ax.set_ylabel(r'$\nu\,F_\nu\;[\mathrm{keV\,cm^{-2}\,s^{-1}}]$', fontsize=16)
        ax.set_title(f"{CONFIG['trigger_name']}  bin {i+1}/{n_bins}  "
                     f"[{bb_starts[i]:.3f}, {bb_stops[i]:.3f}] s", fontsize=13)
        fig_spec.savefig(
            os.path.join(MODEL_PATHS['spectra'], f'band_spec_bin{i+1:03d}.png'),
            dpi=120, bbox_inches='tight')
        plt.close('all')

        # ── Corner plot ───────────────────────────────────────────────────
        fig_corner = bayes.results.corner_plot()
        fig_corner.savefig(
            os.path.join(MODEL_PATHS['corners'], f'band_corner_bin{i+1:03d}.png'),
            dpi=100, bbox_inches='tight')
        plt.close('all')

        # ── Integrated flux ───────────────────────────────────────────────
        flux_df = bayes.results.get_point_source_flux(10 * u.keV, 10 * u.MeV)
        log.info('  10 keV-10 MeV flux = %s', flux_df)

        # Extract the raw floats by stripping the astropy units (.value)
        f_val = flux_df['flux'].iloc[0].value
        f_low = flux_df['low bound'].iloc[0].value
        f_hi  = flux_df['hi bound'].iloc[0].value

        # ── Collect parameters ────────────────────────────────────────────
        df  = bayes.results.get_data_frame('equal tail')
        df  = df.drop('unit', axis=1, errors='ignore')
        param_array = df.to_numpy().flatten()

        stat_df = bayes.results.get_statistic_measure_frame()
        dic_val = stat_df.loc['DIC'].iloc[0] if 'DIC' in stat_df.index else np.nan
        pdic_val = stat_df.loc['PDIC'].iloc[0] if 'PDIC' in stat_df.index else np.nan
        bic_val = stat_df.loc['BIC'].iloc[0] if 'BIC' in stat_df.index else np.nan

        full_row = np.append(param_array, [f_val, f_low, f_hi, dic_val, pdic_val, bic_val])
        fit_rows.append(full_row)
        log.info('  Bin %d / %d complete.', i + 1, n_bins)

    print(f'\nAll {n_bins} bins fitted successfully.')


    # %%
    ## Cell 11 — Save Results to CSV

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    fit_array = np.array(fit_rows)
    band_fit_param = Table(fit_array, names=PARAM_NAMES)

    # ── Attach bin edges as extra columns ────────────────────────────────
    band_fit_param['t_start'] = bb_starts
    band_fit_param['t_stop']  = bb_stops
    band_fit_param['dt']      = bb_stops - bb_starts

    csv_path = os.path.join(MODEL_PATHS['root'], 'fit_params.csv')
    band_fit_param.write(csv_path, overwrite=True)

    print(f'Results saved to: {csv_path}')
    print(f'Shape: {fit_array.shape}  ({n_bins} bins x 16 parameters)')
    band_fit_param


# %% [markdown]
# ## Cell 13 — Time-Resolved Spectral Fitting (Band+BB Model)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    MODEL_PATHS = get_model_paths('bandbb')

    PARAM_NAMES = []
    for p in ['K_1', 'alpha_1', 'xp_1', 'beta_1', 'K_2', 'kT_2']:
        PARAM_NAMES.extend([p, f'{p}_n_err', f'{p}_p_err', f'{p}_avg_err'])
    PARAM_NAMES.extend(['flux', 'flux_low', 'flux_hi', 'DIC', 'PDIC', 'BIC'])

    fit_rows = []

    for i in range(n_bins):
        log.info('==== Time bin %d / %d ====', i + 1, n_bins)

        det_plugins = []
        for det in detectors:
            p = time_resolved_plugins[det][i]
            if det.startswith('b'):
                p.set_active_measurements(CONFIG['bgo_range'])
            else:
                ranges = [r.strip() for r in CONFIG['nai_range'].split(',')]
                p.set_active_measurements(*ranges)
            det_plugins.append(p)

        data_list = DataList(*det_plugins)
        model = build_bandbb_model(CONFIG['ra'], CONFIG['dec'], CONFIG['priors_bandbb'])

        bayes = BayesianAnalysis(model, data_list)
        bayes.set_sampler('multinest')
        bayes.sampler.setup(n_live_points=CONFIG['n_live_points'])
        bayes.sample()

        fig_spec = plot_spectra(
            bayes.results,
            ene_min=1.0, ene_max=1e5, num_ene=200,
            contour_colors='white', contour_style_kwargs=dict(alpha=0.2),
            use_components=True, components_to_use=['total'],
            flux_unit='keV2/(cm2 s keV)',
        )
        ax = plt.gca()
        ax.set_xlabel('Energy (keV)', fontsize=16)
        ax.set_ylabel(r'$\nu\,F_\nu\;[\mathrm{keV\,cm^{-2}\,s^{-1}}]$', fontsize=16)
        ax.set_title(f"{CONFIG['trigger_name']}  bin {i+1}/{n_bins}  [{bb_starts[i]:.3f}, {bb_stops[i]:.3f}] s", fontsize=13)
        fig_spec.savefig(os.path.join(MODEL_PATHS['spectra'], f'bandbb_spec_bin{i+1:03d}.png'), dpi=120, bbox_inches='tight')
        plt.close('all')

        fig_corner = bayes.results.corner_plot()
        fig_corner.savefig(os.path.join(MODEL_PATHS['corners'], f'bandbb_corner_bin{i+1:03d}.png'), dpi=100, bbox_inches='tight')
        plt.close('all')

        flux_df = bayes.results.get_point_source_flux(10 * u.keV, 10 * u.MeV)
        log.info('  10 keV-10 MeV flux = %s', flux_df)
        f_val = flux_df['flux'].iloc[0].value
        f_low = flux_df['low bound'].iloc[0].value
        f_hi  = flux_df['hi bound'].iloc[0].value

        df  = bayes.results.get_data_frame('equal tail')
        df  = df.drop('unit', axis=1, errors='ignore')
        param_array = df.to_numpy().flatten()

        stat_df = bayes.results.get_statistic_measure_frame()
        dic_val = stat_df.loc['DIC'].iloc[0] if 'DIC' in stat_df.index else np.nan
        pdic_val = stat_df.loc['PDIC'].iloc[0] if 'PDIC' in stat_df.index else np.nan
        bic_val = stat_df.loc['BIC'].iloc[0] if 'BIC' in stat_df.index else np.nan

        full_row = np.append(param_array, [f_val, f_low, f_hi, dic_val, pdic_val, bic_val])
        fit_rows.append(full_row)
        log.info('  Bin %d / %d complete.', i + 1, n_bins)

    print(f'\nAll {n_bins} bins fitted successfully.')

    # %%
    ## Cell 14 — Save Results to CSV (Band+BB)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    fit_array = np.array(fit_rows)
    bandbb_fit_param = Table(fit_array, names=PARAM_NAMES)
    bandbb_fit_param['t_start'] = bb_starts
    bandbb_fit_param['t_stop']  = bb_stops
    bandbb_fit_param['dt']      = bb_stops - bb_starts

    csv_path = os.path.join(MODEL_PATHS['root'], 'fit_params.csv')
    bandbb_fit_param.write(csv_path, overwrite=True)
    print(f'Results saved to: {csv_path}')
    print(f'Shape: {fit_array.shape}  ({n_bins} bins x {len(PARAM_NAMES)} parameters)')
    bandbb_fit_param

# %% [markdown]
# ## Cell 15 — Time-Resolved Spectral Fitting (Cutoff+BB Model)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    MODEL_PATHS = get_model_paths('cutoffbb')

    PARAM_NAMES = []
    for p in ['K_1', 'index_1', 'xc_1', 'K_2', 'kT_2']:
        PARAM_NAMES.extend([p, f'{p}_n_err', f'{p}_p_err', f'{p}_avg_err'])
    PARAM_NAMES.extend(['flux', 'flux_low', 'flux_hi', 'DIC', 'PDIC', 'BIC'])

    fit_rows = []

    for i in range(n_bins):
        log.info('==== Time bin %d / %d ====', i + 1, n_bins)

        det_plugins = []
        for det in detectors:
            p = time_resolved_plugins[det][i]
            if det.startswith('b'):
                p.set_active_measurements(CONFIG['bgo_range'])
            else:
                ranges = [r.strip() for r in CONFIG['nai_range'].split(',')]
                p.set_active_measurements(*ranges)
            det_plugins.append(p)

        data_list = DataList(*det_plugins)
        model = build_cutoffbb_model(CONFIG['ra'], CONFIG['dec'], CONFIG['priors_cutoffbb'])

        bayes = BayesianAnalysis(model, data_list)
        bayes.set_sampler('multinest')
        bayes.sampler.setup(n_live_points=CONFIG['n_live_points'])
        bayes.sample()

        fig_spec = plot_spectra(
            bayes.results,
            ene_min=1.0, ene_max=1e5, num_ene=200,
            contour_colors='white', contour_style_kwargs=dict(alpha=0.2),
            use_components=True, components_to_use=['total'],
            flux_unit='keV2/(cm2 s keV)',
        )
        ax = plt.gca()
        ax.set_xlabel('Energy (keV)', fontsize=16)
        ax.set_ylabel(r'$\nu\,F_\nu\;[\mathrm{keV\,cm^{-2}\,s^{-1}}]$', fontsize=16)
        ax.set_title(f"{CONFIG['trigger_name']}  bin {i+1}/{n_bins}  [{bb_starts[i]:.3f}, {bb_stops[i]:.3f}] s", fontsize=13)
        fig_spec.savefig(os.path.join(MODEL_PATHS['spectra'], f'cutoffbb_spec_bin{i+1:03d}.png'), dpi=120, bbox_inches='tight')
        plt.close('all')

        fig_corner = bayes.results.corner_plot()
        fig_corner.savefig(os.path.join(MODEL_PATHS['corners'], f'cutoffbb_corner_bin{i+1:03d}.png'), dpi=100, bbox_inches='tight')
        plt.close('all')

        flux_df = bayes.results.get_point_source_flux(10 * u.keV, 10 * u.MeV)
        log.info('  10 keV-10 MeV flux = %s', flux_df)
        f_val = flux_df['flux'].iloc[0].value
        f_low = flux_df['low bound'].iloc[0].value
        f_hi  = flux_df['hi bound'].iloc[0].value

        df  = bayes.results.get_data_frame('equal tail')
        df  = df.drop('unit', axis=1, errors='ignore')
        param_array = df.to_numpy().flatten()

        stat_df = bayes.results.get_statistic_measure_frame()
        dic_val = stat_df.loc['DIC'].iloc[0] if 'DIC' in stat_df.index else np.nan
        pdic_val = stat_df.loc['PDIC'].iloc[0] if 'PDIC' in stat_df.index else np.nan
        bic_val = stat_df.loc['BIC'].iloc[0] if 'BIC' in stat_df.index else np.nan

        full_row = np.append(param_array, [f_val, f_low, f_hi, dic_val, pdic_val, bic_val])
        fit_rows.append(full_row)
        log.info('  Bin %d / %d complete.', i + 1, n_bins)

    print(f'\nAll {n_bins} bins fitted successfully.')

    # %%
    ## Cell 16 — Save Results to CSV (Cutoff+BB)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    fit_array = np.array(fit_rows)
    cutoffbb_fit_param = Table(fit_array, names=PARAM_NAMES)
    cutoffbb_fit_param['t_start'] = bb_starts
    cutoffbb_fit_param['t_stop']  = bb_stops
    cutoffbb_fit_param['dt']      = bb_stops - bb_starts

    csv_path = os.path.join(MODEL_PATHS['root'], 'fit_params.csv')
    cutoffbb_fit_param.write(csv_path, overwrite=True)
    print(f'Results saved to: {csv_path}')
    print(f'Shape: {fit_array.shape}  ({n_bins} bins x {len(PARAM_NAMES)} parameters)')
    cutoffbb_fit_param

# %% [markdown]
# ## Cell 17 — Time-Resolved Spectral Fitting (Synchrotron Model)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    MODEL_PATHS = get_model_paths('synchrotron')

    PARAM_NAMES = []
    for p in ['K', 'B', 'index', 'gamma_cool']:
        PARAM_NAMES.extend([p, f'{p}_n_err', f'{p}_p_err', f'{p}_avg_err'])
    PARAM_NAMES.extend(['flux', 'flux_low', 'flux_hi', 'DIC', 'PDIC', 'BIC'])

    fit_rows = []

    for i in range(n_bins):
        log.info('==== Time bin %d / %d ====', i + 1, n_bins)

        det_plugins = []
        for det in detectors:
            p = time_resolved_plugins[det][i]
            if det.startswith('b'):
                p.set_active_measurements(CONFIG['bgo_range'])
            else:
                ranges = [r.strip() for r in CONFIG['nai_range'].split(',')]
                p.set_active_measurements(*ranges)
            det_plugins.append(p)

        data_list = DataList(*det_plugins)
        model = build_synchrotron_model(CONFIG['ra'], CONFIG['dec'], CONFIG['priors_synchrotron'])

        bayes = BayesianAnalysis(model, data_list)
        bayes.set_sampler('multinest')
        bayes.sampler.setup(n_live_points=CONFIG['n_live_points'])
        bayes.sample()

        fig_spec = plot_spectra(
            bayes.results,
            ene_min=1.0, ene_max=1e5, num_ene=200,
            contour_colors='white', contour_style_kwargs=dict(alpha=0.2),
            use_components=True, components_to_use=['total'],
            flux_unit='keV2/(cm2 s keV)',
        )
        ax = plt.gca()
        ax.set_xlabel('Energy (keV)', fontsize=16)
        ax.set_ylabel(r'$\nu\,F_\nu\;[\mathrm{keV\,cm^{-2}\,s^{-1}}]$', fontsize=16)
        ax.set_title(f"{CONFIG['trigger_name']}  bin {i+1}/{n_bins}  [{bb_starts[i]:.3f}, {bb_stops[i]:.3f}] s", fontsize=13)
        fig_spec.savefig(os.path.join(MODEL_PATHS['spectra'], f'synchrotron_spec_bin{i+1:03d}.png'), dpi=120, bbox_inches='tight')
        plt.close('all')

        fig_corner = bayes.results.corner_plot()
        fig_corner.savefig(os.path.join(MODEL_PATHS['corners'], f'synchrotron_corner_bin{i+1:03d}.png'), dpi=100, bbox_inches='tight')
        plt.close('all')

        flux_df = bayes.results.get_point_source_flux(10 * u.keV, 10 * u.MeV)
        log.info('  10 keV-10 MeV flux = %s', flux_df)
        f_val = flux_df['flux'].iloc[0].value
        f_low = flux_df['low bound'].iloc[0].value
        f_hi  = flux_df['hi bound'].iloc[0].value

        df  = bayes.results.get_data_frame('equal tail')
        df  = df.drop('unit', axis=1, errors='ignore')
        param_array = df.to_numpy().flatten()

        stat_df = bayes.results.get_statistic_measure_frame()
        dic_val = stat_df.loc['DIC'].iloc[0] if 'DIC' in stat_df.index else np.nan
        pdic_val = stat_df.loc['PDIC'].iloc[0] if 'PDIC' in stat_df.index else np.nan
        bic_val = stat_df.loc['BIC'].iloc[0] if 'BIC' in stat_df.index else np.nan

        full_row = np.append(param_array, [f_val, f_low, f_hi, dic_val, pdic_val, bic_val])
        fit_rows.append(full_row)
        log.info('  Bin %d / %d complete.', i + 1, n_bins)

    print(f'\nAll {n_bins} bins fitted successfully.')

    # %%
    ## Cell 18 — Save Results to CSV (Synchrotron)

# %%
if CONFIG['analysis_type'] in ['time_resolved', 'both']:
    fit_array = np.array(fit_rows)
    synchrotron_fit_param = Table(fit_array, names=PARAM_NAMES)
    synchrotron_fit_param['t_start'] = bb_starts
    synchrotron_fit_param['t_stop']  = bb_stops
    synchrotron_fit_param['dt']      = bb_stops - bb_starts

    csv_path = os.path.join(MODEL_PATHS['root'], 'fit_params.csv')
    synchrotron_fit_param.write(csv_path, overwrite=True)
    print(f'Results saved to: {csv_path}')
    print(f'Shape: {fit_array.shape}  ({n_bins} bins x {len(PARAM_NAMES)} parameters)')
    synchrotron_fit_param

# %% [markdown]
# ## Cell 19 — Time-Integrated Analysis (All Models)

# %%
if CONFIG['analysis_type'] in ['time_integrated', 'both']:
    log.info('==== Starting Time-Integrated Analysis ====')
    
    # Define single bin
    ti_start = np.array([CONFIG['ti_t_start']])
    ti_stop  = np.array([CONFIG['ti_t_stop']])
    ti_dt = ti_stop - ti_start
    
    ti_plugins = {}
    for det in detectors:
        ts = time_series[det]
        ts.create_time_bins(ti_start, ti_stop, method='custom')
        ti_plugins[det] = ts.to_spectrumlike(from_bins=True)[0]
        
    models_to_run = [
        ('cutoffpl', build_cutoffpl_model, CONFIG['priors_cutoffpl']),
        ('band', build_band_model, CONFIG['priors_band']),
        ('bandbb', build_bandbb_model, CONFIG['priors_bandbb']),
        ('cutoffbb', build_cutoffbb_model, CONFIG['priors_cutoffbb']),
        ('synchrotron', build_synchrotron_model, CONFIG['priors_synchrotron'])
    ]
    
    for model_name, builder_func, prior_dict in models_to_run:
        log.info(f'--- Fitting Time-Integrated {model_name.upper()} ---')
        MODEL_PATHS = get_model_paths(f'ti_{model_name}') 
        
        det_plugins = []
        for det in detectors:
            p = ti_plugins[det]
            if det.startswith('b'):
                p.set_active_measurements(CONFIG['bgo_range'])
            else:
                ranges = [r.strip() for r in CONFIG['nai_range'].split(',')]
                p.set_active_measurements(*ranges)
            det_plugins.append(p)
            
        data_list = DataList(*det_plugins)
        model = builder_func(CONFIG['ra'], CONFIG['dec'], prior_dict)
        
        bayes = BayesianAnalysis(model, data_list)
        bayes.set_sampler('multinest')
        bayes.sampler.setup(n_live_points=CONFIG['n_live_points'])
        bayes.sample()
        
        fig_spec = plot_spectra(
            bayes.results,
            ene_min=1.0, ene_max=1e5, num_ene=200,
            contour_colors='white', contour_style_kwargs=dict(alpha=0.2),
            use_components=True, components_to_use=['total'],
            flux_unit='keV2/(cm2 s keV)',
        )
        ax = plt.gca()
        ax.set_xlabel('Energy (keV)', fontsize=16)
        ax.set_ylabel(r'$\nu\,F_\nu\;[\mathrm{keV\,cm^{-2}\,s^{-1}}]$', fontsize=16)
        ax.set_title(f"{CONFIG['trigger_name']} Time-Integrated [{ti_start[0]:.3f}, {ti_stop[0]:.3f}] s", fontsize=13)
        fig_spec.savefig(os.path.join(MODEL_PATHS['spectra'], f'{model_name}_ti_spec.png'), dpi=120, bbox_inches='tight')
        plt.close('all')
        
        fig_corner = bayes.results.corner_plot()
        fig_corner.savefig(os.path.join(MODEL_PATHS['corners'], f'{model_name}_ti_corner.png'), dpi=100, bbox_inches='tight')
        plt.close('all')
        
        flux_df = bayes.results.get_point_source_flux(10 * u.keV, 10 * u.MeV)
        f_val = flux_df['flux'].iloc[0].value
        f_low = flux_df['low bound'].iloc[0].value
        f_hi  = flux_df['hi bound'].iloc[0].value
        
        df = bayes.results.get_data_frame('equal tail')
        df = df.drop('unit', axis=1, errors='ignore')
        param_array = df.to_numpy().flatten()
        
        import pandas as pd
        param_names = df.index.tolist()
        columns = []
        for p in param_names:
            columns.extend([p, f'{p}_n_err', f'{p}_p_err', f'{p}_avg_err'])
        columns.extend(['flux', 'flux_low', 'flux_hi', 't_start', 't_stop', 'dt', 'DIC', 'PDIC', 'BIC'])
        
        stat_df = bayes.results.get_statistic_measure_frame()
        dic_val = stat_df.loc['DIC'].iloc[0] if 'DIC' in stat_df.index else np.nan
        pdic_val = stat_df.loc['PDIC'].iloc[0] if 'PDIC' in stat_df.index else np.nan
        bic_val = stat_df.loc['BIC'].iloc[0] if 'BIC' in stat_df.index else np.nan

        full_row = np.append(param_array, [f_val, f_low, f_hi, ti_start[0], ti_stop[0], ti_dt[0], dic_val, pdic_val, bic_val])
        ti_table = Table(rows=[full_row], names=columns)
        csv_path = os.path.join(MODEL_PATHS['root'], 'ti_fit_params.csv')
        ti_table.write(csv_path, overwrite=True)
        log.info(f'  Saved Time-Integrated {model_name} results to {csv_path}')


# %% [markdown]
# ## Cell 20 — Pipeline Summary (All Models)

# %%
print('=' * 60)
print('  GBM PIPELINE SUMMARY')
print('=' * 60)
print(f'  Trigger          : {CONFIG["trigger_name"]}')
print(f'  Detectors used   : {detectors}')
if 'n_bins' in locals():
    print(f'  Time bins fitted : {n_bins}')
    print(f'  Bin edges (s)    : starts={bb_starts}, stops={bb_stops}')
else:
    print('  Time-resolved    : (skipped)')
if CONFIG['analysis_type'] in ['time_integrated', 'both']:
    print(f'  Time-integrated  : {CONFIG["ti_t_start"]} to {CONFIG["ti_t_stop"]} s')
print('-' * 60)
print('  COMMON OUTPUTS:')
print(f'  Base Directory   : {COMMON_PATHS["root"]}')
print(f'  Lightcurves      : {COMMON_PATHS["lc"]}')
print('-' * 60)
print('  MODEL SPECIFIC OUTPUTS:')

found_models = False
if os.path.exists(COMMON_PATHS['root']):
    for item in sorted(os.listdir(COMMON_PATHS['root'])):
        model_dir = os.path.join(COMMON_PATHS['root'], item)
        if os.path.isdir(model_dir) and item != 'lightcurves':
            found_models = True
            print(f'  Model: [{item.upper()}]')
            print(f'    Spectra plots  : {os.path.join(model_dir, "spectra")}')
            print(f'    Corner plots   : {os.path.join(model_dir, "corners")}')
            csv_file = os.path.join(model_dir, 'fit_params.csv') 
            if os.path.exists(csv_file):
                print(f'    TR CSV results : {csv_file}')
            else:
                print(f'    TR CSV results : (Not generated)')
            
            ti_csv_file = os.path.join(model_dir, 'ti_fit_params.csv')
            if os.path.exists(ti_csv_file):
                print(f'    TI CSV results : {ti_csv_file}')
            print()

if not found_models:
    print('  (No model-specific outputs found yet.)')

print('=' * 60)

# %%
