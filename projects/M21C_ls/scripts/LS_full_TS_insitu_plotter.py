# -*- coding: utf-8 -*-
# Refactor of your script to run for all three networks (tags) and reproduce the same outputs.
# Outputs for each network are written into a separate folder to avoid overwriting.
# Figures keep the SAME names and formatting as your originals.

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

# Optional: maps (cartopy) — required for the three scatter maps at the end
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --------------------
# Configuration
# --------------------
expt_name = 'LS'  # used in the bar-chart filename, kept identical
experiment_names = ['LS_OLv8_M36', 'LS_DAv8_M36']

# Run all three networks
INSITU_TAGS = [
    '_SCAN_SM_1d_c1234smv_25yr',
    '_CalVal_M33_SM_1d__25yr',
    '_USCRN_SM_1d_c1234smv_25yr',
]

# Base directories / file patterns (identical to your originals)
STATS_DIR = '../test_data/M21C_land_sweeper/Evaluation/InSitu/output'
RAW_TS_TEMPLATE = os.path.join(STATS_DIR, '{exp}' + '{tag}' + '_raw_timeseries.mat')
STATS_TEMPLATE  = os.path.join(STATS_DIR, '{exp}' + '{tag}' + '_stats.mat')

# --------------------
# Plot constants (identical to your code)
# --------------------
expt_labels = ["CNTL", "DA"]
title_fontsize = 20
label_fontsize = 20
y_tick_label_fontsize = 18

def _finite_1d(x):
    """Return a clean 1-D array with only finite values."""
    x = np.asarray(x).astype(float).ravel()
    return x[np.isfinite(x)]

def compute_means_cis(Bias, BiasLO, BiasUP,
                      RMSE, RMSELO, RMSEUP,
                      R, RLO, RUP,
                      absBias, absBiasLO, absBiasUP,
                      anomR, anomRLO, anomRUP,
                      ubRMSE, ubRMSELO, ubRMSEUP):
    """Compute means, stds, site counts and confidence interval arrays exactly as in the original."""
    # R
    R_mean = np.around(np.nanmean(R, axis=0), decimals=2)
    R_std  = np.around(np.nanstd(R, axis=0), decimals=3)
    num_sites_r = np.sum(~np.isnan(R), axis=0)  # (depth, exp)
    R_CI_LO = np.around(np.nanmean(RLO, axis=0) / np.sqrt(num_sites_r), decimals=4)
    R_CI_UP = np.around(np.nanmean(RUP, axis=0) / np.sqrt(num_sites_r), decimals=4)
    R_CI = np.array([-R_CI_LO, R_CI_UP])  # shape (2, depth, exp)

    # anomR
    anomR_mean = np.around(np.nanmean(anomR, axis=0), decimals=2)
    anomR_std  = np.around(np.nanstd(anomR, axis=0), decimals=3)
    num_sites_anomr = np.sum(~np.isnan(anomR), axis=0)
    anomR_CI_LO = np.around(np.nanmean(anomRLO, axis=0) / np.sqrt(num_sites_anomr), decimals=4)
    anomR_CI_UP = np.around(np.nanmean(anomRUP, axis=0) / np.sqrt(num_sites_anomr), decimals=4)
    anomR_CI = np.array([-anomR_CI_LO, anomR_CI_UP])

    # Bias
    Bias_mean = np.around(np.nanmean(Bias, axis=0), decimals=3)
    Bias_std  = np.around(np.nanstd(Bias, axis=0), decimals=3)
    num_sites_bias = np.sum(~np.isnan(Bias), axis=0)
    Bias_CI_LO = np.around(np.nanmean(BiasLO, axis=0) / np.sqrt(num_sites_bias), decimals=4)
    Bias_CI_UP = np.around(np.nanmean(BiasUP, axis=0) / np.sqrt(num_sites_bias), decimals=4)
    Bias_CI = np.array([-Bias_CI_LO, Bias_CI_UP])

    # absBias
    absBias_mean = np.around(np.nanmean(absBias, axis=0), decimals=3)
    absBias_std  = np.around(np.nanstd(absBias, axis=0), decimals=3)
    num_sites_abs = np.sum(~np.isnan(absBias), axis=0)
    absBias_CI_LO = np.around(np.nanmean(absBiasLO, axis=0) / np.sqrt(num_sites_abs), decimals=4)
    absBias_CI_UP = np.around(np.nanmean(absBiasUP, axis=0) / np.sqrt(num_sites_abs), decimals=4)
    absBias_CI = np.array([-absBias_CI_LO, absBias_CI_UP])

    # RMSE
    RMSE_mean = np.around(np.nanmean(RMSE, axis=0), decimals=3)
    RMSE_std  = np.around(np.nanstd(RMSE, axis=0), decimals=3)
    num_sites_rmse = np.sum(~np.isnan(RMSE), axis=0)
    RMSE_CI_LO = np.around(np.nanmean(RMSELO, axis=0) / np.sqrt(num_sites_rmse), decimals=4)
    RMSE_CI_UP = np.around(np.nanmean(RMSEUP, axis=0) / np.sqrt(num_sites_rmse), decimals=4)
    RMSE_CI = np.array([-RMSE_CI_LO, RMSE_CI_UP])

    # ubRMSE
    ubRMSE_mean = np.around(np.nanmean(ubRMSE, axis=0), decimals=3)
    ubRMSE_std  = np.around(np.nanstd(ubRMSE, axis=0), decimals=3)
    num_sites_ubrmse = np.sum(~np.isnan(ubRMSE), axis=0)
    ubRMSE_CI_LO = np.around(np.nanmean(ubRMSELO, axis=0) / np.sqrt(num_sites_ubrmse), decimals=4)
    ubRMSE_CI_UP = np.around(np.nanmean(ubRMSEUP, axis=0) / np.sqrt(num_sites_ubrmse), decimals=4)
    ubRMSE_CI = np.array([-ubRMSE_CI_LO, ubRMSE_CI_UP])

    return (R_mean, R_std, R_CI, num_sites_r,
            anomR_mean, anomR_std, anomR_CI, num_sites_anomr,
            Bias_mean, Bias_std, Bias_CI,
            absBias_mean, absBias_std, absBias_CI,
            RMSE_mean, RMSE_std, RMSE_CI,
            ubRMSE_mean, ubRMSE_std, ubRMSE_CI, num_sites_ubrmse)

def plot_bar_grid(output_dir,
                  R_mean, R_CI, num_sites_r,
                  anomR_mean, anomR_CI, num_sites_anomr,
                  ubRMSE_mean, ubRMSE_CI, num_sites_ubrmse,
                  num_expts):
    """Produce the SAME 2x3 bar grid figure and save with identical filename."""
    ind = np.arange(num_expts)
    fig, axs = plt.subplots(2, 3, figsize=(16, 10))

    # Surface R_mean
    axs[0, 0].bar(ind, R_mean[0, :num_expts], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:num_expts])
    axs[0, 0].errorbar(ind, R_mean[0, :num_expts], yerr=R_CI[:, 0, :num_expts], fmt='none', ecolor='grey', capsize=2)
    axs[0, 0].set_ylabel(r'$R$ (-)', fontsize=label_fontsize)
    axs[0, 0].set_ylim(0.5, 0.9)
    axs[0, 0].set_axisbelow(True)
    axs[0, 0].grid(axis='y', color='lightgrey')
    axs[0, 0].set_title(r'Surface $R$ mean (n = {})'.format(num_sites_r[0, 0]), fontsize=title_fontsize)
    axs[0, 0].set_xticks(ind)
    axs[0, 0].set_xticklabels('', fontsize=1)

    # Surface anomR_mean
    axs[0, 1].bar(ind, anomR_mean[0, :num_expts], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:num_expts])
    axs[0, 1].errorbar(ind, anomR_mean[0, :num_expts], yerr=anomR_CI[:, 0, :num_expts], fmt='none', ecolor='grey', capsize=2)
    axs[0, 1].set_ylabel('anomR (-)', fontsize=label_fontsize)
    axs[0, 1].set_ylim(0.5, 0.9)
    axs[0, 1].set_axisbelow(True)
    axs[0, 1].grid(axis='y', color='lightgrey')
    axs[0, 1].set_title(r'Surface anomR mean (n = {})'.format(num_sites_anomr[0, 0]), fontsize=title_fontsize)
    axs[0, 1].set_xticks(ind)
    axs[0, 1].set_xticklabels('', fontsize=1)

    # Surface ubRMSE_mean
    axs[0, 2].bar(ind, ubRMSE_mean[0, :num_expts], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:num_expts])
    axs[0, 2].errorbar(ind, ubRMSE_mean[0, :num_expts], yerr=ubRMSE_CI[:, 0, :num_expts], fmt='none', ecolor='grey', capsize=2)
    axs[0, 2].set_ylabel('ubRMSD ($m^3 \\, m^{-3}$)', fontsize=label_fontsize)
    axs[0, 2].set_ylim(0.02, 0.06)
    axs[0, 2].set_axisbelow(True)
    axs[0, 2].grid(axis='y', color='lightgrey')
    axs[0, 2].set_title(r'Surface ubRMSD mean (n = {})'.format(num_sites_ubrmse[0, 0]), fontsize=title_fontsize)
    axs[0, 2].set_xticks(ind)
    axs[0, 2].set_xticklabels('', fontsize=1)

    # Rootzone R_mean
    axs[1, 0].bar(ind, R_mean[1, :num_expts], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:num_expts])
    axs[1, 0].errorbar(ind, R_mean[1, :num_expts], yerr=R_CI[:, 1, :num_expts], fmt='none', ecolor='grey', capsize=2)
    axs[1, 0].set_ylabel(r'$R$ (-)', fontsize=label_fontsize)
    axs[1, 0].set_ylim(0.5, 0.9)
    axs[1, 0].set_axisbelow(True)
    axs[1, 0].grid(axis='y', color='lightgrey')
    axs[1, 0].set_title(r'Rootzone $R$ mean (n = {})'.format(num_sites_r[1, 0]), fontsize=title_fontsize)
    axs[1, 0].set_xticks(ind)
    axs[1, 0].set_xticklabels(expt_labels[:num_expts], rotation=35, fontsize=label_fontsize)

    # Rootzone anomR_mean
    axs[1, 1].bar(ind, anomR_mean[1, :num_expts], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:num_expts])
    axs[1, 1].errorbar(ind, anomR_mean[1, :num_expts], yerr=anomR_CI[:, 1, :num_expts], fmt='none', ecolor='grey', capsize=2)
    axs[1, 1].set_ylabel('anomR (-)', fontsize=label_fontsize)
    axs[1, 1].set_ylim(0.5, 0.9)
    axs[1, 1].set_axisbelow(True)
    axs[1, 1].grid(axis='y', color='lightgrey')
    axs[1, 1].set_title(r'Rootzone anomR mean (n = {})'.format(num_sites_anomr[1, 0]), fontsize=title_fontsize)
    axs[1, 1].set_xticks(ind)
    axs[1, 1].set_xticklabels(expt_labels[:num_expts], rotation=35, fontsize=label_fontsize)

    # Rootzone ubRMSE_mean
    axs[1, 2].bar(ind, ubRMSE_mean[1, :num_expts], color=plt.rcParams['axes.prop_cycle'].by_key()['color'][:num_expts])
    axs[1, 2].errorbar(ind, ubRMSE_mean[1, :num_expts], yerr=ubRMSE_CI[:, 1, :num_expts], fmt='none', ecolor='grey', capsize=2)
    axs[1, 2].set_ylabel('ubRMSD ($m^3 \\, m^{-3}$)', fontsize=label_fontsize)
    axs[1, 2].set_ylim(0.02, 0.06)
    axs[1, 2].set_axisbelow(True)
    axs[1, 2].grid(axis='y', color='lightgrey')
    axs[1, 2].set_title(r'Rootzone ubRMSD mean (n = {})'.format(num_sites_ubrmse[1, 0]), fontsize=title_fontsize)
    axs[1, 2].set_xticks(ind)
    axs[1, 2].set_xticklabels(expt_labels[:num_expts], rotation=35, fontsize=label_fontsize)

    plt.tight_layout()
    # Save with the SAME filename inside the network-specific folder
    plt.savefig(os.path.join(output_dir, expt_name + '_surf_rz_stats.png'))
    plt.show()

def plot_maps(output_dir, insitu_lat, insitu_lon, R, anomR, ubRMSE, insitu_tag, save=False, outfile=None):
    """
    Make a single 1×3 panel of scatter maps:
      [ΔR, ΔanomR, ΔubRMSD] (surface; depth index 0).
    Extent:
      - SMAP Core ('_CalVal_M33_SM_1d__25yr'): lon [-180, 180], lat [-50, 50]
      - Else (SCAN/USCRN): CONUS
    Colormaps:
      - ΔR, ΔanomR: 'coolwarm_r' (blue=neg, red=pos)
      - ΔubRMSD:    'coolwarm'   (reversed vs above, as in your original)
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    # Extent by tag
    if insitu_tag == '_CalVal_M33_SM_1d__25yr':
        extent = [-180, 180, -50, 50]
    else:
        extent = [-125, -66.5, 24, 49]  # CONUS

    depth_idx = 0  # surface

    # Compute deltas (DA − OL)
    dR      = (R[:, :, 1]      - R[:, :, 0])[:, depth_idx]
    danomR  = (anomR[:, :, 1]  - anomR[:, :, 0])[:, depth_idx]
    dubRMSE = (ubRMSE[:, :, 1] - ubRMSE[:, :, 0])[:, depth_idx]

    # Color levels
    r_levels = np.linspace(-0.06, 0.06, 10)  # fixed
    # symmetric around 0 for anomR; fallback if all NaN/zeros
    maxabs_anom = np.nanmax(np.abs(danomR)) if np.any(np.isfinite(danomR)) else 0.0
    a_levels = np.linspace(-0.01, 0.01, 10) if (maxabs_anom == 0) else np.linspace(-maxabs_anom, maxabs_anom, 10)
    # symmetric for ubRMSE; fallback
    maxabs_u = np.nanmax(np.abs(dubRMSE)) if np.any(np.isfinite(dubRMSE)) else 0.0
    u_levels = np.linspace(-0.01, 0.01, 10) if (maxabs_u == 0) else np.linspace(-maxabs_u, maxabs_u, 10)

    # Colormaps (keep your original choices)
    cmap_R     = plt.cm.get_cmap('coolwarm_r', len(r_levels) - 1)
    cmap_anomR = plt.cm.get_cmap('coolwarm_r', len(a_levels) - 1)
    cmap_u     = plt.cm.get_cmap('coolwarm',    len(u_levels) - 1)  # reversed vs above

    # Build panel
    fig, axes = plt.subplots(
        1, 3, figsize=(18, 6),
        subplot_kw={'projection': ccrs.PlateCarree()}
    )

    panels = [
        ("ΔR (DA − OL)",           dR,      r_levels, cmap_R,     [f"{x:.2f}" for x in r_levels]),
        ("ΔanomR (DA − OL)",       danomR,  a_levels, cmap_anomR, [f"{x:.2f}" for x in a_levels]),
        ("ΔubRMSD (DA − OL) (m³ m⁻³)", dubRMSE, u_levels, cmap_u,     [f"{x:.3f}" for x in u_levels]),
    ]

    for ax, (title, data, levels, cmap, ticklabels) in zip(axes, panels):
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.STATES, edgecolor='black')
        ax.set_extent(extent, crs=ccrs.PlateCarree())

        sc = ax.scatter(
            insitu_lon, insitu_lat,
            c=data, cmap=cmap, s=50, edgecolor='k',
            transform=ccrs.PlateCarree(),
            vmin=levels[0], vmax=levels[-1]
        )

        cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05, shrink=0.8)
        cbar.set_label(title, fontsize=12)
        cbar.set_ticks(levels)
        cbar.ax.set_xticklabels(ticklabels)

        ax.set_xlabel('Longitude', fontsize=12)
        ax.set_ylabel('Latitude', fontsize=12)

    fig.suptitle('Surface-layer deltas (DA − OL)', fontsize=14)
    fig.tight_layout()

    if save:
        # If you want to save, provide outfile or default into the output_dir
        if outfile is None:
            outfile = os.path.join(output_dir, f"maps_panel_surface_{insitu_tag.strip('_')}.png")
        fig.savefig(outfile, dpi=150, bbox_inches='tight')

    plt.show()

def run_for_tag(insitu_tag):
    """One complete run for a given in-situ tag, reproducing your outputs."""
    # Per-tag output directory so filenames can remain identical
    out_dir = os.path.join('outputs_refactor_exact', insitu_tag.strip('_'))
    os.makedirs(out_dir, exist_ok=True)

    # Build stats file list (OL, DA)
    matlab_files = [STATS_TEMPLATE.format(exp=name, tag=insitu_tag) for name in experiment_names]

    # Read first file to get shape
    first = sio.loadmat(matlab_files[0])
    shape = first['Bias'].shape  # (sites, depth)

    num_exp = len(matlab_files)
    # Pre-allocate exactly the same arrays
    Bias      = np.zeros(shape + (num_exp,))
    BiasLO    = np.zeros(shape + (num_exp,))
    BiasUP    = np.zeros(shape + (num_exp,))
    RMSE      = np.zeros(shape + (num_exp,))
    RMSELO    = np.zeros(shape + (num_exp,))
    RMSEUP    = np.zeros(shape + (num_exp,))
    R         = np.zeros(shape + (num_exp,))
    RLO       = np.zeros(shape + (num_exp,))
    RUP       = np.zeros(shape + (num_exp,))
    absBias   = np.zeros(shape + (num_exp,))
    absBiasLO = np.zeros(shape + (num_exp,))
    absBiasUP = np.zeros(shape + (num_exp,))
    anomR     = np.zeros(shape + (num_exp,))
    anomRLO   = np.zeros(shape + (num_exp,))
    anomRUP   = np.zeros(shape + (num_exp,))
    ubRMSE    = np.zeros(shape + (num_exp,))
    ubRMSELO  = np.zeros(shape + (num_exp,))
    ubRMSEUP  = np.zeros(shape + (num_exp,))

    # Fill from MAT files (identical variables)
    for i, f in enumerate(matlab_files):
        mc = sio.loadmat(f)
        Bias[:, :, i]      = mc['Bias']
        BiasLO[:, :, i]    = mc['BiasLO']
        BiasUP[:, :, i]    = mc['BiasUP']
        RMSE[:, :, i]      = mc['RMSE']
        RMSELO[:, :, i]    = mc['RMSELO']
        RMSEUP[:, :, i]    = mc['RMSEUP']
        R[:, :, i]         = mc['R']
        RLO[:, :, i]       = mc['RLO']
        RUP[:, :, i]       = mc['RUP']
        absBias[:, :, i]   = mc['absBias']
        absBiasLO[:, :, i] = mc['absBiasLO']
        absBiasUP[:, :, i] = mc['absBiasUP']
        anomR[:, :, i]     = mc['anomR']
        anomRLO[:, :, i]   = mc['anomRLO']
        anomRUP[:, :, i]   = mc['anomRUP']
        ubRMSE[:, :, i]    = mc['ubRMSE']
        ubRMSELO[:, :, i]  = mc['ubRMSELO']
        ubRMSEUP[:, :, i]  = mc['ubRMSEUP']

    # Compute stats and CIs (identical logic)
    (R_mean, R_std, R_CI, num_sites_r,
     anomR_mean, anomR_std, anomR_CI, num_sites_anomr,
     Bias_mean, Bias_std, Bias_CI,
     absBias_mean, absBias_std, absBias_CI,
     RMSE_mean, RMSE_std, RMSE_CI,
     ubRMSE_mean, ubRMSE_std, ubRMSE_CI, num_sites_ubrmse) = compute_means_cis(
        Bias, BiasLO, BiasUP, RMSE, RMSELO, RMSEUP, R, RLO, RUP,
        absBias, absBiasLO, absBiasUP, anomR, anomRLO, anomRUP, ubRMSE, ubRMSELO, ubRMSEUP
    )

    # Bar-grid figure (saved with same filename inside per-tag folder)
    plot_bar_grid(out_dir,
                  R_mean, R_CI, num_sites_r,
                  anomR_mean, anomR_CI, num_sites_anomr,
                  ubRMSE_mean, ubRMSE_CI, num_sites_ubrmse,
                  num_expts=ubRMSE.shape[2])

    # Raw timeseries file (OL version as in your code)
    m_rs_file = RAW_TS_TEMPLATE.format(exp=experiment_names[0], tag=insitu_tag)
    mat_ts = sio.loadmat(m_rs_file, squeeze_me=True, struct_as_record=False)

    # INSITU_lat / INSITU_lon (flattened as in your plotting)
    INSITU_lat = mat_ts['INSITU_lat'].flatten()
    INSITU_lon = mat_ts['INSITU_lon'].flatten()

    # Reproduce three maps (shown; not saved in original)
    plot_maps(out_dir, INSITU_lat, INSITU_lon, R, anomR, ubRMSE, insitu_tag, save=True)

    dR_surf      = (R[:, :, 1]      - R[:, :, 0])[:, 0]
    dAn_surf     = (anomR[:, :, 1]  - anomR[:, :, 0])[:, 0]
    dUb_surf     = (ubRMSE[:, :, 1] - ubRMSE[:, :, 0])[:, 0]

    dR_rz        = (R[:, :, 1]      - R[:, :, 0])[:, 1]
    dAn_rz       = (anomR[:, :, 1]  - anomR[:, :, 0])[:, 1]
    dUb_rz       = (ubRMSE[:, :, 1] - ubRMSE[:, :, 0])[:, 1]

    # Clean NaNs/infs
    dR_surf_c  = _finite_1d(dR_surf)
    dAn_surf_c = _finite_1d(dAn_surf)
    dUb_surf_c = _finite_1d(dUb_surf)

    dR_rz_c    = _finite_1d(dR_rz)
    dAn_rz_c   = _finite_1d(dAn_rz)
    dUb_rz_c   = _finite_1d(dUb_rz)

    print(f"[{insitu_tag}] valid counts — Surface ΔR:{dR_surf_c.size} ΔanomR:{dAn_surf_c.size} ΔubRMSD:{dUb_surf_c.size} | "
          f"RZ ΔR:{dR_rz_c.size} ΔanomR:{dAn_rz_c.size} ΔubRMSD:{dUb_rz_c.size}")

    return {
        "tag": insitu_tag,
        "surface":  {"delta_R": dR_surf_c,  "delta_anomR": dAn_surf_c,  "delta_ubRMSE": dUb_surf_c},
        "rootzone": {"delta_R": dR_rz_c,    "delta_anomR": dAn_rz_c,    "delta_ubRMSE": dUb_rz_c},
    }

def _zero_ref(ax, lw=1.5, ls='--', color='k'):
    # ensure 0 is inside the y-limits, then draw the line on top
    ymin, ymax = ax.get_ylim()
    if not (ymin <= 0 <= ymax):
        ymin = min(ymin, 0)
        ymax = max(ymax, 0)
        ax.set_ylim(ymin, ymax)
    ax.axhline(0, color=color, lw=lw, ls=ls, zorder=10, alpha=0.7)        

def plot_box_deltas_surface_and_rootzone(results, outfile="boxplot_deltas_surface_rootzone.png"):
    import matplotlib.pyplot as plt
    import numpy as np

    pretty = {
        '_CalVal_M33_SM_1d__25yr': 'SMAP Core',
        '_SCAN_SM_1d_c1234smv_25yr': 'SCAN',
        '_USCRN_SM_1d_c1234smv_25yr': 'USCRN',
    }

    # Keep network order as provided
    labels = [pretty.get(r["tag"], r["tag"]) for r in results]

    # Helper to prep arrays (ensure something is plotted even if empty)
    def _prep(arr):
        arr = np.asarray(arr).astype(float).ravel()
        return arr if arr.size and np.any(np.isfinite(arr)) else np.array([np.nan])

    # Assemble data in the same order for both depths
    surf_dR      = [_prep(r["surface"]["delta_R"])      for r in results]
    surf_dAn     = [_prep(r["surface"]["delta_anomR"])  for r in results]
    surf_dUb     = [_prep(r["surface"]["delta_ubRMSE"]) for r in results]

    rz_dR        = [_prep(r["rootzone"]["delta_R"])      for r in results]
    rz_dAn       = [_prep(r["rootzone"]["delta_anomR"])  for r in results]
    rz_dUb       = [_prep(r["rootzone"]["delta_ubRMSE"]) for r in results]

    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)

    # Row 1: Surface
    axes[0,0].boxplot(surf_dR,  labels=labels, whis=1.5); axes[0,0].set_title(r'Surface  $\Delta R$ (DA − OL)');        axes[0,0].set_ylabel('ΔR (−)'); _zero_ref(axes[0,0]); axes[0,0].grid(axis='y', color='lightgrey')
    axes[0,1].boxplot(surf_dAn, labels=labels, whis=1.5); axes[0,1].set_title(r'Surface  $\Delta$ anomR (DA − OL)');    axes[0,1].set_ylabel('ΔanomR (−)'); _zero_ref(axes[0,1]); axes[0,1].grid(axis='y', color='lightgrey')
    axes[0,2].boxplot(surf_dUb, labels=labels, whis=1.5); axes[0,2].set_title(r'Surface  $\Delta$ ubRMSD (DA − OL)');   axes[0,2].set_ylabel(r'ΔubRMSD ($m^3\,m^{-3}$)'); _zero_ref(axes[0,2]); axes[0,2].grid(axis='y', color='lightgrey')

    # Row 2: Root-zone
    axes[1,0].boxplot(rz_dR,  labels=labels, whis=1.5);  axes[1,0].set_title(r'Root-zone  $\Delta R$ (DA − OL)');        axes[1,0].set_ylabel('ΔR (−)'); _zero_ref(axes[1,0]); axes[1,0].grid(axis='y', color='lightgrey')
    axes[1,1].boxplot(rz_dAn, labels=labels, whis=1.5);  axes[1,1].set_title(r'Root-zone  $\Delta$ anomR (DA − OL)');    axes[1,1].set_ylabel('ΔanomR (−)'); _zero_ref(axes[1,1]); axes[1,1].grid(axis='y', color='lightgrey')
    axes[1,2].boxplot(rz_dUb, labels=labels, whis=1.5);  axes[1,2].set_title(r'Root-zone  $\Delta$ ubRMSD (DA − OL)');   axes[1,2].set_ylabel(r'ΔubRMSD ($m^3\,m^{-3}$)'); _zero_ref(axes[1,2]); axes[1,2].grid(axis='y', color='lightgrey')
    axes[1,2].boxplot(rz_dUb, labels=labels, whis=1.5); axes[1,2].set_title(r'Root-zone  $\Delta$ ubRMSD (DA − OL)');   axes[1,2].set_ylabel(r'ΔubRMSD ($m^3\,m^{-3}$)'); axes[1,2].grid(axis='y', color='lightgrey')

    plt.savefig(outfile, dpi=150)
    plt.show()

if __name__ == "__main__":
    results_all = []  # <<< initialize here
    for tag in INSITU_TAGS:
        print(f"=== Running for tag: {tag} ===")
        res = run_for_tag(tag)
        # Guard in case a run returns None
        if res is not None:
            results_all.append(res)

    print("All done.")

    if results_all:  # only plot if we collected anything
        plot_box_deltas_surface_and_rootzone(
            results_all,
            outfile="boxplot_deltas_all_networks_surface_rootzone.png"
        )
    else:
        print("No results collected; skipping boxplot.")
