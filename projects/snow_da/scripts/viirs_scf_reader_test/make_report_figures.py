#!/usr/bin/env python3
# figures for docs/viirs_scf_reader_test_report.md
#   fig01: obs SCF on one day from ObsFcstAna: MODIS Aqua, VIIRS JPSS-1 with cee4997 QC, VIIRS JPSS-1 with 5d4d636 QC
#   fig02: forecast RMS reduction vs the open loop, per run and per sensor's obs (parsed from compare_vs_ol.py output)
# usage: make_report_figures.py <exp_root> <compare_vs_ol_10d.txt> <out_dir> [day YYYYMMDD]
#   exp_root holds the M21C_test_CF0360_jpss_* experiments (Discover: /discover/nobackup/projects/land_da/M21C_testing)
# python: Discover GEOSpyD py3.13 (netCDF4, matplotlib, cartopy)
import glob, os, re, sys
import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

exp_root, rms_txt, out_dir = sys.argv[1:4]
day = sys.argv[4] if len(sys.argv) > 4 else '20200102'

# reference palette (dataviz skill references/palette.md): chrome, sequential blue, categorical slots 1-3
SURFACE, INK, INK2, MUTED, GRID, AXIS = '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
SEQ = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5', '#2a78d6',
       '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
CAT = ['#2a78d6', '#eb6834', '#1baf7a']
plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 9, 'text.color': INK, 'axes.labelcolor': INK2,
                     'xtick.color': MUTED, 'ytick.color': INK2, 'axes.edgecolor': AXIS,
                     'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE})


def day_obs(exp, descr):
    # obs SCF, lon, lat of species descr for all cycles labelled with `day`
    files = sorted(glob.glob(os.path.join(exp_root, exp, 'output/*/ana/ens_avg/Y*/M*/*.ldas_ObsFcstAna.%s_*.nc4' % day)))
    lon, lat, obs = [], [], []
    for fn in files:
        with nc.Dataset(fn) as d:
            ids = [i for i, s in zip(d['obsparam_species_id'][:], d['obsparam_descr'][:]) if s == descr]
            if not ids or d.dimensions['n_obs'].size == 0:
                continue
            k = d['species'][:] == ids[0]
            lon.append(d['lon'][:][k]); lat.append(d['lat'][:][k]); obs.append(d['obs'][:][k])
    return np.concatenate(lon), np.concatenate(lat), np.concatenate(obs)


# ---------------------------------------------------------------- fig01
import cartopy.crs as ccrs
import cartopy.feature as cfeature

panels = [('M21C_test_CF0360_jpss_cee4997_A', 'MYD10C1', 'MODIS Aqua (species 52)'),
          ('M21C_test_CF0360_jpss_cee4997_B', 'VJ110C1', 'VIIRS JPSS-1 (59), first commit cee4997: Basic_QA = 0 only'),
          ('M21C_test_CF0360_jpss_5d4d636_B', 'VJ110C1', 'VIIRS JPSS-1 (59), commit 5d4d636: Basic_QA 0, 1, 3')]
cmap = LinearSegmentedColormap.from_list('seq_blue', SEQ)
proj = ccrs.PlateCarree()
fig, axs = plt.subplots(3, 1, figsize=(10, 6.6), subplot_kw={'projection': proj}, constrained_layout=True)
for ax, (exp, descr, title) in zip(axs, panels):
    lon, lat, obs = day_obs(exp, descr)
    ax.set_extent([-170, 180, 20, 72], crs=proj)
    ax.add_feature(cfeature.COASTLINE.with_scale('110m'), linewidth=0.4, edgecolor=AXIS)
    sc = ax.scatter(lon, lat, c=obs, s=0.15, marker='s', linewidths=0, cmap=cmap, vmin=0, vmax=1,
                    transform=proj, rasterized=True)
    n60 = (lat >= 60).sum()
    ax.set_title('%s\n%s obs on %s-%s-%s, of which %s north of 60N (none north of ~62N: polar night)' % (title, format(obs.size, ','), day[:4], day[4:6],
                 day[6:], format(n60, ',')), loc='left', fontsize=9, color=INK)
    ax.spines['geo'].set_edgecolor(AXIS); ax.spines['geo'].set_linewidth(0.6)
cb = fig.colorbar(sc, ax=axs, orientation='vertical', shrink=0.6, pad=0.01)
cb.set_label('obs snow cover fraction', color=INK2)
cb.outline.set_edgecolor(AXIS); cb.ax.tick_params(color=AXIS, labelcolor=INK2)
fig.savefig(os.path.join(out_dir, 'fig01_obs_scf_maps.png'), dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- fig02
# parse compare_vs_ol.py output: per sensor block, rows "run mean_fcst omf_mean omf_std rms rms_vs_ol%"
txt = open(rms_txt).read()
blocks = re.split(r'\n(?=MODIS|VIIRS)', txt)
rms = {}                                              # sensor -> run -> % change vs OL
for b in blocks:
    m = re.match(r'(MODIS Aqua|MODIS Terra|VIIRS JPSS-1|VIIRS S-NPP) \(\w+\): N common = (\d+)', b)
    if not m or int(m.group(2)) < 1000:              # skip Terra (5 common pairs; only A2 forecasts at Terra cycles)
        continue
    rms[m.group(1)] = {r: float(p) for r, p in re.findall(r'^\s+(\w+)\s+[\d.]+\s+[+-][\d.]+\s+[\d.]+\s+[\d.]+\s+([+-][\d.]+)%',
                                                       b, flags=re.M)}
sensors = ['MODIS Aqua', 'VIIRS JPSS-1', 'VIIRS S-NPP']
runs = [('A2', 'A2: MODIS Aqua + Terra'), ('B', 'B: VIIRS JPSS-1'), ('C', 'C: VIIRS S-NPP'), ('A', 'A: MODIS Aqua')]
fig, ax = plt.subplots(figsize=(7.5, 3.6), constrained_layout=True)
h = 0.24
for j, s in enumerate(sensors):
    y = np.arange(len(runs)) + (j - 1) * (h + 0.03)
    vals = [-rms[s][r] for r, _ in runs]
    ax.barh(y, vals, height=h, color=CAT[j], edgecolor=SURFACE, linewidth=1.0, label='scored against %s obs' % s)
    # direct label on the top row's bars
    ax.text(vals[0] + 0.3, y[0], s, va='center', ha='left', fontsize=8, color=INK2)
ax.set_yticks(np.arange(len(runs)))
ax.set_yticklabels([l for _, l in runs])
ax.invert_yaxis()
ax.set_xlim(0, 32)
ax.set_xlabel('forecast RMS reduction vs open loop (%)')
ax.xaxis.grid(True, color=GRID, linewidth=0.6); ax.set_axisbelow(True)
for side in ('top', 'right', 'left'):
    ax.spines[side].set_visible(False)
ax.tick_params(axis='y', length=0)
ax.legend(loc='lower right', frameon=False, fontsize=8, labelcolor=INK2)
ax.set_title('Snow cover forecast error, 2020-01-01 to 01-11 (1.33M common tile-cycles)', loc='left', fontsize=9.5)
fig.savefig(os.path.join(out_dir, 'fig02_forecast_rms_vs_ol.png'), dpi=150)
plt.close(fig)
print('wrote', out_dir)
