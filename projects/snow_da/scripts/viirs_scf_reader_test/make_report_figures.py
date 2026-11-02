#!/usr/bin/env python3
# figures for docs/viirs_scf_reader_test_report.md
#   fig01: obs SCF on one day from ObsFcstAna: MODIS Aqua, VIIRS JPSS-1 with cee4997 QC, VIIRS JPSS-1 with 5d4d636 QC
#   fig02: density of matched VIIRS vs MODIS Aqua obs SCF (same tile and cycle), all cycles, with binned mean
#   fig03: forecast RMS reduction vs the open loop, per run and per sensor's obs (parsed from compare_vs_ol.py output)
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
from matplotlib.colors import LogNorm


def cycle_obs(exp, descr):
    # {cycle: (tilenum, obs)} of species descr over all cycles
    out = {}
    for fn in sorted(glob.glob(os.path.join(exp_root, exp, 'output/*/ana/ens_avg/Y*/M*/*.ldas_ObsFcstAna.*.nc4'))):
        with nc.Dataset(fn) as d:
            ids = [i for i, s in zip(d['obsparam_species_id'][:], d['obsparam_descr'][:]) if s == descr]
            if not ids or d.dimensions['n_obs'].size == 0:
                continue
            k = d['species'][:] == ids[0]
            out[fn.split('.')[-2]] = (d['tilenum'][:][k], np.ma.filled(d['obs'][:][k], np.nan))
    return out


modis = cycle_obs('M21C_test_CF0360_jpss_cee4997_A', 'MYD10C1')
edges = np.linspace(0, 1, 21)
edges[-1] = 1.0001                                    # obs SCF = 1 falls in the last bin
centers = 0.5 * (edges[:-1] + np.minimum(edges[1:], 1))
cmap_d = LinearSegmentedColormap.from_list('seq_blue', SEQ)
fig, axs = plt.subplots(1, 2, figsize=(9, 4.3), constrained_layout=True)
for ax, (exp, descr, name) in zip(axs, [('M21C_test_CF0360_jpss_5d4d636_B', 'VJ110C1', 'VIIRS JPSS-1'),
                                        ('M21C_test_CF0360_jpss_5d4d636_C', 'VNP10C1', 'VIIRS S-NPP')]):
    viirs = cycle_obs(exp, descr)
    om, ov = [], []
    for c in sorted(set(modis) & set(viirs)):
        _, im, iv = np.intersect1d(modis[c][0], viirs[c][0], return_indices=True)
        om.append(modis[c][1][im]); ov.append(viirs[c][1][iv])
    om, ov = np.concatenate(om), np.concatenate(ov)
    H, _, _ = np.histogram2d(om, ov, bins=[edges, edges])
    H = np.ma.masked_equal(H, 0)
    pc = ax.pcolormesh(edges[:-1].tolist() + [1.0], edges[:-1].tolist() + [1.0], H.T, cmap=cmap_d,
                       norm=LogNorm(vmin=1, vmax=H.max()), rasterized=True)
    ax.plot([0, 1], [0, 1], color=INK2, linewidth=1.0, label='1:1')
    # binned means both ways: noise in the binning variable flattens either curve (regression dilution), so a
    # pair of curves straddling 1:1 symmetrically means noise, and a shift of both to one side means bias
    ib = np.digitize(om, edges) - 1
    bmean = np.array([ov[ib == k].mean() if (ib == k).any() else np.nan for k in range(len(centers))])
    ax.plot(centers, bmean, color=CAT[1], linewidth=2, marker='o', markersize=4,
            label='mean %s per MODIS bin' % name.split()[-1])
    jb = np.digitize(ov, edges) - 1
    rmean = np.array([om[jb == k].mean() if (jb == k).any() else np.nan for k in range(len(centers))])
    ax.plot(rmean, centers, color=CAT[2], linewidth=2, marker='s', markersize=4,
            label='mean MODIS per %s bin' % name.split()[-1])
    d = ov - om
    ax.text(0.03, 0.97, '%s pairs\ncorrelation %.3f\nmean diff %+.3f' % (format(om.size, ','), np.corrcoef(om, ov)[0, 1],
            d.mean()), transform=ax.transAxes, va='top', ha='left', fontsize=8.5, color=INK2)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect('equal')
    ax.set_xlabel('MODIS Aqua obs SCF'); ax.set_ylabel('%s obs SCF' % name)
    ax.set_title('%s vs MODIS Aqua, same tile and cycle' % name, loc='left', fontsize=9.5)
    ax.legend(loc='lower right', frameon=False, fontsize=8, labelcolor=INK2)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
cb = fig.colorbar(pc, ax=axs, shrink=0.85, pad=0.02)
cb.set_label('matched tile-cycles per 0.05 x 0.05 bin (log scale)', color=INK2)
cb.outline.set_edgecolor(AXIS); cb.ax.tick_params(color=AXIS, labelcolor=INK2)
fig.savefig(os.path.join(out_dir, 'fig02_viirs_vs_modis_density.png'), dpi=150)
plt.close(fig)

# ---------------------------------------------------------------- fig03
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
fig.savefig(os.path.join(out_dir, 'fig03_forecast_rms_vs_ol.png'), dpi=150)
plt.close(fig)
print('wrote', out_dir)
