"""Generate figures for ASCAT Legacy vs H121 CDR integration test report."""

import numpy as np
import glob
import os
from datetime import datetime

import eccodes as ecc
import netCDF4 as nc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ── Config ───────────────────────────────────────────────────────────────────
BUFR_BASE  = '/Users/amfox/Desktop/ASCAT_SSM_CDR/legacy_bufr'
H121_BASE  = '/Users/amfox/Desktop/ASCAT_SSM_CDR/H121'
OFA_BASE   = '/Users/amfox/Desktop/geosldas-analysis/data/hsaf_cdr_test'
OFA_SUFFIX = 'output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/Y2020/M01'
OUT_DIR    = os.path.dirname(os.path.abspath(__file__))

DATE     = datetime(2020, 1, 1)
DATE_STR = '2020-01-01'
GRID_RES = 0.25

PLATFORMS = {
    'Metop-A': {'bufr_prefix': 'M02-ASCA-ASCSMO02-NA-5.0-', 'subdir': 'metop_a'},
    'Metop-B': {'bufr_prefix': 'M01-ASCA-ASCSMO02-NA-5.0-', 'subdir': 'metop_b'},
    'Metop-C': {'bufr_prefix': 'M03-ASCA-ASCSMO02-NA-5.0-', 'subdir': 'metop_c'},
}

EXPTS = {
    'monitor':        'hsaf_cdr_test_DAv8_M36',
    'scale':          'hsaf_cdr_test_DAv8_M36_scale',
    'eumetsat_assim': 'hsaf_cdr_test_DAv8_M36_EUMETSAT_assim',
    'hsaf_assim':     'hsaf_cdr_test_DAv8_M36_HSAF_assim',
}

SPECIES = {
    9:  {'product': 'Legacy', 'platform': 'Metop-A'},
    10: {'product': 'Legacy', 'platform': 'Metop-B'},
    11: {'product': 'Legacy', 'platform': 'Metop-C'},
    15: {'product': 'H121',   'platform': 'Metop-A'},
    16: {'product': 'H121',   'platform': 'Metop-B'},
    17: {'product': 'H121',   'platform': 'Metop-C'},
}
ASCAT_SPECIES = list(SPECIES.keys())
legacy_ids = [9, 10, 11]
h121_ids   = [15, 16, 17]


# ── Readers ──────────────────────────────────────────────────────────────────
def read_bufr(bufr_dir, date, prefix):
    pattern = os.path.join(bufr_dir, f'{prefix}{date.strftime("%Y%m%d")}*.bfr')
    files = sorted(glob.glob(pattern))
    KEYS = {'ssm': 'surfaceSoilMoisture', 'smpf': 'soilMoistureProcessingFlag',
            'smcf': 'soilMoistureCorrectionFlag', 'tpcx': 'topographicComplexity',
            'iwfr': 'inundationAndWetlandFraction', 'lat': 'latitude', 'lon': 'longitude'}
    out = {k: [] for k in KEYS}
    for f in files:
        with open(f, 'rb') as fh:
            bufr = ecc.codes_bufr_new_from_file(fh)
            while bufr is not None:
                ecc.codes_set(bufr, 'unpack', 1)
                arrays = {k: ecc.codes_get_array(bufr, v) for k, v in KEYS.items()}
                n = len(arrays['lat'])
                if any(len(a) not in (1, n) for a in arrays.values()):
                    ecc.codes_release(bufr)
                    bufr = ecc.codes_bufr_new_from_file(fh)
                    continue
                for k, arr in arrays.items():
                    if len(arr) == 1:
                        arrays[k] = np.repeat(arr, n)
                for k, arr in arrays.items():
                    out[k].extend(arr.tolist())
                ecc.codes_release(bufr)
                bufr = ecc.codes_bufr_new_from_file(fh)
    return {k: np.array(v) for k, v in out.items()}

def qc_bufr(d):
    return ((d['ssm'] >= 0) & (d['ssm'] <= 100) & (d['smpf'] == 0) &
            ((d['smcf'] == 0) | (d['smcf'] == 4)) &
            (d['tpcx'] >= 0) & (d['tpcx'] <= 10) &
            (d['iwfr'] >= 0) & (d['iwfr'] <= 10))

def read_h121(h121_dir, date):
    pattern = os.path.join(h121_dir, f'*_{date.strftime("%Y%m%d")}*.nc')
    files = sorted(glob.glob(pattern))
    out = {k: [] for k in ['lat', 'lon', 'ssm', 'pf', 'cf', 'wf', 'tc']}
    for f in files:
        ds = nc.Dataset(f)
        out['lat'].extend(ds.variables['latitude'][:].filled(np.nan).tolist())
        out['lon'].extend(ds.variables['longitude'][:].filled(np.nan).tolist())
        out['ssm'].extend(ds.variables['surface_soil_moisture'][:].filled(np.nan).tolist())
        out['pf'].extend(ds.variables['processing_flag'][:].filled(255).tolist())
        out['cf'].extend(ds.variables['correction_flag'][:].filled(255).tolist())
        out['wf'].extend(ds.variables['wetland_fraction'][:].filled(-128).tolist())
        out['tc'].extend(ds.variables['topographic_complexity'][:].filled(-128).tolist())
        ds.close()
    return {k: np.array(v) for k, v in out.items()}

def qc_h121(d):
    return ((d['ssm'] >= 0) & (d['ssm'] <= 100) & (d['pf'] == 0) &
            ((d['cf'] == 0) | (d['cf'] == 4)) &
            (d['wf'] >= 0) & (d['wf'] <= 10) &
            (d['tc'] >= 0) & (d['tc'] <= 10))

def bin_to_grid(lat, lon, vals, res=0.25):
    n_lat, n_lon = int(180/res), int(360/res)
    i_lat = np.clip(np.floor((lat + 90) / res).astype(int), 0, n_lat-1)
    i_lon = np.clip(np.floor((lon + 180) / res).astype(int), 0, n_lon-1)
    s = np.full((n_lat, n_lon), 0.)
    c = np.full((n_lat, n_lon), 0)
    np.add.at(s, (i_lat, i_lon), vals)
    np.add.at(c, (i_lat, i_lon), 1)
    with np.errstate(invalid='ignore'):
        mean = np.where(c > 0, s/c, np.nan)
    glat = -90 + (np.arange(n_lat) + 0.5) * res
    glon = -180 + (np.arange(n_lon) + 0.5) * res
    return glat, glon, mean

def ofa_dir(key):
    return f'{OFA_BASE}/{EXPTS[key]}/{OFA_SUFFIX}'

def read_ofa(dir_path):
    files = sorted(glob.glob(os.path.join(dir_path, '*.nc4')))
    data = {s: {k: [] for k in ['lat','lon','obs','fcst','ana','innov','incr','assim']}
            for s in ASCAT_SPECIES}
    for f in files:
        ds = nc.Dataset(f)
        sp = ds.variables['species'][:]
        obs = ds.variables['obs'][:]
        fcs = ds.variables['fcst'][:]
        ana = ds.variables['ana'][:]
        lat = ds.variables['lat'][:]
        lon = ds.variables['lon'][:]
        af  = ds.variables['assim_flag'][:]
        ds.close()
        for s in ASCAT_SPECIES:
            m = sp == s
            if m.sum() == 0: continue
            data[s]['lat'].extend(lat[m].tolist())
            data[s]['lon'].extend(lon[m].tolist())
            data[s]['obs'].extend(obs[m].tolist())
            data[s]['fcst'].extend(fcs[m].tolist())
            data[s]['ana'].extend(ana[m].tolist())
            data[s]['innov'].extend((obs[m]-fcs[m]).tolist())
            data[s]['incr'].extend((ana[m]-fcs[m]).tolist())
            data[s]['assim'].extend(af[m].tolist())
    return {s: {k: np.array(v) for k,v in d.items()} for s,d in data.items()}

def map_ax(ax, data, vmin, vmax, cmap, glon, glat):
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=2)
    GLON, GLAT = np.meshgrid(glon, glat)
    sc = ax.pcolormesh(GLON, GLAT, np.ma.masked_invalid(data),
                       vmin=vmin, vmax=vmax, cmap=cmap,
                       transform=ccrs.PlateCarree(), zorder=1)
    return sc

def savefig(name):
    path = os.path.join(OUT_DIR, name)
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved {name}')


# ── Load data ────────────────────────────────────────────────────────────────
print('Reading BUFR (all platforms)...')
bufr_lats, bufr_lons, bufr_ssms = [], [], []
for plat, cfg in PLATFORMS.items():
    d = read_bufr(f'{BUFR_BASE}/{cfg["subdir"]}/2020/01', DATE, cfg['bufr_prefix'])
    q = qc_bufr(d)
    bufr_lats.append(d['lat'][q]); bufr_lons.append(d['lon'][q]); bufr_ssms.append(d['ssm'][q])
    print(f'  {plat}: {q.sum():,} obs after QC')
bufr_lat = np.concatenate(bufr_lats)
bufr_lon = np.concatenate(bufr_lons)
bufr_ssm = np.concatenate(bufr_ssms)

print('Reading H121 (all platforms)...')
h121_lats, h121_lons, h121_ssms = [], [], []
for plat, cfg in PLATFORMS.items():
    d = read_h121(f'{H121_BASE}/{cfg["subdir"]}/2020/01', DATE)
    q = qc_h121(d)
    h121_lats.append(d['lat'][q]); h121_lons.append(d['lon'][q]); h121_ssms.append(d['ssm'][q])
    print(f'  {plat}: {q.sum():,} obs after QC')
h121_lat = np.concatenate(h121_lats)
h121_lon = np.concatenate(h121_lons)
h121_ssm = np.concatenate(h121_ssms)

print('Reading OFA files...')
ofa_all = {key: read_ofa(ofa_dir(key)) for key in EXPTS}

glat, glon, bufr_grid = bin_to_grid(bufr_lat, bufr_lon, bufr_ssm)
_,    _,    h121_grid = bin_to_grid(h121_lat, h121_lon, h121_ssm)
both = np.isfinite(bufr_grid) & np.isfinite(h121_grid)


# ── Fig 1: SSM distributions ─────────────────────────────────────────────────
print('Fig 1...')
fig, ax = plt.subplots(figsize=(7, 4))
bins = np.arange(0, 101, 2)
ax.hist(bufr_ssm, bins=bins, density=True, alpha=0.65, color='steelblue',
        label=f'Legacy BUFR (n={len(bufr_ssm):,})')
ax.hist(h121_ssm, bins=bins, density=True, alpha=0.65, color='darkorange',
        label=f'H121 CDR (n={len(h121_ssm):,})')
ax.set_xlabel('SSM (% saturation)')
ax.set_ylabel('Density')
ax.set_title(f'SSM distributions after QC — all platforms (Metop-A/B/C), {DATE_STR}')
ax.legend()
plt.tight_layout()
savefig('fig1_ssm_distributions.png')

# ── Fig 2: Gridded scatter BUFR vs H121 ──────────────────────────────────────
print('Fig 2...')
b, h = bufr_grid[both], h121_grid[both]
bias = (h-b).mean(); r = np.corrcoef(b,h)[0,1]; rmsd = np.sqrt(((h-b)**2).mean())
fig, ax = plt.subplots(figsize=(5.5, 5.5))
ax.hexbin(b, h, gridsize=60, cmap='viridis', mincnt=1)
ax.plot([0,100],[0,100],'r--',lw=1,label='1:1')
ax.set_xlim(0,100); ax.set_ylim(0,100)
ax.set_xlabel('Legacy BUFR SSM (% sat.)')
ax.set_ylabel('H121 CDR SSM (% sat.)')
ax.set_title(f'0.25° gridded SSM — all platforms (Metop-A/B/C), {DATE_STR}')
ax.legend(fontsize=9)
ax.text(0.05, 0.95, f'n={both.sum():,}\nbias={bias:+.1f}\nRMSD={rmsd:.1f}\nR={r:.3f}',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='white', alpha=0.8))
plt.tight_layout()
savefig('fig2_gridded_ssm_scatter.png')

# ── Fig 3: Innovation distributions (unscaled, all platforms) ────────────────
print('Fig 3...')
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
bins_i = np.linspace(-0.8, 0.8, 80)
mon = ofa_all['monitor']
for ax, ids, title in [
        (axes[0], [11,17], 'Metop-C only'),
        (axes[1], legacy_ids+h121_ids, 'All platforms')]:
    for id_list, label, color in [(legacy_ids if ids==legacy_ids+h121_ids else [11],
                                    'Legacy', 'steelblue'),
                                   (h121_ids if ids==legacy_ids+h121_ids else [17],
                                    'H121', 'darkorange')]:
        if ids == legacy_ids+h121_ids:
            inn = np.concatenate([mon[s]['innov'] for s in id_list])
        else:
            inn = mon[ids[0] if label=='Legacy' else ids[1]]['innov']
        ax.hist(inn, bins=bins_i, density=True, alpha=0.65, color=color,
                label=f'{label}  bias={inn.mean():+.3f}, std={inn.std():.3f}')
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.set_xlabel('O−F innovation (unscaled, deg. sat.)')
    ax.set_ylabel('Density')
    ax.set_title(f'Innovation distributions — {title}')
    ax.legend(fontsize=8)

# redo properly
plt.close()
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
bins_i = np.linspace(-0.8, 0.8, 80)
for ax, (s_leg, s_h121), title in [
        (axes[0], (11, 17), 'Metop-C only'),
        (axes[1], (None, None), 'All platforms')]:
    for ids, label, color in [(legacy_ids, 'Legacy', 'steelblue'),
                               (h121_ids,   'H121',   'darkorange')]:
        inn = np.concatenate([mon[s]['innov'] for s in ids]) if s_leg is None else \
              mon[11]['innov'] if label=='Legacy' else mon[17]['innov']
        ax.hist(inn, bins=bins_i, density=True, alpha=0.65, color=color,
                label=f'{label}  bias={inn.mean():+.3f}  std={inn.std():.3f}')
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.set_xlabel('O−F innovation (unscaled, deg. sat.)')
    ax.set_ylabel('Density')
    ax.set_title(f'Unscaled innovations — {title}, {DATE_STR}')
    ax.legend(fontsize=8)
plt.tight_layout()
savefig('fig3_innovation_distributions.png')

# ── Fig 4: Effect of scaling ──────────────────────────────────────────────────
print('Fig 4...')
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle(f'Effect of z-score scaling — {DATE_STR}', fontsize=11)
bins_obs = np.linspace(0, 0.6, 80)
bins_inn = np.linspace(-0.4, 0.4, 80)
for row, (ids, label, color) in enumerate([
        (legacy_ids, 'Legacy ASCAT', 'steelblue'),
        (h121_ids,   'H121 CDR',     'darkorange')]):
    obs_raw  = np.concatenate([ofa_all['monitor'][s]['obs']   for s in ids])
    obs_scl  = np.concatenate([ofa_all['scale'][s]['obs']     for s in ids])
    fcst_scl = np.concatenate([ofa_all['scale'][s]['fcst']    for s in ids])
    inn_raw  = np.concatenate([ofa_all['monitor'][s]['innov'] for s in ids])
    inn_scl  = np.concatenate([ofa_all['scale'][s]['innov']   for s in ids])
    ax = axes[row, 0]
    ax.hist(obs_raw, bins=bins_obs, density=True, alpha=0.6, color=color,
            label=f'obs unscaled  mean={obs_raw.mean():.3f}')
    ax.hist(obs_scl, bins=bins_obs, density=True, alpha=0.6, color='gray',
            label=f'obs scaled    mean={obs_scl.mean():.3f}')
    ax.hist(fcst_scl, bins=bins_obs, density=True, alpha=0.35, color='k',
            label=f'fcst          mean={fcst_scl.mean():.3f}')
    ax.set_xlabel('deg. sat. (0-1)'); ax.set_ylabel('Density')
    ax.set_title(f'{label} — obs before/after scaling')
    ax.legend(fontsize=8)
    ax = axes[row, 1]
    ax.hist(inn_raw, bins=bins_inn, density=True, alpha=0.6, color=color,
            label=f'unscaled  bias={inn_raw.mean():+.3f}  std={inn_raw.std():.3f}')
    ax.hist(inn_scl, bins=bins_inn, density=True, alpha=0.6, color='gray',
            label=f'scaled    bias={inn_scl.mean():+.3f}  std={inn_scl.std():.3f}')
    ax.axvline(0, color='k', lw=0.8, ls='--')
    ax.set_xlabel('O−F innovation (deg. sat.)'); ax.set_ylabel('Density')
    ax.set_title(f'{label} — innovations before/after scaling')
    ax.legend(fontsize=8)
plt.tight_layout()
savefig('fig4_scaling_effect.png')

# ── Fig 5: Post-scaling gridded innovation scatter ───────────────────────────
print('Fig 5...')
scl = ofa_all['scale']
_, _, leg_scl_grid  = bin_to_grid(
    np.concatenate([scl[s]['lat'] for s in legacy_ids]),
    np.concatenate([scl[s]['lon'] for s in legacy_ids]),
    np.concatenate([scl[s]['innov'] for s in legacy_ids]))
_, _, h121_scl_grid = bin_to_grid(
    np.concatenate([scl[s]['lat'] for s in h121_ids]),
    np.concatenate([scl[s]['lon'] for s in h121_ids]),
    np.concatenate([scl[s]['innov'] for s in h121_ids]))
both_scl = np.isfinite(leg_scl_grid) & np.isfinite(h121_scl_grid)
li_s, hi_s = leg_scl_grid[both_scl], h121_scl_grid[both_scl]
bias_s = (hi_s-li_s).mean(); rmsd_s = np.sqrt(((hi_s-li_s)**2).mean())
r_s = np.corrcoef(li_s, hi_s)[0,1]
fig, ax = plt.subplots(figsize=(5.5, 5.5))
lim = 0.15
ax.hexbin(li_s, hi_s, gridsize=60, cmap='viridis', mincnt=1, extent=[-lim,lim,-lim,lim])
ax.plot([-lim,lim],[-lim,lim],'r--',lw=1,label='1:1')
ax.axhline(0,color='k',lw=0.5,ls=':'); ax.axvline(0,color='k',lw=0.5,ls=':')
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel('Legacy O−F (scaled, deg. sat.)')
ax.set_ylabel('H121 O−F (scaled, deg. sat.)')
ax.set_title(f'Post-scaling innovations — 0.25° grid, {DATE_STR}')
ax.legend(fontsize=9)
ax.text(0.05,0.95, f'n={both_scl.sum():,}\nbias={bias_s:+.4f}\nRMSD={rmsd_s:.4f}\nR={r_s:.4f}',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round',fc='white',alpha=0.8))
plt.tight_layout()
savefig('fig5_postscaling_scatter.png')

# ── Fig 6: Analysis increment maps ───────────────────────────────────────────
print('Fig 6...')
def grid_incr(expt_key, ids):
    d_all = ofa_all[expt_key]
    lats, lons, incrs = [], [], []
    for s in ids:
        d = d_all[s]; mask = d['assim'] == 1
        if mask.sum() == 0: continue
        lats.append(d['lat'][mask]); lons.append(d['lon'][mask])
        incrs.append(d['incr'][mask])
    if not lats:
        return np.full((int(180/0.25), int(360/0.25)), np.nan)
    _, _, g = bin_to_grid(np.concatenate(lats), np.concatenate(lons), np.concatenate(incrs))
    return g

incr_eu   = grid_incr('eumetsat_assim', legacy_ids)
incr_hsaf = grid_incr('hsaf_assim',     h121_ids)
both_i    = np.isfinite(incr_eu) & np.isfinite(incr_hsaf)
diff_i    = np.where(both_i, incr_hsaf - incr_eu, np.nan)

fig, axes = plt.subplots(3, 1, figsize=(14, 11),
                          subplot_kw={'projection': ccrs.PlateCarree()})
fig.suptitle(f'0.25° gridded analysis increments (A−F) — {DATE_STR}', fontsize=11)
plots = [(incr_eu,   'Legacy ASCAT (EUMETSAT_assim)', -0.05, 0.05, 'RdBu_r'),
         (incr_hsaf, 'H121 CDR (HSAF_assim)',          -0.05, 0.05, 'RdBu_r'),
         (diff_i,    'H121 − Legacy increment',         -0.03, 0.03, 'PuOr_r')]
for ax, (data, title, vmin, vmax, cmap) in zip(axes, plots):
    sc = map_ax(ax, data, vmin, vmax, cmap, glon, glat)
    ax.set_title(title, fontsize=10)
    plt.colorbar(sc, ax=ax, orientation='vertical', fraction=0.02, pad=0.02,
                 label='A−F (deg. sat.)')
plt.tight_layout()
savefig('fig6_analysis_increments.png')

# ── Fig 7: Increment scatter ──────────────────────────────────────────────────
print('Fig 7...')
li_i = incr_eu[both_i]; hi_i = incr_hsaf[both_i]
diff_iv = hi_i - li_i
bias_i = diff_iv.mean(); rmsd_i = np.sqrt((diff_iv**2).mean()); r_i = np.corrcoef(li_i,hi_i)[0,1]
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
ax = axes[0]
lim = 0.06
ax.hexbin(li_i, hi_i, gridsize=60, cmap='viridis', mincnt=1, extent=[-lim,lim,-lim,lim])
ax.plot([-lim,lim],[-lim,lim],'r--',lw=1,label='1:1')
ax.axhline(0,color='k',lw=0.5,ls=':'); ax.axvline(0,color='k',lw=0.5,ls=':')
ax.set_xlim(-lim,lim); ax.set_ylim(-lim,lim)
ax.set_xlabel('Legacy A−F increment (deg. sat.)')
ax.set_ylabel('H121 A−F increment (deg. sat.)')
ax.set_title(f'Gridded increments — {DATE_STR}')
ax.legend(fontsize=9)
ax.text(0.05,0.95, f'n={both_i.sum():,}\nbias={bias_i:+.5f}\nRMSD={rmsd_i:.5f}\nR={r_i:.4f}',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round',fc='white',alpha=0.8))
ax = axes[1]
bins_inc = np.linspace(-0.15, 0.15, 80)
leg_inc  = np.concatenate([ofa_all['eumetsat_assim'][s]['incr']
                            [ofa_all['eumetsat_assim'][s]['assim']==1] for s in legacy_ids])
h121_inc = np.concatenate([ofa_all['hsaf_assim'][s]['incr']
                            [ofa_all['hsaf_assim'][s]['assim']==1]     for s in h121_ids])
ax.hist(leg_inc,  bins=bins_inc, density=True, alpha=0.65, color='steelblue',
        label=f'Legacy  mean={leg_inc.mean():+.4f}  std={leg_inc.std():.4f}')
ax.hist(h121_inc, bins=bins_inc, density=True, alpha=0.65, color='darkorange',
        label=f'H121    mean={h121_inc.mean():+.4f}  std={h121_inc.std():.4f}')
ax.axvline(0,color='k',lw=0.8,ls='--')
ax.set_xlabel('A−F increment (deg. sat., assimilated obs only)')
ax.set_ylabel('Density')
ax.set_title('Analysis increment distributions')
ax.legend(fontsize=8)
plt.tight_layout()
savefig('fig7_increment_scatter.png')

print('All figures saved.')
