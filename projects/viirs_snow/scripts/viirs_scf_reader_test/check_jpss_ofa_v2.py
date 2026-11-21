#!/usr/bin/env python3
# checks 5-8 of the VIIRS SCF reader rerun: NC4 ObsFcstAna of several runs, side by side
# usage: check_jpss_ofa_v2.py <map_day YYYYMMDD> <out_png> <label>=<exp_dir>:<descr>[,<descr>] [...]
#   first run is the MODIS reference, e.g. A=M21C_test_CF0360_jpss_cee4997_A:MYD10C1
import glob, os, sys
import numpy as np
import netCDF4 as nc

map_day, out_png = sys.argv[1], sys.argv[2]
runs = []                                               # (label, exp_dir, [descr])
for arg in sys.argv[3:]:
    label, rest = arg.split('=', 1)
    exp_dir, descrs = rest.rsplit(':', 1)
    runs.append((label, exp_dir, descrs.split(',')))

bands = [('60-90N', 60, 90), ('40-60N', 40, 60), ('20-40N', 20, 40), ('<20N', -90, 20)]
scf_classes = [('0', 0, 0.01), ('partial', 0.01, 0.99), ('1', 0.99, 2)]


def read_run(exp_dir):
    # returns {descr: {cycle: dict(tilenum, lon, lat, obs, fcst, ana)}}
    files = sorted(glob.glob(os.path.join(exp_dir, 'output/*/ana/ens_avg/Y*/M*/*.ldas_ObsFcstAna.*.nc4')))
    out = {}
    for fn in files:
        cyc = fn.split('.')[-2]
        with nc.Dataset(fn) as d:
            if d.dimensions['n_obs'].size == 0:
                continue
            descr = dict(zip(d['obsparam_species_id'][:].tolist(), d['obsparam_descr'][:]))
            sp = d['species'][:]
            v = {k: np.ma.filled(d[k][:], np.nan) for k in ('tilenum', 'lon', 'lat', 'obs', 'fcst', 'ana')}
        for s in np.unique(sp):
            i = sp == s
            out.setdefault(descr[int(s)], {})[cyc] = {k: a[i] for k, a in v.items()}
    return out, len(files)


data = {}
for label, exp_dir, descrs in runs:
    data[label], nfiles = read_run(exp_dir)
    print('%s: %s  %d ObsFcstAna files, species with obs: %s' % (label, exp_dir, nfiles, sorted(data[label])))
cols = [(label, s) for label, _, descrs in runs for s in descrs]


def cat(label, s, k):
    d = data[label].get(s, {})
    return np.concatenate([x[k] for x in d.values()]) if d else np.array([])


print('\n== check 5: obs counts (total over all cycles) by latitude band')
print('%-22s %10s' % ('run:species', 'total') + ''.join('%10s' % b[0] for b in bands) + '   share 60-90N  40-60N')
for label, s in cols:
    lat = cat(label, s, 'lat')
    n = [((lat >= lo) & (lat < hi)).sum() for _, lo, hi in bands]
    sh = [k / max(lat.size, 1) for k in n[:2]]
    print('%-22s %10d' % ('%s:%s' % (label, s), lat.size) + ''.join('%10d' % k for k in n) +
          '   %11.3f %7.3f' % tuple(sh))

print('\n== check 6: value ranges and means')
for label, s in cols:
    o, f, a = cat(label, s, 'obs'), cat(label, s, 'fcst'), cat(label, s, 'ana')
    if o.size == 0:
        print('%-22s no obs' % ('%s:%s' % (label, s)))
        continue
    print('%-22s N=%8d  obs min/mean/max %.4f/%.4f/%.4f  N obs>1 %6d | fcst %.4f/%.4f/%.4f | ana %.4f/%.4f/%.4f | NaN fcst/ana %d/%d' %
          ('%s:%s' % (label, s), o.size, np.nanmin(o), np.nanmean(o), np.nanmax(o), (o > 1).sum(),
           np.nanmin(f), np.nanmean(f), np.nanmax(f), np.nanmin(a), np.nanmean(a), np.nanmax(a),
           np.isnan(f).sum(), np.isnan(a).sum()))

print('\n== check 7: VIIRS vs MODIS (same tile and cycle), overall and by MODIS SCF class')
ref_label, _, ref_descrs = runs[0]
ref = data[ref_label].get(ref_descrs[0], {})
per_day = []                                            # (name, day, obs VIIRS, obs MODIS) of matched pairs
for label, s in cols[1:]:
    d = data[label].get(s, {})
    if not d:
        continue
    o_v, o_m = [], []
    for c in sorted(set(d) & set(ref)):
        _, iv, im = np.intersect1d(d[c]['tilenum'], ref[c]['tilenum'], return_indices=True)
        o_v.append(d[c]['obs'][iv]); o_m.append(ref[c]['obs'][im])
    days = np.concatenate([np.full(len(x), c[:8]) for c, x in zip(sorted(set(d) & set(ref)), o_v)])
    o_v, o_m = np.concatenate(o_v), np.concatenate(o_m)
    per_day.append(('%s:%s' % (label, s), days, o_v, o_m))
    print('%s:%s vs %s:%s' % (label, s, ref_label, ref_descrs[0]))
    for name, lo, hi in [('all', -1, 9)] + scf_classes:
        i = (o_m >= lo) & (o_m < hi)
        diff = o_v[i] - o_m[i]
        corr = np.corrcoef(o_v[i], o_m[i])[0, 1] if i.sum() > 2 and o_m[i].std() > 0 else np.nan
        flag = '  <-- |bias| > 0.05' if abs(diff.mean()) > 0.05 else ''
        print('    %-8s N pairs = %8d  corr = %7.4f  mean diff = %+.4f  std diff = %.4f%s' %
              (name, i.sum(), corr, diff.mean(), diff.std(), flag))

print('\n== check 7b: VIIRS vs MODIS per day (analysis cycles 03z-00z labelled by cycle date)')
for name, days, o_v, o_m in per_day:
    print(name)
    print('    day        N pairs  corr all  corr partial  bias all  bias 0   bias partial  bias 1')
    for day in np.unique(days):
        j = days == day
        v, m = o_v[j], o_m[j]
        b = [(v[i] - m[i]).mean() if i.any() else np.nan for i in
             [m < 0.01, (m >= 0.01) & (m < 0.99), m >= 0.99]]
        p = (m >= 0.01) & (m < 0.99)
        print('    %s %8d  %8.4f  %12.4f  %+8.4f  %+7.4f  %+12.4f  %+7.4f' %
              (day, j.sum(), np.corrcoef(v, m)[0, 1], np.corrcoef(v[p], m[p])[0, 1], (v - m).mean(), *b))

print('\n== check 8: map of obs SCF on', map_day)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, len(cols), figsize=(6 * len(cols), 4.2), constrained_layout=True)
for ax, (label, s) in zip(np.atleast_1d(axs), cols):
    d = data[label].get(s, {})
    x = [d[c] for c in sorted(d) if c.startswith(map_day)]
    if not x:
        ax.set_title('%s:%s  no obs' % (label, s)); continue
    lon = np.concatenate([a['lon'] for a in x]); lat = np.concatenate([a['lat'] for a in x])
    obs = np.concatenate([a['obs'] for a in x])
    sc = ax.scatter(lon, lat, c=obs, s=0.05, vmin=0, vmax=1, cmap='Blues_r', rasterized=True)
    ax.set_title('%s:%s  obs SCF %s  (N=%d)' % (label, s, map_day, obs.size))
    ax.set_xlim(-180, 180); ax.set_ylim(-60, 90); ax.set_facecolor('0.85')
fig.colorbar(sc, ax=axs, shrink=0.9, label='obs SCF')
fig.savefig(out_png, dpi=130)
print('wrote', out_png)
