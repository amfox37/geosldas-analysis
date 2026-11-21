#!/usr/bin/env python3
# score the forecast SCF of the OL and DA runs against the obs of each sensor (from the OL innovations run)
#  - with FOV=0 the obs operator is the tile SCF, so a run's fcst at (tile, cycle) is its model SCF there,
#    whichever of its own species it is reported under
#  - per sensor, stats use only the (tile, cycle) pairs common to ALL runs, so columns are comparable
# usage: compare_vs_ol.py <label>=<exp_dir> [...]     (first label must be the OL run)
import glob, os, sys
import numpy as np
import netCDF4 as nc

runs = [a.split('=', 1) for a in sys.argv[1:]]
sensors = [('MYD10C1', 'MODIS Aqua'), ('MOD10C1', 'MODIS Terra'), ('VJ110C1', 'VIIRS JPSS-1'), ('VNP10C1', 'VIIRS S-NPP')]


def read_run(exp_dir):
    # returns {cycle: {descr: (tilenum, obs, fcst)}}
    out = {}
    for fn in sorted(glob.glob(os.path.join(exp_dir, 'output/*/ana/ens_avg/Y*/M*/*.ldas_ObsFcstAna.*.nc4'))):
        cyc = fn.split('.')[-2]
        with nc.Dataset(fn) as d:
            if d.dimensions['n_obs'].size == 0:
                continue
            descr = dict(zip(d['obsparam_species_id'][:].tolist(), d['obsparam_descr'][:]))
            sp, t = d['species'][:], d['tilenum'][:]
            o, f = (np.ma.filled(d[k][:], np.nan) for k in ('obs', 'fcst'))
        out[cyc] = {descr[int(s)]: (t[sp == s], o[sp == s], f[sp == s]) for s in np.unique(sp)}
    return out


data = {label: read_run(e) for label, e in runs}
ol_label = runs[0][0]
labels = [label for label, _ in runs]


def model_fcst(label, cyc):
    # tile -> fcst SCF of run `label` at cycle cyc (from any species it reports)
    m = {}
    for t, _, f in data[label].get(cyc, {}).values():
        m.update(zip(t.tolist(), f.tolist()))
    return m


print('O-F of each run\'s forecast SCF vs each sensor\'s obs, on (tile, cycle) pairs common to all runs')
print('runs: ' + ', '.join('%s=%s' % (l, e) for l, e in runs))
for descr, name in sensors:
    obs_all, f_all = [], {l: [] for l in labels}
    for cyc in sorted(data[ol_label]):
        if descr not in data[ol_label][cyc]:
            continue
        t, o, _ = data[ol_label][cyc][descr]
        maps = {l: (model_fcst(l, cyc) if l != ol_label else dict(zip(t.tolist(), data[ol_label][cyc][descr][2].tolist())))
                for l in labels}
        keep = np.array([all(tt in maps[l] for l in labels) for tt in t.tolist()], dtype=bool)
        if not keep.any():
            continue
        obs_all.append(o[keep])
        for l in labels:
            f_all[l].append(np.array([maps[l][tt] for tt in t[keep].tolist()]))
    if not obs_all:
        print('\n%s (%s): no common pairs' % (name, descr)); continue
    o = np.concatenate(obs_all)
    print('\n%s (%s): N common = %d   mean obs = %.4f' % (name, descr, o.size, o.mean()))
    print('    %-4s %9s %9s %9s %9s %10s' % ('run', 'mean fcst', 'O-F mean', 'O-F std', 'RMS', 'RMS vs OL'))
    rms_ol = None
    for l in labels:
        f = np.concatenate(f_all[l]); d = o - f
        rms = np.sqrt((d ** 2).mean())
        rms_ol = rms if l == ol_label else rms_ol
        print('    %-4s %9.4f %+9.4f %9.4f %9.4f %+9.1f%%' % (l, f.mean(), d.mean(), d.std(), rms, 100 * (rms / rms_ol - 1)))
