#!/usr/bin/env python3
"""CYGNSS L1 DA update-health diagnostics for one or more GEOSldas experiments over one month.

Per run: wrong-way analyses (sign(A-F) != sign(O-F) for |O-F|>0.5 dB), |O-A|>|O-F| share, median |A-F|/|O-F|,
ens-mean RZEXC/SRFEXC increment p99/p99.9/max (catch_progn_incr HISTORY collection must be enabled),
local obs-error covariance conditioning cond(R) for obs within +-xcompact of each obs (Gaussian R as in
assemble_obs_cov), and q = d^T (R + diag Pf)^-1 d / n (about 1 if R and Pf are consistent with the innovations).

Usage:
  da_l1_update_diagnostics.py --run DA_L1_full_fixedop:0.625 --run DA_L1_full_xc015_fixedop:0.15 --month 202001
"""
import argparse
import glob
import numpy as np
import netCDF4 as nc

p = argparse.ArgumentParser()
p.add_argument('--exp-path', default='/discover/nobackup/projects/land_da/cygl1_operator_test')
p.add_argument('--run', action='append', required=True, help='EXP_ID:xcorr (xcorr in deg, as in the run nml)')
p.add_argument('--month', required=True, help='yyyymm')
p.add_argument('--xcompact', type=float, default=1.25)
p.add_argument('--species', type=int, default=13, help='ObsFcstAna species index of CYGNSS L1')
p.add_argument('--stride', type=int, default=3, help='use every n-th obs as a local-set centre (speed)')
a = p.parse_args()
Y, M = a.month[:4], a.month[4:]

print('run                            N_assim wrongway% |O-A|>|O-F|% med|A-F|/|O-F| RZEXC p99  p99.9    max  n>25mm | SRFEXC p99   max | cond(R) med      p90      max  q_med')
for spec in a.run:
    exp, L = spec.rsplit(':', 1)
    L = float(L)
    omf, amf, cond, q = [], [], [], []
    for f in sorted(glob.glob(f'{a.exp_path}/{exp}/output/*/ana/ens_avg/Y{Y}/M{M}/*ObsFcstAna*.nc4')):
        d = nc.Dataset(f)
        if d.dimensions['n_obs'].size:
            m = (d['species'][:] == a.species) & (d['assim_flag'][:] == 1)
            o, fc, an, v, fv, lon, lat = [np.asarray(d[x][:][m], 'f8') for x in ('obs', 'fcst', 'ana', 'obsvar', 'fcstvar', 'lon', 'lat')]
            omf += list(o - fc)
            amf += list(an - fc)
            for i in range(0, len(lon), a.stride):
                s = np.where((np.abs(lon - lon[i]) < a.xcompact) & (np.abs(lat - lat[i]) < a.xcompact))[0]
                if len(s) < 3:
                    continue
                R = np.exp(-0.5 * ((lon[s, None] - lon[None, s]) ** 2 + (lat[s, None] - lat[None, s]) ** 2) / L ** 2) * np.sqrt(v[s, None] * v[None, s])
                ev = np.linalg.eigvalsh(R)
                cond.append(ev.max() / max(ev.min(), 1e-300))
                dd = (o - fc)[s]
                q.append(dd @ np.linalg.solve(R + np.diag(fv[s]), dd) / len(s))
        d.close()
    omf, amf = np.array(omf), np.array(amf)
    oma = omf - amf
    big = np.abs(omf) > 0.5
    inc = {'RZEXC_INCR': [], 'SRFEXC_INCR': []}
    for f in glob.glob(f'{a.exp_path}/{exp}/output/*/cat/ens_avg/Y{Y}/M{M}/*catch_progn_incr*.nc4'):
        d = nc.Dataset(f)
        for k in inc:
            x = np.abs(d[k][:].filled(0).ravel())
            inc[k].append(x[x > 0])
        d.close()
    r = np.concatenate(inc['RZEXC_INCR']) if inc['RZEXC_INCR'] else np.array([np.nan])
    s_ = np.concatenate(inc['SRFEXC_INCR']) if inc['SRFEXC_INCR'] else np.array([np.nan])
    cond, q = np.array(cond), np.array(q)
    print(f'{exp[:30]:30s} {len(omf):7d}   {100*np.mean(np.sign(amf[big]) != np.sign(omf[big])):6.1f}     {100*np.mean(np.abs(oma) > np.abs(omf)):6.1f}        '
          f'{np.median(np.abs(amf[big]) / np.abs(omf[big])):.3f}    {np.nanpercentile(r, 99):7.2f} {np.nanpercentile(r, 99.9):7.1f} {np.nanmax(r):7.1f} {int((r > 25).sum()):6d} | '
          f'{np.nanpercentile(s_, 99):6.2f} {np.nanmax(s_):6.1f} | {np.median(cond):8.1e} {np.percentile(cond, 90):8.1e} {cond.max():8.1e}  {np.median(q):.2f}')
