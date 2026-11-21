#!/usr/bin/env python3
# compare MODIS Aqua+Terra DA (A2) with VIIRS-only DA (B, C) and MODIS Aqua-only DA (A)
#  1. per run and species: N obs, mean obs/fcst/ana, O-F mean/std, mean increment (ana-fcst) in obs space
#  2. cross-check on (tile, cycle) pairs that have both MODIS Aqua (A2 sp. 52) and VIIRS obs (B or C):
#     forecast of each run vs MODIS Aqua obs and vs VIIRS obs
#  3. model snow state in the end restarts (total SWE, snow-covered tiles) vs the initial restart
# usage: compare_a2.py <exp_A> <exp_A2> <exp_B> <exp_C> <end YYYYMMDD>
import glob, os, sys
import numpy as np
import netCDF4 as nc

exps = dict(zip(['A', 'A2', 'B', 'C'], sys.argv[1:5]))
end = sys.argv[5]


def read_run(exp_dir):
    files = sorted(glob.glob(os.path.join(exp_dir, 'output/*/ana/ens_avg/Y*/M*/*.ldas_ObsFcstAna.*.nc4')))
    out = {}
    for fn in files:
        cyc = fn.split('.')[-2]
        with nc.Dataset(fn) as d:
            if d.dimensions['n_obs'].size == 0:
                continue
            descr = dict(zip(d['obsparam_species_id'][:].tolist(), d['obsparam_descr'][:]))
            sp = d['species'][:]
            v = {k: np.ma.filled(d[k][:], np.nan) for k in ('tilenum', 'lat', 'obs', 'fcst', 'ana')}
        for s in np.unique(sp):
            i = sp == s
            out.setdefault(descr[int(s)], {})[cyc] = {k: a[i] for k, a in v.items()}
    return out


data = {r: read_run(e) for r, e in exps.items()}

print('== 1. obs-space summary per run and species (all cycles)')
print('%-4s %-8s %9s %7s %7s %7s %8s %7s %8s' % ('run', 'species', 'N', 'obs', 'fcst', 'ana', 'O-F mean', 'O-F std', 'ana-fcst'))
for r in ['A', 'A2', 'B', 'C']:
    for s in sorted(data[r]):
        c = {k: np.concatenate([x[k] for x in data[r][s].values()]) for k in ('obs', 'fcst', 'ana')}
        omf = c['obs'] - c['fcst']
        print('%-4s %-8s %9d %7.4f %7.4f %7.4f %+8.4f %7.4f %+8.4f' %
              (r, s, c['obs'].size, c['obs'].mean(), c['fcst'].mean(), c['ana'].mean(), omf.mean(), omf.std(),
               (c['ana'] - c['fcst']).mean()))

print('\n== 2. cross-check on common (tile, cycle): forecast of each run vs MODIS Aqua obs and vs VIIRS obs')
aqua = data['A2']['MYD10C1']
for rv, sv in (('B', 'VJ110C1'), ('C', 'VNP10C1')):
    vi = data[rv][sv]
    o_m, o_v, f_a2, f_v = [], [], [], []
    for c in sorted(set(aqua) & set(vi)):
        _, ia, iv = np.intersect1d(aqua[c]['tilenum'], vi[c]['tilenum'], return_indices=True)
        o_m.append(aqua[c]['obs'][ia]); f_a2.append(aqua[c]['fcst'][ia])
        o_v.append(vi[c]['obs'][iv]); f_v.append(vi[c]['fcst'][iv])
    o_m, o_v, f_a2, f_v = (np.concatenate(x) for x in (o_m, o_v, f_a2, f_v))
    print('%s (%s DA) vs A2 (MODIS Aqua+Terra DA): N common = %d' % (rv, sv, o_m.size))
    for oname, o in (('MODIS Aqua obs', o_m), ('%s obs' % sv, o_v)):
        for fname, f in (('A2 fcst', f_a2), ('%s fcst' % rv, f_v)):
            d = o - f
            print('    %-15s - %-8s  mean %+.4f  std %.4f  RMS %.4f  corr %.4f' %
                  (oname, fname, d.mean(), d.std(), np.sqrt((d ** 2).mean()), np.corrcoef(o, f)[0, 1]))
    print('    fcst A2 - %s: mean %+.4f  mean |diff| %.4f' % (rv, (f_a2 - f_v).mean(), np.abs(f_a2 - f_v).mean()))


def swe(exp_dir, date):
    f = glob.glob(os.path.join(exp_dir, 'output/*/rs/ens0000/Y%s/M%s/*.catch_internal_rst.%s_0000' %
                               (date[:4], date[4:6], date)))[0]
    with nc.Dataset(f) as d:
        return sum(np.asarray(d['WESNN%d' % i][:], dtype=np.float64) for i in (1, 2, 3))


print('\n== 3. model snow state (catch_internal_rst): total SWE = sum WESNN1-3 [kg m-2]')
w0 = swe(exps['A'], '20200101')
print('initial 20200101: tiles with SWE>0 = %d  mean SWE = %.3f' % ((w0 > 0).sum(), w0.mean()))
w_end = {r: swe(e, end) for r, e in exps.items()}
for r in ['A', 'A2', 'B', 'C']:
    w = w_end[r]
    dA2 = w - w_end['A2']
    print('%-3s %s: tiles SWE>0 = %7d  mean SWE = %.3f  | vs A2: tiles differ (>0.1 kg m-2) = %7d  mean diff = %+.4f  '
          'tiles snow in run but not A2 = %6d  in A2 but not run = %6d' %
          (r, end, (w > 0).sum(), w.mean(), (np.abs(dA2) > 0.1).sum(), dA2.mean(),
           ((w > 0) & (w_end['A2'] == 0)).sum(), ((w == 0) & (w_end['A2'] > 0)).sum()))
