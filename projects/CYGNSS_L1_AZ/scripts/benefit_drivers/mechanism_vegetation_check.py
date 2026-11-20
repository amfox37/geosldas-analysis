#!/usr/bin/env python3
"""PRELIMINARY checks on the analysis table (mask tiles).
1. Operator variance budget: SM-term vs geometry-term std, sensitivity slope = sd_smterm / sfmc_std (dB per m3/m3), and their
   rank correlations with SM-term skill and covariates. Model-derived (OL SM run through the operator), not measurements.
2. Vegetation: Spearman and partial Spearman (controls: OL SM std, OL SM mean, elevation) of LAI/opacity/canopy covariates with the
   main outcomes, on all mask tiles and within open shrub (IGBP 7); 2x2 deg block-bootstrap CIs (300 draws) for the partial rho.
Writes output/benefit_drivers/{mechanism_check.csv,vegetation_check.csv}.
"""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata, spearmanr

OUT = Path(__file__).resolve().parents[2] / 'output' / 'benefit_drivers'
d = pd.read_csv(OUT / 'analysis_table.csv')
m = d[d.in_l1_mask.astype(bool)].copy()

# 1. variance budget
m['slope'] = m.op_sd_smterm / m.op_sfmc_std
m['sm_share'] = m.op_sd_smterm ** 2 / (m.op_sd_smterm ** 2 + m.op_sd_geom ** 2)
q = lambda c: m[c].quantile([.25, .5, .75]).round(3).tolist()
rows = [('sd_smterm_dB', *q('op_sd_smterm')), ('sd_geom_dB', *q('op_sd_geom')), ('sens_dB_per_m3m3', *q('slope')), ('sm_share_of_var', *q('sm_share'))]
pd.DataFrame(rows, columns=['quantity', 'q25', 'median', 'q75']).to_csv(OUT / 'mechanism_check.csv', index=False)
def rho(a, b):
    k = m[[a, b]].dropna(); return round(spearmanr(k[a], k[b])[0], 2)
print('variance budget (q25, median, q75):'); [print(' ', r) for r in rows]
for a, b in [('op_rs_smterm', 'op_sfmc_std'), ('op_rs_smterm', 'op_sfmc_mean'), ('op_sd_smterm', 'op_sfmc_std'), ('slope', 'op_sfmc_mean'),
             ('slope', 'vg_lai_clim_mean'), ('slope', 'st_canopy_height_m'), ('slope', 'op_opac_med')]:
    print(f'  rho({a}, {b}) = {rho(a, b)}')

# 2. vegetation
veg = ['vg_lai_clim_mean', 'vg_ndvi_clim_mean', 'vg_vegopacity_clim_mean', 'vg_lai_MAM_mean', 'vg_ndvi_MAM_mean', 'op_opac_med', 'st_canopy_height_m',
       'ig_frac_shrub', 'ig_frac_grass', 'ig_frac_barren']
veg = [v for v in veg if v in m]
out = ['out_omf_indep', 'out_iv_smapL3', 'out_omf_ascat', 'out_iv_ascatH121', 'out_tc_dR2']
ctrl = ['sm_sfmc_std', 'sm_sfmc_mean', 'st_bcs_mean_elevation_m']
m['block'] = np.floor(m.com_lon / 2).astype(int).astype(str) + '_' + np.floor(m.com_lat / 2).astype(int).astype(str)
rng = np.random.default_rng(0)
def resid(y, Z):
    Z = np.column_stack([np.ones(len(y))] + [rankdata(z) for z in Z.T]); ry = rankdata(y)
    return ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
def stat(df, x, y, partial):
    k = df[[x, y] + ctrl].dropna()
    if len(k) < 25: return np.nan, len(k)
    if partial: return np.corrcoef(resid(k[x].values, k[ctrl].values), resid(k[y].values, k[ctrl].values))[0, 1], len(k)
    return spearmanr(k[x], k[y])[0], len(k)
CI = {'vg_lai_clim_mean', 'vg_vegopacity_clim_mean', 'vg_lai_MAM_mean', 'st_canopy_height_m'}
rows = []
for sub, name in [(m, 'all'), (m[m.ig_dominant_class == 7], 'open shrub')]:
    gb = {b: g for b, g in sub.groupby('block')}; blocks = np.array(list(gb))
    for y in out:
        for x in veg:
            r0, n = stat(sub, x, y, False); r1, _ = stat(sub, x, y, True); lo = hi = np.nan
            if x in CI:
                v = [stat(pd.concat([gb[b] for b in rng.choice(blocks, len(blocks))]), x, y, True)[0] for _ in range(300)]
                lo, hi = np.nanpercentile(v, [2.5, 97.5])
            rows.append((name, y, x, n, round(r0, 2), round(r1, 2), round(lo, 2), round(hi, 2)))
r = pd.DataFrame(rows, columns=['subset', 'outcome', 'veg', 'n', 'rho_raw', 'rho_partial', 'lo', 'hi'])
r.to_csv(OUT / 'vegetation_check.csv', index=False)
print(r[r.veg.isin(CI)].to_string(index=False))
