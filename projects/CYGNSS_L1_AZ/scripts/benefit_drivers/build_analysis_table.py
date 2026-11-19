#!/usr/bin/env python3
"""Join the benefit-drivers tables into ONE tile-level analysis table (output/benefit_drivers/analysis_table.csv).

Outcomes (out_*, positive = L1 DA better than OL) and covariates (cov_*), keyed on tile_axis_index.
Combined score = mean of z(TC dR2_mod) and z(independent-monitor O-F gain), z-scored over the L1 mask.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parents[2] / 'output' / 'benefit_drivers'
D = OUT / 'discover'
rd = lambda p: pd.read_csv(p)

b = rd(OUT / 'benefit_table.csv')
cv = rd(OUT / 'covariate_table.csv')
ig = rd(OUT / 'igbp_landcover_per_tile.csv')
tc = rd(D / 'tc_iv_per_tile.csv')
op = rd(D / 'operator_skill_per_tile.csv')
st = rd(D / 'static_covariates.csv')
vg = rd(D / 'vegetation_seasonal.csv')
cl = rd(D / 'climate_per_tile.csv')
sm = rd(D / 'ol_sm_climate.csv')

key = ['tile_axis_index']
base = b[['tile_axis_index', 'tile_id', 'idx0', 'com_lon', 'com_lat', 'n_l1_obs', 'in_l1_mask']].copy()
for name, d in dict(cv=cv, ig=ig, tc=tc, op=op, st=st, vg=vg, cl=cl, sm=sm).items():   # key sanity
    assert d.tile_axis_index.is_unique and len(d) == 909, name
    if 'com_lon' in d:
        assert np.nanmax(np.abs(d.set_index('tile_axis_index').com_lon - base.set_index('tile_axis_index').com_lon)) < 1e-3, name
    if 'tile_id' in d:
        assert (d.set_index('tile_axis_index').tile_id == base.set_index('tile_axis_index').tile_id).all(), name

def pref(d, p, drop=()):
    skip = set(key) | {'tile_id', 'idx0', 'ig', 'jg', 'com_lon', 'com_lat', 'lon', 'lat', 'n_l1_obs', 'in_l1_mask', 'l1area', 'block'} | set(drop)
    return d[[key[0]] + [c for c in d if c not in skip]].rename(columns=lambda c: c if c == key[0] else f'{p}_{c}')

# outcomes
o = base[key].copy()
o['out_tc_dR2'] = tc.set_index(key[0]).d_tc_R2mod_L1_minus_OL.reindex(o.tile_axis_index).values
o['out_tc_dR2_L3'] = tc.set_index(key[0]).d_tc_R2mod_L3_minus_OL.reindex(o.tile_axis_index).values
o['out_tc_dR2_L1_minus_L3'] = tc.set_index(key[0]).d_tc_R2mod_L1_minus_L3.reindex(o.tile_axis_index).values
o['out_tc_dSig2pct'] = tc.set_index(key[0]).d_tc_pctSigma2inv_L1_minus_OL.reindex(o.tile_axis_index).values
for s in ('smapL3', 'ascatH121'):
    o[f'out_iv_{s}'] = tc.set_index(key[0])[f'd_iv_R_{s}_L1_minus_OL'].reindex(o.tile_axis_index).values
    o[f'out_iv_{s}_L3'] = tc.set_index(key[0])[f'd_iv_R_{s}_L3_minus_OL'].reindex(o.tile_axis_index).values
bb = b.set_index(key[0])
for m in ('smos', 'smap', 'ascat', 'l3obs', 'l1obs'):
    o[f'out_omf_{m}'] = bb[f'omf_{m}_pct_better_L1'].reindex(o.tile_axis_index).values
    o[f'out_omf_{m}_L3'] = bb[f'omf_{m}_pct_better_L3'].reindex(o.tile_axis_index).values
    for y in (2020, 2021, 2022):
        o[f'out_omf_{m}_{y}'] = bb[f'omf_{m}_pct_better_L1_{y}'].reindex(o.tile_axis_index).values
o['out_omf_indep'] = o[['out_omf_smap', 'out_omf_smos']].mean(axis=1, skipna=False)   # independent Tb monitors (SMAP, SMOS)
o['out_omf_indep_L3'] = o[['out_omf_smap_L3', 'out_omf_smos_L3']].mean(axis=1, skipna=False)

mask = base.in_l1_mask.values.astype(bool)
def z(x):
    x = pd.Series(x)
    return (x - x[mask].mean()) / x[mask].std()
o['out_combined'] = pd.concat([z(o.out_tc_dR2), z(o.out_omf_indep)], axis=1).mean(axis=1, skipna=False)
s1, s2 = np.sign(o.out_tc_dR2), np.sign(o.out_omf_indep)
o['out_agree'] = np.where(o.out_tc_dR2.isna() | o.out_omf_indep.isna(), 'na',
                          np.where(s1 > 0, np.where(s2 > 0, 'both_better', 'tc_better_only'),
                                   np.where(s2 > 0, 'omf_better_only', 'both_worse')))

cov = base.copy()
for p, d, drop in [('cv', cv, [c for c in cv if c.startswith('y_')] + ['geo_frac_cell']),
                   ('ig', ig, ['valid_frac', 'pixel_pct_sum_mean', 'igbp_year', 'dominant_name']),
                   ('op', op, []), ('st', st, []), ('vg', vg, []), ('cl', cl, []), ('sm', sm, [])]:
    cov = cov.merge(pref(d, p, drop), on=key, how='left')
cov = cov.merge(o, on=key, how='left')
cov.to_csv(OUT / 'analysis_table.csv', index=False)
print(cov.shape, 'L1-mask', mask.sum())
print(cov.loc[mask, ['out_tc_dR2', 'out_omf_indep', 'out_omf_ascat', 'out_combined']].describe().round(4).to_string())
print(cov.loc[mask, 'out_agree'].value_counts().to_string())
print(cov.loc[mask, ['out_tc_dR2', 'out_omf_indep', 'out_omf_smap', 'out_omf_smos', 'out_iv_smapL3', 'out_iv_ascatH121', 'out_omf_ascat']].corr(method='spearman').round(2).to_string())
