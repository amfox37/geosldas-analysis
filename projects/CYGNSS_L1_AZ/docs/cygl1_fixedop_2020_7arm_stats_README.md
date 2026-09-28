# CYGNSS L1 operator test: fixed-operator 2020 arms, O-F stats bundle

*Packaged 2026-09-27. Arizona limited domain (lon −118..−106, lat 29..40), EASEv2 M36 (909 land tiles), 2020-01-01 → 2020-12-31,
GEOSldas fixed-operator build (GEOSldas.x md5 6a70b496). All L1 numbers in this bundle are post-fix and valid; everything packaged
before 2026-09-25 is pre-fix.*

Results and interpretation: `cygl1_fixedop_2020_arm_comparison_report.md` (in this bundle; the source copy is in the geosldas-analysis
repo, commit 4540a50).

## Contents

```
README.md
cygl1_fixedop_2020_arm_comparison_report.md    full write-up: summary, month x arm tables, noise-vs-gain split, next steps
species_index.csv                              species axis (13) -> species id, descr, varname
tile_coords.csv                                tile axis (909) -> tile_id, center lon/lat
stats/                                         28 files = 7 arms x {DA, OL cross-mask} x {spatial pkl, temporal nc4}
  spatial_stats_DA_paired_<tag>_202001_202012.pkl
  spatial_stats_OL_paired_monitor_xmask_<tag>_202001_202012.pkl
  temporal_stats_DA_paired_<tag>_20200101_20201231.nc4
  temporal_stats_OL_paired_monitor_xmask_<tag>_20200101_20201231.nc4
summary/
  month_arm_omf_stdv_pct_vs_OL.{csv,md}          species-group-mean % (DA-OL)/OL of O-F stdv, per month + full year
  noise_gain_<tag>_202001_202012.csv             noise-vs-gain split of the MSE change, per month + pooled year, per monitor group
```

## Arms (`<tag>`)

| tag | GEOSldas EXP_ID | what is assimilated |
|---|---|---|
| `full_xc015_fixedop` | DA_L1_full_xc015_fixedop | full CYGNSS L1 stream, xcorr/ycorr 0.15°, errstd 2.75 dB (**benchmark**) |
| `full_xc015_err39_fixedop` | DA_L1_full_xc015_err39_fixedop | benchmark, errstd 3.9 |
| `full_xc015_coh040216_fixedop` | DA_L1_full_xc015_coh040216_fixedop | benchmark + coherency_ratio filter 0.40–2.16 + L1 z-score clim rebuilt on the screened obs |
| `full_xc015_coh040216_err39_fixedop` | DA_L1_full_xc015_coh040216_err39_fixedop | filter + rebuilt clim + errstd 3.9 (**best L1 arm**) |
| `dense075_coh05_fixedop` | DA_L1_dense075_coh05_fixedop | thinned L1 (0.75° min separation, coherency ≥ 0.5), xcorr 0.625°, errstd 2.75 |
| `smap_fixedop` | DA_SMAP_fixedop | SMAP L1C Tb only (L1 monitor-only) |
| `l3_fixedop` | DA_L3_fixedop | CYGNSS L3 SM (`CYGNSS_SM_6hr`) only (L1 monitor-only) |

The other 12 species are monitor-only in every arm (all 13 in the SMAP- and L3-only arms, apart from the one assimilated), and all obs
are z-score scaled.

## DA vs OL cross-mask: which file to compare against which

There is a single open loop, `OLv8_M36_AZ_fixedop` (unscaled, full L1 stream). For each arm the toolkit
(`postproc_ObsFcstAna`) produced two sets of files:

- `DA_paired_<tag>`: the DA arm's own O, F, A.
- `OL_paired_monitor_xmask_<tag>`: the OL's F (and A = F) scored against **that arm's** scaled obs O, on exactly the
  same obs population (`use_obs=True` cross-mask).

**Always compare `DA_paired_<tag>` with `OL_paired_monitor_xmask_<tag>` of the same tag.** Each OL file is specific to its arm:
the arms have different obs populations (the coherency filter drops about 20% of L1, and the thinned arm keeps about 18%) and
different scaling clims. Report the change as `% (DA − OL)/OL`, the same sign convention for every species (negative = DA better).

## File formats

The species axis (13) and the tile axis (909) are the same in every file and every arm (verified). See `species_index.csv`:

| axis index | species |
|---|---|
| 0–3 | SMOS Tbh_A, Tbh_D, Tbv_A, Tbv_D (K) |
| 4–7 | SMAP L1C Tbh_A, Tbh_D, Tbv_A, Tbv_D (K) |
| 8–10 | ASCAT H SAF MetOp-A/B/C surface SM degree of saturation |
| 11 | CYGNSS L3 `CYGNSS_SM_6hr` (m3/m3) |
| 12 | CYGNSS L1 `CYGNSS_L1_DDM3X5_CROP_SCALAR` (dB) |

All values are in scaled obs space (the obs were rescaled to the model climatology, so the units are those of the model-equivalent
forecast).

**`spatial_stats_*.pkl`**: a Python dict of domain-wide statistics, one row per month:

- `date_vec`: list of 12 `'yyyymm'`
- `O_mean, O_stdv, F_mean, F_stdv, A_mean, A_stdv, OmF_mean, OmF_stdv, OmA_mean, OmA_stdv, N_data`: ndarray `(12 months, 13 species)`

To get a full-year value, pool the months N-weighted:
`var = Σ N(stdv² + mean²)/ΣN − pooled_mean²`. For example, SMAP species mean for the best arm: −1.51% pooled from this pkl,
vs −1.52% from the toolkit's full-year run in the summary table (rounding).

```python
import pickle, numpy as np
da = pickle.load(open('stats/spatial_stats_DA_paired_full_xc015_coh040216_err39_fixedop_202001_202012.pkl', 'rb'))
ol = pickle.load(open('stats/spatial_stats_OL_paired_monitor_xmask_full_xc015_coh040216_err39_fixedop_202001_202012.pkl', 'rb'))
pct = 100 * (da['OmF_stdv'] - ol['OmF_stdv']) / ol['OmF_stdv']     # (12, 13) monthly % vs OL
print(pct[:, 4:8].mean(axis=1))                                     # SMAP species mean per month
```

**`temporal_stats_*.nc4`**: per-tile statistics over the whole year, dims `(tile=909, species=13)`, fill value −9999:
`O_mean/stdv, F_mean/stdv, A_mean/stdv, OmF_mean/stdv, OmA_mean/stdv, OmF_norm_mean/stdv` (normalized by the ensemble-derived
expected std) and `N_data`. The files contain no lon/lat, so use `tile_coords.csv` (row i = tile axis index i).

## Summary tables

- `month_arm_omf_stdv_pct_vs_OL.*`: columns `arm, period (yyyymm or 2020), SMOS, SMAP, ASCAT, L3, L1, N_L1`. SMOS/SMAP/ASCAT
  are the means of their species' % changes. The `2020` row is one pooled full-year statistic, not the mean of the months. L1 is
  the arm's own CygL1 O-F change on its own obs population.
- `noise_gain_<tag>_*.csv`: columns `month (yyyy-mm or all), grp (SMOS/SMAP/ASCAT/L3/L1), N, dMSE, noise, gain, rms_dF, corr,
  alpha_opt, dMSE_alpha`. It is computed directly from the ens_avg ObsFcstAna files (DA matched to the OL on time/species/tile/lon/lat,
  O from DA): with inn = O − F_OL and dF = F_DA − F_OL,
  - `dMSE = noise + gain`, all as % of the OL MSE
  - `noise = Σ dF²`
  - `gain = −2 Σ dF·inn`
  - `corr = corr(dF, inn)`
  - `alpha_opt = Σ dF·inn / Σ dF²` (the increment scaling that would minimise MSE)

  MSE includes the O-F bias, so it differs from the stdv table.

## Headline (full year 2020, % change of O-F stdv vs OL; negative = better)

| arm | SMOS | SMAP | ASCAT | L3 | L1 own |
|---|---:|---:|---:|---:|---:|
| full_xc015 (benchmark) | +0.12 | +0.17 | +3.23 | −1.65 | −0.36 |
| full_xc015_err39 | −0.42 | −0.49 | +0.79 | −2.29 | −0.34 |
| full_xc015_coh040216 | −1.29 | −1.51 | +1.99 | −3.01 | −0.81 |
| **full_xc015_coh040216_err39** | **−1.35** | **−1.52** | **+0.23** | **−2.99** | −0.71 |
| dense075_coh05 | −0.06 | +0.09 | +0.78 | −0.34 | −0.26 |
| smap (SMAP-only) | −14.41 | −11.81 | +8.29 | +1.80 | −0.67 |
| l3 (L3-only) | −1.67 | −1.75 | −4.89 | −7.61 | −0.51 |

In all L1 arms, May–June Tb gets worse: corr(dF, inn) for SMAP/SMOS drops to 0.02–0.05 and alpha_opt ≈ 0. That is a signal problem,
not an R-amplitude problem (see report §2).

## Provenance

- Scoring: `geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/postproc_drivers/score_cygl1_arm.py` (wraps the NC4-capable
  `postproc_ObsFcstAna.py` in the hsaf_cdr_test `obsfcstana-nc4-postproc` worktree). Month × arm table: `build_month_arm_table.py`
  (same dir). Noise-gain: `scripts/cygl1_noise_gain_split.py`.
- Source files on Discover: `geosldas-analysis/projects/CYGNSS_L1_AZ/output/{postproc_paired_density/stats_output,
  month_arm_table_2020, noise_gain_split}/`. Single-month pkl/nc4 files and the monthly sums are there too; they are not included
  here.
- Runs: `/discover/nobackup/projects/land_da/cygl1_operator_test/<EXP_ID>/`. Every run reached cap_restart 20210101 with no errors.
