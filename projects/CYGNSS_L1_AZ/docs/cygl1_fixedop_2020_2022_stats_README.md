# CYGNSS L1 operator test: fixed-operator arms 2020–2022, O-F stats bundle

*Packaged 2026-09-28. Arizona limited domain (lon −118..−106, lat 29..40), EASEv2 M36 (909 land tiles), GEOSldas fixed-operator build
(GEOSldas.x md5 6a70b496). Four arms run 2020-01-01 → 2022-12-31; three more stop at 2020-12-31. All L1 numbers are post-fix and valid.
This bundle supersedes `cygl1_fixedop_2020_7arm_stats.tar.gz`: it contains all of that bundle's arms and adds 2021–2022.*

Results and interpretation: `cygl1_fixedop_2020_arm_comparison_report.md` (in this bundle; §3 covers 2021–2022). The source copy is in
the geosldas-analysis repo, commit 8c95c73.

## Contents

```
README.md
cygl1_fixedop_2020_arm_comparison_report.md    write-up: 2020 (§1-2), 2021-2022 extension (§3), next steps (§4)
species_index.csv                              species axis (13) -> species id, descr, varname
tile_coords.csv                                tile axis (909) -> tile_id, center lon/lat
stats/                                         WHOLE PERIOD 2020-01..2022-12, 4 arms x {DA, OL cross-mask} x {pkl, nc4} = 16 files
  spatial_stats_DA_paired_<tag>_202001_202212.pkl                    36 monthly rows
  spatial_stats_OL_paired_monitor_xmask_<tag>_202001_202212.pkl
  temporal_stats_DA_paired_<tag>_20200101_20221231.nc4                per-tile stats over all 36 months
  temporal_stats_OL_paired_monitor_xmask_<tag>_20200101_20221231.nc4
stats/per_year/                                per-tile nc4 for each calendar year (2020, 2021, 2022), same 4 arms = 24 files
  temporal_stats_{DA_paired,OL_paired_monitor_xmask}_<tag>_<yyyy>0101_<yyyy>1231.nc4
stats/2020_only_arms/                          3 arms that stop at 2020 (same files as the 2020 bundle) = 12 files
  spatial_stats_*_<tag>_202001_202012.pkl, temporal_stats_*_<tag>_20200101_20201231.nc4
summary/
  month_arm_omf_stdv_pct_vs_OL_2020.{csv,md}       2020, all 7 arms: species-group-mean % (DA-OL)/OL of O-F stdv, per month + year
  month_arm_omf_stdv_pct_vs_OL_2021_2022.{csv,md}  2021-2022, 4 arms: per month, per year (2021, 2022) and pooled 2021-22 (ALL)
  noise_gain_<tag>_202001_202212.csv               noise-vs-gain split, 4 arms, per month + pooled 36 months ("all")
  noise_gain_<tag>_202001_202012.csv               the same for the 3 2020-only arms
```

Per-year pkl files are not included because the whole-period pkl already has one row per month; rows 0–11 / 12–23 / 24–35 are
identical to the per-year pkl files (verified). To get a year or any other sub-period, pool those rows (formula below).

## Arms (`<tag>`)

| tag | GEOSldas EXP_ID | period | what is assimilated |
|---|---|---|---|
| `full_xc015_coh040216_fixedop` | DA_L1_full_xc015_coh040216_fixedop | 2020–2022 | full CYGNSS L1 stream, xcorr 0.15°, coherency_ratio filter 0.40–2.16, L1 z-score clim rebuilt on the screened obs, errstd 2.75 dB |
| `full_xc015_coh040216_err39_fixedop` | DA_L1_full_xc015_coh040216_err39_fixedop | 2020–2022 | the same with errstd 3.9 (**best L1 arm**) |
| `l3_fixedop` | DA_L3_fixedop | 2020–2022 | CYGNSS L3 SM (`CYGNSS_SM_6hr`) only (L1 monitor-only) |
| `smap_fixedop` | DA_SMAP_fixedop | 2020–2022 | SMAP L1C Tb only (L1 monitor-only) |
| `full_xc015_fixedop` | DA_L1_full_xc015_fixedop | 2020 | full L1 stream, xcorr 0.15°, errstd 2.75 dB (2020 benchmark) |
| `full_xc015_err39_fixedop` | DA_L1_full_xc015_err39_fixedop | 2020 | benchmark, errstd 3.9 |
| `dense075_coh05_fixedop` | DA_L1_dense075_coh05_fixedop | 2020 | thinned L1 (0.75° min separation, coherency ≥ 0.5), xcorr 0.625°, errstd 2.75 |

The 2021–2022 segments are straight continuations of the 2020 runs: same build, same exeinp, CAP.rc END_DATE moved to 20230101. The
other 12 species are monitor-only in every arm, and all obs are z-score scaled.

## DA vs OL cross-mask: which file to compare against which

There is a single open loop, `OLv8_M36_AZ_fixedop` (unscaled, full L1 stream). For each arm the toolkit (`postproc_ObsFcstAna`)
produced two sets of files:

- `DA_paired_<tag>`: the DA arm's own O, F, A.
- `OL_paired_monitor_xmask_<tag>`: the OL's F (and A = F) scored against **that arm's** scaled obs O, on exactly the same obs
  population (`use_obs=True` cross-mask).

**Always compare `DA_paired_<tag>` with `OL_paired_monitor_xmask_<tag>` of the same tag and the same period.** Each OL file is
specific to its arm: the arms have different obs populations and scaling clims. Report the change as `% (DA − OL)/OL` for every
species (negative = DA better).

## File formats

The species axis (13) and tile axis (909) are the same in every file. See `species_index.csv`:

| axis index | species |
|---|---|
| 0–3 | SMOS Tbh_A, Tbh_D, Tbv_A, Tbv_D (K) |
| 4–7 | SMAP L1C Tbh_A, Tbh_D, Tbv_A, Tbv_D (K) |
| 8–10 | ASCAT H SAF MetOp-A/B/C surface SM degree of saturation |
| 11 | CYGNSS L3 `CYGNSS_SM_6hr` (m3/m3) |
| 12 | CYGNSS L1 `CYGNSS_L1_DDM3X5_CROP_SCALAR` (dB) |

All values are in scaled obs space (obs rescaled to the model climatology).

**`spatial_stats_*.pkl`**: a Python dict of domain-wide statistics, one row per month:

- `date_vec`: list of `'yyyymm'` (36 for the whole-period files, 12 for the 2020-only arms)
- `O_mean, O_stdv, F_mean, F_stdv, A_mean, A_stdv, OmF_mean, OmF_stdv, OmA_mean, OmA_stdv, N_data`: ndarray `(n_months, 13)`

To get a value for a year or any sub-period, pool the months N-weighted:
`var = Σ N(stdv² + mean²)/ΣN − pooled_mean²`. Then average the per-species % changes within a sensor group.

```python
import pickle, numpy as np
t = 'full_xc015_coh040216_err39_fixedop'
da = pickle.load(open(f'stats/spatial_stats_DA_paired_{t}_202001_202212.pkl', 'rb'))
ol = pickle.load(open(f'stats/spatial_stats_OL_paired_monitor_xmask_{t}_202001_202212.pkl', 'rb'))

def pooled_stdv(d, rows, sp):
    N, m, s = (np.asarray(d[k])[rows, sp] for k in ('N_data', 'OmF_mean', 'OmF_stdv'))
    ok = (N > 0) & ~np.isnan(m); N, m, s = N[ok], m[ok], s[ok]
    pm = (N * m).sum() / N.sum()
    return np.sqrt((N * (s**2 + m**2)).sum() / N.sum() - pm**2)

rows = slice(12, 24)                                   # 2021; slice(None) = whole period
smap = np.mean([100 * (pooled_stdv(da, rows, i) - pooled_stdv(ol, rows, i)) / pooled_stdv(ol, rows, i) for i in range(4, 8)])
print(f'SMAP 2021: {smap:+.2f}%')                       # -0.70 (summary table: -0.70)
```

**`temporal_stats_*.nc4`**: per-tile statistics over the file's period, dims `(tile=909, species=13)`, fill value −9999:
`O_mean/stdv, F_mean/stdv, A_mean/stdv, OmF_mean/stdv, OmA_mean/stdv, OmF_norm_mean/stdv` (normalized by the ensemble-derived
expected std) and `N_data`. The files contain no lon/lat, so use `tile_coords.csv` (row i = tile axis index i). Per-tile stats
cannot be pooled across years from the whole-period file, which is why `stats/per_year/` is included.

## Summary tables

- `month_arm_omf_stdv_pct_vs_OL_*.csv`: columns `arm, period, L1, N_L1, SMOS, SMAP, ASCAT, L3` (the same in both files).
  `period` is `yyyymm`, a year (`2020`, `2021`, `2022`), or `ALL` (pooled 2021-01..2022-12). Year and ALL rows are one
  pooled statistic, not the mean of the months. SMOS/SMAP/ASCAT are the means of their species' % changes. L1 is the arm's own CygL1
  O-F change on its own obs population.
- `noise_gain_<tag>_*.csv`: columns `month (yyyy-mm or all), grp (SMOS/SMAP/ASCAT/L3/L1), N, dMSE, noise, gain, rms_dF, corr,
  alpha_opt, dMSE_alpha`. It is computed directly from the ens_avg ObsFcstAna files (DA matched to the OL on
  time/species/tile/lon/lat, O from DA). With inn = O − F_OL and dF = F_DA − F_OL:
  - `dMSE = noise + gain`, all as % of the OL MSE
  - `noise = Σ dF²`
  - `gain = −2 Σ dF·inn`
  - `corr = corr(dF, inn)`
  - `alpha_opt = Σ dF·inn / Σ dF²` (the increment scaling that would minimise MSE; > 1 = increments too small, < 1 = too big,
    ≈ 0 = no information)

  MSE includes the O-F bias, so dMSE can have the opposite sign to the stdv tables in a given month.

## Headline (% change of O-F stdv vs OL; negative = better)

| arm | period | SMOS | SMAP | ASCAT | L3 | L1 own |
|---|---|---:|---:|---:|---:|---:|
| full_xc015_coh040216 | 2020 | −1.29 | −1.51 | +1.99 | −3.01 | −0.81 |
| | 2021–22 | −0.77 | −1.02 | +2.10 | −1.15 | −0.71 |
| | 2020–22 | −1.02 | −1.29 | +2.18 | −1.86 | −0.76 |
| **full_xc015_coh040216_err39** | 2020 | −1.35 | −1.52 | +0.23 | −2.99 | −0.71 |
| | 2021–22 | −0.71 | −0.93 | +0.61 | −1.51 | −0.64 |
| | **2020–22** | **−0.97** | **−1.20** | **+0.55** | **−2.08** | −0.68 |
| l3 (L3-only) | 2020 | −1.67 | −1.75 | −4.89 | −7.61 | −0.51 |
| | 2021–22 | −0.87 | −1.20 | −2.57 | −5.09 | −0.55 |
| | 2020–22 | −1.14 | −1.48 | −3.61 | −6.01 | −0.54 |
| smap (SMAP-only) | 2020 | −14.41 | −11.81 | +8.29 | +1.80 | −0.67 |
| | 2021–22 | −13.02 | −10.26 | +12.70 | +3.88 | −0.65 |
| | 2020–22 | −14.12 | −11.48 | +12.59 | +2.94 | −0.67 |

The 2020–22 rows are pooled from the whole-period pkl files, and the other rows come from the toolkit's own period runs. 2020-only arms:
see `summary/month_arm_omf_stdv_pct_vs_OL_2020.md`.

**The spring degradation recurs every year in both L1 arms, with a moving window:** May–Jun 2020, Apr–Jun 2021, Feb–May 2022. It has two
mechanisms:
- 2020 and 2022: signal collapse (SMAP corr(dF, inn) < 0.1, α_opt ≈ 0).
- 2021: overshoot (corr normal, α_opt ≈ 0.45–0.75, noise 3–5× winter).

See report §3.

## Provenance

- Scoring: `geosldas-analysis/projects/CYGNSS_L1_AZ/scripts/postproc_drivers/score_cygl1_arm.py` (wraps the NC4-capable
  `postproc_ObsFcstAna.py` in the hsaf_cdr_test `obsfcstana-nc4-postproc` worktree). Month × arm tables: `build_month_arm_table.py`.
  Noise-gain: `scripts/cygl1_noise_gain_split.py`.
- Jobs: whole-period stats + noise-gain `output/full_period_2020_2022/run_full_period.sh` (58621001); 2021–22 tables
  `output/month_arm_table_2021_2022/score_arm.sh` (58620093–96); 2020 as in the 2020 bundle.
- Source files on Discover: `geosldas-analysis/projects/CYGNSS_L1_AZ/output/{postproc_paired_density/stats_output,
  month_arm_table_2020, month_arm_table_2021_2022, full_period_2020_2022, noise_gain_split}/`.
- Runs: `/discover/nobackup/projects/land_da/cygl1_operator_test/<EXP_ID>/`. The four 2020–2022 arms reached cap_restart 20230101 and
  the 2020-only arms reached 20210101, all with no `LDAS ERROR`/`forrtl`.
