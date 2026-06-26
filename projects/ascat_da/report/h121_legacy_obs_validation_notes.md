# ASCAT Legacy vs H121 observation validation notes

Last updated: 2026-06-17

## Purpose

We are checking whether raw ASCAT observations can be processed locally into
GEOSldas-like super-observations before comparing the legacy EUMETSAT BUFR ASCAT
stream with the H SAF H121 CDR stream. The intent is to prove the mechanics
first, then use the notebook to investigate product differences.

Main notebook:

- `projects/ascat_da/notebooks/legacy_vs_h121_obs.ipynb`

## Local sample data

Observation sample root:

- `/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample`

Raw observations:

- Legacy BUFR:
  - `legacy_bufr/metop_a/Y2020/M06`
  - `legacy_bufr/metop_b/Y2020/M06`
  - `legacy_bufr/metop_c/Y2020/M06`
- H121:
  - `H121/metop_a/Y2020/M06`
  - `H121/metop_b/Y2020/M06`
  - `H121/metop_c/Y2020/M06`

GEOSldas innovation/OFA files:

- `/Users/amfox/Desktop/ASCAT_SSM_CDR/hsaf_cdr_test_DAv8_M36_202006_innov/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/Y2020/M06`

GEOS tile metadata copied from the experiment `rc_out`:

- `/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample/tilecoord/hsaf_cdr_test_DAv8_M36_202006_innov.ldas_tilecoord.bin`
- `/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample/tilecoord/hsaf_cdr_test_DAv8_M36_202006_innov.ldas_tilegrids.bin`

## Code touched/added

ASCAT helpers:

- `projects/ascat_da/lib/readers.py`
- `projects/ascat_da/lib/superob.py`

Diagnostics:

- `projects/ascat_da/scripts/run_legacy_vs_h121_smoke.py`
- `projects/ascat_da/scripts/check_raw_superobs_vs_ofa.py`
- `projects/ascat_da/scripts/build_global_superob_cache.py`
- `projects/ascat_da/scripts/check_global_superobs_vs_ofa.py`

The notebook now starts with an OFA validation table before moving into
Legacy-vs-H121 differences.

## GEOSldas conventions copied from source

GEOSldas source checked:

- `/Users/amfox/Desktop/GEOSldas_develop/GEOSldas/src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90`

Important points:

- File discovery is wider than the assimilation window so files that overlap the
  window are not missed.
- The actual obs selection happens later, at the per-observation level.
- ASCAT obs are assigned to model tiles with `get_tile_num_for_obs`, which calls
  `get_tile_num_from_latlon`.
- Super-obs are arithmetic means of all accepted obs in the same tile for the
  same species/cycle.
- H121/HSAF and legacy BUFR use slightly different interval-boundary conventions:
  - Legacy BUFR: keep approximately `[date_time - dt/2, date_time + dt/2)`
  - H121/HSAF: keep approximately `(date_time - dt/2, date_time + dt/2]`
- The analysis cycles are centered on GEOSldas `date_time`, i.e. for 3-hourly
  analysis:
  - `0000z`, `0300z`, `0600z`, ..., `2100z`
  - not simple midnight-anchored `00-03`, `03-06`, etc. bins.

The Python readers now return both:

- `window`: cycle label modulo 8 (`0=0000z`, `1=0300z`, ..., `7=2100z`)
- `cycle`: day-aware cycle index so next-day `0000z` is not mixed with same-day
  `0000z`

For H121/HSAF QC, the GEOSldas reader version mirrored by the
`geos_cycle_global_v2_geos_hsaf_qc` cache applied:

- reject if `surface_flag` open-water bit `0x01` is set,
- reject if `processing_flag` bad bits `0x01|0x02` are set,
- reject if `wetland_fraction >= 10%`,
- reject if `topographic_complexity >= 10%`,
- reject if `subsurface_scattering_probability >= 10%`,
- do **not** screen `correction_flag`.

This matters. The first Python cache screened `correction_flag` and did not
mirror the GEOS HSAF QC exactly. Rebuilding with the looser GEOS-style H121 QC
substantially improved the global H121/OFA match fraction and removed most of
the global mean bias.

Update, 2026-06-26: the active H121/H139 development branch
`GEOSldas_GridComp:feature/amfox/ascat-hsaf-v8` now tightens H121 QC beyond
this v2 cache baseline: `subsurface_scattering_probability >= 5%`,
`surface_soil_moisture_sensitivity <= 1 dB`, and `backscatter40_flag` bit 4
are rejected. Those changes match the current Python `QC_DEFAULT_H121`; v2
cache results should be read as the pre-tightening GEOS-style baseline.

## OFA species mapping

Legacy BUFR ASCAT:

- Metop-A: species 9, `ASCAT_META_SM`
- Metop-B: species 10, `ASCAT_METB_SM`
- Metop-C: species 11, `ASCAT_METC_SM`

H121/HSAF ASCAT:

- Metop-A: species 14, `ASCAT_HSAF_META_SM`
- Metop-B: species 15, `ASCAT_HSAF_METB_SM`
- Metop-C: species 16, `ASCAT_HSAF_METC_SM`

## Derived caches

Regional notebook cache:

- `projects/ascat_da/.cache/legacy_vs_h121_obs/`
- Current notebook cache tag: `geos_cycle_v4_geos_hsaf_qc`
- This is zoom-domain only (`35--38 N`, `-99---96 E`) and is rebuilt by the
  notebook when needed.

Global super-ob cache:

- `projects/ascat_da/.cache/global_superobs/`
- Current global cache tag: `geos_cycle_global_v2_geos_hsaf_qc`
- One pickle per day:
  - `ascat_global_superobs_YYYYMMDD_geos_cycle_global_v2_geos_hsaf_qc.pkl`
- Summary:
  - `ascat_global_superobs_summary_geos_cycle_global_v2_geos_hsaf_qc.csv`

Each cache row is one global GEOS tile/cycle super-ob:

- `date`
- `product`
- `platform`
- `tilenum`
- `cycle`
- `window`
- `i_indg`, `j_indg`
- `lat`, `lon`
- `ssm_pct`
- `n_obs`
- `ssm_std_pct`
- `ssm_min_pct`
- `ssm_max_pct`

Build command:

```bash
/Users/amfox/mamba/envs/regrid/bin/python \
  projects/ascat_da/scripts/build_global_superob_cache.py \
  --start-date 2020-06-01 --end-date 2020-06-10 \
  --version geos_cycle_global_v2_geos_hsaf_qc
```

The v2 global cache took about 12 minutes to build and is about 230 MB for
10 daily pickle files.

Ten-day global cache totals:

| product/platform | raw obs read | global tile/cycle superobs | obs contributing to superobs |
|---|---:|---:|---:|
| H121 Metop-A | 6,562,210 | 1,041,946 | 6,536,095 |
| H121 Metop-B | 6,613,565 | 1,048,003 | 6,587,268 |
| H121 Metop-C | 6,541,687 | 1,035,975 | 6,515,690 |
| Legacy Metop-A | 665,522 | 436,408 | 659,089 |
| Legacy Metop-B | 597,234 | 394,943 | 591,109 |
| Legacy Metop-C | 656,960 | 431,721 | 650,691 |

## Current OFA validation result

Global validation command:

```bash
/Users/amfox/mamba/envs/regrid/bin/python \
  projects/ascat_da/scripts/check_global_superobs_vs_ofa.py \
  --start-date 2020-06-01 --end-date 2020-06-10 \
  --version geos_cycle_global_v2_geos_hsaf_qc
```

Outputs:

- `global_superobs_vs_ofa_summary_20200601_20200610_geos_cycle_global_v2_geos_hsaf_qc.csv`
- `global_superobs_vs_ofa_pairs_20200601_20200610_geos_cycle_global_v2_geos_hsaf_qc.csv`
- `global_superobs_vs_ofa_daily_20200601_20200610_geos_cycle_global_v2_geos_hsaf_qc.csv`

Ten-day global aggregate:

| product | matched OFA tile-cycles | OFA total | match fraction | weighted bias | weighted RMSE | weighted MAE | weighted median abs |
|---|---:|---:|---:|---:|---:|---:|---:|
| H121 | 2,206,338 | 2,423,366 | 0.910 | -0.108 | 7.786 | 5.313 | 3.616 |
| Legacy | 548,376 | 785,345 | 0.698 | -0.006 | 6.301 | 4.067 | 2.600 |

Platform breakdown:

| product/platform | matched OFA tile-cycles | OFA total | match fraction | weighted bias | weighted RMSE | weighted MAE | weighted median abs |
|---|---:|---:|---:|---:|---:|---:|---:|
| H121 Metop-A | 727,379 | 798,974 | 0.910 | -0.071 | 7.790 | 5.328 | 3.633 |
| H121 Metop-B | 742,828 | 815,687 | 0.911 | -0.118 | 7.768 | 5.291 | 3.598 |
| H121 Metop-C | 736,131 | 808,705 | 0.910 | -0.133 | 7.798 | 5.321 | 3.616 |
| Legacy Metop-A | 186,201 | 268,060 | 0.695 | -0.050 | 6.313 | 4.069 | 2.600 |
| Legacy Metop-B | 170,053 | 245,879 | 0.692 | 0.054 | 6.313 | 4.076 | 2.600 |
| Legacy Metop-C | 192,122 | 271,406 | 0.708 | -0.020 | 6.277 | 4.057 | 2.600 |

Interpretation:

- H121 global match fraction improved from about 0.73 to about 0.91 when the
  Python cache was rebuilt with GEOS-style HSAF QC.
- H121 global mean bias relative to OFA is now close to zero for all platforms.
- H121 RMSE remains around 7.8% saturation, so the match is not exact. Remaining
  differences likely involve model-based QC, edge cases, or exact file/window
  handling.
- Legacy global match fraction remains around 0.70. This is expected to be
  lower because the legacy GEOS reader also uses an external ASCAT mask and
  additional model-based QC.

## Legacy vs H121 global super-ob comparison

The notebook now compares Legacy and H121 directly in global GEOS tile/cycle
space over all 10 days. Matching is by:

- analysis date,
- platform,
- `tilenum`,
- GEOS cycle.

Current 10-day global matched counts:

| platform | matched Legacy-H121 tile/cycles | bias H121-Legacy | RMSE | MAE | median abs | r | Legacy median n_obs | H121 median n_obs |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Metop-A | 344,783 | -8.066 | 19.658 | 15.385 | 12.191 | 0.795 | 1 | 8 |
| Metop-B | 316,573 | -5.621 | 18.316 | 14.344 | 11.414 | 0.806 | 1 | 8 |
| Metop-C | 345,224 | -7.725 | 19.445 | 15.185 | 11.940 | 0.794 | 1 | 8 |
| All | 1,006,580 | -7.180 | 19.172 | 14.989 | 11.840 | 0.797 | 1 | 8 |

Interpretation:

- Legacy and H121 are correlated but not interchangeable in tile/cycle space.
- H121 is substantially denser in this sample: median `n_obs` is about 8 versus
  1 for Legacy.
- The negative H121-Legacy bias and large RMSE should be interpreted alongside
  the sampling/super-ob structure difference, not just as a retrieval-value
  difference.

## Current notebook status

The notebook has been updated to:

- use the local Discover-style sample paths,
- read GEOS `tilecoord`/`tilegrids`,
- assign raw observations to GEOS `tilenum` and centered `cycle`,
- validate raw super-obs against ObsFcstAna globally before product comparison,
- compare Legacy and H121 globally on matched `date + platform + tilenum + cycle`
  pairs,
- make zoom-domain spatial plots only for date/cycle windows with enough
  observations instead of plotting empty cycles.

Spatial-map selection:

- The notebook scans all 10 days and all 8 GEOS cycles per day.
- It plots a zoom-map figure only when at least one platform/product has
  `MIN_SPATIAL_MAP_OBS` raw observations in the zoom box.
- Current default: `MIN_SPATIAL_MAP_OBS = 250`.
- Optional draft cap: `MAX_SPATIAL_MAP_FIGURES = None`.

Notebook execution with the `regrid` kernel succeeds via:

```bash
MPLCONFIGDIR=/private/tmp/matplotlib XDG_CACHE_HOME=/private/tmp \
/Users/amfox/mamba/envs/regrid/bin/python -m jupyter nbconvert \
  --to notebook --execute projects/ascat_da/notebooks/legacy_vs_h121_obs.ipynb \
  --output /private/tmp/legacy_vs_h121_obs.executed.ipynb \
  --ExecutePreprocessor.timeout=300
```

Note: full notebook execution can still be slow because the regional cache is
rebuilt when the regional cache tag changes and because many spatial figures may
be generated. For quick debugging, reduce `MAX_SPATIAL_MAP_FIGURES`.

## Next planned work

- Decide whether GEOSldas should add stricter H121 QC, especially whether
  `correction_flag` should be screened. The Python-side looser QC mirrors
  current GEOS; it is not necessarily the scientifically preferred final QC.
- Diagnose remaining H121/OFA scatter after v2 QC:
  - outliers by `n_obs`,
  - outliers by cycle,
  - spatial clusters,
  - model-based QC effects.
- Add compact global maps for:
  - H121-only coverage,
  - Legacy-only coverage,
  - matched H121-Legacy mean difference,
  - observation count differences.
- Consider replacing notebook-heavy figure production with a script once the
  figure set stabilizes.
