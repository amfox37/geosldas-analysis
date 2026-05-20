# IMS snow-cover validation methods

This note documents the IMS table and Figure 8-style map workflow.

Current workflow:

1. Download/convert daily IMS ASCII files:
   `projects/IMS/scripts/download_ims_ascii_to_nc.py`
2. Regrid IMS categories to the GEOS-LDAS M36 grid:
   `projects/IMS/scripts/regrid_ims_to_m36_nearest.py`
3. Compute OL/DA contingency metrics:
   `projects/IMS/scripts/run_ims_ol_da_cell_metrics.py`
4. Make manuscript tables and maps from precomputed outputs:
   `projects/IMS/notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`

The optional notebook
`projects/IMS/notebooks/ims_daily_seasonal_ol_da_scf.ipynb` contains the same
direct daily-validation logic, but the current recommended path is script-first.

## Data provenance

The validation uses the NOAA/National Ice Center Interactive Multisensor Snow
and Ice Mapping System (IMS) Daily Northern Hemisphere Snow and Ice Analysis
from NSIDC/NOAA data set G02156:

- NSIDC landing page:
  <https://nsidc.org/data/g02156/>  
- DOI:
  <https://doi.org/10.7265/N52R3PMC>
- Access root used by the downloader:
  `https://noaadata.apps.nsidc.org/NOAA/G02156`
- Product used in this workflow:
  24 km daily ASCII grids named like
  `imsYYYYDDD_00UTC_24km_v1.3.asc.gz`

Important terminology: NSIDC cites the data set as G02156 Version 1, while the
24 km ASCII filenames carry the IMS file-version tag `v1.3`.

IMS provides daily Northern Hemisphere snow/ice categorical maps. The workflow
uses the 24 km product and treats the IMS categorical field as the reference
binary snow/no-snow observation after the code mapping described below.

## Download and annual NetCDF conversion

`download_ims_ascii_to_nc.py` downloads one daily `.asc.gz` file per day and
writes an annual NetCDF file, for example:

- `projects/IMS/output/ims_snowcover_24km_YYYY.nc4`

The script reads the IMS 24 km latitude/longitude file
`projects/IMS/scripts/ims_24km_latlon.nc4`, parses the daily ASCII grids without
assuming a fixed header length, validates each grid against the expected
1024 x 1024 24 km IMS shape, and stores the raw IMS categorical code as
`ims_snowcover`. Missing/unreadable days can be written as `NaN` if the
`--allow-missing-days` option is used.

## Regridding to M36

`regrid_ims_to_m36_nearest.py` maps IMS categorical data to the
`SMAP_EASEv2_M36_GLOBAL` GEOS-LDAS grid. The method is deliberately
categorical-safe:

- nearest-neighbor mapping only;
- no interpolation, smoothing, or averaging of IMS categories;
- one deterministic representative GEOS-LDAS tile per M36 cell;
- representative-tile rule: maximum `frac_cell`, tie-break by minimum `tile_id`;
- nearest IMS source point found with a unit-sphere `cKDTree`;
- default maximum allowed source-target distance: 60 km;
- unmapped or rejected cells use fill value `-32768`.

The regridded annual products are named like:

- `projects/IMS/output/ims_snowcover_24km_YYYY_on_m36_nearest.nc4`

The output variable used by the validation script is `ims_category`. The
regridded files also store mapping diagnostics such as `within_max_distance`,
`nearest_distance_km`, source row/column indices, representative tile ID,
representative `frac_cell`, representative elevation, and M36 cell lat/lon.

## Model data

The IMS validation compares IMS against daily GEOS-LDAS OL and DA snow-cover
fraction fields from:

- OL experiment: `LS_OLv8_M36`
- DA experiment: `LS_DAv8_M36`
- domain: `SMAP_EASEv2_M36_GLOBAL`
- daily file pattern:
  `<exp>.tavg24_1d_lnd_Nt.YYYYMMDD_1200z.nc4`

The validation script searches for the first available model snow-cover fraction
variable from:

- `FRLANDSNO`
- `FRSNO`
- `SNCOVFR`
- `SNOWCOVERFR`
- `SCF`

Model snow is defined as snow-cover fraction greater than or equal to 0.5.
Model fill values greater than `1e14` and physically invalid negative values are
treated as missing.

## IMS binary coding

IMS categories are converted to binary snow/no-snow before scoring. The default
settings are:

- IMS snow code: `4`
- IMS no-snow codes: `0, 1, 2, 3`
- IMS fill value: `-32768`

The validation script can auto-infer the code sets from the observed categories.
For the current workflow, this recovers the same convention: code `4` is snow
and all other finite non-fill categories are no snow.

## Collocation and fair OL/DA mask

The comparison is performed daily on the representative M36 land-cell vector.
For each day, the script reads IMS and both model experiments and builds a
common valid mask:

`valid = finite(IMS binary snow/no-snow) & finite(OL SCF) & finite(DA SCF) & land/mapping-valid cell`

The same daily collocated cell set is therefore used for OL and DA. This avoids
giving either experiment a different sample population on a given day.

Daily records with successful OL and DA model reads are retained in:

- `ims_ol_da_daily_counts_*.csv/.parquet`
- `ims_ol_da_pair_daily_*.csv/.parquet`

## Domain filtering

For the map/table outputs used by the manuscript workflow, cells are filtered to
retain only locations with at least 10 IMS observed-snow days over the all-period
record. In the script this is:

- `--min-ims-snow-days 10`
- eligibility based on `A + C`, i.e. IMS observed-snow days

Cells failing this threshold are zeroed before map-ready metrics and summary
tables are regenerated. This keeps the domain focused on locations with at
least minimal IMS snow occurrence.

## Contingency counts and metrics

For each collocated cell/day, the script accumulates binary contingency counts:

- `A`: model snow and IMS snow (hit)
- `B`: model snow and IMS no snow (false alarm)
- `C`: model no snow and IMS snow (miss)
- `D`: model no snow and IMS no snow (correct rejection)
- `N = A + B + C + D`

Metrics are computed from accumulated counts:

- `accuracy = (A + D) / N`
- `hit_rate = A / (A + C)`
- `miss_rate = C / (A + C)`
- `false_alarm_ratio = B / (A + B)`
- `correct_rejection_rate = D / (B + D)`

The script writes metrics for all-period, season-all-years, year, and
year-season scopes.

## Confidence intervals

The default full daily workflow uses a paired day-block bootstrap:

- resample days with replacement;
- keep OL and DA counts paired for each resampled day;
- aggregate `A/B/C/D` for OL and DA;
- recompute metrics;
- report percentile 95% confidence intervals from the 2.5th and 97.5th
  percentiles;
- default bootstrap count: 2000;
- default seed: 42.

When `--reuse-cell-counts-nc` fast mode is used, the script regenerates tables
from existing per-cell counts and uses a cell bootstrap because day-level rows
are not reconstructed.

## Current Figure 8 / table products

The plotting notebook is configured for:

- years: 2000-2024
- SCF threshold: 0.5
- minimum IMS snow days: 10
- cache tag:
  `SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10`

Expected precomputed inputs on Discover:

- `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_comparison_table_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.csv`
- `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_pair_daily_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.csv`
- `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_scope_metadata_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.csv`
- `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.nc4`

The all-period table is written by the plotting notebook to:

- `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/tables_ims_maps_and_tables/ims_all_period_ol_da_ci_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10.csv`

The Figure-8-style all-period map panel is written to:

- `/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/figures_ims_maps_and_tables/ims_all_period_delta_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50_imsSnowDaysGe10_nh_robinson_2x3.png`

The screenshot table values correspond to:

| Metric | OL (95% CI) | DA (95% CI) |
| --- | --- | --- |
| accuracy | 0.901 [0.900, 0.902] | 0.912 [0.911, 0.913] |
| hit_rate | 0.824 [0.822, 0.826] | 0.863 [0.861, 0.864] |
| miss_rate | 0.176 [0.174, 0.178] | 0.137 [0.136, 0.139] |
| false_alarm_ratio | 0.116 [0.113, 0.119] | 0.118 [0.116, 0.121] |
| correct_rejection_rate | 0.942 [0.941, 0.944] | 0.938 [0.937, 0.940] |

## Draft methods text

Daily snow-cover validation used the NOAA/National Ice Center IMS Daily Northern
Hemisphere Snow and Ice Analysis distributed by NSIDC as data set G02156
Version 1. We used the 24 km daily ASCII product with filename version tag
`v1.3`, downloaded from the NOAA@NSIDC archive and converted to annual NetCDF
files. The raw IMS categories were regridded to the GEOS-LDAS
`SMAP_EASEv2_M36_GLOBAL` grid using nearest-neighbor categorical mapping. No
interpolation or averaging was applied. For each M36 grid cell, one
representative GEOS-LDAS tile was selected using the largest tile
`frac_cell`, with ties broken by smallest tile ID; the nearest IMS source pixel
to that representative tile was used if the source-target distance was no more
than 60 km.

IMS categorical values were converted to binary snow/no-snow observations before
verification. In the current processing, IMS code `4` was treated as snow and
finite non-fill codes `0, 1, 2, 3` were treated as no snow; `-32768` and
non-finite values were treated as missing. Model snow cover was taken from daily
GEOS-LDAS OL and DA `cat/ens_avg` files. The script used the first available
snow-cover fraction variable among `FRLANDSNO`, `FRSNO`, `SNCOVFR`,
`SNOWCOVERFR`, and `SCF`, and classified a model grid cell as snow-covered when
snow-cover fraction was at least 0.5.

All scores were computed from daily collocated samples using a common valid mask
for OL and DA: finite IMS binary snow/no-snow, finite OL snow-cover fraction,
finite DA snow-cover fraction, and a valid land/mapping cell. Thus the OL and DA
statistics used the same grid cells on each day. Cells with fewer than 10 IMS
observed-snow days over the full analysis period were excluded from the
map/table domain.

For each experiment, contingency counts were accumulated as hits (`A`: model
snow and IMS snow), false alarms (`B`: model snow and IMS no snow), misses
(`C`: model no snow and IMS snow), and correct rejections (`D`: model no snow
and IMS no snow). From these counts we computed accuracy, hit rate, miss rate,
false alarm ratio, and correct rejection rate. Confidence intervals for the
summary table were estimated with a paired day-block bootstrap: days were
resampled with replacement, OL and DA daily counts were kept paired, contingency
counts were reaggregated, and metrics were recomputed. The reported intervals
are the 2.5th and 97.5th percentiles from 2000 bootstrap replicates.

Figure 8 shows DA minus OL spatial differences in the all-period IMS
contingency metrics over the Northern Hemisphere using the filtered M36 domain.
Positive values indicate larger DA than OL metric values; for miss rate and
false alarm ratio, lower values indicate better performance, so sign
interpretation differs from accuracy, hit rate, and correct rejection rate.

## Items still worth verifying

- Confirm whether the final manuscript figure uses the all-period 2x3 map above
  or one of the seasonal `ims_seasonal_delta_*_nh_robinson_2x2.png` products.
- Confirm whether the all-period table in the manuscript was generated from the
  full daily paired bootstrap or from the fast-mode cell bootstrap. The file
  names are the same style, but the script notes the difference in NetCDF
  metadata when fast mode is used.
