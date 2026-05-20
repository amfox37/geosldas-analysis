# GHCN snow-depth validation methods

This note documents the GHCN-Daily snow-depth workflow behind the manuscript
Figure 10-style Northern Hemisphere SNWD panel:

- `projects/GHCN_snwd/notebooks/ghcn_snwd_global_1998_present_build.ipynb`
- `projects/GHCN_snwd/notebooks/ghcn_snwd_daily_seasonal_ol_da_snwd_baseline_basic.ipynb`
- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/figures/ghcn_baseline_core_snodpland_2x3_bars_maps_nh_da_minus_ol_ALL_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.png`

The figure is the baseline-core GHCN SNWD 2 x 3 summary: top-row seasonal
OL/DA station-mean bars for RMSE, unbiased RMSE, and absolute bias; bottom-row
Northern Hemisphere maps of station-level DA - OL differences for the same
metrics.

## Data provenance

The observations are NOAA/NCEI GHCN-Daily snow-depth records for element
`SNWD`. The build notebook uses the GHCN-Daily `by_year` CSV files rather than
per-station `.dly` files:

- source root: `https://www.ncei.noaa.gov/pub/data/ghcn/daily`
- yearly files: `https://www.ncei.noaa.gov/pub/data/ghcn/daily/by_year/YYYY.csv.gz`
- station metadata: `ghcnd-stations.txt`
- station inventory: `ghcnd-inventory.txt`
- readme: `readme.txt`

The local build manifest is:

- `ghcn_snwd_build/output/ghcn_snwd_manifest.json`

The manifest records:

- dataset: `NOAA/NCEI GHCN-Daily SNWD subset`
- element: `SNWD`
- units: `mm`
- source years: `2000` to `2025`
- `drop_qflagged = true`
- `drop_missing = true`
- build time: `2026-03-19T23:04:43.869553+00:00`

GHCN-Daily is a living archive rather than a fixed versioned product in this
workflow, so the safest manuscript wording is to report it as a NOAA/NCEI
GHCN-Daily archive snapshot built on 19 March 2026 from the GHCN-Daily
`by_year` files.

## Observation build

The build notebook streams one yearly CSV at a time, keeps only
`ELEMENT == "SNWD"`, and writes a partitioned parquet data set:

- `ghcn_snwd_build/output/ghcn_snwd_parquet/year=YYYY/part-*.parquet`

The NOAA `by_year` rows contain:

- `ID`
- `YYYYMMDD`
- `ELEMENT`
- `DATA_VALUE`
- `MFLAG`
- `QFLAG`
- `SFLAG`
- `OBS_TIME`

The build applies the following observation QC:

- retain only `SNWD`
- drop `DATA_VALUE == -9999`
- drop rows with nonblank `QFLAG`
- parse daily dates from `YYYYMMDD`
- store `SNWD` as `snwd_mm`

The build also writes:

- `ghcn_snwd_station_metadata.parquet`
- `ghcn_snwd_inventory.parquet`
- `ghcn_snwd_station_coverage.parquet`
- `ghcn_snwd_year_summary.csv`
- `ghcn_snwd_build_log.csv`

For the current local build, `ghcn_snwd_year_summary.csv` contains 26 yearly
rows from 2000 to 2025 and `90,271,710` retained SNWD rows after the missing
and QFLAG filters.

## Model data

The validation compares GHCN-Daily SNWD against daily GEOS-LDAS OL and DA
output on the `SMAP_EASEv2_M36_GLOBAL` grid:

- OL experiment: `LS_OLv8_M36`
- DA experiment: `LS_DAv8_M36`
- daily file pattern:
  `<exp>.tavg24_1d_lnd_Nt.YYYYMMDD_1200z.nc4`

The model snow-depth quantity used for this workflow is grid-cell-average snow
depth:

- `model_snwd_mm = SNODPLAND * SCF * 1000`

`SNODPLAND` is read in meters. The notebook searches for the first available
snow-cover fraction variable from:

- `FRLANDSNO`
- `FRSNO`
- `SNCOVFR`
- `SNOWCOVERFR`
- `SCF`

If the snow-cover fraction is stored as percent, values are divided by 100.
Snow-cover fraction is then clipped to `[0, 1]`.

## Analysis window and temporal collocation

The validation window used for the figure is:

- `2000-01-01` to `2024-12-31`, inclusive

Temporal matching is daily and exact by date. GHCN daily observations are
pivoted to station-by-date matrices and reindexed to the GEOS-LDAS daily file
dates. No subdaily interpolation or nearest-time search is applied.

The raw collocated cache for the current run is:

- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_snwd_raw_timeseries_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.nc`

The notebook output reports:

- raw cache dimensions: `time=9132`, `station=10119`, `exp=2`
- experiments in the cache: `OL`, `DA`
- domain: `SMAP_EASEv2_M36_GLOBAL`

## Station selection and spatial collocation

Before model extraction, stations must have finite latitude and longitude and
must pass observation-coverage filters:

- at least `1500` valid GHCN SNWD days in 2000-2024
- at least `5.0` observed-snow days per year on average, where observed snow is
  `SNWD > 0`

Each retained station is mapped to the nearest GEOS-LDAS M36 tile using
great-circle distance (`haversine_km`). The current constraints are:

- maximum station-to-tile distance: `40 km`
- maximum absolute tile minus station elevation difference: `500 m`
- common station membership across OL and DA

The station metadata and map diagnostics are written to:

- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_station_selection_summary_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.csv`
- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_station_tile_map_OL_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.csv`
- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_station_tile_map_DA_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.csv`

The current selection summary contains `10,119` common OL/DA mapped stations.
All selected stations are in the Northern Hemisphere in this run
(`station_lat` range `32.9542` to `81.1667` degrees).

## Baseline-core metric sample

The manuscript figure uses the `baseline_core` station metrics, not the
`obspos` or `filled_obs_snowmonths` diagnostic variants. In the baseline-core
workflow:

- only reported, QC-passing GHCN SNWD observations are used;
- missing GHCN observation days are not filled as zero;
- statistics use finite paired observation-model samples only.

Before computing station metrics, the notebook applies a common OL/DA paired
sample filter using the ALL-season mask:

- at least `1500` common paired OL/DA days;
- at least `5.0` common observed-snow days per year on average.

This pre-stat filter reduces the plotted baseline-core metric set from `10,119`
mapped stations to `9,451` stations. The station-metric table is:

- `projects/GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_station_metrics_baseline_core_SMAP_EASEv2_M36_GLOBAL_20000101_20241231.csv`

It contains `94,510` rows, corresponding to `9,451` stations, 2 experiments,
and 5 seasons.

## Quality control and missing data

Observation QC begins in the GHCN build:

- non-SNWD elements are discarded;
- missing value `-9999` is discarded;
- nonblank GHCN `QFLAG` rows are discarded.

Validation then applies physical range checks to both observations and model
values:

- minimum valid snow depth: `0 mm`
- maximum valid snow depth: `< 5000 mm`
- finite observation and model values are required for each paired sample.

Model `SNODPLAND` values greater than `1e14` are treated as fill and set to
missing. Negative model `SNODPLAND` is set to missing. Snow-cover fraction fill
values greater than `1e14` are set to missing before clipping.

No independent frozen-ground filter is applied. Because the target variable is
snow depth, frozen conditions are part of the validation sample. No snow-on
threshold is used when computing the baseline metrics beyond the station
selection requirements described above.

## Seasons

The notebook computes metrics for:

- `ALL`
- `DJF`
- `MAM`
- `JJA`
- `SON`

For the current Northern Hemisphere figure:

- plotted bar seasons are `ALL`, `DJF`, `MAM`, and `SON`
- JJA is excluded from the bar panels
- bottom-row maps use `ALL`
- `ALL` excludes JJA for Northern Hemisphere stations

The notebook also contains Southern Hemisphere logic that excludes DJF from ALL
for Southern Hemisphere stations. It does not affect the current figure because
all selected stations are Northern Hemisphere.

## Metrics

Metrics are computed from finite paired samples for each station, experiment,
and season:

- `Bias = mean(model - observation)`
- `abs_bias = abs(Bias)`
- `RMSE = sqrt(mean((model - observation)^2))`
- `ubRMSE = sqrt(mean((error - Bias)^2))`
- `R`: Pearson correlation
- `anomR`: anomaly correlation after removing month-of-year climatology from
  paired valid samples

The manuscript 2 x 3 SNWD figure uses:

- RMSE
- unbiased RMSE
- absolute bias

For all three plotted metrics, smaller values are better.

## Confidence intervals and station counts

The top-row bars are station-level means. For each season, metric, and
experiment, the notebook bootstraps the finite station metric values:

- bootstrap replicates: `1000`
- confidence interval: percentile 95 percent interval
- random seed: `42`

For the current figure, finite station counts are:

- `9,451` stations for `ALL`, `DJF`, `MAM`, and `SON`
- `8,735` stations for `JJA`, which is not plotted in the 2 x 3 figure

The top-row station-mean values for the plotted baseline-core bars are:

| Metric | Season | OL | DA |
| --- | --- | ---: | ---: |
| RMSE (mm) | ALL | 126.253 | 124.375 |
| RMSE (mm) | DJF | 148.926 | 146.164 |
| RMSE (mm) | MAM | 144.530 | 142.731 |
| RMSE (mm) | SON | 34.390 | 34.253 |
| ubRMSE (mm) | ALL | 100.633 | 99.128 |
| ubRMSE (mm) | DJF | 102.077 | 100.514 |
| ubRMSE (mm) | MAM | 95.424 | 94.150 |
| ubRMSE (mm) | SON | 31.110 | 30.503 |
| abs(Bias) (mm) | ALL | 66.802 | 65.651 |
| abs(Bias) (mm) | DJF | 92.789 | 90.810 |
| abs(Bias) (mm) | MAM | 94.329 | 93.044 |
| abs(Bias) (mm) | SON | 12.075 | 12.842 |

## Figure-specific DA - OL convention

The bottom-row maps use:

- mapped value = `DA - OL`
- negative values indicate lower DA error and therefore DA improvement

For the current `ALL`-season map panels:

- RMSE map: `9,451` stations, mean DA - OL = `-1.878 mm`
- ubRMSE map: `9,451` stations, mean DA - OL = `-1.505 mm`
- absolute-bias map: `9,451` stations, mean DA - OL = `-1.151 mm`

The current Figure 10-style PNG was last modified locally on
`2026-04-12 13:42:12`. The baseline-core station metrics CSV was last modified
locally on `2026-04-08 11:56:16`.

## Draft methods text

Daily snow-depth observations were taken from the NOAA/NCEI GHCN-Daily archive
for element `SNWD`. We built a local GHCN-Daily SNWD snapshot from the
GHCN-Daily `by_year` files for 2000-2025, retaining only rows with
`ELEMENT == SNWD`, non-missing `DATA_VALUE`, and blank `QFLAG`. The resulting
station-level parquet archive stores snow depth in millimeters and was built on
19 March 2026. The validation used observations from 1 January 2000 through
31 December 2024.

Stations were required to have finite coordinates, at least 1500 valid SNWD
observation days, and at least 5 observed-snow days per year on average during
the validation period. Each station was collocated to the nearest GEOS-LDAS M36
tile using great-circle distance. Stations were retained only if the nearest
tile was within 40 km and the absolute tile minus station elevation difference
was no greater than 500 m, with common station membership required for OL and
DA. The raw collocated cache contains 10,119 Northern Hemisphere stations.

GHCN daily SNWD observations were matched to GEOS-LDAS daily output by exact
date. Model snow depth was computed as `SNODPLAND * SCF * 1000`, where
`SNODPLAND` is snow depth in meters and `SCF` is the first available
snow-cover-fraction variable among `FRLANDSNO`, `FRSNO`, `SNCOVFR`,
`SNOWCOVERFR`, and `SCF`. Model fill values and negative model snow depths were
set to missing; snow-cover fraction was converted from percent if needed and
clipped to [0, 1]. Metrics were computed from finite paired samples with both
observed and modeled snow depth in the range 0 to less than 5000 mm.

The baseline-core figure uses only reported GHCN SNWD days; missing GHCN days
were not filled as zero. A common OL/DA pre-statistical filter required at least
1500 paired days and at least 5 observed-snow days per year on average,
retaining 9,451 stations for the plotted baseline metrics. The all-season
Northern Hemisphere statistic excludes JJA; the plotted seasons are ALL, DJF,
MAM, and SON. Station-level RMSE, unbiased RMSE, and absolute bias were
averaged across stations, with 95 percent confidence intervals estimated from
1000 bootstrap resamples of stations. Map panels show station-level DA - OL
differences for the ALL season, so negative values indicate lower DA error.

## Items still worth verifying

- If the manuscript needs a formal GHCN-Daily citation, add the exact NOAA/NCEI
  citation text and archive access date alongside the local build date.
- The figure shown here is the baseline-core version. Do not substitute the
  `obspos` or `filled_obs_snowmonths` outputs unless the methods text is also
  changed.
