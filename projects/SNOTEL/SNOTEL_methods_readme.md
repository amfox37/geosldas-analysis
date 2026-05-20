# SNOTEL snow validation methods

This note documents the SNOTEL SWE and snow-depth validation workflow behind
the manuscript Figure 9-style SWE panel:

- `projects/SNOTEL/notebooks/snow_daily_seasonal_ol_da_swe_snwd.ipynb`
- `projects/SNOTEL/outputs_snotel_ol_da_validation/figures/snotel_swe_2x3_bars_allsites_maps_da_minus_ol_elevfilt500_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.png`

The figure shown in the manuscript is the SWE-only 2 x 3 summary: top-row
seasonal OL/DA bars for RMSE, unbiased RMSE, and absolute bias; bottom-row maps
of station-level DA - OL differences for the same metrics.

## Data provenance

The observations are daily SNOTEL data downloaded from USDA NRCS AWDB using:

- `projects/SNOTEL/scripts/download_snotel_swe_snwd.py`
- AWDB SOAP/WSDL station discovery and metadata endpoint:
  `https://wcc.sc.egov.usda.gov/awdbWebService/services?WSDL`
- AWDB REST v1 daily data endpoint:
  `https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1/data`

The downloader requests:

- network code: `SNTL`
- duration: `DAILY`
- elements: `WTEQ` and `SNWD`
- default download window: `1999-10-01` to `2025-09-30`

The local combined source file used by the validation workflow is:

- `projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet`

That parquet was last modified locally on `2026-03-01 12:57:21` and contains
daily `WTEQ`, `SNWD`, `stationTriplet`, station name, latitude, longitude, and
elevation. The downloader preserves AWDB values as returned; unit conversion is
done only in the validation notebook.

## Model data

The validation compares SNOTEL observations with daily GEOS-LDAS OL and DA
output on the `SMAP_EASEv2_M36_GLOBAL` grid:

- OL experiment: `LS_OLv8_M36`
- DA experiment: `LS_DAv8_M36`
- daily file pattern:
  `<exp>.tavg24_1d_lnd_Nt.YYYYMMDD_1200z.nc4`

The snow variables read from the model are:

- `SNOMASLAND`: model SWE, in `kg m-2`
- `SNODPLAND`: model snow depth, in `m`

The Figure 9-style panel uses only `SNOMASLAND` against SNOTEL `WTEQ`.

## Observation preprocessing and units

The notebook requires `date`, `stationTriplet`, `latitude`, `longitude`,
`WTEQ`, and `SNWD`. Rows missing date, station ID, latitude, or longitude are
dropped before collocation.

SNOTEL variables are converted as:

- `WTEQ` inches of water equivalent to SWE: `WTEQ * 25.4 kg m-2`
- `SNWD` inches to snow depth: `SNWD * 0.0254 m`
- station elevation from AWDB feet to meters: `elevation * 0.3048`

If duplicate station/day records are present, the notebook aggregates to one
daily record per station by taking the mean observation value for that day and
the first station metadata values.

## Temporal collocation

The analysis window is:

- `2000-06-01` to `2024-06-01`, inclusive

Temporal matching is daily and exact by date. SNOTEL daily observations are
aligned to the corresponding GEOS-LDAS daily output date. No subdaily
interpolation or nearest-time search is applied in this workflow.

The raw collocated cache for the current run is:

- `projects/SNOTEL/outputs_snotel_ol_da_validation/snotel_raw_timeseries_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.nc`

The notebook output reports:

- daily observation rows in the analysis window: `6,927,973`
- unique stations with observations in the analysis window: `911`
- raw cache dimensions: `time=8767`, `station=910`, `exp=2`
- experiments in the cache: `OL`, `DA`

## Spatial collocation and representativeness

Each SNOTEL station is mapped to the nearest GEOS-LDAS M36 tile using
great-circle distance (`haversine_km`). The configured maximum station-to-tile
distance is `40 km`. A squared-degree threshold of `0.10` is present in the
notebook as a fallback setting, but the current cache uses the haversine method.

The cache stores station and tile coordinates, tile index, station elevation,
tile elevation, horizontal distance, and elevation mismatch:

- `distance_km`
- `elev_diff_m = tile_elev_m - station_elev_m`

The top-row bars use all stations with finite station metrics for each season,
metric, and experiment. The bottom-row maps apply an additional
representativeness filter:

- retain stations with finite OL elevation mismatch and `abs(elev_diff_m) < 500 m`

This elevation filter is used only for the mapped station deltas in the
manuscript-style SWE panel. The all-site map extent is computed from all
elevation-filtered SWE sites, and each bottom-row metric then requires finite
OL and DA station metrics. The current plotted map panels report `611` valid
stations for RMSE, unbiased RMSE, and absolute bias.

## Quality control and missing data

Observation QC in this workflow is structural and finite-value based:

- drop records missing station ID, date, latitude, or longitude;
- convert `WTEQ`, `SNWD`, and elevations with numeric coercion;
- aggregate duplicate station/day rows by mean;
- use only finite paired observation-model samples when computing metrics.

Model values greater than `1e14` are treated as fill and set to missing. Negative
model snow values are treated as physically invalid and set to missing. Missing
observations, missing model values, and invalid model values are omitted
pairwise for each station, experiment, variable, and season.

No additional independent SNOTEL quality-flag screen, frozen-ground filter, or
snow-on/snow-off threshold is applied. Because the validation target is snow
itself, frozen conditions are not removed. Instead, the seasonal aggregation
handles snow-season representativeness explicitly: the "ALL" statistic excludes
JJA, and JJA is shown as not applicable in the Figure 9-style bars.

## Seasons

The notebook computes station and domain metrics for:

- `ALL`
- `DJF`
- `MAM`
- `JJA`
- `SON`

For the manuscript-style SWE panel:

- `EXCLUDE_JJA_FROM_ALL = True`
- `ALL` includes non-JJA days only
- the figure shows `ALL`, `SON`, `DJF`, and `MAM`
- JJA is omitted/shown as not applicable for the bar summary

The notebook reports `6558` ALL-mask days out of `8767` total daily dates when
JJA is excluded.

## Metrics

Metrics are computed from paired finite samples for each station, experiment,
variable, and season:

- `Bias = mean(model - observation)`
- `RMSE = sqrt(mean((observation - model)^2))`
- `ubRMSE = sqrt(mean(((model - mean(model)) - (observation - mean(observation)))^2))`
- `NSE = 1 - sum((observation - model)^2) / sum((observation - mean(observation))^2)`

The SWE figure uses:

- RMSE
- unbiased RMSE
- absolute bias, `abs(Bias)`

For all three plotted SWE metrics, smaller values are better.

## Confidence intervals and station counts

The top-row bars are station-level means, not pooled daily-sample metrics. For
each season, metric, and experiment, the notebook takes the finite station
metric values and estimates the mean and 95 percent confidence interval using
a station bootstrap:

- bootstrap replicates: `1000`
- confidence interval: percentile 95 percent interval
- random seed: `42`

The machine-readable table used for the top-row panels is:

- `projects/SNOTEL/outputs_snotel_ol_da_validation/tables/snotel_swe_toprow_bar_values_ci_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.csv`

Current top-row station counts are:

- `878` finite SWE stations for `ALL` and `SON`
- `876` finite SWE stations for `DJF` and `MAM`

The station-metric table for the full validation is:

- `projects/SNOTEL/outputs_snotel_ol_da_validation/snotel_station_metrics_SMAP_EASEv2_M36_GLOBAL_20000601_20240601.csv`

It contains `18,200` rows, corresponding to `910` stations, 2 experiments, 2
variables, and 5 seasons.

## Figure-specific DA - OL convention

The bottom-row maps use the convention:

- mapped value = `DA - OL`
- for bias maps, both OL and DA are converted to absolute bias before differencing
- negative values therefore indicate DA improvement for RMSE, unbiased RMSE,
  and absolute bias

The current Figure 9-style PNG and top-row CSV were last modified locally on
`2026-04-12 13:34:09`.

## Draft methods text

Daily snow observations from USDA NRCS SNOTEL were obtained through AWDB for
network `SNTL`, using daily `WTEQ` and `SNWD` records. Station metadata were
obtained from the AWDB SOAP service and daily time series from the AWDB REST
service. The downloaded archive covered 1 October 1999 to 30 September 2025;
the validation used 1 June 2000 to 1 June 2024. AWDB values were preserved
during download, then converted in the validation workflow: `WTEQ` from inches
of water equivalent to `kg m-2` by multiplying by 25.4, `SNWD` from inches to
meters by multiplying by 0.0254, and station elevation from feet to meters.

For each SNOTEL station, daily observations were collocated to the nearest
GEOS-LDAS M36 tile by great-circle distance, with a maximum allowed station-tile
distance of 40 km. The cache stores station-tile distance and tile minus station
elevation difference. To reduce representativeness errors in the station-delta
maps, mapped SWE differences were shown only for stations with
`abs(tile elevation - station elevation) < 500 m`; the map panels require finite
paired OL and DA station metrics and include 611 stations. The seasonal bar
panels use all stations with finite station-level SWE metrics for the relevant
season and experiment.

GEOS-LDAS daily OL and DA fields were compared to SNOTEL by exact daily date.
SWE validation used model `SNOMASLAND` against SNOTEL `WTEQ`; snow-depth
validation in the same workflow used model `SNODPLAND` against SNOTEL `SNWD`.
Rows missing station ID, date, latitude, or longitude were removed, duplicate
station/day observations were averaged, model fill values greater than `1e14`
and negative model snow values were set to missing, and statistics were computed
from finite paired observation-model samples only. No additional frozen-ground
or snow-threshold filter was applied because the target variable is snow; to
avoid summer no-snow dominance in the aggregate, the all-season SWE statistic
excluded JJA, and JJA was not plotted in the seasonal bar figure.

For each station, experiment, and season we computed bias as
`mean(model - observation)`, RMSE, and unbiased RMSE. The manuscript SWE figure
plots RMSE, unbiased RMSE, and absolute bias. Bar values are station-mean
statistics with 95 percent confidence intervals estimated by 1000-replicate
station bootstrapping. Map colors show DA - OL station differences; for all
three plotted metrics, negative values indicate improvement in the DA run.

## Items still worth verifying

- Add the exact AWDB/SNOTEL access date if it is needed in the manuscript; the
  local parquet timestamp is `2026-03-01 12:57:21`.
- Confirm whether the manuscript should cite the figure as the all-site extent
  version or one of the alternate tight-site/western-CONUS versions.
