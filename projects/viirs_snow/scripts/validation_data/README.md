# validation_data: SNOTEL, GHCN and IMS for the VIIRS transition analysis (analysis C)

Fetches the independent snow data needed to score the five VIIRS_testing arms for 2020-08-01 to 2026-07-31
(`docs/viirs_transition_plan.md`, analysis C). Data are written to
`/discover/nobackup/projects/land_da/VIIRS_testing/validation_data/` (land_da project space, not personal nobackup);
nothing here is committed. Built 2026-11-21. Run on a Discover login node with the GEOSpyD py3.13 python
(`/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3`, use `-I`).

| Script | What it does |
|---|---|
| `get_snotel_extension.py` | SNOTEL daily WTEQ (SWE) and SNWD (snow depth) 2025-10-01..2026-07-31 for all AWDB SNTL stations (batched REST calls). Checks a 7-day overlap against the existing consolidated parquet. |
| `get_ghcn_snwd.py` | GHCN-Daily SNWD, one parquet per calendar year, 2020-08-01..2026-07-31 (same rules as `GHCN_snwd/notebooks/ghcn_snwd_global_1998_present_build.ipynb`). |
| `get_ims_2025_2026.py` | IMS 24 km daily snow cover 2025-01-01..2026-07-31 (one NetCDF per year, raw `.gz` kept only until verified) and a byte-level spot check against four existing 2024 raw files. Uses the functions of `projects/IMS/scripts/download_ims_ascii_to_nc.py`. |
| `verify_validation_data.py` | Coverage, counts, ranges and spot checks against the existing products; writes `validation_report.json`; `--delete-ims-gz` removes the raw IMS cache after verification. |

The scripts import or copy logic from `projects/SNOTEL`, `projects/GHCN_snwd` and `projects/IMS`; none of those files
were modified.

## Products (all under `validation_data/`)

### SNOTEL (`snotel/`)
- `snotel_daily_wteq_snwd_20251001_20260731.parquet`: columns `date` (datetime64), `stationTriplet` (e.g. `301:CA:SNTL`),
  `WTEQ` (SWE), `SNWD` (snow depth), `name`, `latitude`, `longitude`, `elevation`. **Units: inches** (AWDB native; same as
  the existing `projects/SNOTEL/output/all_stations_daily_wteq_snwd.parquet`, which covers 1999-10-01 to 2025-09-30).
  Concatenate the two files for the full record; convert to mm / kg m-2 yourself (1 in = 25.4 mm). NaN = not reported.
  There is no QC-flag column (as in the existing product). Values for recent dates are NRCS provisional and may be revised.
- `snotel_stations_metadata.csv`: AWDB station metadata (934 stations listed; lat/lon/elevation in degrees and feet).
- `snotel_overlap_differences.csv`, `snotel_extension_report.json`: overlap check against the existing parquet.

### GHCN-Daily SNWD (`ghcn/`)
- `ghcn_snwd_<YYYY>.parquet` (2020 .. 2026): columns `station_id`, `date`, `snwd_mm` (float32, **mm**), `mflag`, `qflag`,
  `sflag`, `obs_time`. 2020 starts 2020-08-01 and 2026 ends 2026-07-31. Rules: SNWD only; missing dropped; **any value with a
  non-blank QFLAG dropped** (so `qflag` is always blank); no other QC.
- `ghcn_snwd_station_coverage.parquet`: per station `first_date`, `last_date`, `n_days`, `max_snwd_mm` within the window,
  plus `lat`, `lon`, `elev` (m), `state`, `name`, `gsn_flag`, `hcn_crn_flag`, `wmo_id`. Station coordinates live here, not in
  the yearly files.
- `ghcn_snwd_station_metadata.parquet` (NOAA inventory rows for SNWD with `last_year >= 2020`), `ghcnd-stations.txt`,
  `ghcnd-inventory.txt`, `ghcnd-readme.txt`, `ghcn_snwd_build_log.csv`, `ghcn_snwd_manifest.json`.
- Station selection for scoring (distance to tile, elevation difference, minimum valid days, snow days per year) is NOT
  part of the build; the existing OL/DA validation used max 40 km, max 500 m elevation difference, at least 1500 valid
  days and at least 5 snow days per year (attributes of `GHCN_snwd/outputs_ghcn_snwd_ol_da_validation/ghcn_snwd_raw_timeseries_*.nc`).
  Those thresholds assume a 25-year record; with a 6-year window `min_valid_days` has to be rescaled.

### IMS (`ims/`)
- `ims_snowcover_24km_2025.nc4` (365 days), `ims_snowcover_24km_2026.nc4` (212 days, to 2026-07-31): variable `ims_snowcover`
  (time, y, x), 1024 x 1024 polar stereographic, `lat`/`lon` (y, x) included; `time` is days since `<year>-01-01`. Codes:
  0 outside, 1 open water, 2 land without snow, 3 sea/lake ice, 4 snow-covered land.
- For comparison with the M36 tiles, regrid with the existing script (see below); the IMS-on-M36 files for 2020-2024 already
  exist in `projects/IMS/output/ims_snowcover_24km_<YYYY>_on_m36_nearest.nc4` (variable `ims_category`, 379 x 964 M36 grid,
  nearest neighbour, max distance 60 km), so only 2025 and 2026 were fetched.
- `ims_snowcover_24km_2025_on_m36_nearest.nc4`, `..._2026_on_m36_nearest.nc4`: same layout as the existing files.

```
python projects/IMS/scripts/regrid_ims_to_m36_nearest.py --ims-nc <ims_snowcover_24km_YYYY.nc4> \
  --out-nc <ims_snowcover_24km_YYYY_on_m36_nearest.nc4> \
  --tilecoord /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rc_out/LS_OLv8_M36.ldas_tilecoord.bin
```

## Caveats
- IMS is not independent of MODIS/VIIRS imagery (see the plan); SNOTEL and GHCN are.
- SNOTEL covers the western US (and Alaska); GHCN SNWD is sparse and mostly mid-latitude.
- Recent SNOTEL values are provisional (the overlap check found 10 revised values in 6,347 rows, all small).
