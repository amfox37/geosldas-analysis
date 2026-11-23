# IMS Validation Workflow

This project validates GEOS-LDAS OL/DA snow cover against IMS daily snow-cover products on the M36 grid.

## Latest workflow run order

1. Build annual IMS NetCDF from daily ASCII:
   - `scripts/download_ims_ascii_to_nc.py`

2. Regrid IMS annual files to M36 EASE (nearest-neighbor categorical mapping):
   - `scripts/regrid_ims_to_m36_nearest.py`

3. Build OL/DA vs IMS validation metrics (recommended scripted pipeline):
   - `scripts/run_ims_ol_da_cell_metrics.py`

4. Make maps/tables from precomputed outputs:
   - `notebooks/ims_maps_and_tables_from_precomputed_outputs.ipynb`

5. Optional direct notebook workflow (no precompute step):
   - `notebooks/ims_daily_seasonal_ol_da_scf.ipynb`

6. Optional regridding QC map checks:
   - `notebooks/ims_regrid_qc_colorado.ipynb`

## Inputs

- IMS daily source data (`imsYYYYDDD_00UTC_24km_v1.3.asc.gz`) and IMS 24-km lat/lon file.
- M36 tilecoord (used by regridding script).
- GEOS-LDAS daily model outputs for:
  - `LS_OLv8_M36`
  - `LS_DAv8_M36`

## Outputs

- IMS annual files (example):
  - `projects/IMS/output/ims_snowcover_24km_YYYY.nc4`
- Regridded IMS annual files:
  - `projects/IMS/output/ims_snowcover_24km_YYYY_on_m36_nearest.nc4`
- Validation outputs:
  - `projects/IMS/output/outputs_ims_ol_da_validation/`
  - includes metric caches/tables such as:
    - `ims_ol_da_cell_counts_metrics_*.nc4`
    - `ims_ol_da_scope_metadata_*.csv`
    - `ims_ol_da_comparison_table_*.parquet/.csv`
    - `ims_ol_da_pair_daily_*.parquet/.csv`
    - `ims_ol_da_daily_counts_*.parquet/.csv`

## Notes

- Current default rerun path is script-first (`run_ims_ol_da_cell_metrics.py`) plus plotting notebook.
- `run_ims_ol_da_cell_metrics.py --fast-mode` reuses cached counts and skips model re-read.
- `--min-ims-snow-days` controls analysis domain filtering (default `10`).

## 2024+ orientation fix (2026-11-22/23)

**Problem.** The 24 km IMS ASCII files from 2024-01-01 onward are stored bottom-up (header: "(1,1) starts at: lower left
corner"); the 2000-2023 arrays were built top-down. `download_ims_ascii_to_nc.py` does not flip, so 2024 came out upside
down. After the nearest-neighbour regrid to M36, the share of land cells labelled IMS code 1 ("sea") jumped from about 6%
(2023) to about 43% (2024), and 2024 skill was biased low (OL accuracy 0.736, DA-OL accuracy -0.0006).

**Fix.** Flip the 2024 field vertically (`raw[::-1]`, lat/lon untouched), regrid, recompute the 2024 per-cell counts, swap
the 2024 contribution into the 2000-2024 counts, and regenerate the tables with `--reuse-cell-counts-nc`. The 10-snow-day
mask is rebuilt from the corrected counts. Scripts (Discover-only paths, kept for provenance):
`scripts/flip2024.py`, `merge_counts.py`, `run_metrics.sh`, `run2024.sh`, `regen.j`, `regen2.j`, `extract2024.j`,
`p24_new.j`, `p24_old.j`. Work area on Discover: `/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/`.

**Effect on the 2000-2024 products (ALL_PERIOD, OL to DA).** Accuracy 0.901 to 0.912 became 0.921 to 0.934; hit rate
0.824 to 0.863 became 0.825 to 0.864; false alarm ratio 0.116 to 0.118 became 0.044 to 0.045; correct rejection rate 0.942 to
0.938 became 0.978 to 0.976. DA-OL hit-rate improvement is unchanged (+0.0389). 2024 DA-OL accuracy is now +0.0170. Eligible
cells dropped from 48,721 to 46,076 because the old mask let wrong cells pass. Counts for 2000-2023 on cells eligible in both
versions are identical. The 2000-2007 Terra/Aqua scope products predate 2024 and are unchanged.

**Checked locally.** Staged table has 650 rows and agrees with the parquet; ALL_PERIOD counts equal the sum of the YEAR counts
and match the table; the corrected 2024 file is the original reversed north-south for the 7 days the old local copy held.

**Not verified / still open.**
- `merge_counts.py check` (old-IMS 2024 rerun reproduces the original counts) was not re-run locally; it needs the
  intermediate files on Discover.
- The old local `ims_snowcover_24km_2024.nc4` held only 7 valid days, so the flip of the other 359 days rests on the
  Discover run and the sea-share check (4.6 to 6.1% of land cells).
- 2025 and 2026 IMS files are not corrected. Flip them at extraction before use, or teach `download_ims_ascii_to_nc.py` to
  honour the header.
- `ims_ol_da_pair_daily_*` and `ims_ol_da_daily_counts_*` are empty fast-mode placeholders (written whenever
  `--reuse-cell-counts-nc` is used); nothing in this repo reads them.
- `ims_daily_seasonal_ol_da_scf.ipynb` (YEAR_END = 2023) reads `/gpfsm/dnb06/projects/p163/IMS/`, a different directory from
  `projects/IMS/output/`; check that copy before extending it past 2023. `ims_regrid_qc_colorado.ipynb` is fixed to 2020 but
  would show a flipped year if `YEAR` is set to 2024 or later.

