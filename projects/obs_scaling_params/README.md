# Observation Scaling Parameters

Python scripts that compute the climatological mean and standard deviation of
observations and model forecasts from the nc4 `ldas_ObsFcstAna` files of a
GEOSldas experiment (usually an open loop). GEOSldas uses these statistics for
the z-score scaling of observations before assimilation:

```text
obs_scaled = m_mean + (obs - o_mean) * m_std / o_std
```

These are the Python versions of the MATLAB scripts in
`GEOSldas_GridComp/GEOSldas_App/util/inputs/obs_scaling_params/` and are written
to move there. With the MATLAB fixes on GridComp branch
`feature/amfox/matlab-obsfcstana-nc4`, MATLAB and Python give the same results
(checked on the CF0360 OL1 experiment, Oct 2026).

Only the files in this directory move to GridComp. `cygnss_l1/` (the CYGNSS L1
owner-tile method) stays in geosldas-analysis; it imports `obs_scaling_utils.py`
from here, and nothing here imports from `cygnss_l1/`.

## Files

| File | Contents |
|---|---|
| `run_obs_scaling_params.py` | command-line driver |
| `get_obs_scaling_latlon.py` | soil moisture obs (e.g. ASCAT, CYGNSS) on a global lat/lon grid; nc4 output for `scale_obs_sfmc_zscore` |
| `get_obs_scaling_tb.py` | Tb obs (e.g. SMAP, SMOS) in tile space, one orbit per run; binary output for `scale_obs_Tb_zscore` |
| `obs_scaling_utils.py` | readers (obsparam, tilecoord, tilegrids, nc4 ObsFcstAna) and time helpers |
| `jobs/run_obs_scaling_params.sbatch` | Discover batch job; passes all arguments to the driver |
| `tests/` | `python -m pytest -q tests` |
| `cygnss_l1/` | CYGNSS L1 owner-tile method (geosldas-analysis only, see below); tests: `python -m pytest -q cygnss_l1/tests` |

Other scripts (not part of the scaling run): `scripts/gapfill_scaling_params.py`
(fill gaps in a lat/lon scaling file), `scripts/make_superob_obsfcstana.py`,
`scripts/verify_scaling_omf.py`.

## Method

For each day of the climatological year, the statistics are computed from all
obs/forecast pairs within a moving window of `--window-days` days (75 by
default), over all years of the period, and stored at the center of the window.
Cells/tiles with fewer than `--ndata-min` obs (20 by default) get no-data.

Output goes to `<exp-path>/<exp-run>/output/<domain>/stats/<out-dir>/`
(`--out-dir` default: `z_score_clim`). The `ldas_obsparam` file is read from
`rc_out/Y<start>/M<start>/` unless `--obsparam` is given.

### `--method latlon`

The grid cell of each obs is computed from its lat/lon with `ceil()` in single
precision, as in GEOSldas (`scale_obs_sfmc_zscore`). By default all species
given with `--species` are combined (`_sp_ALL_` in the file name); use
`--no-combine-species` for one file per species. Output: one nc4 file with all
73 pentads (`--no-all-pentads` to skip), and optionally one file per pentad
(`--each-pentad`) or per DOY (`--each-doy`). Variables `o_mean`, `o_std`,
`m_mean`, `m_std`, `n_data` are (pentad, lon, lat); `m_min`, `m_max` are
(lon, lat).

`--exclude-tiles-csv` drops obs at tiles flagged in a CSV file (columns
`tilenum` and, by default, `reject`), e.g. to apply a QC that postdates the
ObsFcstAna files.

```bash
python run_obs_scaling_params.py --method latlon \
  --exp-path /discover/nobackup/projects/land_da/M21C_testing --exp-run OL1 \
  --domain CF0360x6C_GLOBAL --start 2019-02 --end 2020-08 \
  --species ASCAT_HSAF_META_SM,ASCAT_HSAF_METB_SM,ASCAT_HSAF_METC_SM \
  --prefix CF0360_OL1_HSAF_zscore_stats_
```

### `--method tb`

Selects one H-pol and one V-pol Tb species per angle (`--angles`, default 40)
for the orbit (`--orbit A|D`) whose descr contains `--descr-contains` (default
`SMAP_L1C`). Statistics are computed on EASEv2 M36 "administering" tiles, or on
all model tiles with `--no-convert-grid`. One file per pentad (and optionally
per DOY) in the big-endian Fortran sequential format read by GEOSldas
(`read_Tb_scale_param.py` in GridComp reads it too).

```bash
python run_obs_scaling_params.py --method tb \
  --exp-path /path/to/exps --exp-run OL1 --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2015-04 --end 2021-03 --descr-contains SMAP_L1C --orbit D \
  --prefix OL1_SMAP_zscore_stats_
```

### CYGNSS L1 owner-tile method (`cygnss_l1/`, not for GridComp)

`cygnss_l1/run_cygnss_l1_scaling_params.py` computes the statistics for the
CYGNSS L1 species (`--descr`, default `CYGNSS_L1_DDM3X5_CROP_SCALAR`, dB) on
the model tiles of an EASEv2 M36 experiment, from an unscaled (monitoring) run
with at most one CYGNSS L1 obs per tile and ObsFcstAna file. Pentads also need
the obs and model stdv above `--std-epsilon` (1e-6). Output: one nc4 file with
all 73 pentads, with the tile statistics on the dense global EASEv2 M36 grid
(`x`, `y` = 0-based global column and row); default `--out-dir`
`cygnss_l1_z_score_clim`. `--keep-list` uses only the obs that pass the
coherency screen in a parquet file from
`CYGNSS_L1_AZ/scripts/cygl1_coherency_filter.py annotate`.

```bash
python cygnss_l1/run_cygnss_l1_scaling_params.py \
  --exp-path /discover/nobackup/projects/land_da/cygl1_operator_test \
  --exp-run OLv8_M36_AZ_fixedop --start 2020-01 --end 2022-12 \
  --obsparam <exp>/output/SMAP_EASEv2_M36_GLOBAL/rc_out/Y2020/M01/OLv8_M36_AZ_fixedop.ldas_obsparam.20200101_0000z.txt
```

## History

- Ported from MATLAB in 2025 (`projects/matlab2python`), Tb added Jul 2026.
- Oct 2026: refactored into plain scripts in GridComp style; nc4 ObsFcstAna
  only; removed the binary ObsFcstAna reader, the legacy de-duplication option,
  the MATLAB reference copy and fixtures. The CYGNSS L1 owner-tile method moved
  to `cygnss_l1/` (not ready for GridComp). Output is unchanged (bit-identical
  to the previous code on CF0360 OL1 test cases, and for CYGNSS L1 on
  `OLv8_M36_AZ_fixedop`).
- `docs/tb_discover_validation.md` documents the Jul 2026 validation of the Tb
  method against MATLAB and GEOSldas; it refers to the code before the Oct 2026
  refactor (commit `cd3eb01`). The scripts in `docs/tb_discover_validation_scripts/`
  are kept as a record: they import the removed `obs_scaling` package and do not
  run with the current code (check out `cd3eb01` to rerun them).
