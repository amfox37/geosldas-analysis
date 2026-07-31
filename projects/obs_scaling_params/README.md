# Observation Scaling Parameters

Tools for generating GEOS-LDAS observation z-score scaling climatologies from
`ldas_ObsFcstAna` files.

This project was promoted out of `projects/matlab2python/obs_scaling_params`
once the Python port became the maintained workflow. The MATLAB files are kept
under `matlab_reference/` as provenance and parity references.

Two output paths are supported:

- `scripts/run_scaling_params.py` writes the regular lat/lon NetCDF format used
  by the maintained ASCAT workflow.
- `scripts/run_tb_scaling_params.py` begins the tile-space SMAP Tb port and
  writes the big-endian sequential-binary format consumed directly by
  GEOS-LDAS `scale_obs_Tb_zscore`.

The Tb path preserves the MATLAB field order, H/V and incidence-angle
dimensions, EASEv2 M36 administering-tile selection, moving-window statistics,
and pentad naming. Its initial production scope is `hscale=0`, one orbit per
run, and EASEv2 M09/M36 model tiles converted to M36 observation support.

## What This Computes

For each moving time window, observation species group, and 0.25-degree
lat/lon grid cell, the workflow accumulates matched observation and model
forecast statistics:

- `o_mean`: observation mean
- `o_std`: observation standard deviation
- `m_mean`: model forecast mean
- `m_std`: model forecast standard deviation
- `n_data`: number of matched samples
- `m_min`: model forecast minimum
- `m_max`: model forecast maximum

The scaling implied by these fields is:

```text
obs_scaled = m_mean + (obs_raw - o_mean) * (m_std / o_std)
```

## Layout

- `obs_scaling/`: importable Python package.
- `scripts/run_scaling_params.py`: CLI driver with land-sweeper ASCAT defaults.
- `scripts/prepare_fixture_case.py`: builds a tiny GEOS-LDAS-like fixture tree
  from `test_data/inputs`.
- `jobs/run_scaling_params.sbatch`: Discover batch-job template.
- `jobs/run_python_fixture.sbatch`: Python one-timestamp fixture run.
- `jobs/run_matlab_fixture.sbatch`: MATLAB reference one-timestamp fixture run.
- `matlab_reference/`: frozen MATLAB reference implementation.
- `docs/matlab_python_parity.md`: local MATLAB/Python parity check notes.

## Default Run

From this project directory:

```bash
python scripts/run_scaling_params.py
```

The defaults target the current production ASCAT scaling configuration:

```text
Experiment: /discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36
Domain:     SMAP_EASEv2_M36_GLOBAL
Period:     2007-06 to 2024-05
Species:    ASCAT_META_SM, ASCAT_METB_SM, ASCAT_METC_SM
Window:     75 days
Nmin:       20
Grid:       0.25 degree lat/lon
Output:     output/<domain>/stats/python_z_score_clim_quarter_degree/
```

By default, production runs write the all-pentads NetCDF file only. Pass
`--each-doy` when you need one NetCDF file per output DOY.

For another experiment/species set, pass explicit arguments:

```bash
python scripts/run_scaling_params.py \
  --exp-path /discover/nobackup/path/to/archive \
  --exp-run LS_OLv8_M36 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2007-06 \
  --end 2024-05 \
  --species ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM
```

ObsFcstAna input format can be selected with `--obsfcstana-format auto|bin|nc4`.
The default is `nc4` for current production archives. Use
`--obsfcstana-format bin` for legacy binary ObsFcstAna archives, or `auto` to
try `.bin` first and then `.nc4`.

Legacy observation de-duplication is available with `--dedup`, but is off by
default. The de-dup path matches the MATLAB reference behavior by using a
rounded tile/species/lon/lat/obs/year key and keeping the later recent
occurrence.

Daily DOY output is available with `--each-doy`, but is off by default to avoid
large numbers of intermediate NetCDF writes.

On Discover:

```bash
cd /path/to/geosldas-analysis/projects/obs_scaling_params
sbatch jobs/run_scaling_params.sbatch
```

For legacy binary ObsFcstAna archives:

```bash
sbatch jobs/run_scaling_params.sbatch --obsfcstana-format bin
```

## SMAP Tb Tile-Space Run

The SMAP driver requires an explicit experiment and output prefix:

```bash
python scripts/run_tb_scaling_params.py \
  --exp-path /discover/nobackup/path/to/archive \
  --exp-run SPL4SM_OL7000 \
  --domain SMAP_EASEv2_M09_GLOBAL \
  --start 2015-04 \
  --end 2021-03 \
  --orbit D \
  --prefix L4SM_OL7000_SMAPL1CR17000_zscore_stats_
```

Defaults match the MATLAB SMAP setup: `SMAP_L1C` Tb, 40-degree incidence,
75-day windows, 20 samples minimum, three-hourly assimilation cycles, and
EASEv2 M36 administering tiles. Use `--obsparam` when the observation table is
not in the standard dated `rc_out` location.

The generated files can be inspected independently of GEOS-LDAS with
`obs_scaling.seqbin.read_tb_scaling_file()`.

For MATLAB byte-contract parity, the second integer in the file header remains
`0` (the legacy writer's version slot). The requested `--ndata-min` is still
applied when screening every H/V statistic and remains encoded in the filename.

On Discover, the same arguments can be supplied to the batch template through
environment variables:

```bash
EXP_PATH=/discover/nobackup/path/to/archive \
EXP_RUN=SPL4SM_OL7000 \
START_MONTH=2015-04 END_MONTH=2021-03 \
PREFIX=L4SM_OL7000_SMAPL1CR17000_zscore_stats_ \
sbatch jobs/run_python_SMAP_Tb.sbatch
```

## CYGNSS L1 Owner-Tile Run

CYGNSS L1 scaling uses the preprocessor-selected owner tile rather than
reconstructing a grid cell from the specular-point longitude and latitude:

```bash
python scripts/run_cygnss_l1_scaling_params.py \
  --exp-path /discover/nobackup/path/to/archive \
  --exp-run LS_OLv8_M36_CYGNSS \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2020-01 \
  --end 2022-12
```

The driver resolves `CYGNSS_L1_DDM3X5_CROP_SCALAR` through the experiment's
`ldas_obsparam` file and requires NetCDF ObsFcstAna input with
`obsparam_scale=0`. Observation and forecast values are accumulated only as
finite pairs. Multiple records for the same owner tile and cycle are rejected
because the current CYGNSS preprocessor and GEOS-LDAS reader contract permits
at most one selected observation per tile and species.

Output is one compressed NetCDF file with dimensions `pentad=73`, `y=406`, and
`x=964`. The spatial indices are normalized from the tile files as
`i_indg-i_offg-ind_base` and `j_indg-j_offg-ind_base`. Statistics are in dB and
include `o_mean`, `o_std`, `m_mean`, `m_std`, `n_data`, `m_min`, and `m_max`.
Use `--obsparam`, `--tilecoord`, or `--tilegrids` when those files are outside
their normal `rc_out` locations.

## MATLAB/Python Fixture on Discover

The repository includes a compact fixture under `test_data/inputs`. These jobs
copy it into a GEOS-LDAS-like tree and run the same one-timestamp ASCAT scaling
case through Python and the MATLAB reference implementation:

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params
sbatch jobs/run_python_fixture.sbatch
sbatch jobs/run_matlab_fixture.sbatch
```

Both jobs honor these optional environment variables:

```bash
REPO_ROOT=/path/to/geosldas-analysis
FIXTURE_ROOT=/discover/nobackup/$USER/obs_scaling_fixture
OFA_FORMAT=nc4
MATLAB_MODULE=matlab
```

Expected outputs are written under:

```text
$FIXTURE_ROOT/LS_DAv8_M36_as_test2/output/SMAP_EASEv2_M36_GLOBAL/stats/
  python_fixture_zscore_1deg/
  matlab_fixture_zscore_1deg/
```

After both jobs finish, compare the NetCDF outputs:

```bash
python scripts/compare_fixture_outputs.py \
  --fixture-root "$FIXTURE_ROOT" \
  --key-contains DOY121
```

The comparison reports dimensions, finite-cell counts, finite-mask mismatches,
and max/mean absolute differences for `o_mean`, `o_std`, `m_mean`, `m_std`,
`n_data`, `m_min`, and `m_max`.
