# Observation Scaling Parameters

Tools for generating GEOS-LDAS observation z-score scaling climatologies from
`ldas_ObsFcstAna` files.

This project was promoted out of `projects/matlab2python/obs_scaling_params`
once the Python port became the maintained workflow. The MATLAB files are kept
under `matlab_reference/` as provenance and parity references.

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

The defaults match the local parity-test configuration:

```text
Experiment: /discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36
Domain:     SMAP_EASEv2_M36_GLOBAL
Period:     2007-06 to 2024-05
Species:    ASCAT_META_SM, ASCAT_METB_SM, ASCAT_METC_SM
Window:     75 days
Nmin:       20
Grid:       0.25 degree lat/lon
Output:     output/<domain>/stats/python_z_score_dedup_clim_quarter_degree/
```

For another experiment/species set, pass explicit arguments:

```bash
python scripts/run_scaling_params.py \
  --exp-path /discover/nobackup/path/to/archive \
  --exp-run LS_OLv8_M36 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2007-06 \
  --end 2024-05 \
  --species ASCAT_META_SM,ASCAT_METB_SM,ASCAT_METC_SM \
  --obsfcstana-format auto
```

ObsFcstAna input format can be selected with `--obsfcstana-format auto|bin|nc4`.
The default `auto` preserves legacy behavior by trying `.bin` first and then
`.nc4`. Use `--obsfcstana-format nc4` to force NetCDF ObsFcstAna input.

On Discover:

```bash
cd /path/to/geosldas-analysis/projects/obs_scaling_params
sbatch jobs/run_scaling_params.sbatch
sbatch jobs/run_scaling_params.sbatch --obsfcstana-format nc4
```

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
OFA_FORMAT=auto
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
