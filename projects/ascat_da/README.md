# ascat_da

Diagnostics and utilities for ASCAT-based land data assimilation experiments. The notebooks here generate diagnostic figures and statistics for different ASCAT filtering/footprint configurations, while the scripts automate data preparation.

For the older ASCAT/SMAP paper figure workflow, see `ASCAT_SMAP_paper_provenance.md`. That workflow predates this project layout, so the main publication figure notebook lives under `projects/utils`.

## Current H121 DA project

- `report/h121_da_current_goals.md` - current goals, status, and evidence for
  deciding whether GEOSldas should switch ASCAT assimilation from legacy
  EUMETSAT BUFR to H SAF H121 CDR. This supersedes the older one-day
  integration-test framing for the current production DA comparison.
- `notebooks/h121_legacy_omf_figures.ipynb` - rerunnable O-F summary figure
  notebook comparing SMAP, legacy ASCAT, and H121 ASCAT metrics, with
  species-level diagnostics underneath. Outputs go to
  `output/omf_h121_legacy_figures`.
- `notebooks/h121_iv_tc_skill_figures.ipynb` - full-span IV/TC figure notebook
  for `data/step4_h121_cdr_test_20260725`, using Robinson maps clipped at
  60S and tile-area-weighted spatial means. Outputs go to
  `output/h121_iv_tc_skill_figures`.
- `notebooks/h121_ismn_skill_figures.ipynb` - ISMN station-skill figure
  notebook for `data/ismn_ol_da_skill_bundle`, including paired OL deltas,
  H121-vs-legacy station maps, and tile-area-weighted station/network means.
  Outputs go to `output/h121_ismn_skill_figures`.

## Key notebooks
- `notebooks/h121_iv_tc_skill_figures.ipynb` - full-span IV/TC summary and
  map figures from the H121 step-4 bundle.
- `notebooks/h121_ismn_skill_figures.ipynb` - ISMN in-situ OL/DA skill figures
  from the local ISMN bundle.
- `notebooks/h121_legacy_omf_figures.ipynb` - current H121 vs legacy DA O-F summary figures from `data/omf_compare_sums`, including combined-family and species-level O-F stddev improvement diagnostics.
- `notebooks/legacy_vs_h121_obs.ipynb` - validates raw Legacy BUFR and H121/H139 observations after QC and GEOS M36 tile/cycle super-ob formation, including comparisons with `ObsFcstAna` (OFA). The notebook uses the corrected `geos_cycle_global_v7_gridoriginfix_bsflag_noise` cache described below.
- `notebooks/compare_legacy_bufr_vs_H121.ipynb` – compares legacy BUFR ASCAT (ASCSMR02) against H121 CDR, including observation-space diagnostics from a GEOSldas diagnostic run (per-platform stats, O−F innovations, obs and innovation maps for Jan 2020).
- `notebooks/compare_comb_fp_043024.ipynb` – latest comparison of combined FP experiments with OL, including map plots and PDF summaries.
- `notebooks/ASCAT_masking_021324.ipynb` – builds and inspects ASCAT observation masks for assimilation tests.
- `notebooks/regrid_ASCAT_mask_022724.ipynb` – regrids ASCAT masks to the analysis grid for quick QA.

## Supporting scripts
- `scripts/check_ascat_duplicates.py` – helper to spot duplicate ASCAT files before ingest.
- `scripts/filename_lister_v2.py` – utility for listing/renaming raw observation files.
- `scripts/run_ismn_ol_da_skill.py` – ISMN in-situ soil-moisture skill for OL vs the DA runs (see below).

## ISMN in-situ validation

`scripts/run_ismn_ol_da_skill.py` scores OL and the DA experiments against
International Soil Moisture Network stations over one window: the full
experiment span, `2015-04-01` to `2021-03-31`. Submit it with
`jobs/run_ismn_ol_da_skill.sbatch`.

The method follows `projects/M21C_ls/ISMN_methods_readme.md`, with three
deliberate departures for this project:

- **No predetermined network list.** Every network in the local archive is
  considered; 60 of the 88 networks (2,104 of 3,279 stations) have
  soil-moisture sensors overlapping the window. `--network` can still restrict
  the run, but the default is all of them.
- **One window, so no cross-window common-site intersection.** A station is
  kept when it clears `--nmin` (default 1000) paired obs/model days, and only
  when *every* run clears it for that domain, so run-to-run differences never
  reflect a shifting site population.
- **Generic root zone instead of per-network strict layer rules.** The
  M21C_ls `matlab_strict` composites are hand-tuned for six networks and do
  not generalize. Root zone here is the profile-weighted 0–1 m average, which
  requires at least `--rz-min-sensors` (3) finite layers spanning both a
  shallow (≤0.20 m) and a deep (≥0.50 m) sensor at each timestamp. Surface
  skill is therefore available at far more stations than root-zone skill, and
  the two site counts are reported separately.

Surface observations come from the sensor depth nearest 0.05 m among sensors
no deeper than `--surface-max-depth-m` (0.10), so a station whose shallowest
sensor is 1 m deep yields no surface series rather than a mislabeled one. Only
ISMN `G`-flagged records are retained; observations are daily-averaged and
shifted 12 h to align with the model day.

Model values for **all four runs** are `sm_surface`/`sm_rootzone` from the
`SMAP_L4_SM_gph` collection, averaged over each file's eight 3-hourly
instants. That collection is the analysis state and is the only daily
collection `DAv7_M36_SMAP_type_13_comb_fp_scaled` carries, so using it
everywhere keeps the four runs directly comparable.

Stations are matched to the nearest M36 tile by squared lat/lon distance with
a 0.1 deg² cutoff. Per-station statistics are Pearson `R`, anomaly `R`
(day-of-year climatology, 31-day circular window), `bias`, `RMSE`, and
`ubRMSE`, computed by the shared
`projects/matlab2python/scripts/sm_skill_vs_insitu.py` helpers.

Outputs land in `--output-dir` (default
`/discover/nobackup/projects/land_da/Evaluation/ISMN/ascat_da_ol_da_skill`):

| File | Contents |
| --- | --- |
| `cache_obs_daily.nc` | daily ISMN surface/root-zone series per station |
| `cache_model_daily_<run>.nc` | daily model series + station→tile mapping per run |
| `ismn_station_inventory.csv` | per-station metadata and obs-day counts |
| `ismn_skill_stations.csv` | per station/domain/run statistics |
| `ismn_skill_network_summary.csv` | per network/domain means and deltas vs `--reference-run` |

Both caches are reused on rerun; pass `--overwrite-obs` / `--overwrite-model`
to rebuild. Network-summary deltas are signed so **positive always means the
run beat the reference** (`R`/`anomR` gain, `RMSE`/`ubRMSE` shrink, `bias`
compared in magnitude).

The `ismn` Python package is not installed on Discover, so `lib/ismn_io.py`
parses the archive's `.stm` files and `python_metadata/ISMN_data.csv`
directly.

Full method, provenance, selection-funnel counts, and known limitations:
`report/ismn_insitu_validation_methods.md`.

## M36 tile-assignment correction (July 2026)

Commit `68f7799` corrected the EASEv2 M36 grid origin used by
`lib/superob.py:latlon_to_ij`. The old code treated the geographic pole as
the grid edge, producing a systematic row-index error for Python-derived
Legacy and H121/H139 super-obs. The GEOSldas `ObsFcstAna` files were not
produced by this Python function and are not affected.

Rebuilding the 2020-06-01 through 2020-06-10 global caches with the same
current H121 QC increased raw-to-OFA agreement substantially:

| Product | OFA tile/cycle match, old -> corrected | RMSE, old -> corrected |
| --- | --- | --- |
| Legacy (three platforms) | 69-71% -> 96-99% | about 6.3% -> 1.1-1.2% |
| H121 (three platforms) | 78% -> 84-85% | about 7.9% -> 1.5-1.6% |

Use cache version `geos_cycle_global_v7_gridoriginfix_bsflag_noise` for the
corrected validation. Earlier Python global super-ob caches, including
`geos_cycle_global_v5_bsflag_noise`, retain the incorrect tile assignment and
must not be used as current validation evidence. Derived products that call
`form_super_obs` or `latlon_to_ij` must be regenerated, including Python H121
IV/TC daily pairs. Production scaling parameters that were not generated
through this Python super-ob path are unaffected.
