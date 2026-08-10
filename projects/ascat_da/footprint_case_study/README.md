# Peat-QC footprint case study — instructions for a local AI assistant

## Context

This directory is a self-contained data bundle from a GEOSldas land data
assimilation investigation on NASA's Discover HPC cluster. It answers one
question: **why does a peat-quality-control (peat-QC) filter on ASCAT soil
moisture observations leave a narrow band of very-low-but-nonzero
observation counts near peatland boundaries in Canada/Alaska?**

The short answer, already established and verified: a per-observation
"footprint peat fraction" (the peat-classified area within each
observation's ~25 km search radius, weighted by a Gaussian distance kernel
and model-tile area) predicts, with **100% agreement** across 2,720
candidate observations near 6 named boundary tiles, whether an observation
was kept or rejected by the QC. Rejected observations have footprint peat
fraction >= 0.10; kept observations sit at ~0%. This is not a correlation —
it is the actual QC threshold logic from the GEOSldas Fortran source,
reimplemented in Python and cross-checked against the model's own porosity
BCS file, with every claim backed by a hard assertion (see "Validation"
below), not just a printed number someone eyeballed.

## What's in this directory

- **`plot_footprint_case_study.py`** — run this. Self-contained: only needs
  `numpy`, `pandas`, `matplotlib`. No cluster paths, no netCDF4, no network
  access. Reads the CSVs below, recomputes the footprint peat fraction per
  observation, runs three hard validation assertions, and regenerates all
  figures plus the results CSV.
- `data/local_tile_geometry_porosity.csv` — model (catchment) tile centers, areas,
  and porosity (`POROS`) for the neighborhoods around the 6 boundary tiles.
- `data/obs_candidates.csv` — every candidate ASCAT observation (species 8/9,
  H SAF Metop-A/B soil moisture) near those 6 tiles, with its baseline
  lon/lat/obs/fcst/assim_flag and whether it was `kept` or lost by the
  peat-QC run. **Built directly from raw OFA files by `extract_footprint_data.py`
  — this is not sourced from an external/undocumented CSV.**
- `data/constants.csv` — every physical constant used (FOV, search-radius factor,
  porosity threshold, max peat fraction, Earth radius), each with the exact
  source file/line it came from in the GEOSldas Fortran source.
- `validation_report.txt` — the log of every assertion that passed during
  bundle construction (stage 1, on Discover): key uniqueness, N_data
  closure. Read this to see the evidence, not just take the "kept" column
  on faith.
- `extract_footprint_data.py` — the stage-1 script that built this bundle
  from the live Discover filesystem (raw OFA files + BCS porosity + tile
  geometry + N_data closure check). It will NOT run locally (references
  `/discover/...` paths and needs `netCDF4`) — included for provenance/audit
  so you can see exactly how every CSV was derived, line by line.
- `output/footprint_case_study_results.csv` — the reproducible per-observation
  reconstruction, including the calculated peat fraction and predicted QC
  decision. This generated table is intentionally not committed.
- `figures/footprint_case_study_maps.png`,
  `figures/footprint_case_study_threshold.png`, and
  `figures/footprint_case_study_timeseries_3823.png` — overview diagnostics.
- `figures/footprint_ellipse_maps/tile_<tilenum>_footprint_ellipses.png` — one paired
  map per target tile, comparing the geographically nearest actual retained
  and rejected OFA observations over a pcolormesh of local POROS. The solid
  ellipse is the 12.5 km Gaussian weighting scale; the dashed ellipse is the
  25 km contribution cutoff. The labelled magenta rings identify contributing
  model-tile centers and their normalized area-and-distance weights. The
  plotted rectangles are display cells inferred from tile centers, not exact
  catchment polygons.

## What "kept" means (exact predicate)

A baseline candidate observation is `kept=True` if a **valid** row exists
for the same `(timestamp, species, tilenum)` key in the peat-QC output,
where "valid" = the row exists in the OFA file AND its `obs` value is
finite (not the netCDF fill value). Otherwise it is `kept=False` ("lost").

This was not assumed — four candidate predicates were tested against
production `N_data` counts (row-exists / finite-obs / finite-obs-and-fcst /
`assim_flag==1`), and **all four turned out to be numerically identical**
for these two ASCAT HSAF species on these 6 tiles (12/12 tile-species pairs
exact match, both experiments — see `validation_report.txt`, Stage 1b). So
there is no ambiguity in what "kept" means here: an OFA row for these
assimilated species is apparently never written except when finite and
`assim_flag==1`. If you ever extend this to a species/tile set where that
stops being true, the predicate would need to be revisited explicitly.

## What to do

1. `cd` into this directory.
2. `python3 plot_footprint_case_study.py`
3. It will print, and hard-assert on, three things:
   - all 2,720 candidates have a reconstructable footprint (>=1 tile in the
     search ellipse with valid porosity);
   - the reconstructed QC decision agrees with the actual kept/lost outcome
     for every single candidate (100% agreement — the confusion matrix has
     zero off-diagonal entries);
   - no observation lies on the wrong side of the 0.10 threshold relative to
     its actual outcome.
   **If any of these assertions fails, the script raises and stops** — that
   means either the local floating-point environment behaves differently
   than the reference run, or the bundle was edited. Investigate; don't
   relax the assertion to make it pass.
4. Open `figures/footprint_case_study_threshold.png` — kept obs should sit at ~0.0
   peat fraction, lost obs scattered from ~0.10 to ~0.49, with a clean gap
   straddling the threshold.
5. Open `figures/footprint_case_study_maps.png` — six panels, one per boundary tile:
   local model-tile porosity (color), peat tiles (porosity >= 0.90) ringed
   in red, observation footprint centers marked green (kept) vs black (lost).
6. Open `figures/footprint_case_study_timeseries_3823.png` — tile 3823's full actual
   observation time series (soil moisture vs date, April-September 2015):
   446 candidates in grey/red (441 rejected), 5 retained in green. This is
   the plain-language version of "446 actual observations reduced to five."
7. Open the six maps under `figures/footprint_ellipse_maps/`. These show how moving
   an observation by only a few kilometres changes the contributing model
   tiles and can move the weighted peat fraction across the 0.10 threshold.

## Scope and honest limitations

- This covers only 6 named boundary tiles (3823, 4338, 6548, 9449, 4813,
  5230) and 6 months of data (April-September 2015), not the full
  Canada/Alaska domain (~179 similarly affected tiles exist per the earlier
  domain-wide count). The causal mechanism is proven for these 6 tiles; a
  full-domain reconstruction is future work, not required to establish that
  the effect is real and caused by the footprint peat threshold.
- Do not substitute a different/generic peat map or threshold value if
  extending this analysis — the whole point of this exercise was avoiding
  exactly that shortcut. Extending requires the same BCS porosity file and
  tile geometry pulled fresh from Discover via `extract_footprint_data.py`
  (or an equivalent), not an approximation from public land-cover data.
- If you find an actual disagreement between a local rerun and the
  100%-agreement result described above, treat that as the interesting
  finding — investigate rather than dismiss it.
