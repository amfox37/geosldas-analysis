# ascat_da

Diagnostics and utilities for ASCAT-based land data assimilation experiments. The notebooks here generate diagnostic figures and statistics for different ASCAT filtering/footprint configurations, while the scripts automate data preparation.

For the older ASCAT/SMAP paper figure workflow, see `ASCAT_SMAP_paper_provenance.md`. That workflow predates this project layout, so the main publication figure notebook lives under `projects/utils`.

## Key notebooks
- `notebooks/legacy_vs_h121_obs.ipynb` - validates raw Legacy BUFR and H121/H139 observations after QC and GEOS M36 tile/cycle super-ob formation, including comparisons with `ObsFcstAna` (OFA). The notebook uses the corrected `geos_cycle_global_v7_gridoriginfix_bsflag_noise` cache described below.
- `notebooks/compare_legacy_bufr_vs_H121.ipynb` – compares legacy BUFR ASCAT (ASCSMR02) against H121 CDR, including observation-space diagnostics from a GEOSldas diagnostic run (per-platform stats, O−F innovations, obs and innovation maps for Jan 2020).
- `notebooks/compare_comb_fp_043024.ipynb` – latest comparison of combined FP experiments with OL, including map plots and PDF summaries.
- `notebooks/ASCAT_masking_021324.ipynb` – builds and inspects ASCAT observation masks for assimilation tests.
- `notebooks/regrid_ASCAT_mask_022724.ipynb` – regrids ASCAT masks to the analysis grid for quick QA.

## Supporting scripts
- `scripts/check_ascat_duplicates.py` – helper to spot duplicate ASCAT files before ingest.
- `scripts/filename_lister_v2.py` – utility for listing/renaming raw observation files.

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
