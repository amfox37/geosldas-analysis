# snow_da

Snow data assimilation diagnostics and figure notebooks. Focuses on LS snow DA experiments, QC of new snow observations, and experiment-to-experiment comparisons.

## Key notebooks
- `notebooks/analyze_new_snow_qc_051525.ipynb` – evaluates new snow QC logic with ObsFcstAna tiles.
- `notebooks/analyze_snowDA_061925.ipynb` – most recent snow DA analysis package (global maps and stats).
- `notebooks/plot_LS_snowda_021825.ipynb` – plotting utility for LS snow DA experiments across multiple regions.

## VIIRS snow cover reader test (2026-11)
- `docs/viirs_scf_reader_test_report.md` – test of the GEOSldas VIIRS SCF reader (`read_obs_JPSS_SCF`, species 59-61) on
  CF0360: reader checks, QC fix, VIIRS vs MODIS Aqua obs, and DA vs open-loop forecast skill (Discover, Jan 2020).
- `scripts/viirs_scf_reader_test/` – obs download, log and ObsFcstAna checks, comparisons, report figures (GEOSpyD py3.13).
- `runs/viirs_scf_reader_test/` – exeinp, namelists and build scripts (`configs/`) and check outputs (`results/`).
