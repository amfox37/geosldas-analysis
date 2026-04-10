# ISMN Workflow

This project area documents local ISMN exploration and network-coverage checks used to support model-vs-insitu analysis.

## Latest notebook run order

1. `notebooks/interface.ipynb`
   - Primary ISMN reader walkthrough using `ismn.interface.ISMN_Interface`.
   - Initializes metadata cache (`python_metadata/ISMN_data.csv`) if missing.
   - Used to inspect networks, stations, sensors, and basic data access/filtering.

2. `notebooks/ismn_network_listing.ipynb`
   - Quick local inventory of available networks/stations.
   - Includes focused checks for selected networks/stations and coverage-style summaries.

3. `notebooks/ismn_network_date_ranges_from_metadata.ipynb`
   - Reads metadata CSV and computes earliest/latest data availability by network.
   - Writes:
     - `projects/ISMN/notebooks/ismn_network_timerange_summary.csv`

## Inputs

- Local ISMN dataset root (currently set in notebook cells), typically:
  - `/Users/amfox/Desktop/ISMN_data`
- Metadata CSV produced by `ISMN_Interface`:
  - `/Users/amfox/Desktop/ISMN_data/python_metadata/ISMN_data.csv`

## Notes

- These notebooks currently use local absolute paths in top cells; update those paths for other environments before running.
- If the local ISMN data content changes, delete/regenerate `python_metadata` to avoid stale metadata cache issues.
