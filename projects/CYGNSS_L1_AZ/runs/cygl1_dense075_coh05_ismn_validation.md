# CYGNSS L1 AZ: ISMN in-situ validation (OL vs dense075_coh05)

Purpose: check the project's headline `dense075_coh05` result (own-fit skill improvement
documented in `runs/cygl1_coh05_density_spectrum_and_24mo_result.md`) against independent
in-situ soil moisture, rather than only other GEOSldas observing systems' O-F statistics.

## Station coverage

AZ domain box: lat 29-40, lon -118 to -106 (matches the GEOSldas exeinp
`MINLON`/`MAXLON`/`MINLAT`/`MAXLAT`). ISMN archive:
`/discover/nobackup/projects/land_da/ISMN_data`
(`python_metadata/ISMN_data.csv`, snapshot through 2026-02-15).

24-month window (2020-01-01 to 2021-12-31, matching the `dense075_coh05` headline period):
**189 unique stations** across 5 networks:

| network | stations |
| --- | --- |
| SNOTEL | 126 |
| SCAN | 38 |
| USCRN | 13 |
| iRON | 7 |
| COSMOS | 5 |

No SOILSCAPE stations. 47 SOILSCAPE stations fall inside this box (Jornada/Kendall/Lucky
Hills clusters near the NM border), but the network has a coverage gap there: an early
"node" deployment ran 2015-2017, the current deployment didn't start until 2022 — nothing
was live during 2020-2021.

## Model side

Reused `projects/ascat_da/scripts/run_ismn_ol_da_skill.py` unmodified: its ISMN reader is
network-agnostic and the model-side reader is fully CLI-parameterized (domain, collection,
variable names, tile-matching tolerance, `--nmin`), so no new Python was needed. Differences
from its `ascat_da`/`hsaf_cdr_test` usage:

- Collection: `.tavg24_1d_lnd_Nt.` (single already-daily-averaged file) instead of
  `.SMAP_L4_SM_gph.` — confirmed present for all 365/366 days of both 2020 and 2021 in both
  experiment directories.
- Variables: `SFMC`/`RZMC` instead of `sm_surface`/`sm_rootzone`.
- `--nmin 90` (down from the 1000-day default, which is tuned for multi-year `hsaf_cdr_test`
  runs) since this is a 24-month window.
- Domain is the AZ box, not global. The GEOSldas grid/tilecoord is itself restricted to the
  box, so nearest-tile matching (existing `--max-distance-deg2` logic, unchanged) naturally
  drops any station outside it.

Runs compared:
- `OL` = `OLv8_M36_all_sensors_AZ`
- `DA_dense075_coh05` = `DAv8_M36_AZ_paired_cygl1_dense075_coh05`

## Job

`projects/CYGNSS_L1_AZ/jobs/run_ismn_ol_da_skill_az.sbatch`, account `s3208`. Outputs to
`projects/CYGNSS_L1_AZ/output/ismn_dense075_coh05_ol_da_skill/` (gitignored, per project
convention):

- `cache_obs_daily.nc`, `cache_model_daily_OL.nc`, `cache_model_daily_DA_dense075_coh05.nc`
- `ismn_station_inventory.csv`
- `ismn_skill_stations.csv`
- `ismn_skill_network_summary.csv`

## Results

Job `58298184` completed 2026-09-09 (14 min walltime, clean exit). 1146 ISMN stations had
usable observations in the window; 177 mapped to an AZ-domain tile within the default
tile-matching tolerance (`--max-distance-deg2 0.1`, slightly fewer than the 189 stations
inside the raw lat/lon box, since it also requires proximity to an actual model tile).
After the `--nmin 90` day threshold, SCAN/SNOTEL/USCRN/iRON all cleared it; COSMOS's 5
stations did not (too few paired obs/model days in this window).

**Site counts:** surface 165 sites / 4 networks, root zone 149 sites / 4 networks.

**Global mean skill (reference OL):**

| domain | run | R | anomR | ubRMSE |
| --- | --- | --- | --- | --- |
| surface | OL | 0.5469 | 0.4547 | 0.0682 |
| surface | DA_dense075_coh05 | 0.5358 | 0.4476 | 0.0687 |
| rz | OL | 0.6341 | 0.5571 | 0.0453 |
| rz | DA_dense075_coh05 | 0.6333 | 0.5760 | 0.0456 |

DA is essentially flat-to-slightly-worse than OL against independent ISMN in-situ: surface
R/ubRMSE both move slightly the wrong way (ΔR -0.011, ΔubRMSE +0.0005); root zone is mixed —
R and ubRMSE also very slightly worse, but anomaly correlation improves (ΔanomR +0.019).
Per-network detail (`ismn_skill_network_summary.csv`) shows the same pattern site-group by
site-group, no single network driving it.

This is a much smaller, more neutral signal than the own-fit O-F improvement reported in
`cygl1_coh05_density_spectrum_and_24mo_result.md`, and is consistent with the CygL1 R-sweep
finding (`cygl1_assim_R_sweep.md`) that this observation correlates only weakly with the
model background (r ~ 0.35) — its increments improve the fit to itself without necessarily
improving the fit to the true land state that independent in-situ sensors see.

## Status

Done. Job `58298184`, 2026-09-09.
