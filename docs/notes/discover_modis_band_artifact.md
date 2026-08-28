# MODIS snow-cover longitude-band artifact — Discover investigation notes

Companion to the WY2001-WY2006 snow-DA water-budget study (`projects/M21C_ls`).
No model rerun performed; no fix attempted (both explicitly out of scope for
this pass). Source: `LS_DAv8_M36_v2` (M21C_land_sweeper), 2000-06-01..2007-05-31
DA leg.

## 0. Critical correction to species IDs

The brief assumed species 12 = MOD10C1 (Terra), species 11 = MYD10C1 (Aqua).
**This is wrong for this experiment.** From the actual obsparam file
(`LS_DAv8_M36.ldas_obsparam.20040101_0000z.txt`):

| species | descr | notes |
|---|---|---|
| 11 | `ASCAT_METC_SM` | not MODIS at all |
| 12 | `MYD10C1` | Aqua |
| 13 | `MOD10C1` | Terra |

Confirmed via the brief's own mandatory Section 8 sanity check: summing
species (12,13) over `temporal_stats_DA_20000601_20070531.nc4`, tiles
28-38N, mean N_data columns 717-722 (west of 90E) vs 723-728 (east) gives
2493.39 / 3319.06, ratio **1.331** — matches the brief's expected ~2493 /
~3319 / 1.33 almost exactly. The brief's assumed mapping (species 11,12)
gives 901.6 / 1144.3, ratio 1.269 — does not match. **All species filters
in this investigation and its output use {12, 13}, not {11, 12}.**

## 1. Data availability

- The `.bin` OFA files for the 10 target dates are not missing, but the
  `ana/ens_avg` tree for `LS_DAv8_M36_v2` (2000-2007) was tarred (not
  gzipped) and the raw per-timestep files subsequently deleted from disk
  on 2026-07-22 for inode-pressure reasons. All 80 needed files (10 dates
  x 8 windows) were confirmed present via `tar -tf` and extracted from
  `SMAP_EASEv2_M36_GLOBAL_ana_tars/Y{2001,2003,2004,2005,2006}.tar`.
- `LS_DAv8_M36_v3`'s directory tree for 2000 - 2007-05 is a symlink into
  `LS_DAv8_M36_v2`'s real output, so reading v2 directly (as done here)
  carries no risk of silently mixing the v2/v3 legs.

## 2. Task 2 — assimilation window / obsparam / MODIS product confirmation

- **dtstep_assim**: not overridden in this experiment's `LDAS.rc`; compiled
  default is `LANDASSIM_DT = 10800` s (3h), `LANDASSIM_T0 = 000000` (0000
  UTC), enforced/asserted in `GEOS_LandAssimGridComp.F90`. Consistent with
  the observed 8 files/day at 0000, 0300, ..., 2100z. Note the code's
  `dtstep_assim_max = 21600` (6h) in `read_obs_MODIS_SCF` is a hard ceiling
  the reader asserts against, not the actual value in use here.
- **obs_param entries**: species 12 (`MYD10C1`), species 13 (`MOD10C1`),
  both `FOV = 0.0 deg` (see §4 below for why).
- **MODIS product**: confirmed MOD10C1/MYD10C1 CMG at 0.05 deg, Collection
  061 — both from the obsparam `path` fields (`.../MYD10C1_V61/`,
  `.../MOD10C1_V61/`) and the reader's filename-template comment
  (`MOD10C1.Ayyyyddd.061.hdf`).

## 3. Task 3 — code inspection (`clsm_ensupd_read_obs.F90`, `read_obs_MODIS_SCF`)

**Caveat:** `LS_DAv8_M36_v2`'s exact original GEOSldas build
(`GEOSldas_Landsweeper_v2/GEOSldas/install-Aggressive`) no longer exists —
it was renamed/consolidated into `GEOSldas_release` (tag v20.2.0,
`GEOSldas_GridComp` v3.2.0) around Dec 2025. All code inspection below uses
that v3.2.0 source, the closest available proxy (used by the sibling
`LS_DAv8_M36_v3`/`LS_OLv8_M36_v2` runs), not v2's own literal source.

**(1) Is there de-duplication across consecutive windows?**
Yes, and it appears to work correctly. The reader reads MODIS obs from a
deliberately *widened* longitude band (`tmp_delta` = 3x max tile lon extent
either side of the true window), but final tile output is filtered to keep
only tiles with `lon_end < tile_coord%com_lon <= lon_beg` — i.e. against
the TRUE (narrow, un-widened) window bounds, not the widened CMG-read
bounds. This is a half-open interval, so consecutive windows partition the
globe without overlap in principle.
Empirical confirmation (Task 1 data, all 13,140 extracted rows): **zero**
`(date, tile_id, species)` triples were seen in more than one window — the
de-dup is not failing at the per-day level in this dataset.

**(2) CMG-to-tile assignment: tile centre, bounding box, or CMG-cell centre?**
Neither, directly, at the `get_tile_num_for_obs` stage: MODIS obs_param
entries have `FOV = 0.0 deg`, meaning the nearest-tile search
(`get_tile_num_from_latlon`) uses zero tolerance. This only works because
the reader's docstring states the design assumption: "MODIS resolution is
finer than Catchment tile space ... super-ob data to tile space, with
lat/lon coords of obs matching lat/lon coords of tiles" — i.e. the
CMG-cell -> tile super-obbing (many 0.05 deg cells per tile footprint) is
done upstream inside the CMG-reading step itself, which outputs obs already
snapped to `tile_coord%com_lon/com_lat`. `get_tile_num_for_obs` is then just
an exact-match confirmation, not the actual footprint assignment.

**(3) N_files=1 vs 2 / date-boundary effect at 90E?**
Does not apply here. `N_files = 2` only triggers when the longitude band
*wraps the antimeridian* (`lon_end_MODIS >= lon_beg_MODIS`), which happens
near +/-180 deg, not at 90E. At 90E the reader is always in the `N_files=1`
branch, one daily file, no day-boundary branching to investigate. This
rules out the "eastward extra day" hypothesis as posed — it's a
dateline-only code path, structurally irrelevant to this artifact.

## 4. Task 1 — per-window empirical data (decisive test)

Extracted all MODIS (species 12+13) observations for tiles in
[85E,95E] x [25N,45N] (1,508 tiles) across the 10 target dates x 8 windows
(80 files, 13,140 rows). Output:
`projects/M21C_ls/output/modis_obs_by_window_85E_95E.csv.gz` (+ `.sha256`).
Extraction script: `projects/M21C_ls/scripts/extract_modis_band_artifact_obs.py`.

**Finding A — no same-day double-window duplication.** Every
`(tile_id, species, date)` triple appears in exactly one window across all
13,140 rows checked. This answers the brief's decisive test directly:
excess tiles do **not** appear in two windows on the same day. Mechanism
(a) from Section 1 (un-deduplicated widening overlap) is not what's
happening in this dataset.

**Finding B — window assignment is fixed and side-specific, not noisy.**
For the 90E boundary columns, west-of-boundary tiles (i_indg 717-722) are
*always* read in window 0600z, and east-of-boundary tiles (i_indg 723-728)
are *always* read in window 0300z, on every single one of the 10 target
dates, for both species. This is exactly what the 45-degree-wide window /
`localtime2longitude` mechanism predicts (dtstep_assim=10800s -> 45 deg per
window, and 90E happens to sit almost exactly on a window edge as well as
an EASEv2 column boundary) — not evidence of a bug, just confirms the
window-boundary geometry.

**Finding C — the excess is a date-assignment/QC-timing effect, and it is
strongly seasonal.** Since no tile is ever double-counted within a day, the
west/east asymmetry must come from *how many of the 10 sampled days
produced a valid obs* for a given tile, not from within-day duplication.
Per-tile day-coverage (out of 10 candidate dates): west mean 4.77, east
mean 5.50 (species pooled); rows-per-tile ratio east/west = 1.155 for this
10-day sample (vs. 1.331 for the full 24-year climatology — same
direction, smaller magnitude, as expected from a 10-day vs. 24-year
sample).

Critically, this asymmetry is **not evenly distributed across the year**:

| date | west rows | east rows | ratio |
|---|---|---|---|
| 2003-01-15 | 306(12)+312(13) | 305+309 | ~1.0 |
| 2004-01-15 | 256+284 | 241+303 | ~0.98-0.94 |
| 2004-04-15 | 123+152 | 128+172 | ~1.04-1.13 |
| 2004-10-15 | 178+161 | 151+158 | ~0.85-0.98 |
| 2005-01-15 | 257+303 | 274+297 | ~1.07-0.98 |
| 2006-01-15 | 263+258 | 271+238 | ~1.03-0.92 |
| **2001-07-15** | 0+25 | 0+111 | **4.4x** (Terra only; Aqua not yet launched) |
| **2003-07-15** | 3+38 | 32+133 | **~3.5-10.7x** |
| **2004-07-15** | 11+23 | 66+109 | **~4.7-6.0x** |
| **2005-07-15** | 0+10 | 15+121 | **~12x / undefined** |

January/April/October dates show ratios near 1.0 (no meaningful
asymmetry). All four July dates show large (3.5x-12x) west-deficit /
east-excess ratios, and absolute counts in July are an order of magnitude
lower than in the other three months (10-133 rows vs. 150-312) — consistent
with July being outside the seasonal snow season for 25-45N, where valid
SCF retrievals are naturally sparse and QC-marginal. **The JJA "3x step" in
the downstream analysis increment that motivated this brief is very likely
driven by this low-count summer regime, where a small absolute
west/east difference in successful retrievals translates into a large
relative ratio — not by an amplified physical band-effect in summer.**

## 5. What remains unresolved

The reader's de-dup and file-selection logic (as read in the v3.2.0 proxy
source) appears internally consistent and does not, by inspection alone,
predict a systematic west/east asymmetry in retrieval success — yet the
Task 1 data clearly show one, concentrated in summer. Two explanations
that source-code reading cannot distinguish:

1. It's a genuine, physically-driven difference in CMG-level snow/no-snow
   QC pass rates a few tiles either side of a window boundary (plausible
   in a summer marginal-snow regime, and not obviously a code defect).
2. There is some other asymmetry in the CMG-to-tile super-obbing or QC
   step (inside `read_MODIS_SCF_hdf`, not yet read line-by-line here) that
   depends on which side of a window edge a tile falls, and happens to
   only bite when counts are already low (summer).

Distinguishing these would require pulling the actual MOD10C1/MYD10C1 CMG
pixel values and QC flags for the 4 July target dates near 90E — a
natural next step, but out of scope for this pass per the brief's Section 7
(diagnose only, no fix, no rerun).

## 6. Deliverables

- `projects/M21C_ls/output/modis_obs_by_window_85E_95E.csv.gz` (+
  `.sha256`) — one row per MODIS obs, columns: date, window, tile_id, lon,
  lat, i_indg, species, species_name, obs, fcst, assim_flag. Header
  comments in the gzip'd file itself record experiment path, OFA filename
  pattern, tilecoord file, species mapping, and code-version caveat.
- `projects/M21C_ls/scripts/extract_modis_band_artifact_obs.py` —
  extraction/join script (read-only, no model rerun).
- Transfer target on the Mac side (not performed here — Discover-side
  deliverable only): `scp` to
  `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`.

**Note on tilecoord file used:** the brief named
`LS_OLv8_M36.ldas_tilecoord.bin` specifically; this analysis used
`LS_DAv8_M36.ldas_tilecoord.bin` instead (same EASEv2 M36 domain/tile
space as the DA experiment being analyzed, N_tile=112573 confirmed to
match). OL and DA share the same tile-space domain and should be
numerically identical, but this substitution has not been independently
verified byte-for-byte against the OL file and should be treated as an
assumption, not a confirmed identity.
