# ISMN in-situ soil-moisture validation — methods and provenance

Validation of the ASCAT DA experiments against International Soil Moisture
Network stations, covering OL and three DA runs over one window: the full
experiment span, **2015-04-01 to 2021-03-31**.

Built 2026-07-24. First production run: SLURM job `57326217`.

Code:

| Path | Role |
| --- | --- |
| `projects/ascat_da/lib/ismn_io.py` | direct ISMN `.stm` archive reader |
| `projects/ascat_da/lib/test_ismn_io.py` | tests pinning QC / depth / root-zone rules |
| `projects/ascat_da/scripts/run_ismn_ol_da_skill.py` | driver: observations, model extraction, statistics |
| `projects/ascat_da/jobs/run_ismn_ol_da_skill.sbatch` | production submission |

The workflow derives from `projects/M21C_ls/ISMN_methods_readme.md`; the
section "Departures from the M21C_ls workflow" below records every place it
deliberately differs and why.

---

## 1. Data provenance

### ISMN archive

Root: `/discover/nobackup/projects/land_da/ISMN_data`

| Property | Value |
| --- | --- |
| Snapshot covers data through | **2026-02-15** (from the `.stm` filename end-date, uniform across networks) |
| Archive extracted on Discover | 2026-02-20 (88 of 93 network directory mtimes) |
| `python_metadata/ISMN_data.csv` built | 2026-02-21 16:41 |
| Networks in archive | 93 directories, **88** with soil-moisture sensors |
| Soil-moisture sensor records | **17,606** |
| Distinct stations with soil moisture | **3,279** |

This resolves the open item carried in `projects/M21C_ls/ISMN_methods_readme.md`
("Recover the ISMN archive download/access date, if possible") **for this
archive copy**. Note the M21C_ls figures were produced on a local macOS archive;
this is the Discover copy, so treat the date as authoritative only for work
run from this path.

### Model experiments

All four runs are on `SMAP_EASEv2_M36_GLOBAL`, 112,573 tiles.

| Run label | Root |
| --- | --- |
| `OL` | `/discover/nobackup/projects/land_da/hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121` |
| `DA_H121` | `/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_H121` |
| `DA_legacy` | `/discover/nobackup/projects/land_da/hsaf_cdr_test/DAv7_M36_ASCAT_type_13_legacy` |
| `DA_SMAP_comb_fp_scaled` | `/discover/nobackup/amfox/Experiments/DAv7_M36_SMAP_type_13_comb_fp_scaled/DAv7_M36_SMAP_type_13_comb_fp_scaled` |

Daily `SMAP_L4_SM_gph` file counts over the span: OL 2,192; DA_H121 2,688;
DA_legacy 2,440; DA_SMAP_comb_fp_scaled 2,192. The window itself is 2,191 days,
so every run has full coverage (the larger counts extend past the window).

---

## 2. Reading ISMN without the `ismn` package

The `ismn` Python package is **not installed in any GEOSpyD environment on
Discover** (checked across all `/usr/local/other/GEOSpyD/*/*/envs/*`), so the
M21C_ls approach of driving `ISMN_Interface` is not available here.

`lib/ismn_io.py` therefore parses the archive directly:

- **Inventory** from `python_metadata/ISMN_data.csv`. It has a two-row header;
  columns are `(name, sub)` pairs. The workflow uses `network`, `station`,
  `latitude`, `longitude`, `variable` (+ `depth_from` / `depth_to`),
  `timerange_from` / `timerange_to`, and `file_path` — the last giving the
  archive-relative path to each sensor's `.stm` file, which sidesteps the
  station-name mangling in the `.stm` headers (spaces become underscores).
- **Sensor series** from the `.stm` "header_values" format: one header line,
  then `YYYY/MM/DD HH:MM value ismn_flag provider_flag`.

Verified: all 8,046 window-overlapping sensor files resolve on disk (0 missing).

This is a dependency reduction, not a workaround — it also runs in batch with
no `TkAgg` backend hazard, which the M21C_ls notebook has to defend against
with a lazy import.

---

## 3. Station and sensor selection

No predetermined network list. The funnel, with counts as of the first
production run:

| Stage | Result |
| --- | --- |
| Soil-moisture sensors in archive | 17,606 across 88 networks / 3,279 stations |
| Sensor time range overlaps 2015-04-01..2021-03-31 | 8,046 sensors → **2,104 stations across 60 networks** |
| Has ≥1 sensor yielding QC-passing data inside the window | station-dependent |
| Nearest M36 tile within 0.1 deg² | station-dependent |
| ≥ `--nmin` (1000) paired obs/model days, **in every run**, per domain | final site set |

Sensors with a missing metadata time bound are **kept**, not dropped — the
`.stm` file is treated as the authority and the metadata only as a pre-filter.

The final gate is applied per domain (surface / root zone) and requires all
four runs to clear it, so a station never contributes to one run's mean and not
another's. This replaces M21C_ls's cross-window common-site intersection, which
has no meaning in a single-window design.

---

## 4. Observation construction

**Quality control.** Only records flagged `G` ("good") are retained
(`--required-flag`). As in M21C_ls, this is what excludes frozen-soil and other
geophysically flagged records; **no independent snow, frozen-soil, or
snow-cover mask is applied.**

**Replicate sensors.** Multiple sensors at the same nominal depth (e.g. SCAN's
Hydraprobe A/B/C) are averaged into one series per depth. Duplicate timestamps
within a series are averaged first.

**Surface.** The sensor depth nearest `--target-surface-depth-m` (0.05) **among
sensors no deeper than `--surface-max-depth-m` (0.10)**. A station whose
shallowest sensor is deeper than 0.10 m yields *no* surface series.

> This depth cap is the single most important departure from M21C_ls, which
> took the nearest depth unconditionally. That was safe across its six curated
> networks; across all 88 it would silently promote a 0.5 m or 1.0 m sensor to
> "surface soil moisture." `test_build_station_series_rejects_deep_only_station`
> pins this behavior.

**Root zone.** Profile-weighted average over `[--rz-top-m, --rz-bottom-m]` =
0–1 m. Each depth is weighted by the layer thickness it represents (edges
midway between neighbouring sensors, clamped to the profile bounds). A
timestamp produces a value only when it has:

- at least `--rz-min-sensors` (3) finite layers, **and**
- a finite shallow sensor ≤ `--rz-shallow-max-m` (0.20), **and**
- a finite deep sensor ≥ `--rz-deep-min-m` (0.50).

so a partial profile is never passed off as root zone. Root-zone site counts
are consequently much lower than surface and are reported separately — never
combine the two counts.

**Daily alignment.** Sub-daily records are trimmed to the window as they are
read (ISMN histories run decades longer than any experiment), averaged to daily
means, then shifted `--daily-shift-hours` (12) to sit at local noon, matching
the model day.

---

## 5. Model extraction

**Collection: `SMAP_L4_SM_gph`, variables `sm_surface` / `sm_rootzone`,
averaged over each file's eight 3-hourly instants — for all four runs.**

Rationale:

1. `DAv7_M36_SMAP_type_13_comb_fp_scaled` has **no `tavg24_1d_lnd_Nt`
   collection at all**; its only daily collections are `SMAP_L4_SM_gph`,
   `inst3_1d_lndfcstana_Nt`, and `catch_progn_incr`.
2. `SMAP_L4_SM_gph` is the analysis state (post-update), which is the intended
   quantity for DA evaluation.
3. All four runs carry it, so using it everywhere keeps them directly
   comparable. Mixing `tavg24` for three runs and `gph` for the fourth would
   confound a collection change with the experiment difference.

`inst3_1d_lndfcstana_Nt` was considered and rejected: it splits into explicit
`SFMC_FCST` / `SFMC_ANA`, so it *could* supply the analysis, but only three of
the four runs would then match on collection semantics.

**Fill handling.** Values `> 1e14` become NaN before averaging. Exact `0.0`
after tile subsetting is treated as missing, matching the legacy skill scripts
(`sm_skill_vs_insitu.py` lineage).

**Station → tile matching.** Nearest tile by squared lat/lon degree distance,
rejected beyond `--max-distance-deg2` (0.1). Computed per run so a run with a
different tile space cannot silently misalign; the resulting `tile_index` and
`tile_distance_deg2` are stored in each run's cache. A stale cache built
against a different station set is caught by an explicit station-axis check.

### Tilecoord note — filesystem change outside the repo

`DAv7_M36_SMAP_type_13_comb_fp_scaled`'s extracted output tree contains only
`ana/` and `cat/` — **no `rc_out/`, therefore no `.ldas_tilecoord.bin`**, which
every tile-space workflow needs. (The 256 GB source tarball does contain a
nested per-restart `rc_out/Y*/M*/`, but not the flat top-level tilecoord.)

On 2026-07-24, with the user confirming the tile spaces are identical, the
tilecoord was copied in:

```
cp .../DAv7_M36_ASCAT_type_13_legacy/output/SMAP_EASEv2_M36_GLOBAL/rc_out/DAv7_M36_ASCAT_type_13_legacy.ldas_tilecoord.bin \
   .../DAv7_M36_SMAP_type_13_comb_fp_scaled/output/SMAP_EASEv2_M36_GLOBAL/rc_out/DAv7_M36_SMAP_type_13_comb_fp_scaled.ldas_tilecoord.bin
```

Supporting evidence: the three `hsaf_cdr_test` DAv7/OLv7 tilecoords are
**byte-identical** (all md5 `928999c17503eee32685c470e8c09aba`), and the SMAP
run's own output carries the same 112,573-tile dimension. This same copied file
is also what the `iv_tc` step-2 pairs for that run depend on.

---

## 6. Statistics

Computed by the shared helpers in
`projects/matlab2python/scripts/sm_skill_vs_insitu.py`, so this workflow stays
consistent with the rest of the repo's in-situ skill lineage.

Per station × domain × run, from pairwise-complete obs/model samples:

| Metric | Definition |
| --- | --- |
| `N_pairs` | paired finite obs/model days |
| `R` | Pearson correlation, CI via lag-1-autocorrelation-adjusted effective sample size |
| `anomR` | correlation of anomalies |
| `bias` | mean(model − obs) |
| `RMSE`, `ubRMSE` | `sqrt(MSE)`, `sqrt(MSE − bias²)` |

**Anomalies** remove a day-of-year climatology built with a 31-day circular
window, requiring `--nmin-day` (30) samples per day-of-year. Over six years a
31-day window holds up to ~186 samples, so this threshold is not binding here.
Anomalies are computed on each run's own paired record, matching M21C_ls.

**Network summary** (`ismn_skill_network_summary.csv`) reports per-network,
per-domain means and paired deltas against `--reference-run` (default `OL`),
with standard error. Deltas are signed so **positive always means the run beat
the reference**:

- `R`, `anomR`: `run − reference`
- `RMSE`, `ubRMSE`: `reference − run`
- `bias`: `|reference| − |run|` (magnitude, since either sign is a miss)

---

## 7. Outputs

Default `--output-dir`:
`/discover/nobackup/projects/land_da/Evaluation/ISMN/ascat_da_ol_da_skill`

| File | Contents |
| --- | --- |
| `cache_obs_daily.nc` | daily ISMN surface/root-zone series, `(time, station)` |
| `cache_model_daily_<run>.nc` | daily model series + station→tile mapping per run |
| `ismn_station_inventory.csv` | per-station metadata, chosen surface depth, obs-day counts |
| `ismn_skill_stations.csv` | per station/domain/run statistics |
| `ismn_skill_network_summary.csv` | per network/domain means and deltas vs the reference run |

Both caches are reused on rerun; `--overwrite-obs` / `--overwrite-model` force
a rebuild. Re-scoring with a different `--nmin` is therefore cheap — the
expensive stages are cached.

---

## 8. Running it

```bash
cd /gpfsm/dnb06/projects/p284/geosldas-analysis/projects/ascat_da/jobs
sbatch run_ismn_ol_da_skill.sbatch
```

Useful flags for narrower runs: `--network SCAN` (repeatable), `--max-stations N`
(debug cap), `--nmin`, `--surface-max-depth-m`, `--reference-run`.

Observed timings (job `57326217`): observation stage ~250 stations/2 min, so
~17 min for 2,104; model extraction ~0.3 s/day/run, so ~45 min for four runs ×
2,191 days; statistics ~15 min. Roughly 1.5 h total against a 6 h walltime.
Memory stays well under the requested 16 GB.

---

## 9. Departures from the M21C_ls workflow

| Topic | M21C_ls | Here | Why |
| --- | --- | --- | --- |
| Networks | 6 curated (SNOTEL, SCAN, USCRN, SMOSMANIA, OZNET, ARM) | every network in the archive | requested; nothing in the pipeline required a fixed list |
| Windows | 3 (pre-ASCAT / pre-SMAP / SMAP-era) | 1 (full experiment span) | matches the experiment design |
| Common-site rule | intersection across windows | intersection across **runs** | single window; still prevents a shifting site population |
| Surface depth | nearest to 0.05 m, unbounded | nearest to 0.05 m among sensors ≤0.10 m | unbounded search mislabels deep-only stations once all networks are in scope |
| Root zone | `matlab_strict` per-network layer composites | generic profile-weighted 0–1 m | strict rules are hand-tuned for 6 networks and do not generalize |
| ISMN access | `ismn` package v1.5.2 (`ISMN_Interface`) | direct `.stm` + metadata CSV parsing | package unavailable on Discover |
| Model source | `tavg24_1d_lnd_Nt` `SFMC`/`RZMC` | `SMAP_L4_SM_gph` `sm_surface`/`sm_rootzone`, 8-instant mean | the SMAP run has no `tavg24`; gph is the analysis and is common to all four runs |
| Execution | notebook | script + SLURM, cached stages | 2,104 stations × 4 runs is not notebook work |

---

## 10. Known limitations and open items

- **Root-zone comparability across networks.** The generic 0–1 m composite
  weights whatever depths a station has. Two stations with different sensor
  layouts produce root-zone series representing slightly different effective
  profiles. This is acceptable for OL-vs-DA *differences* at a fixed station
  (the layout is constant across runs) but means absolute root-zone skill is
  not strictly comparable between networks.
- **Not MATLAB-parity.** Unlike the `iv_tc` engines, nothing here has been
  checked against a live MATLAB reference; the M21C_ls notebook it derives from
  was itself Python. Treat the numbers as internally consistent, not as
  reproducing a legacy MATLAB product.
- **No snow / frozen-soil mask** beyond the ISMN `G` flag, same caveat as
  M21C_ls. Relevant if reviewers ask about cold-season skill.
- **Networks are unweighted in the global mean.** `report()` prints a plain
  mean over stations, so station-dense networks (SNOTEL, SCAN) dominate.
  Per-network means in `ismn_skill_network_summary.csv` are the safer basis for
  conclusions; a network-weighted global aggregate is not yet implemented.
- **Surface depth cap of 0.10 m is a judgement call.** It admits ISMN's usual
  0–0.05 and 0–0.10 surface layers and excludes 0.20 m and deeper. Stations
  whose only shallow sensor sits at, say, 0.15 m are dropped from surface
  scoring rather than included with a caveat.
