# CYGNSS L1 AZ: ISMN in-situ validation of the fixed-operator arms, 2020–2022

*2026-09-28, rewritten the same day.* The first version averaged over every ISMN station in the GEOSldas box. Most of those stations
sit where the CYGNSS L1 stream never reaches, so it understated L1 and was not a fair L1-vs-L3 comparison. This version scores only the
**area where L1 obs exist**. The all-station numbers are kept in the appendix for reference. This is the post-fix follow-up to
`cygl1_dense075_coh05_ismn_validation.md`, which was pre-fix and is invalid.

## Setup

- **Driver and job:** `projects/ascat_da/scripts/run_ismn_ol_da_skill.py`, run by `jobs/run_ismn_ol_da_skill_az.sbatch` (job 58621216,
  22 min).
  - Model data: daily `tavg24_1d_lnd_Nt` SFMC/RZMC, all 1096 days found for every run.
  - Settings: `--nmin 90`, `--max-distance-deg2 0.1`.
  - Station pool: every ISMN station in the GEOSldas box (118–106° W, 29–40° N) matched to a model tile. That gives 195 surface and
    167 rz stations, with station-to-tile-center offsets of 0.15° median and 0.31° max.
  - Skill is computed per station, so restricting to a subset of stations needs no rerun.
  - Output: `output/ismn_fixedop_20200101_20221231/` (gitignored).
- **Runs:**

  | label | experiment |
  |---|---|
  | OL (reference) | `OLv8_M36_AZ_fixedop` |
  | L1coh | `DA_L1_full_xc015_coh040216_fixedop` |
  | L1coh_err39 | `DA_L1_full_xc015_coh040216_err39_fixedop` |
  | L3 | `DA_L3_fixedop` |
  | SMAP | `DA_SMAP_fixedop` |

- **L1 area = tiles with ≥ 100 CYGNSS L1 obs over 2020–2022 in the unfiltered L1 stream.** The counts are the `N_data` of species 12
  in the L3-only arm's 2020–2022 temporal stats, where L1 is monitor-only and unfiltered. The same mask is applied to every run.
  - 574 of the 909 domain tiles have any L1 obs, with a median of 785 per tile.
  - The L1 stream is limited to specular points within 200 km of Arizona (the CYGNSS_operator preprocessing region) and south of about
    37.9° N (the CYGNSS orbit). It also runs thin above about 37.4° N.
- **Why this matters:** only 43 of the 195 surface stations are on a tile with any L1 obs. The median station latitude is 37.8° N; half
  of the stations are Utah/Colorado mountain sites (mostly SNOTEL) beyond CYGNSS coverage. Breakdown of the 195 surface stations:

  | location | stations | on a tile with L1 obs |
  |---|---:|---:|
  | Arizona + 200 km, south of 37.4° N | 45 | 41 |
  | Arizona + 200 km, 37.4° N or further north | 54 | 0 |
  | outside Arizona + 200 km | 96 | 2 |

- **Stations scored:**
  - Surface: 42 stations on 25 tiles (SNOTEL 14, SCAN 12, SOILSCAPE 9, USCRN 7).
  - Root zone: 33 stations on 19 tiles (SNOTEL 13, SCAN 11, SOILSCAPE 6, USCRN 3).
  - The median is about 600 L1 obs per tile.
- **Statistics:**
  - Each Δ is the mean over stations of the paired difference.
  - 95% CIs are from a **tile-cluster bootstrap**: tiles are resampled with all their stations, because stations on the same tile share
    one model time series.
  - \* marks a CI that excludes 0. "better" is the share of stations that improve.

## Mean skill (L1 area)

| domain | run | R | anomR | ubRMSE |
|---|---|---:|---:|---:|
| surface | OL | 0.521 | 0.444 | 0.0515 |
| | L1coh | 0.528 | 0.450 | 0.0512 |
| | L1coh_err39 | 0.529 | 0.449 | 0.0511 |
| | L3 | 0.525 | 0.437 | 0.0512 |
| | SMAP | **0.588** | **0.532** | **0.0490** |
| rz | OL | 0.468 | 0.491 | 0.0354 |
| | L1coh | 0.511 | 0.518 | 0.0349 |
| | L1coh_err39 | 0.497 | 0.510 | 0.0351 |
| | L3 | 0.479 | 0.495 | 0.0354 |
| | SMAP | **0.534** | **0.568** | **0.0343** |

## Change vs OL (L1 area)

| domain | metric | L1coh | L1coh_err39 | L3 | SMAP |
|---|---|---|---|---|---|
| surface | ΔR | +.007 [−.008, +.020] 62% | +.008 [−.002, +.016] 69% | +.005 [−.010, +.018] 55% | **+.067 [+.033, +.097]\* 86%** |
| surface | ΔanomR | +.006 [−.010, +.022] 63% | +.005 [−.006, +.017] 68% | −.007 [−.018, +.002] 46% | **+.088 [+.060, +.110]\* 90%** |
| surface | ΔubRMSE ×10⁻³ | −.27 [−.61, +.12] 50% | **−.31 [−.54, −.02]\* 67%** | −.25 [−.67, +.24] 60% | **−2.42 [−3.62, −.88]\* 81%** |
| rz | ΔR | **+.044 [+.019, +.067]\* 88%** | **+.030 [+.013, +.044]\* 91%** | +.012 [−.018, +.037] 55% | **+.066 [+.034, +.112]\* 76%** |
| rz | ΔanomR | **+.027 [+.005, +.063]\* 71%** | **+.019 [+.004, +.040]\* 71%** | +.003 [−.021, +.027] 61% | **+.077 [+.028, +.141]\* 74%** |
| rz | ΔubRMSE ×10⁻³ | **−.53 [−.79, −.22]\* 88%** | **−.36 [−.51, −.18]\* 85%** | −.07 [−.36, +.25] 61% | **−1.18 [−1.69, −.70]\* 82%** |

## L1 vs L3 directly (L1 area)

| domain | metric | L1coh − L3 | L1coh_err39 − L3 |
|---|---|---|---|
| surface | ΔR | +.002 [−.009, +.020] 36% | +.003 [−.008, +.019] 45% |
| surface | ΔanomR | **+.012 [+.002, +.027]\* 61%** | **+.012 [+.003, +.024]\* 63%** |
| surface | ΔubRMSE ×10⁻³ | −.03 [−.54, +.33] 38% | −.06 [−.52, +.29] 48% |
| rz | ΔR | **+.032 [+.011, +.065]\* 70%** | +.018 [−.001, +.046] 64% |
| rz | ΔanomR | +.024 [−.005, +.068] 52% | +.015 [−.010, +.047] 52% |
| rz | ΔubRMSE ×10⁻³ | **−.46 [−.73, −.26]\* 82%** | **−.29 [−.56, −.09]\* 70%** |

## Reading

1. **Where L1 obs exist, L1 assimilation clearly improves the root zone against in-situ data.** rz R rises by +0.03 to +0.04 over the OL
   and rz anomR by +0.02 to +0.03, at 70–90% of stations, with tile-cluster CIs excluding 0. The surface gain is smaller (+0.005 to
   +0.008) and not significant.
2. **L3 assimilation does not significantly improve anything against in-situ data in the same area,** and its surface anomR is slightly
   negative (−0.007).
3. **L1 beats L3 in situ:** surface anomR by +0.012 (both L1 arms, significant), and rz R / ubRMSE for the L1coh arm (+0.032, significant).
   **This is the opposite of the O-F ranking,** where L3-only beat the best L1 arm on every monitor:
   - The L3-only arm's L3-monitor score is own-fit, since L3 is assimilated there.
   - On the fully independent Tb monitors, L3's O-F lead is only about 0.2–0.3 points.
   - The large ASCAT O-F advantage of L3 (−3.6% vs +0.6%) is not reflected in situ. One possibility, not tested: L3 and ASCAT are both
     SM retrievals with shared error structure.
4. **Within the L1 arms, errstd 2.75 (L1coh) does better in situ than errstd 3.9** (rz R +0.044 vs +0.030), while O-F slightly preferred
   errstd 3.9 because of ASCAT. The larger R costs in-situ root-zone skill.
5. **SMAP-only is still clearly the strongest,** at both surface and root zone.
6. **Caveats:**
   - The sample is small (19–25 independent tiles).
   - The station mix is point-scale against 36 km tiles.
   - 2022 SOILSCAPE is one year on 3 tiles (below).
   - Next checks: O-F scored over the same L1-area tiles (does L3's O-F lead survive there?), and a seasonal split of the in-situ skill
     (does the spring O-F degradation appear in situ?).

## SOILSCAPE

Of particular interest because the SOILSCAPE instruments are run by the group developing the L1 operator. All SOILSCAPE data in the
window are from **2022 only** (the current deployment), with 115–242 paired days per station.

**The sites are clustered into very few model tiles, so they are not independent samples.** Stations in one M36 tile are all compared
with the same model time series, so the effective sample size is the number of tiles:

| tile | tile center | sites | surface / rz stations | direct CYGNSS L1 obs on tile, 2020–22 (unfiltered stream) |
|---|---|---|---|---|
| 591 | 31.62 N, 109.98 W | Walnut Gulch: Kendall, Lucky Hills | 9 / 6 | 1,303 (filtered arms: 437 / 434 / 323 per year) |
| 867 | 32.62 N, 106.62 W | Jornada: JR-1, JR-2, JR-3 | 10 / 0 | **0** |
| 885 | 37.08 N, 106.24 W | CO-Z1 | 4 / 0 | **0** |

Station-bootstrap CIs over SOILSCAPE (as used in the sections above) are therefore not meaningful. The surface is effectively n = 3, and
the root zone is n = 1 (tile 591).

Mean over stations:

| domain | run | R | anomR | ubRMSE |
|---|---|---:|---:|---:|
| surface (23 stations, 3 tiles) | OL | 0.649 | 0.530 | 0.0303 |
| | L1coh | 0.661 | 0.539 | 0.0299 |
| | L1coh_err39 | 0.658 | 0.535 | 0.0300 |
| | L3 | 0.666 | 0.528 | 0.0296 |
| | SMAP | **0.716** | **0.590** | **0.0283** |
| rz (6 stations, tile 591) | OL | 0.669 | 0.623 | 0.0270 |
| | L1coh | **0.766** | 0.630 | **0.0259** |
| | L1coh_err39 | 0.733 | 0.626 | 0.0263 |
| | L3 | 0.749 | 0.640 | 0.0263 |
| | SMAP | 0.708 | **0.652** | **0.0259** |

Change vs OL by tile (mean over the tile's stations):

| tile | metric | L1coh | L1coh_err39 | L3 | SMAP |
|---|---|---:|---:|---:|---:|
| 591 Walnut Gulch | sfc ΔR | +.022 | +.019 | +.040 | +.136 |
| | sfc ΔanomR | +.007 | +.004 | +.006 | +.123 |
| | rz ΔR | **+.097** | +.064 | +.080 | +.039 |
| | rz ΔanomR | +.008 | +.003 | +.017 | +.029 |
| 867 Jornada | sfc ΔR | +.008 | +.005 | +.007 | +.025 |
| | sfc ΔanomR | **+.013** | +.007 | −.005 | +.010 |
| 885 CO-Z1 | sfc ΔR | −.001 | −.000 | −.011 | +.014 |
| | sfc ΔanomR (1 station) | .000 | .000 | −.034 | +.006 |

Reading:
- **Walnut Gulch is the only place in the domain where an L1 arm beats SMAP-only on any in-situ metric (rz R +.097 vs +.039).** All
  runs improve rz R there by a lot, though, and the gain is almost entirely in raw R, not anomaly R (L1coh rz ΔanomR +.008). So DA mostly
  corrects the 2022 seasonal cycle at this tile rather than the day-to-day variability. At the surface, SMAP-only dominates as elsewhere.
- **At Jornada the L1 filter arm has the best surface anomaly R of all runs** (+.013, vs SMAP +.010 and L3 −.005). But the tile receives
  **no CYGNSS L1 obs at all**, so this is spill-over from neighbouring-tile updates (xcorr 0.15°), not direct L1 information.
- **Why Jornada and CO-Z1 get no L1 obs (resolved 2026-09-28): they are outside the preprocessing region.** The CYGNSS_operator
  pipeline (`scripts/sbatch_cygnss_m36_best_obs_daily_workflow.sh`: `REGION=arizona`, `REGION_BUFFER_KM=200`) keeps only specular
  points inside the Arizona polygon or within 200 km of it (`point_in_region`). Jornada is **215 km** from the Arizona border and CO-Z1
  **271 km**; Walnut Gulch is inside Arizona. In the 598k preprocessed obs for 2020–2022 there is not a single specular point within 0.25°
  of either tile. The preprocessed stream stops at about 106.8° W, while the GEOSldas domain extends to 106° W.
  Breakdown of the 335 of 909 domain tiles with no L1 obs:
  - 238 are outside Arizona + 200 km. This is 28% of the GEOSldas box, including Jornada, CO-Z1, and the NW/SW/SE corners.
  - 90 are in-region but at the CYGNSS latitude limit: every tile at 37.4–37.8° N, plus the rows north of that. The specular points are
    too sparse there for the scaling climatology's N_min.
  - 7 are isolated in-region tiles: the Gulf of California / Baja coast, the western edge, and the Grand Canyon (tile 391).
  To get L1 at Jornada the obs would have to be re-preprocessed with a larger buffer (≥ 250 km) or a `bounds` region matching the
  GEOSldas box, and the L1 scaling climatology rebuilt.
- With one year and effectively 1–3 independent tiles, this is a case study, not evidence of skill.

## Appendix: all stations in the GEOSldas box (misleading for L1, kept for reference)

195 surface / 167 rz stations; 78% of them are on tiles with no L1 obs. Mean Δ vs OL (station bootstrap):

| domain | metric | L1coh | L1coh_err39 | L3 | SMAP |
|---|---|---:|---:|---:|---:|
| surface | ΔR | +.002 | +.002 | −.000 | +.043 |
| surface | ΔanomR | +.002 | +.002 | −.004 | +.079 |
| rz | ΔR | +.009 | +.006 | +.003 | +.004 |
| rz | ΔanomR | +.005 | +.004 | +.001 | +.045 |

The L1 effect in this table is diluted about 4–5× by stations that L1 can reach only through spill-over from neighbouring tiles.
