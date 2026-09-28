# CYGNSS L1 AZ: ISMN in-situ validation of the fixed-operator arms, 2020–2022

*2026-09-28.* This is the post-fix follow-up to `cygl1_dense075_coh05_ismn_validation.md`. That run was pre-fix, so its DA result is
invalid. Same driver and settings (`projects/ascat_da/scripts/run_ismn_ol_da_skill.py`, daily `tavg24_1d_lnd_Nt` SFMC/RZMC, `--nmin 90`,
`--max-distance-deg2 0.1`). Job `jobs/run_ismn_ol_da_skill_az.sbatch`, whose runs are now passed in the `RUNS` variable (job 58621216,
22 min, all 1096 days found for every run). Output: `output/ismn_fixedop_20200101_20221231/` (gitignored).

Runs: OL = `OLv8_M36_AZ_fixedop` (reference); L1coh = `DA_L1_full_xc015_coh040216_fixedop`; L1coh_err39 = `…_err39_fixedop`;
L3 = `DA_L3_fixedop`; SMAP = `DA_SMAP_fixedop`. Sites: 195 surface and 167 root zone, from SNOTEL 116/114, SCAN 36/34, SOILSCAPE 20/4,
USCRN 13/7 and iRON 6/6. SOILSCAPE enters because its 2022 deployment is inside the window.

## Mean skill

| domain | run | R | anomR | ubRMSE |
|---|---|---:|---:|---:|
| surface | OL | 0.5383 | 0.4599 | 0.0636 |
| | L1coh | 0.5401 | 0.4617 | 0.0636 |
| | L1coh_err39 | 0.5403 | 0.4615 | 0.0636 |
| | L3 | 0.5380 | 0.4556 | 0.0636 |
| | SMAP | **0.5809** | **0.5386** | **0.0616** |
| rz | OL | 0.5841 | 0.5289 | 0.0440 |
| | L1coh | 0.5933 | 0.5343 | 0.0439 |
| | L1coh_err39 | 0.5905 | 0.5330 | 0.0439 |
| | L3 | 0.5873 | 0.5299 | 0.0440 |
| | SMAP | 0.5880 | **0.5737** | 0.0437 |

## Paired per-station change vs OL (mean Δ, 95% bootstrap CI over stations, % of stations improved)

| domain | metric | L1coh | L1coh_err39 | L3 | SMAP |
|---|---|---|---|---|---|
| surface | ΔR | +.0018 [−.0007, +.0043] 58% | **+.0020 [+.0004, +.0036] 60%** | −.0003 [−.0022, +.0017] 56% | **+.043 [+.036, +.049] 85%** |
| surface | ΔanomR | +.0017 [−.0009, +.0044] 56% | +.0016 [−.0003, +.0037] 59% | **−.0043 [−.0066, −.0022] 46%** | **+.079 [+.071, +.086] 92%** |
| rz | ΔR | **+.0091 [+.0051, +.0139] 70%** | **+.0063 [+.0035, +.0094] 72%** | +.0031 [−.0003, +.0068] 63% | +.004 [−.008, +.018] 45% |
| rz | ΔanomR | **+.0054 [+.0010, +.0114] 60%** | **+.0041 [+.0013, +.0078] 62%** | +.0010 [−.0028, +.0048] 56% | **+.045 [+.030, +.060] 72%** |

Excluding SNOTEL (mountain snow sites, 60% of the sample), n = 78 surface / 53 rz:
- The L1 root-zone gain roughly doubles: ΔR L1coh +.020 [+.009, +.034], L1coh_err39 +.013 [+.006, +.021].
- Surface ΔR for L1coh_err39 is +.004 [+.001, +.008].
- L3 surface ΔanomR stays negative: −.005 [−.009, −.001].

## Excluding SNOTEL

78 surface sites (SCAN 36, SOILSCAPE 23, USCRN 13, iRON 6) and 53 root-zone sites (SCAN 34, USCRN 7, SOILSCAPE 6, iRON 6).

| domain | run | R | anomR | ubRMSE |
|---|---|---:|---:|---:|
| surface | OL | 0.6180 | 0.4786 | 0.0367 |
| | L1coh | 0.6228 | 0.4802 | 0.0366 |
| | L1coh_err39 | 0.6222 | 0.4796 | 0.0366 |
| | L3 | 0.6205 | 0.4736 | 0.0365 |
| | SMAP | **0.6610** | **0.5459** | **0.0353** |
| rz | OL | 0.5544 | 0.5248 | 0.0253 |
| | L1coh | 0.5746 | 0.5350 | 0.0251 |
| | L1coh_err39 | 0.5672 | 0.5313 | 0.0252 |
| | L3 | 0.5615 | 0.5269 | 0.0252 |
| | SMAP | 0.5781 | 0.5481 | 0.0253 |

Paired change vs OL (mean Δ [95% bootstrap CI], % of stations improved; bold = CI excludes 0):

| domain | metric | L1coh | L1coh_err39 | L3 | SMAP |
|---|---|---|---|---|---|
| surface | ΔR | +.0049 [−.0002, +.0104] 59% | **+.0043 [+.0010, +.0080] 64%** | +.0026 [−.0016, +.0068] 54% | **+.043 [+.032, +.056] 85%** |
| surface | ΔanomR | +.0016 [−.0039, +.0074] 49% | +.0010 [−.0030, +.0048] 53% | **−.0050 [−.0091, −.0014] 37%** | **+.067 [+.056, +.079] 89%** |
| surface | ΔubRMSE ×10⁻³ | −.11 [−.25, +.02] 55% | **−.12 [−.22, −.03] 65%** | **−.24 [−.40, −.08] 68%** | **−1.35 [−1.89, −.84] 71%** |
| rz | ΔR | **+.020 [+.009, +.034] 70%** | **+.013 [+.006, +.021] 70%** | +.007 [−.002, +.017] 64% | +.024 [−.003, +.053] 57% |
| rz | ΔanomR | +.010 [−.001, +.028] 57% | **+.007 [+.000, +.016] 57%** | +.002 [−.006, +.009] 59% | +.023 [−.009, +.060] 65% |
| rz | ΔubRMSE ×10⁻³ | **−.24 [−.37, −.12] 75%** | **−.16 [−.25, −.08] 72%** | **−.14 [−.25, −.03] 70%** | −.00 [−.45, +.45] 47% |

Mean ΔR by network (surface / rz):

| network | L1coh | L1coh_err39 | L3 | SMAP |
|---|---|---|---|---|
| SCAN | +.005 / +.014 | +.004 / +.008 | −.004 / −.005 | +.038 / +.032 |
| SOILSCAPE | +.012 / +.097 | +.009 / +.064 | +.017 / +.080 | +.067 / +.039 |
| USCRN | −.006 / +.003 | −.003 / +.003 | −.004 / +.007 | +.038 / +.009 |
| iRON | +.001 / +.001 | +.001 / +.001 | +.001 / +.001 | −.006 / −.022 |

**Caveat: SOILSCAPE is 2022 only** (115–242 paired days per station) and drives much of the rz mean: its 6 rz sites contribute
+.06 to +.10. Excluding both SNOTEL and SOILSCAPE (55 sfc / 47 rz):
- rz ΔR stays significant for both L1 arms: L1coh +.010 [+.002, +.023], L1coh_err39 +.006 [+.001, +.014].
- L1 surface is not significant: ΔR +.002, ΔanomR ≈ 0.
- L3 becomes negative: sfc ΔR −.003 [−.007, −.000], sfc ΔanomR −.006 [−.012, −.002], rz ΔR −.002.
- SMAP is unchanged: sfc ΔanomR +.070.

## SOILSCAPE only

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
- **Open question for the operator team:** Jornada and CO-Z1 are among the 335 of 909 domain tiles with zero L1 obs over 2020–2022 in
  every arm, including the unfiltered monitor stream. Jornada is flat, sparsely vegetated desert, which should be a good CYGNSS target,
  so the reason it gets no L1 obs is worth checking: preprocessing/owner-grid assignment, a QC or topography mask, or the scaling-clim
  N_min. SMAP Tb (1,623 obs) and CYGNSS L3 (1,098 obs) are both available on that tile.
- With one year and effectively 1–3 independent tiles, this is a case study, not evidence of skill.

## Reading

- **In situ, SMAP-only is the only arm with a large gain.** Surface anomR is +0.08, and it improves at 92% of stations. This matches its
  −10 to −14% Tb O-F. The ASCAT (+12.7%) and L3 (+3.9%) O-F degradations seen in the O-F comparison do not show up as in-situ harm.
- **The L1 arms give a small but statistically robust in-situ gain, mostly in the root zone.** rz R is +0.006 to +0.009 (CI excludes 0,
  improves at 70% of stations), and larger away from SNOTEL. The medians are ≈ 0, so the gain comes from a minority of stations with
  real changes, while most stations barely move. That is consistent with the ~1% Tb O-F gains.
- **L3-only, the best non-SMAP arm in O-F space (L3 O-F −5 to −6%, ASCAT −2.6 to −3.6%), is neutral to slightly negative in situ.** Its
  surface anomR is −0.004 (CI excludes 0). The L3-vs-L1 ranking from O-F statistics does not hold up against in-situ data; if anything,
  L1 ≥ L3 in situ.
- The effect sizes for L1/L3 (±0.002–0.009 in R) are small compared with the station-to-station spread. Treat them as "not harmful,
  slightly positive for L1", not as a strong skill result.

Possible follow-ups: seasonal split (does the spring O-F degradation show in situ?), and per-station maps of ΔR against CYGNSS
L1 obs density.
