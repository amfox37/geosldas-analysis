# Legacy BUFR vs H121 QC flags: what we believed, what we found, what we changed

Last updated: 2026-06-26

## Purpose

This note documents the QC flag investigation that explains why GEOSldas
ingests far more H121 observations than legacy BUFR observations for the same
period/satellite, and what (if anything) should change in the readers as a
result. It supersedes an earlier mistaken interpretation of the legacy BUFR
flag tables (see "False lead" below) — keep this version as the reference.

Sources used:

- `H121 PUM v0.2` (HSAF/CDOP4/PUM, 2025-04-09) — H SAF's own product user
  manual. Its internal table cross-references for the BUFR flags turned out
  to be ambiguous/inconsistent (see below) — do not trust Table 4.7's "see
  Table 4.8 / 4.9" pointers in isolation.
- `Hahn et al. (2026)`, "Next-generation Metop ASCAT Surface Soil Moisture
  datasets", ESSD discussion paper (essd-2025-746) — describes the H121/H139
  algorithm and is the source for the `surface_soil_moisture_sensitivity`
  and `subsurface_scattering_probability` recommendations used below.
- `EUM/OPS-EPS/SPE/07/0343` v4A (27 Apr 2015), "ASCAT Level 2 Soil Moisture:
  Product Format Specification" — EUMETSAT's own native L2 SM format spec.
  This is the authoritative source for `CORRECTION_FLAGS` (Table 9) and
  `PROCESSING_FLAGS` (Table 10) bit definitions; it resolved the PUM
  ambiguity.
- GEOSldas Fortran source:
  `src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90`,
  subroutines `read_obs_sm_ASCAT_EUMET` (legacy) and `read_obs_sm_ASCAT_HSAF`
  (H121/H139).
- Local sample data: one day (2020-06-01), Metop-A, at
  `/Users/amfox/Desktop/ASCAT_SSM_CDR/discover_sample/`, decoded directly
  with `eccodes`/`netCDF4` and run through
  `projects/ascat_da/lib/qc.py` / `lib/readers.py`.

---

## 1. The headline number this explains

Ten-day global validation (2020-06-01 to 2020-06-10, all three Metop
satellites; see `h121_legacy_obs_validation_notes.md`):

| product | OFA total (10 days, 3 platforms) |
|---|---:|
| Legacy BUFR | 785,345 |
| H121        | 2,423,366 |

Ratio: **~3.09x**. For a single platform/day (Metop-A, 2020-06-01), the raw
post-QC obs ratio is closer to **~10x** (65,967 legacy vs 657,903 H121) before
GEOSldas superobbing onto 36 km tiles compresses it back down. This note is
about the QC half of that gap (resolution explains roughly 4x on its own —
see prior conversation; QC strictness explains the rest).

---

## 2. Field name -> bit table mapping (corrected)

The BUFR mnemonics read by `read_obs_sm_ASCAT_EUMET` are:

- `SMPF` = `soilMoistureProcessingFlag` (WMO descriptor `040006`)
- `SMCF` = `soilMoistureCorrectionFlag` (WMO descriptor `040005`)

The **PUM** (Table 4.7) cross-references these as "see Table 4.8" (SMPF) and
"see Table 4.9" (SMCF) — but Table 4.8 is captioned "correction flag" and
Table 4.9 is captioned "processing flag" in the same document, i.e. the
caption-to-field-name mapping looks swapped relative to the "see Table X"
pointers. This is almost certainly a LaTeX float-placement artifact (Table
4.8 prints on page 13, ahead of the citing text on page 14), not a real
ambiguity in the underlying product — but it is enough to cause a
misreading if you only have the PUM.

The **native EUMETSAT format spec** (EUM/OPS-EPS/SPE/07/0343, Tables 9 & 10)
removes the ambiguity, because its field names and table captions agree
directly:

### `PROCESSING_FLAGS` (16-bit native; maps to BUFR `SMPF`)

| bit | value | meaning |
|---|---:|---|
| 1 | 1 | Not meaningful soil measurement (insufficient valid Hamming-window neighbours, or invalid > valid) |
| 2 | 2 | Sensitivity to soil moisture ≤ 2 dB |
| 3 | 4 | Azimuthal noise ≥ 1 dB |
| 4 | 8 | Backscatter Fore-Aft beam out of range |
| 5 | 16 | Slope Mid-Fore beam out of range (> 6x noise of the slope) |
| 6 | 32 | Slope Mid-Aft beam out of range (> 6x noise of the slope) |
| 7 | 64 | Surface soil moisture below -20% |
| 8 | 128 | Surface soil moisture above 120% |
| 9-16 | — | Reserved |

All 16 bits set = flags not available.

### `CORRECTION_FLAGS` (8-bit native; maps to BUFR `SMCF`)

| bit | value | meaning |
|---|---:|---|
| 1 | 1 | Soil moisture between -20% and 0% |
| 2 | 2 | Soil moisture between 100% and 120% |
| 3 | 4 | Correction of wet backscatter reference applied |
| 4 | 8 | Correction of dry backscatter reference applied |
| 5 | 16 | Correction of volume scattering in sand applied |
| 6-8 | — | Reserved |

All 8 bits set = flags not available.

**Conclusion: `SMPF` carries genuine per-pass retrieval-quality information
(slope/backscatter geometry, sensitivity, azimuthal noise). `SMCF` carries
"a correction was applied" bookkeeping plus the SM-clip indicators.**

### False lead (do not repeat)

Earlier in this investigation we initially assumed the PUM's "see Table 4.8"
pointer for `SMPF` was correct as written, decoded the BUFR sample using
Table 4.8's bit content, and concluded that `SMPF`'s dominant rejection
causes were "dry backscatter reference correction" (bit 16) and "wet
backscatter reference correction" (bit 32) — i.e. that GEOSldas's
`smpf == 0` rule was throwing away routine corrections, analogous to the
correction_flag over-screening we'd already fixed for H121. **This was
wrong.** The native EUMETSAT spec shows `SMPF` is the 16-bit
`PROCESSING_FLAGS` table, not the 8-bit `CORRECTION_FLAGS` table, and the
empirical bit pattern in real data (next section) only makes sense under the
corrected mapping.

---

## 3. What the legacy BUFR sample actually shows

One day, Metop-A, 2020-06-01, decoded directly with `eccodes` (233,385 obs
with `0 <= ssm <= 100`):

### `SMPF` (`PROCESSING_FLAGS`) — fraction of valid-range obs with each bit set

| bit value | meaning | % of obs |
|---:|---|---:|
| 1 | insufficient neighbours | 0.0% |
| 2 | sensitivity ≤ 2 dB | 5.7% |
| 4 | azimuthal noise ≥ 1 dB | 0.4% |
| 8 | backscatter Fore-Aft OOR | 0.9% |
| **16** | **slope Mid-Fore OOR** | **40.5%** |
| **32** | **slope Mid-Aft OOR** | **38.5%** |
| 64 | SM < -20% | 0.9% |
| 128 | SM > 120% | 1.3% |

Current GEOSldas rule (`smpf == 0`, exact match) passes **40.1%** of these
obs. The dominant rejection cause is genuine retrieval noise (beam-slope
estimate out of range), not a "correction was applied" flag. **This QC rule
is doing legitimate, intentional quality control** — it is a structural
property of the legacy algorithm (per-pass local slope estimation is noisy;
see Hahn et al. 2026 Sect. 3.3.3 for the contrast with H121's smoothed
climatological approach) and is not something we are recommending changing.

### `SMCF` (`CORRECTION_FLAGS`) — fraction of valid-range obs with each bit set

| bit value | meaning | % of obs |
|---:|---|---:|
| 1 | SM between -20% and 0% | 6.1% |
| 2 | SM between 100% and 120% | 1.9% |
| 4 | wet backscatter ref correction applied | 51.4% |
| 8 | dry backscatter ref correction applied | 0.0% |
| 16 | volume-scattering-in-sand correction applied | 0.0% |

Current GEOSldas rule (`smcf in {0, 4}`) passes **92.0%** of these obs, and
explicitly **allows "wet correction applied" through** (bit 4 is in the
accepted set). Dry correction and volume-scattering correction never fire
in this sample, so there is no live asymmetry to fix here either.

**Net conclusion on legacy QC: nothing needs to change.** The two filters
that matter (`smpf == 0`, `smcf in {0,4}`) are both doing what they were
designed to do, and the dominant attrition (~60% of valid-range obs) comes
from a real, structural noise property of the per-pass beam-slope estimate
in the 25 km legacy algorithm — not from over-aggressive screening of
routine corrections.

The Fortran reader (`read_obs_sm_ASCAT_EUMET`) also applies two QC steps that
our Python mirror does **not** currently replicate:

- `ALFR >= 0.9` (land fraction must be ≥90%)
- an external static ASCAT obs mask (netCDF, regular lat/lon grid)

These make the real legacy yield somewhat lower than our Python-side
estimate; see "Legacy global match fraction remains around 0.70" in
`h121_legacy_obs_validation_notes.md`.

---

## 4. What H121's QC did at the start, and the two real gaps we found

The initial `read_obs_sm_ASCAT_HSAF` draft on the H121/H139 development branch
screened:

- `surface_flag` bit 0x01 (open water) → reject
- `processing_flag` bits 0x01 | 0x02 (model_parameter_not_usable,
  backscatter40_not_usable) → reject
- `wetland_fraction >= 10%` → reject
- `topographic_complexity >= 10%` → reject
- `subsurface_scattering_probability >= 10%` → reject
- `correction_flag` — **not screened at all** (deliberate; screening it
  previously hurt OFA match fraction and introduced bias — see
  `h121_legacy_obs_validation_notes.md`)

**Update, 2026-06-26:** on `GEOSldas_GridComp` branch
`feature/amfox/ascat-hsaf-v8`, the Fortran reader now also screens
`surface_soil_moisture_sensitivity <= 1 dB`, tightens
`subsurface_scattering_probability` to `>= 5%`, and rejects
`backscatter40_flag` bit 4 (`noise_out_of_limits`). That branch matches the
current Python `QC_DEFAULT_H121` on these H121-specific QC additions.

We compared this against two explicit, quantified recommendations in
Hahn et al. (2026):

1. **Sect. 4.1.1**: *"[sensitivity] serves as an indicator of retrieval
   uncertainty with values below 1 dB typically pointing to densely
   vegetated areas with a low backscatter signal variation."* — H121
   exposes `surface_soil_moisture_sensitivity` as a continuous variable
   (dB) with **no corresponding flag bit**; the initial GEOSldas draft did
   not screen it.
2. **Sect. 3.5** (the paper's own validation methodology): *"internal
   quality flags ... are applied to mask observations affected by
   subsurface scattering (subsurface scattering probability > 5%)."* — the
   initial GEOSldas draft used `>= 10%`, twice as permissive as what the
   dataset's authors validated against.

### Why `slope40`/`curvature40` do *not* need an analogous "out of range" flag

Legacy's `SMPF` bits 5/6 ("slope out of range") exist because the legacy
algorithm estimates slope **locally, per swath pass** — inherently noisy,
so each pass needs its own sanity check. H121 instead derives slope and
curvature from an **offline, kernel-smoothed multi-year climatology per
day-of-year** (Hahn et al. 2026, Sect. 3.3.3: Epanechnikov kernel, ±21-day
bandwidth). The "is this slope/curvature estimate trustworthy for this
location" question is answered once, upstream, when the climatology is
built — and the answer is already exposed via `processing_flag` bit 1
(`model_parameter_not_usable`), which GEOSldas **already screens**. Adding
a separate slope/curvature threshold would be redundant.

`backscatter40_flag` (sigma0_usable / slightly_degraded / noise_out_of_limits)
is the closest remaining analogue to legacy's "Backscatter Fore-Aft out of
range" bit, but that legacy bit only fired on 0.9% of obs in our sample —
minor. H121's `backscatter40_flag` bit 4 is now screened on
assimilation-design grounds because GEOSldas otherwise assigns the same fixed
observation error to noisy and clean H121 observations.

---

## 5. Empirical impact of the two new H121 thresholds

Same sample (Metop-A, 2020-06-01), 657,903 obs currently pass H121 QC.

### Sensitivity (`surface_soil_moisture_sensitivity`)

| reject if sensitivity ≤ | obs lost |
|---:|---:|
| 1.0 dB | 5.0% |
| 1.5 dB | 9.7% |
| 2.0 dB (legacy's literal threshold, before we knew it didn't apply here) | 26.4% |
| 2.5 dB | 42.4% |
| 3.0 dB | 54.5% |

Median sensitivity among currently-passing obs is 2.80 dB; 25th percentile
is 1.97 dB — so a 2 dB cutoff would remove roughly a quarter of all
currently-accepted H121 obs. We are using **1 dB**, per the paper's own
recommendation, not 2 dB (legacy's number doesn't transfer cleanly anyway,
since legacy's threshold was calibrated for a noisier per-pass estimate —
see Sect. 4 above).

### Subsurface scattering probability

| `subsfc_max` | obs lost from currently-passing |
|---:|---:|
| 10 (current default) | — (baseline) |
| 5 (paper's validation threshold) | 10.5% |

(70.5% of currently-passing obs have `subsfc == 0`; 21.8% are in `(0,5]`;
7.7% are in `(5,10]` — that last band is what the tighter threshold removes.)

---

## 6. What we changed

### Python mirror (`projects/ascat_da/lib/`)

- `qc.py`: added `sens_min` key and a corresponding `sens > qc['sens_min']`
  check in `apply_h121_qc`. Initially landed disabled (`sens_min=None`) so
  existing behaviour was unchanged unless explicitly overridden.
- `readers.py`: `read_h121` now also loads `surface_soil_moisture_sensitivity`
  and passes it into the QC arrays dict as `sens`.
- A 10-day global superob cache was rebuilt with `{'sens_min': 1.0,
  'subsfc_max': 5}` (version tag `geos_cycle_global_v3_sens1_subsfc5`) and
  validated — see Sect. 7 and 8 below for results.
- **Update:** after reviewing those results, `QC_DEFAULT_H121` itself was
  changed to make `sens_min=1.0` and `subsfc_max=5` (was `10`) the actual
  defaults, matching Hahn et al.'s recommendations. This is a behavioural
  change for every caller that doesn't explicitly override these keys
  (`read_h121`, `build_global_superob_cache.py`,
  `check_raw_superobs_vs_ofa.py`, `run_legacy_vs_h121_smoke.py`, and the
  notebook's `QC_H121`, which already set the same values explicitly and is
  now a redundant-but-harmless override). The old default (`subsfc_max=10`,
  no sensitivity screen) is what `geos_cycle_global_v2_geos_hsaf_qc` was
  built with; any future rebuild under those names would need an explicit
  override to reproduce the old behaviour.
- Added `bsflag_bad_bits` key and a corresponding bitwise check in
  `apply_h121_qc`; `readers.py` now also loads `backscatter40_flag` into
  the QC arrays dict as `bsflag`. After testing `bsflag_bad_bits=4`
  (reject `noise_out_of_limits`) vs `=6` (also reject `slightly_degraded`)
  against legacy and finding `=4` strictly more efficient (Sect. 9),
  **`QC_DEFAULT_H121` was updated again to set `bsflag_bad_bits=4` as the
  default**, decided on assimilation-design grounds rather than the
  (inconclusive) legacy cross-check — see Sect. 9 for the full reasoning.

### GEOSldas Fortran (`feature/amfox/ascat-hsaf-v8`)

Implemented changes to `read_obs_sm_ASCAT_HSAF` in `clsm_ensupd_read_obs.F90`:

- Add `real, parameter :: thr_sens = 1.0` and screen
  `surface_soil_moisture_sensitivity <= thr_sens` (NC_INT, scale factor
  `1e-7`, units dB; fill value is large-negative so it fails the threshold
  naturally, no separate fill check needed).
- Change `thr_subsfc` from `10.` to `5.`.
- Reject `backscatter40_flag` bit 4 (`noise_out_of_limits`).
- Update the subroutine's header QC-summary comment to match.

Recommendation: validate both the Python-side and Fortran-side changes
against OFA bias/match-fraction (the same way the original H121
`correction_flag` decision was validated) before promoting the feature branch.

---

## 7. Validation: does the new H121 QC improve agreement with legacy?

### A wrong first attempt

Our first validation attempt rebuilt the global superob cache with the new
H121 QC and compared it against the *existing* OFA files via
`check_global_superobs_vs_ofa.py`. That comparison's `bias_pct`/`rmse_pct`
columns measure `ofa_obs_pct - ssm_pct`, i.e. how closely a locally-rebuilt
superob reproduces the `obs` value GEOSldas's *production* OFA already
contains — a code-correctness check against the **existing** (old-QC)
Fortran reader, not an accuracy check of a **new** QC rule that doesn't
exist in production yet. Comparing a new QC regime against output generated
by the old regime can only look like drift, regardless of whether the new
rules are actually better. We discarded this comparison rather than report
it as evidence either way.

### The right test: Legacy vs H121, matched tile/cycle space

The actual goal (matching what `legacy_vs_h121_obs.ipynb` Sect. 8 already
does) is to compare H121 directly against legacy as two independent
retrievals, merged on `(date, platform, tilenum, cycle)`, using the global
superob cache directly — no OFA involved. Ten-day window, 2020-06-01 to
2020-06-10, all three Metop satellites:

| version | platform | matched tile/cycles | bias (H121−Legacy) | RMSE | MAE | median abs | r | H121 median n_obs |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| v2 (old QC) | Metop-A | 344,783 | -8.066 | 19.658 | 15.385 | 12.191 | 0.795 | 8 |
| v3 (sens≥1dB, paper's recommendation) | Metop-A | 322,092 | -8.733 | 19.902 | 15.621 | 12.493 | 0.788 | 7 |
| v4 (sens≥2dB, legacy's literal threshold) | Metop-A | 268,980 | -6.170 | 17.434 | 13.604 | 10.855 | 0.827 | 7 |
| v2 | Metop-B | 316,573 | -5.621 | 18.316 | 14.344 | 11.414 | 0.806 | 8 |
| v3 | Metop-B | 294,747 | -6.280 | 18.464 | 14.493 | 11.585 | 0.801 | 7 |
| v4 | Metop-B | 249,634 | -4.291 | 16.536 | 12.917 | 10.390 | 0.833 | 7 |
| v2 | Metop-C | 345,224 | -7.725 | 19.445 | 15.185 | 11.940 | 0.794 | 8 |
| v3 | Metop-C | 320,198 | -8.462 | 19.704 | 15.435 | 12.258 | 0.786 | 7 |
| v4 | Metop-C | 268,156 | -6.007 | 17.320 | 13.493 | 10.720 | 0.825 | 7 |
| v2 | All | 1,006,580 | -7.180 | 19.172 | 14.989 | 11.840 | 0.797 | 8 |
| v3 | All | 937,037 | -7.869 | 19.392 | 15.203 | 12.109 | 0.791 | 7 |
| v4 | All | 786,770 | -5.518 | 17.115 | 13.348 | 10.655 | 0.828 | 7 |

**Result at 1 dB (v3): the new QC does not improve agreement with legacy —
every metric is marginally worse** (bias more negative, RMSE/MAE up
slightly, r down slightly, ~7% fewer matched tile/cycles, H121 median obs
per tile down 8→7).

**Result at 2 dB (v4): every metric improves, and by a meaningful margin**
— bias shrinks (-7.18 → -5.52), RMSE drops (19.17 → 17.12), MAE drops
(14.99 → 13.35), r rises (0.797 → 0.828) — at the cost of ~22% fewer
matched tile/cycles than the v2 baseline.

**Why the 2 dB result isn't necessarily vindication of 2 dB as "more
correct," though:** per Sect. 8 below, a 2 dB cutoff removes ~92% of
rainforest obs (vs ~62% at 1 dB) — i.e. it strips out almost the entire
population where both legacy and H121 already struggle and disagree most.
Removing the worst-agreeing tail of any two noisy datasets will almost
always improve aggregate bias/RMSE/r on what remains, largely as a
mechanical consequence of the exclusion, not necessarily because the
surviving H121 retrievals became more accurate. It doesn't distinguish "the
QC made H121 better" from "we threw out the hard cases that disagree with
everything."

**Caveat on interpretation:** legacy isn't ground truth either — it's a
different retrieval with its own (likely larger) error characteristics, so
"agrees more/less with legacy" doesn't directly translate to "more/less
accurate." Taken together, this test gives **mixed, threshold-dependent
signal**: no support at 1 dB, an apparent (but likely partly mechanical)
improvement at 2 dB. The actual independent-truth validation behind the
1 dB and 5% recommendations is Hahn et al.'s own comparison against ISMN,
GLDAS, and ESA CCI (Sect. 4.2-4.3 of the paper), which remains the
better-grounded basis for choosing a threshold than either result here.

---

## 8. Regional impact: tropical rainforest coverage

The sensitivity screen specifically targets dense vegetation (Hahn et al.
2026, Sect. 4.1.1), so we checked whether it disproportionately removes
coverage over tropical rainforest. Same sample (Metop-A, 2020-06-01), full
new QC (`sens_min=1.0` + `subsfc_max=5`) vs. the old QC:

| region | obs before | sens≥1dB obs after | % lost (1 dB) |
|---|---:|---:|---:|
| Global | 657,903 | 555,882 | **15.5%** |
| Amazon (lat -10 to 5, lon -75 to -50) | 23,626 | 9,704 | **58.9%** |
| Congo basin (lat -5 to 5, lon 10 to 30) | 14,990 | 4,841 | **67.7%** |
| Indonesia/Maritime SE Asia (lat -10 to 5, lon 95 to 140) | 7,073 | 2,998 | **57.6%** |
| Rainforest (combined) | 45,689 | 17,543 | **61.6%** |
| Rest of globe | 612,214 | 538,339 | **12.1%** |

**Rainforest coverage drops ~4-5x more than the global average (~62% vs
~12%).** This matches Hahn et al.'s own framing of the sensitivity flag and
their independent finding (Sect. 4.2.1) that dense canopy in the
Amazon/Congo/Indonesia already suppresses retrieval skill before this
screen is even applied — so the effect is consistent with the flag doing
what it's designed to do (removing low-confidence retrievals), not an
artifact of a bad threshold.

### At 2 dB (legacy's literal threshold), rainforest coverage is essentially eliminated

We also checked the 2 dB cutoff (legacy's literal `SMPF` bit-2 threshold —
not recommended for H121, see Sect. 5, but checked here specifically to see
whether it converges on the same near-total rainforest dropout that legacy
shows):

| region | obs before | sens≥2dB obs after | % lost (2 dB) |
|---|---:|---:|---:|
| Global | 657,903 | 425,877 | **35.3%** |
| Amazon | 23,626 | 1,503 | **93.6%** |
| Congo basin | 14,990 | 1,911 | **87.3%** |
| Indonesia | 7,073 | 225 | **96.8%** |
| Rainforest (combined) | 45,689 | 3,639 | **92.0%** |
| Rest of globe | 612,214 | 422,238 | **31.0%** |

At 2 dB, rainforest coverage drops by ~92% combined (Indonesia loses 96.8%,
Amazon 93.6%), while the rest of the globe loses "only" 31%. This converges
toward the same near-total dropout legacy already exhibits in these
regions, consistent with dense canopy genuinely defeating C-band retrieval
regardless of which algorithm processes the backscatter — at a strict
enough sensitivity floor, both products converge on "essentially no usable
signal in rainforest," they just need different threshold values to get
there because the underlying sensitivity distributions aren't on the same
footing (per Sect. 5's caveat about legacy's noisier per-pass estimate vs.
H121's smoothed climatology).

**Operational implication, independent of retrieval-accuracy arguments:**
adopting this threshold in GEOSldas would mean the analysis over the
Amazon/Congo/Indonesia leans much more heavily on the model background
during assimilation than it does today, simply from reduced observation
density. Worth weighing before adopting the threshold, separately from
whether the screened-out obs are actually less accurate.

---

## 9. `backscatter40_flag`: a third QC candidate, decided on different grounds

`backscatter40_flag` (PUM Table 4.15) carries three independent bits:
bit 1 = `sigma0_usable`, bit 2 = `sigma0_slightly_degraded`, bit 4 =
`sigma0_noise_out_of_limits`. None were screened. We tested two options
(both layered on top of the v3 default QC) over the full 10-day global
window:

| version | total H121 obs | rainforest obs | bias (All) | RMSE (All) | r (All) |
|---|---:|---:|---:|---:|---:|
| v3 (baseline) | 16,663,093 | 497,287 (3.0%) | -7.869 | 19.392 | 0.791 |
| v5: reject bit 4 (`noise_out_of_limits`) only | 15,805,936 (-5.1%) | 482,922 (-2.9%) | -7.971 | 19.525 | 0.793 |
| v6: reject bits 2 or 4 (`slightly_degraded`\|`noise_out_of_limits`) | 15,263,511 (-8.4%) | 473,959 (-4.7%) | -8.118 | 19.568 | 0.793 |

Unlike sensitivity, **this flag is geographically neutral-to-favorable**:
both options remove proportionally *less* from rainforest than globally
(v6 removes 8.4% globally but only 4.7% from rainforest) — the opposite of
the sensitivity screen's pattern. But the legacy cross-check alone is
**inconclusive**: bias/RMSE on the matched "All" population get marginally
*worse*, not better, for both options — split as a modest improvement
within rainforest (bias -2.57 → -2.00/-2.01; RMSE 15.84 → 15.56/15.54) and
a modest degradation everywhere else, which dominates the aggregate since
rainforest is only ~1.5% of the matched sample. `v6` buys nothing over
`v5` — the extra `slightly_degraded` bit discards ~3 more percentage
points of data for statistically indistinguishable bias/RMSE/r.

**Decided on a different basis than the legacy cross-check.** Checking
`read_obs_sm_ASCAT_HSAF` directly (line ~2618):
`ASCAT_sm_std(ii) = this_obs_param%errstd / 100.` — every accepted H121
observation is assigned the **same fixed, namelist-configured observation
error std**, regardless of flag status. There is no per-observation noise
down-weighting anywhere in this reader (this is the same static-`errstd`
pattern used for every other obs type in the file — legacy ASCAT, SMOS,
SMAP, CYGNSS, etc.). That means a `noise_out_of_limits` observation
currently gets *exactly* the same weight in the Kalman update as a clean
one. With no fallback down-weighting mechanism, hard rejection is the only
available lever against elevated backscatter noise — so we are adopting
**`bsflag_bad_bits=4`** (reject `noise_out_of_limits` only, i.e. `v5`,
not `v6`) as a conservative default, despite the inconclusive local
cross-check. The cross-check's limitations (legacy not ground truth, small
effect size, one-month/10-day sample, rainforest signal diluted by the
much larger rest-of-globe population) mean it shouldn't be read as
evidence *against* the change either.

---

## 10. Bottom line

- Legacy QC (`smpf == 0`, `smcf in {0,4}`): **no change recommended** —
  confirmed to be legitimate, intentional quality control, not
  over-aggressive correction-screening.
- H121 sensitivity screen + tightened subsurface threshold
  (`subsfc_max=5`), mechanically implemented and tested at two sensitivity
  cutoffs (Sect. 5, 7, 8):
  - **1 dB (Hahn et al.'s recommendation, now the Python default):** no
    improvement vs legacy (every metric marginally worse), and loses ~62%
    of rainforest obs vs ~12% globally.
  - **2 dB (legacy's literal threshold, not paper-recommended for H121):**
    clear improvement vs legacy on every metric, but likely largely
    mechanical — it loses ~92% of rainforest obs, i.e. it mostly works by
    removing the population where both products already disagree most.
  - Neither local result should override Hahn et al.'s own
    ISMN/GLDAS/CCI-based validation, which is the actual independent-truth
    basis for the 1 dB recommendation. The legacy cross-check here is, at
    best, a sanity check, not a substitute.
  - Whichever threshold (if any) is adopted, the rainforest coverage
    trade-off in Sect. 8 should be weighed explicitly — it is a much larger
    and more concentrated effect than the global average suggests.
- H121 `backscatter40_flag` screen (`bsflag_bad_bits=4`, now the Python
  default): adopted on assimilation-design grounds (no per-obs noise
  down-weighting exists in `read_obs_sm_ASCAT_HSAF`), not because the
  legacy cross-check showed a clear benefit — it didn't, and shouldn't be
  read as showing harm either. See Sect. 9.
- **All three changes are now baked into `QC_DEFAULT_H121`** in
  `projects/ascat_da/lib/qc.py` (`sens_min=1.0`, `subsfc_max=5`,
  `bsflag_bad_bits=4`) and implemented in the Fortran HSAF reader on
  `GEOSldas_GridComp` branch `feature/amfox/ascat-hsaf-v8`. The companion
  implementation/rationale note is `geosldas_fortran_handoff_h121_qc.md`.

---

## 11. Summary table: legacy vs H121 QC, side by side

| concept | legacy BUFR (`SMPF`/`SMCF`) | H121 NetCDF | currently screened? | action |
|---|---|---|---|---|
| open water | external static mask + land fraction `ALFR>=0.9` | `surface_flag` bit 1 | yes (both) | none |
| model/backscatter not usable | `SMPF` bits 1,4 (insufficient neighbours; backscatter Fore-Aft OOR) | `processing_flag` bits 1,2 | yes (both) | none |
| local slope/curvature unreliable | `SMPF` bits 5,6 (slope Mid-Fore/Mid-Aft OOR) — per-pass check | n/a — handled upstream by climatology; falls under `processing_flag` bit 1 if model param unusable | yes, via different mechanism | none — not a gap |
| sensitivity too low | `SMPF` bit 2 (≤2 dB) | `surface_soil_moisture_sensitivity` (continuous, dB) | legacy: yes; H121: **no** | **added** (`sens_min=1.0` in Python; Fortran instructions drafted) |
| azimuthal noise | `SMPF` bit 3 (≥1 dB) | handled structurally by the 12-direction empirical azimuth correction (Hahn et al. 2026 Sect. 3.3.1) | n/a | none — algorithm change removed the need |
| SM out of physical range | `SMPF` bits 7,8 / `SMCF` bits 1,2 | `processing_flag` bits 3,4 / `correction_flag` bits 2,3 | legacy: yes; H121: implicitly via `ssm_min`/`ssm_max` | none |
| correction applied (wet/dry/vol-scatter) | `SMCF` bits 3,4,5 | `correction_flag` bits 1 (wet only — no dry/vol-scatter equivalent) | legacy: wet allowed, dry/vol-scatter would be rejected (neither ever fires); H121: not screened at all (deliberate) | none — confirmed not an asymmetry |
| wetland fraction | n/a (BUFR has no separate field beyond IWFR) | `wetland_fraction >= 10%` | yes | none |
| inundation/wetland (legacy) | `IWFR <= 10%` | (see above) | yes | none |
| topographic complexity | `TPCX <= 10%` | `topographic_complexity >= 10%` | yes (both, same 10% threshold) | none |
| subsurface scattering | n/a (not in legacy product) | `subsurface_scattering_probability` | H121: yes, but at 10% not the paper's 5% | **tightened to 5%** |
| backscatter noise | `SMPF` bit 8 (backscatter Fore-Aft OOR, 0.9% of obs — minor) | `backscatter40_flag` bit 4 (`noise_out_of_limits`) | legacy: yes (minor); H121: **no** | **added** (`bsflag_bad_bits=4` in Python; decided on assimilation-design grounds — see Sect. 9 — not the legacy cross-check; Fortran instructions drafted) |
| frozen soil / snow | not in legacy BUFR fields used here | static probability flags exist but GEOSldas uses ERA5 model state instead (`qc_model_based_for_sat_sfmc`) | yes, via model state (both) | none |
| RFI | not flagged in either product | not flagged in either product | no | can't fix — no flag exists |
