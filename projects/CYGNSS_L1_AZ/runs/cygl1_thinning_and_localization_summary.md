# CYGNSS L1 paired-density experiment: thinning discrepancy + why localization can't substitute

## Context
Testing whether "DA-intermediate" (a thinned CYGNSS L1 obs stream, designed to sit between
fully-isolated single-obs and the full/dense obs stream) shows clean single-obs-like EnKF
update behavior, vs. "DA-dense" (full stream) which shows a sharp update-mechanism breakdown.
Two questions came up: (1) is intermediate's obs stream actually as spatially isolated as its
own design claims, and (2) could tightening the localization radius (xcompact/ycompact) achieve
the same isolation instead of thinning obs.

## Finding 1: thinning script's min-separation guarantee doesn't hold in the actual output data

`thin_cygl1_nested_density_6mo.py` builds "DA-intermediate" via a greedy admission algorithm
(`build_intermediate_candidate()`): starting from the "DA-sparse" force-kept set, it walks
remaining candidate observations in fixed order and admits one only if it is `>= min_sep_deg`
(locked at 2.40°) from every observation already kept in that same 3-hour assimilation window.
By construction this should make it mathematically impossible for two kept obs in the same
window to end up closer than 2.40° apart.

Checking the actual written obs files directly (not the design intent, not memory — the live
netCDF files under `CYGNSS_L1_thinned_intermediate_6mo/`), across the full Jan 2020–Dec 2021
record (20,587 kept obs, 731 daily files):

- Aggregate stats look right on the surface: median same-window nearest-neighbor separation
  ≈ 2.2–2.4°, p90 ≈ 2.7–2.75° (~220–270 km) — consistent with the 2.40° design target and the
  2.5° Gaspari-Cohn localization-overlap radius (`2×xcompact`).
- But the full distribution tells a different story: **45–57% of same-window nearest-neighbor
  pairs are closer than the nominal 2.40° minimum** (checked separately for the original
  Jan–Oct-2020 build and the later Nov-2020–Dec-2021 extension build — both show this, not just
  one run).
- The single closest pair found: two obs 0.029° apart (~3 km), essentially simultaneous
  (timestamps 0.5s apart, both round to the identical 3-hr window anchor) — i.e., genuinely two
  near-duplicate specular-point observations that both survived the admission filter.

This should not be possible given the algorithm as written — the admission check compares
every remaining candidate against *all* currently-kept members of its window, in a single
linear pass over a globally-shared `window_members` dict, regardless of which source-day file
or satellite an obs came from. Root cause is not yet identified. Candidates not yet ruled out:
a stale/mixed-version output file, an ordering/day-boundary edge case in how `remaining` is
processed, or an interaction between the two separate script invocations (original build vs.
extension build) that somehow isn't as isolated as intended.

**Practical implication:** the interpretation that "DA-intermediate stays clean because its obs
are properly isolated" is currently on shakier ground than previously stated — roughly half of
its same-window pairs are not actually isolated at the level the design assumed. This needs to
be resolved (or at minimum quantified: how much of the "clean" event-level result is coming
from the genuinely-isolated half vs. the accidentally-crowded half) before leaning further on
the density-experiment conclusions.

**Related, smaller housekeeping bug:** the thinning script writes its log CSV
(`cygl1_intermediate_6mo_thinning_log.csv`) to a fixed filename regardless of the `--beg-date`/
`--end-date` CLI overrides used to build separate date ranges. The later Nov-2020–Dec-2021
extension run overwrote the original Jan–Oct-2020 log. The actual obs netCDF data files were
not affected (they're written into per-year/month subdirectories), only the human-readable log
of which obs got kept and why.

## Finding 2: reducing xcompact/ycompact (localization radius) instead of thinning is a dead end

Already tested directly (2026-08-22) and re-verified against the current code/config
(2026-08-25, unchanged):

- GEOSldas hard-aborts if `xcompact/ycompact < 2 × (largest relevant correlation length)`
  (`check_compact()` in `clsm_ensupd_upd_routines.F90`, `min_compact_div_corr = 2.0`).
- In this project's active nml (`nml_paired_intermediate/LDASsa_SPECIAL_inputs_ensupd.nml`),
  the binding constraint is **not** CYGNSS L1's own obs-error correlation length
  (`xcorr=ycorr=0.25°`, would only require a 0.5° floor) — it's the **forcing-perturbation**
  correlation length (`xcorr_force_pert%pcp/sw/lw = 0.5°`), which sets the floor at
  **2 × 0.5° = 1.00°**. Current setting is 1.25°, leaving only ~0.25° of legal headroom before
  the run aborts.
- A sweep of xcompact=ycompact ∈ {1.25, 1.20, 1.10, 1.05, 1.00} against the full unthinned CYGNSS
  L1 obs stream showed the local-interaction-count metric (how many other obs fall inside a
  given obs's localization ellipse per cycle) drops only **median 60 → 42** (~30% reduction)
  across the *entire legal range*.
- DA-intermediate's actual thinning-achieved target was median=1 (a **~60x reduction** from the
  full-stream baseline) — localization tightening alone reaches nowhere near that, even at the
  code-enforced legal floor.

**Why it can't be pushed further:** the floor is set by the forcing-perturbation ensemble's
spatial correlation length, not by anything specific to CYGNSS L1. Lowering it further would
mean reducing `xcorr_force_pert`/`ycorr_force_pert` themselves — but that changes the spatial
structure of the meteorological-forcing perturbation ensemble used across *all* assimilated
species project-wide, not a scoped, CYGNSS-L1-specific change. It's a materially different and
much higher-blast-radius experiment than "just tighten the localization radius," with its own
implications for ensemble spread and error representation everywhere, not a drop-in substitute
for obs thinning.

**Bottom line:** localization-radius reduction cannot substitute for obs thinning as a way to
achieve DA-intermediate-level isolation — the legal range is far too shallow. Thinning remains
the only practical lever, which is why Finding 1's discrepancy (thinning not actually enforcing
its own isolation guarantee) is the higher-priority thing to resolve next, not localization.
