# H SAF ASCAT Reader, QC & Peatland Footprint QC — Development Log

A chronological engineering record of building H SAF H121/H139 ASCAT support into
GEOSldas, tuning its quality control, running the production DA experiments that
validated it, and merging it to `develop`.

- **Repo:** `GEOSldas_GridComp`
- **Branches:** `feature/amfox/ascat-hsaf-v8` → `feature/amfox/ascat-peatland-qc`
- **Merged as:** PR #186
- **Project dir:** `hsaf_cdr_test`
- **Span:** 2026-06-12 → 2026-08-12

**At a glance:** 21 commits across 2 branches · 1 merged PR (#186) · 2 six-year DA
production runs · 6 major HPC incidents diagnosed.

---

## Phase 0 — Reader inception (2026-06-12)

Added `read_obs_sm_ASCAT_HSAF()` to `clsm_ensupd_read_obs.F90` on a new branch,
`feature/amfox/ascat-hsaf-v8`: three new observation species
(`ASCAT_HSAF_META_SM`, `METB_SM`, `METC_SM`) reading H SAF H121 CDR / H139 ICDR
NetCDF swath files on the native 12.5 km Fibonacci grid, in place of the legacy
25 km EUMETSAT BUFR product. `N_obs_species_nml` bumped 55→58. QC at this point
used only the flags embedded in the netCDF files themselves — no external mask
file needed, unlike the legacy reader.

**Landed:** commit `4f97d0d` + two same-day follow-ups (`2e7034e` obs_dir_hier
fix, `4e52094` drop redundant snow/frozen QC, deferring to the existing
model-based check).

## Phase 1 — QC criteria tuning (2026-06-15 → 2026-06-20)

Iterated the QC chain over five commits, each targeting a specific known ASCAT
retrieval failure mode:

- **Fill-value handling** (`8e21406`) — wetland/topography/subsurface-scattering
  probability bytes use `-128` as fill; fixed misinterpretation and updated
  MetOp-B/C data paths.
- **Subsurface fill semantics** (`1889ba3`) — special-cased a fill value on
  `subsurface_scattering_probability` to mean "no scattering detected" (accept),
  not missing data. Wetland/topo fill stayed a reject (static-database coverage
  unknown there).
- **Subsurface threshold + SSM sensitivity** (`dda64f0`) — tightened the
  subsurface-scattering reject threshold 10%→5%, and added a reject on
  `surface_soil_moisture_sensitivity ≤ 1 dB`, following Hahn et al. (2026, ESSD).
- **Backscatter noise** (`a9d0ba8`) — hard reject on `backscatter40_flag` bit 4
  (noise-out-of-limits) only, deliberately excluding the "slightly degraded" bit
  since the reader assigns a fixed `errstd` with no per-obs down-weighting.

> **Finding.** A 10-day innovation comparison (baseline vs. tuned QC,
> monitor-only so O−F = O−A exactly) showed the new QC removed ~15.5% of obs but
> slightly *increased* aggregate bias/RMSE — counter to expectation. Attribution
> by replicating each check against 72 raw H121 granules found the
> subsurface-threshold tightening (10%→5%) alone responsible for ~74% of new
> rejections, and it was preferentially removing *well-fitting* observations,
> concentrated 10–50°N. The literature-motivated additions (SSM sensitivity,
> backscatter noise) were each ~10% of rejections and well-targeted. Flagged for
> reconsideration ahead of full reprocessing; the fix ultimately shipped anyway
> pending a longer validation window.

## Phase 2 — Open-loop validation & a real legacy bug (2026-06-21 → 2026-06-24)

Rebuilt `OLv7_M36_MULTI_type_13_H121` (a 2024 sensitivity-sweep cell) to add
H121 obs alongside SMAP and legacy ASCAT, then validated at 1-day and 10-day
scale against the original 2024 run. Forecast values matched at floating-point
precision (correlation 1.000000) — the rebuild faithfully reproduces model
physics. Every obs-count discrepancy was individually traced to a specific,
non-code cause:

- Three config-drift bugs from schema changes since 2024 (`SNOW_ALBEDO_INFO`
  silently defaulting on, `out_ObsFcstAna` switching from logical to integer
  type, a dead SMAP input path) — all fixed.
- An ASCAT observation mask revision and an SMAP L1C_TB reprocessing-release
  mismatch (`SPL4SM_Vv7032` vs. `OL8000`) — both external-data provenance
  issues, resolved by pointing at the exact release the 2024 run used.

> **Bug found (pre-existing, in old code).** Confirmed a real historical
> defect: the 2024-vintage codebase double-counted legacy EUMETSAT ASCAT
> observations across *consecutive* 3-hour assimilation windows (a broken
> `.and.` exclusion check that never actually excluded anything, later fixed
> upstream as `.or.`). Quantified at ~35–37% obs inflation in the old data.
> Verified **zero** such duplicates in the current codebase. Any historical
> analysis of that archived run's ASCAT obs counts should be treated as
> inflated.

## Phase 3 — DA production launch: H121 vs. legacy (2026-07-02)

Built fresh z-score scaling climatologies from the validated OL run's O−F
statistics, then launched two parallel 6-year (2015-04-01 → 2021-04-01)
24-member ensemble DA production experiments sharing identical
forcing/BCS/restart lineage and scaling — one assimilating H121 HSAF ASCAT, one
assimilating legacy EUMETSAT ASCAT — to make the observation-product swap the
only variable.

> **Gotcha (documented, recurred later).** Two setup traps found and fixed
> before launch, then written down for reuse: Fortran's
> `obs_param_nml%scalename` field is `character(80)` and silently truncates
> longer scaling-file names (fix: short symlinks); and a custom `HISTRC_FILE`
> copied from another experiment keeps that experiment's old `EXPID:` line,
> which `ldas_setup` does not rewrite.

## Phase 4 — The checkpoint-write stall (2026-07-04 → 2026-07-06)

Both production jobs began timing out mid-run. An initial `NUM_SGMT`/walltime
sizing mistake was fixed quickly, but the jobs kept dying anyway — not during
physics (fast, ~2h/simulated month) but during the ensemble
restart-checkpoint write phase, which would silently stall for 8–10 hours
until SLURM killed the job. This recurred three times running.

- A live `gdb` attach to the hung MPI ranks (a first for this project, via
  `srun --overlap`) caught it mid-deadlock: rank 0 sat in `MPI_Barrier`
  (already in ESMF teardown) while worker ranks waited in `MPI_Gatherv` for a
  partner that had already moved on — a genuine MPI collective call-count
  mismatch.
- A code hypothesis (unsynchronized `npert` used as a collective loop bound in
  `LandPertGridComp::Finalize`) was proposed, then disproven on closer reading
  — `npert` is actually broadcast-synchronized everywhere it matters.

> **Root cause, confirmed.** Stale `*_internal_checkpoint` files left in
> `scratch/` by earlier walltime-killed jobs, combined with
> `overwrite_checkpoint:F` (NetCDF NOCLOBBER), caused a silent `NC_EEXIST`
> (status = −35) failure on the writer rank for whichever ensemble member's
> file already existed. That rank fell through into ESMF finalize while the
> other 23 ranks waited on it inside a collective — permanent deadlock. Not a
> science-code bug; a restart-hygiene gap in the resubmit procedure. Fixed by
> clearing scratch checkpoint debris before every resubmit; the chain then ran
> cleanly for 7+ consecutive segments.

## Phase 5 — Walltime/timeout chase to completion (2026-07-10 → 2026-07-22)

With the stall resolved, both jobs still hit a slower-burning problem: a
genuine, gradual ~25–30% per-simulated-month slowdown from 2018 onward,
timing-consistent with MetOp-C (launched Nov 2018) ramping into the
operational H SAF record and swelling obs volume. Walltime was raised in steps
(10:30→11:55) for each experiment as it individually started missing its
segment, until both sat within 5 minutes of the compute partition's hard
12:00:00 ceiling. A project-wide inode quota exhaustion (798,844/800,000
files) briefly blocked resubmission and required an unrelated cleanup pass
before the runs could resume. Both experiments reached 92–96% of their 6-year
span by 2026-07-22, completing shortly after.

## Phase 6 — Comparison pipeline, then a FOV/history bug (2026-07-14 → 2026-07-29)

Fixed a pre-existing field-misalignment bug in the Python ObsFcstAna
postprocessor's `read_obs_param()` (two fields inserted upstream between
`units` and `path` in the current obsparam text format had shifted everything
after them) and used the corrected pipeline to build the O−F comparison shown
in the science report.

Separately, a new `DAv7_M36_ASCAT_type_13_H121_FOV12p5` experiment (tightening
the H SAF footprint radius to 12.5 km) crashed at startup with `Cannot Find
ENSAVG`: an upstream commit had split one History averaging component into
four (`METFORCEAVG`/`LANDAVG`/`LANDICEAVG`/`ROUTEAVG`), and this experiment's
`HISTORY.rc` template hadn't been updated. Fixed by re-mapping ~150
field/source pairs to the correct component, verified against each GridComp's
actual export list.

## Phase 7 — Peatland footprint QC (2026-07-28 → 2026-08-04)

Opened a new branch, `feature/amfox/ascat-peatland-qc`, adding a footprint-based
peat/organic-soil screen on top of the existing QC chain: for each ASCAT
retrieval, a Gaussian-weighted peat fraction is computed over all model tiles
inside the sensor footprint, and the observation is rejected if that fraction
exceeds 10%. The peat indicator itself was iterated three times before landing
— first from a soil-class flag, then a composite soil class, finally from
**PEATCLSM porosity** directly (commit `95f0317`), which is what the validated
results use.

Getting a clean 6-month test run out of this branch took four attempts:

> **Incident 1 — EXPID gotcha, recurred.** Despite already being documented
> from the July 2 launch, the same `HISTRC_FILE`/`EXPID` mismatch was hit
> again: the run "completed" (SLURM exit 0, restarts correct) but 5 of 6
> months of diagnostic history were written under the wrong EXPID and
> silently never bundled, then wiped by the next job's scratch cleanup. ~16
> hours of node time lost. This forced a change in process, not just a fix:
> `EXPID:` is now grepped and confirmed immediately after every `ldas_setup`
> that uses a custom `HISTRC_FILE`, and a GEOSldas run is never called "done"
> from SLURM exit code alone — `output/.../cat/` is checked for real
> populated files first.

> **Incident 2 — shared-checkout binary collision.** The relaunch crashed
> with `obs_param_nml%fcstvarname must be NULL on input` after ~4h47m. Root
> cause: an unrelated CYGNSS branch was checked out and rebuilt into the
> *same shared* `install/bin/GEOSldas.x` this job was running against — live,
> mid-run — changing `N_obs_species_nml` 58→59 under a namelist frozen at the
> old 58-species template. Fixed by reverting the branch and rebuilding;
> formalized as a standing rule: never switch branches or rebuild in a shared
> checkout while a job depending on it is running.

The 6-month test (2015-04 → 2015-10) finally completed cleanly against commit
`95f0317` (PEATCLSM porosity), the version used for every result reported.

## Phase 8 — Comparison results (2026-08-05)

Compared the peatland-QC run against its pre-QC baseline and, separately,
investigated a specific recollection of DA underperforming around Hudson Bay.
Full numbers and interpretation are in the companion science report; in
short: the QC removes ~6.4% of ASCAT obs, concentrated almost entirely in
high-latitude peatlands, and tightens ASCAT self-fit there by ~4.4% with no
measurable side effect on the independent SMAP check — and confirmed a real,
localized exception at Hudson Bay where DA degrades the independent check,
which the peat QC shrinks but does not eliminate.

## Phase 9 — Footprint mechanism validation (2026-08-10)

To explain a narrow band of very-low-but-nonzero obs counts right at
peatland boundaries, reconstructed the exact Fortran QC algorithm (ellipse
search, Gaussian footprint weighting, porosity threshold) in Python from
source and tested it against 2,720 candidate observations across 6 named
boundary tiles, using the actual BCS porosity field and tile geometry.

> **Validated.** 100% agreement — 2,690/2,690 correctly-predicted rejections
> and 30/30 correctly-predicted acceptances, zero disagreements. Closes the
> causal loop: footprint peat fraction crossing 10% is demonstrably why these
> observations are rejected, driven by fixed high-porosity peat tiles sitting
> immediately adjacent to mineral owner tiles. Not yet extended past these 6
> tiles to the full Canada/Alaska domain.

## Phase 10 — Extension & merge into develop (2026-08-11 → 2026-08-12)

Two final commits (`bfd293b` "Apply peatland QC to sfmc observations",
`5a807cd` "Document peatland QC in changelog") extended the QC to a second
observation type beyond what the 6-month validation run had tested, and
closed out the changelog. **Neither commit carries a Claude co-authorship
trailer** — they were made directly, outside any logged session.

**PR #186**, "Add H SAF ASCAT H121/H139 soil moisture reader, incl. QC of
retrievals," opened from `feature/amfox/ascat-hsaf-v8` (by then carrying the
full peatland-QC chain) against `develop`, was reviewed and retested by a
co-developer, then merged 2026-08-12. Reviewer note on the PR: *"0-diff for
existing tests (incl. SMAP Tb assimilation), but not 0-diff for assimilation
of any observations of 'sfmc' or 'sfds' (owing to addition of peatland QC for
those observation types)."*

---

## Current state

*As verified directly against the live git history and GitHub, 2026-09-15.*

| Status | Note |
|---|---|
| ✅ Merged | The full H SAF ASCAT reader + all QC (including peatland QC and the sfmc extension) is in `develop` via PR #186, merge commit `e7b1647`. |
| ⚠️ Stale | The shared `GEOSldas_develop` checkout on Discover is still on the pre-merge `feature/amfox/ascat-hsaf-v8` tip — needs `git checkout develop && git pull` and a rebuild before any new experiment picks up the merged version. |
| 🔲 Open | Footprint-QC mechanism validation covers only 6 named boundary tiles, not the full Canada/Alaska domain. |
| 🔲 Open | The sfmc-observation extension (`bfd293b`) has never been exercised in any experiment tracked here — the 6-month validation run predates it. |
| 🔲 Open | Hudson Bay residual degradation (peat QC shrinks it ~41% but doesn't close it) has no identified further cause yet. |

## Appendix A — Commit ledger

Every commit in the ASCAT H SAF / peatland-QC chain, oldest first, as it stood
immediately pre-merge.

| Hash | Date | Author | Message |
|---|---|---|---|
| `4f97d0d` | 2026-06-12 | amfox37 | Add H SAF ASCAT SSM reader for H121 CDR and H139 ICDR (v8, 12.5 km) |
| `2e7034e` | 2026-06-12 | amfox37 | Use obs_dir_hier=1 in H SAF ASCAT reader flist calls |
| `4e52094` | 2026-06-12 | amfox37 | Remove snow/frozen QC from H SAF reader; defer to model-based QC |
| `8e21406` | 2026-06-15 | amfox37 | Fix fill-value QC in H SAF ASCAT reader; update MetOp-B/C data paths |
| `1889ba3` | 2026-06-16 | amfox37 | Accept fill (−128) for subsurface_scattering_probability in H SAF ASCAT QC |
| `dda64f0` | 2026-06-20 | amfox37 | Tighten subsfc QC threshold and add SSM sensitivity QC for H SAF ASCAT reader |
| `a9d0ba8` | 2026-06-20 | amfox37 | Reject noisy backscatter (backscatter40_flag bit 4) in H SAF ASCAT reader |
| `41b06ac` | 2026-07-20 | amfox37 | Merge branch 'develop' into feature/amfox/ascat-hsaf-v8 |
| `d11860f` | 2026-07-20 | amfox37 | Document H SAF ASCAT support |
| `3776e08` | 2026-07-21 | Rolf Reichle | Merge branch 'develop' into feature/amfox/ascat-hsaf-v8 |
| `2063791` | 2026-07-24 | Rolf Reichle | added/edited comments; white-space changes |
| `9d05ba5` | 2026-07-28 | amfox37 | Limit H SAF ASCAT flist entries to hourly granules |
| `7f488e2` | 2026-07-28 | amfox37 | Add ASCAT peatland footprint QC |
| `823a839` | 2026-07-29 | amfox37 | Fix ASCAT peatland QC mwRTM handling |
| `2216974` | 2026-08-03 | amfox37 | Bundle ASCAT peat QC halo exchange |
| `f9a7e69` | 2026-08-03 | amfox37 | Use catchment soil class for ASCAT peat QC |
| `e093f09` | 2026-08-03 | amfox37 | Use composite soil class for ASCAT peat QC |
| `5effa7e` | 2026-08-03 | amfox37 | Use PEATCLSM porosity for ASCAT peat QC |
| `2aead1a` | 2026-08-03 | amfox37 (+ Claude Opus 5) | Scope ASCAT peat QC to peat QC; fix stale comments |
| `bfd293b` | 2026-08-11 | amfox37 | Apply peatland QC to sfmc observations |
| `5a807cd` | 2026-08-11 | amfox37 | Document peatland QC in changelog |
| `e7b1647` | 2026-08-12 | — | **Merge PR #186** into develop |

## Appendix B — Key SLURM jobs

Milestones and incidents, not every resubmit in every retry chain.

| Job | Experiment | Date | Outcome |
|---|---|---|---|
| 56750413 | OLv7_M36_MULTI_type_13_H121 (6yr OL) | 2026-06-23 | completed |
| 56909732 / 56909815 | DAv7 H121 / legacy launch | 2026-07-02 | NUM_SGMT timeout |
| 56934953 / 56934954 | DAv7 H121 / legacy resubmit | 2026-07-04 | checkpoint stall |
| 56979409 / 56979410 | DAv7 H121 / legacy (live-captured) | 2026-07-06 | deadlock, root-caused live |
| 56996819 / 56996820 | DAv7 H121 / legacy, clean restart | 2026-07-06 | chain unblocked |
| 57274390 / 57274391 | DAv7 H121 / legacy, post-quota | 2026-07-22 | 92–96% complete |
| 57563557 | DAv7 H121 FOV12p5 | 2026-07-29 | Cannot Find ENSAVG |
| 57596551 | DAv7 H121 FOV12p5 (baseline) | 2026-07-30 | completed — pre-QC baseline |
| 57605282 / 57615209 | peatlandqc 6mo, 1st attempt | 2026-07-30/31 | EXPID gotcha — unusable diagnostics |
| 57618217 | peatlandqc 6mo, relaunch | 2026-07-31 | shared-checkout binary collision |
| 57640668 | peatlandqc, post-collision | 2026-08-01 | completed (246b103) |
| 57695424 | peatlandqc, final rebuild | 2026-08-04 | completed (95f0317) — used for all results |

---

*Reconstructed from project memory, live git history (`GEOSldas_GridComp`), and
the GitHub API. Companion document: `ascat_reprocessing_impact_report.md`
(science report).*
