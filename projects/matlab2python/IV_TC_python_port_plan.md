# IV / TC Python port — working plan

Status as of 2026-07-15. Written for whoever (human or AI) picks this up next,
possibly in a different Claude Code session on a different machine — this
file is the durable record, not chat history.

## Why

Three new production-scale HSAF CDR comparison experiments are running on
Discover (`/discover/nobackup/projects/land_da/hsaf_cdr_test/`), all
2015-04-01 to 2021-04-01, all using the per-species z-score ASCAT scaling
stats validated in `projects/obs_scaling_params/` (see that project's O-F
verification results):

| Role | EXP_ID | Status as of 2026-07-15 |
|---|---|---|
| OL baseline | `OLv7_M36_MULTI_type_13_H121` | finished (it's also the source run for the scaling stats) |
| DA, H121 CDR ASCAT assimilated | `DAv7_M36_ASCAT_type_13_H121` | RUNNING (job 57140802) |
| DA, legacy EUMETSAT ASCAT assimilated | `DAv7_M36_ASCAT_type_13_legacy` | RUNNING (job 57140819) |

Once these finish we'll want independent validation (IV) and triple
collocation (TC) skill metrics for a **3-way comparison** (OL vs DA-H121 vs
DA-legacy), the same kind of analysis `common/matlab/IVs/` and
`common/matlab/TC/` already do for the CYGNSS OLv8/DAv8 experiments. The
MATLAB pipeline is 2-way (`D1` vs `D2`) and hardcoded per dataset; doing this
comparison well means finishing the Python port instead of hand-adapting
more one-off MATLAB scripts.

## Inventory: what exists today

**MATLAB (`common/matlab/`), the "current" pipeline** (see
`common/matlab/TC/Compute_TC_two_input_files.m` /
`Compute_TC_three_input_files.m` and `common/matlab/IVs/scripts/step{2,3,4,5}/`):
- Step2: daily obs/model pair extraction, one script per obs source
  (`step2_asc` ASCAT L4/BUFR, `step2_smosic` SMOS-IC, `step2_cyg` CYGNSS,
  `step2/Save_SMPL3_LDAS_tavg24_nc4_daily.m` SMAP L3).
- Step3: pentad climatology (`Compute_clim*.m`).
- Step4: IV skill (IVD/IVS lag-2-day anomaly correlation).
- Step5: R-diff (skill difference between two experiments).
- TC: 2- and 3-input triple collocation error-variance solvers
  (`Compute_TC_two_input_files.m`, `Compute_TC_three_input_files.m`) —
  genuinely generic over which datasets are supplied.
- SMAP L3 product in use: **SPL3SMP v009 (R19240)**, the current NSIDC
  release — confirmed live on disk at
  `/discover/nobackup/projects/land_da/Evaluation/IVs/data/SPL3SMP_v009/`
  and read by `step2/Save_SMPL3_LDAS_tavg24_nc4_daily.m`. An older
  v008/R18290 path still exists (`step2/Save_SMPL3_LDAS_gph_nc4_daily.m`,
  and everything under `matlab_postprocess/TC/` in this staging dir) but is
  disconnected from the live step2→step3→step4→TC chain — treat it as dead
  reference code, not something to port.

**Python (this repo), what actually runs today:**
- `projects/SMOS_IC/notebooks/ivs_tc_ascat_smosic_python_workflow.ipynb` —
  a working, *notebook-only* port of steps 2–5 + TC, but:
  - Hardcoded to exactly two experiments (`OLv8_M36_cd`, `DAv8_M36_cd`,
    the CYGNSS-era run names) and exactly the ASCAT+SMOS-IC+model triplet.
  - TC solver (`run_tc()`, cell 7) is the 3-input case only, inlined,
    not reusable.
  - No SMAP L3 reader at all — SMAP only exists on the MATLAB side.
  - Already executed once for real (2018-08-01 to 2024-06-30, OLv8/DAv8
    CYGNSS-era comparison) — outputs cleared, code kept in
    `projects/SMOS_IC/notebooks/`, committed `4c270a9`.
- `common/python/stats/` has no IV/TC code at all (trend analysis,
  Theil-Sen, mask-maker only).
- `projects/matlab2python/matlab_postprocess/TC|IVs/` is reference MATLAB
  copied here for translation, not Python.

## What "tidy up and port" means, concretely

Priority order — do these roughly in sequence, each should be independently
useful:

1. **Extract, don't rewrite.** Pull the working functions out of
   `ivs_tc_ascat_smosic_python_workflow.ipynb` verbatim into
   `common/python/stats/` (new files, e.g. `iv_skill.py` for step2-5,
   `triple_collocation.py` for TC). The notebook becomes a thin driver that
   imports and calls them — same pattern as `common/python/io/read_GEOSldas.py`
   being imported by notebooks today.
2. **Generalize `Config`/run_roots to N experiments**, not just a
   `{OL, DA}` pair. The 3-way OL/DA-H121/DA-legacy comparison needs step5
   R-diff computed pairwise (DA-H121 minus OL, DA-legacy minus OL, and
   probably DA-H121 minus DA-legacy) — don't hardcode a single diff like the
   current `RUN_OL`/`RUN_DA` globals do.
3. **Generalize the TC solver to N inputs**, mirroring the MATLAB split
   between `Compute_TC_two_input_files.m` and `Compute_TC_three_input_files.m`
   rather than hardcoding the ASCAT+SMOS-IC+model triplet. Three
   experiments each need TC against the same obs pair, so the solver should
   take a list of (obs_a, obs_b, model) and not care which model run it came
   from.
4. **Add a SMAP L3 (SPL3SMP v009/R19) reader**, porting
   `step2/Save_SMPL3_LDAS_tavg24_nc4_daily.m`. This is the one real gap
   versus MATLAB — everything else (ASCAT, SMOS-IC) already has a Python
   equivalent. Read from
   `/discover/nobackup/projects/land_da/Evaluation/IVs/data/SPL3SMP_v009/`;
   don't touch the disconnected v008/R18290 path.
5. **Point the generalized workflow at the new hsaf_cdr_test experiments**
   once `DAv7_M36_ASCAT_type_13_H121` / `_legacy` finish — this is the
   actual deliverable, everything above is groundwork for it. Output paths
   for these runs follow the standard GEOSldas layout
   (`<exp>/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg/...`), same as the CYGNSS
   runs the notebook already reads.

## Explicit non-goals for now

- Don't port CYGNSS L3 TC support (`Compute_TC_three_input_files.m`'s
  ASCL4_SMPL3_CYGL3 case) — not needed for the OL/DA-H121/DA-legacy
  comparison.
- Don't touch the dead v008/R18290 MATLAB paths — no reason to preserve
  parity with code nothing downstream reads.
- Follow the repo's existing notebook hygiene rules (clear outputs before
  commit, watch for the f-string `\n` / no-cell-id gotchas) — see root
  `CLAUDE.md`.
