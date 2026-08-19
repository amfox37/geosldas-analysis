# TWLAND water-year endpoint reconstruction (Discover brief Task 2)

Prepared 2026-08-19. Companion to `discover_peatclsm_budget_brief.md`
(Task 2: instantaneous `TWLAND` for the storage-endpoint term). This
resolves Task 2 in full for WY2001-WY2006 — no rerun, no approximation.

---

## What this is

Figure 14's water budget residual decomposes into two terms (see the main
brief, Section 1). Term B is a storage-endpoint sampling artefact: the
manuscript's `dStorage` used the difference of *monthly-mean* `TWLAND`
around each water-year boundary, when what the budget actually needs is the
*instantaneous* value at 00Z Oct 1. The brief proposed a fallback — average
the September and October monthly means — that was shown to remove 78% of
the resulting bias, and authorized adopting it **only** because a model
rerun to get true instantaneous states wasn't considered worth the cost for
a term worth ~7% of the total residual.

That trade-off turned out to be unnecessary. True instantaneous states
already exist on disk for this exact window, in `catch_internal_rst`
restarts — nobody needs to choose between "rerun" and "approximate."

## Where the restarts were found

Both experiments' restarts for 2000-2006 live under **predecessor runs**,
not the directories used for the monthly diagnostic archive:

| Experiment | Restart source | Coverage |
|---|---|---|
| OL | `/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rs/` | monthly, 2000-06 → 2018-06 |
| DA | `/discover/nobackup/projects/land_da/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/rs/` | monthly, 2000-06 → 2007-05 |

`M21C_land_sweeper_OLv8_M36` is the run `LS_OLv8_M36_v2` was itself
restarted from (`RESTART_PATH` in `LS_OLv8_M36_v2/LS_OLv8_M36.txt` points
directly at it) — the same v2→v3-style leg structure DA already has. This
was easy to miss because the monthly `tavg24_1d_lnd_Nt` diagnostic archive
under `LS_OLv8_M36_v2` looks continuous back to 2000 (it's been stitched
across legs), while `LS_OLv8_M36_v2`'s own `rs/` directory only holds
restarts from 2018 onward — the early restarts were never copied forward
into the `_v2` directory, only the diagnostics were.

Restarts are per-ensemble-member (24 members, `ens0000`-`ens0023`), no
`ens_avg`. Each restart file already contains every prognostic the
reconstruction formula needs — no separate BCS lookup required:

```
TWLAND = CDCR2/(1-WPWET) - CATDEF + RZEXC + SRFEXC + CAPAC + WESNN1 + WESNN2 + WESNN3
```

`CATDEF`, `RZEXC`, `SRFEXC`, `CAPAC`, `WESNN1-3` are all `kg m-2`; `CDCR2`
is `kg m-2`; `WPWET` is a dimensionless degree-of-saturation fraction — same
formula the brief specified, taken from `lsm_routines.F90`'s
`catch_calc_wtotl`.

## What was computed

For each of OL and DA, `TWLAND` was reconstructed at 00Z Oct 1 for every
year 2000-2006 (7 dates, bounding water years WY2001-WY2006), as the
straight mean of the formula applied to all 24 ensemble members
individually (not applied to an already-ensemble-averaged restart, since no
such restart exists — the mean is taken *after* computing `TWLAND` per
member).

## Files delivered

```
OLv8_twland_wy_endpoints_2000_2006.nc   sha256: cd9974f097f70e96a91631e67dacb6811fe9f3e0f0f36a6e5c8f572d773bb39f
DAv8_twland_wy_endpoints_2000_2006.nc   sha256: a5958befeef3958323e66592171a255d64abd3f854a7d4ac02a0e74122cd3b21
```

Location: `projects/M21C_ls/output/monthly_flux_states/` (repo-local,
Discover/gpfsm filesystem — same transfer caveat as the peat-FSW files:
copy to `/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/`
via `scp` when ready, no direct write access to that path from here).

Structure: `TWLAND(wy_boundary=7, tile=112573)`, float32, zlib complevel 4,
`lat`/`lon` copied verbatim from the audited `OLv8_water_budget_...`
product (same tile order guarantee as the peat-FSW files), `time` as CF
`days since 2000-06-01` giving Oct-1-of-each-year, global attributes record
the exact restart source path, ensemble size, and formula used.

## Domain-mean TWLAND at each boundary (kg m⁻², all 112,573 land tiles)

| Date | OL | DA |
|---|---:|---:|
| 2000-10-01 | 622.379 | 625.096 |
| 2001-10-01 | 615.195 | 619.651 |
| 2002-10-01 | 609.087 | 616.887 |
| 2003-10-01 | 616.696 | 625.891 |
| 2004-10-01 | 619.112 | 629.417 |
| 2005-10-01 | 616.740 | 627.053 |
| 2006-10-01 | 615.297 | 625.771 |

## Accuracy vs. the Sep/Oct-mean approximation — aggregate is fine, annual isn't

Comparing the true restart-based (DA−OL) annual `dStorage` against the
brief's authorized fallback (mean of Sep and Oct monthly-mean `TWLAND`),
per water year:

| WY | true Δstorage (kg m⁻²) | approx Δstorage (kg m⁻²) | diff |
|---|---:|---:|---:|
| WY2001 | 1.739 | 1.812 | −0.074 |
| WY2002 | 3.344 | 3.135 | +0.209 |
| WY2003 | 1.395 | 1.647 | −0.252 |
| WY2004 | 1.110 | 1.082 | +0.028 |
| WY2005 | 0.009 | −0.076 | +0.086 |
| WY2006 | 0.160 | −0.006 | +0.166 |
| **6-yr aggregate** | **7.756** | **7.594** | **+0.163** |

Annualized 6-year mean: true 1.293 vs approx 1.266 kg m⁻² yr⁻¹ — only ~2%
apart, consistent with the brief's "78% bias reduction" claim holding up in
aggregate.

**But individual years don't fit nearly as well.** WY2003 is off by 0.25 on
a true value of 1.40 (~18% relative error); WY2005 and WY2006 are small
enough in magnitude that the approximation gets the wrong sign entirely.
The approximation's year-to-year errors are effectively uncorrelated noise
that partially cancels over 6 years — which is exactly why the 6-year
number looked fine while no single year would have.

## Recommendation

**Use the restart-based `TWLAND` values in these two files directly as
Term B's storage endpoint for WY2001-WY2006, in both the 6-year aggregate
and any water-year-resolved figure/table.** Don't fall back to the Sep/Oct
approximation for this window — it's no longer needed, and it's measurably
worse per-year even though it's fine in aggregate. The approximation
remains the right (and only available) choice outside this window, where
no restarts exist: DA post-2007 (`LS_DAv8_M36_v3` has no restarts at all —
its `rs/` output is tarred but empty of loose files under the checked
paths) and OL post-2018 in the sense that `LS_OLv8_M36_v2`'s own restarts
only start there (though `M21C_land_sweeper_OLv8_M36` restarts do extend to
mid-2018, so there is limited room to extend this exact method a bit
further forward if a later figure ever needs it).
