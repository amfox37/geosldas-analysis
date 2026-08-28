# Discover task brief: PEATCLSM free-standing water and the snow-DA budget residual

Prepared 2026-08-19. Audience: an assistant with access to NASA Discover and the
M21C land-sweeper archive. Everything below was established locally from the
compressed monthly products and the GEOSldas source; the tasks are the parts
that require Discover.

---

## 1. Background: what we found and why it matters

Manuscript Figure 14 presents a differential water budget for MODIS snow-cover
assimilation over WY2001–WY2006:

    I_snow = dRunoff + dET + dStorage + residual        (all DA − OL)

The residual is persistently negative, −2.67 kg m⁻² yr⁻¹, or −4.6% of the
58.2 kg m⁻² yr⁻¹ snow-DA input, and is currently unexplained in the caption.

We traced it. It decomposes exactly (to 1.7e-8) into two terms:

| term | six-year mean | share |
|---|---:|---:|
| **A** = dET + dRunoff + ∫dWCHANGELAND | +1.69 kg m⁻² yr⁻¹ | 63% |
| **B** = dStorage − ∫dWCHANGELAND − I_snow | +0.98 kg m⁻² yr⁻¹ | 37% |

residual = −(A + B).

**Term A is PEATCLSM free-standing water.** In `catchment.F90` (~line 1305):

```fortran
FSW_CHANGE(N) = 0.
IF(POROS(N) >= PEATCLSM_POROS_THRESHOLD) THEN
   pr = trainc(n)+trainl(n)+tsnow(n)+tice(n)
   FSW_CHANGE(N) = PR - EVAP(N) - RUNOFF(N) - WCHANGE(N)
ENDIF
```

`PEATCLSM_POROS_THRESHOLD = 0.90` (`catch_constants.F90:152`). The free-standing
water store is deliberately excluded from `WTOT`, which `catch_calc_wtotl`
(`lsm_routines.F90`) builds as

    wtotl = cdcr2/(1-wpwet) - catdef + rzexc + srfexc + capac + sum(wesnn)

so on peat tiles `P - ET - runoff - WCHANGE` is non-zero **by construction**.

We verified this against porosity from
`GEOSldas_diagnostics/test_data/clsm/soil_param.dat` (112,573 rows, column 7 =
POROS, matching the M36 tile set):

- 4,189 tiles globally have POROS >= 0.90 (3.72%)
- of the 3,359 tiles with a water-budget non-closure worse than −1000 kg m⁻²
  over 24 years, **100.0% are peat**
- peat tiles carry **99.2%** of the total area-weighted non-closure
- mean non-closure: peat **−137.4**, non-peat **−0.0** kg m⁻² yr⁻¹
- excluding peat, the OL and DA budgets both close to **−0.27 kg m⁻² over six
  years** (7e-5 of precipitation, i.e. float32 rounding)

Effect on Figure 14 (domain-mean recomputation, not a pipeline rerun):

| | residual | % of input |
|---|---:|---:|
| as published | −2.673 | −4.6% |
| peat excluded | −0.186 | −0.3% |
| peat excluded + endpoint fix | −0.041 | −0.07% |

**Term B is a storage-endpoint problem.** `dStorage` is taken as the difference
of *September monthly means* of `TWLAND`, while fluxes are integrated
1 Oct – 30 Sep. A monthly mean approximates mid-month, so the estimate is offset
by ~2 weeks at each end.

**Correction, 19 Aug 2026.** An earlier version of this brief stated that using
the mean of the September and October monthly means reduces the bias by 78%.
**That figure was wrong.** It was computed on a peat-excluded domain and does
not generalise. On the full seasonal-snow domain all monthly-mean estimators
perform about equally, and the Sep/Oct midpoint is marginally *worse* than the
September mean:

| estimator | bias | RMS |
|---|---:|---:|
| September mean (current) | +0.983 | 1.252 |
| Sep/Oct midpoint | +1.143 | 1.196 |
| cubic spline at 1 Oct | +1.126 | 1.159 |

kg m⁻² yr⁻¹, DA−OL, against `∫WCHANGELAND + I_snow`. Because the choice of
sampling point barely changes the answer, Term B is not a timing artefact and
**cannot be fixed by a better monthly-mean estimator**. It needs true
instantaneous endpoints — see Task 2.

---

## 2. Task 1 (highest priority): extract `PEATCLSM_FSWCHANGE`

If this field is in the daily archive, the budget can be closed properly with
peatlands retained, instead of excluding 6.7% of the seasonal-snow domain.

**It is requested in the HISTORY collection we use.** `GEOSldas_HIST.rc:142`
lists `'PEATCLSM_FSWCHANGE' , 'GridComp'` inside `tavg24_1d_lnd_Nt.fields:` —
the same collection the compressed monthly products were built from. Its
long name elsewhere in that file is `free_surface_water_on_peat_flux`.

### Steps

1. **Confirm the field is present** in the daily files:

   ```
   /discover/nobackup/projects/land_da/M21C_land_sweeper/LS_OLv8_M36_v2/
       LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/ens_avg/Y<YYYY>/M<MM>/
       *tavg24_1d_lnd_Nt*.nc4
   ```

   and the DA equivalent (`LS_DAv8_M36_v2` / `LS_DAv8_M36`). Report the exact
   variable short name, long name, units and sign convention. **If absent, stop
   and report** — the field was never written and recovering it needs a model
   rerun, which is a separate decision.

2. **Build monthly means** for OL and DA, one value per tile per month, on the
   112,573-tile M36 set. Minimum required period is Oct 2000 – Sep 2006
   (WY2001–WY2006, the Figure 14 window). Preferred is the full record,
   June 2000 – May 2024, to match the other monthly products.

   The local analogue of the processing code is
   `projects/utils/notebooks/process_tavg24_1d_lnd_Nt_files_monthly_040125.ipynb`
   in the `geosldas-analysis` repo: it globs the daily files per Y/M directory,
   selects variables by long name, means over `time`, and writes one file per
   month. The production version that made the existing archive lives on
   Discover — its provenance is recorded in the global attributes of
   `OLv8_water_budget_2000_2024_compressed.nc`:
   `Date = 2024-12-21`, `file_prefix = LS_OLv8_M36`,
   `note = 'Monthly means from tavg24_1d_lnd_Nt.*; lat/lon attached on tile.'`
   Reuse that code path so the output is byte-compatible with the existing
   products.

3. **Deliver** as a standalone pair of files, one per experiment, matching the
   existing products exactly. Full specification in Section 2a below.

4. **Verify before delivering**: on peat tiles (POROS >= 0.90), confirm

       PRECTOTCORRLAND - EVLAND - RUNSURFLAND - BASEFLOWLAND - WCHANGELAND
           ≈ PEATCLSM_FSWCHANGE

   to within rounding, and that it is ~0 on non-peat tiles. Report the residual
   of that identity. This is the decisive test of the whole diagnosis: if it
   holds, the budget can be closed exactly with a fifth term.

---

## 2a. Packaging specification for local analysis

The local analysis validates every input file against a machine-readable
contract (`config/trend_breakpoint_inputs.json`) before any variable is read.
A file that does not match is rejected outright, so please follow this exactly.

### File naming and placement

Two new files, one per experiment, alongside the existing products:

    OLv8_peat_fsw_2000_2024_compressed.nc
    DAv8_peat_fsw_2000_2024_compressed.nc

**Do not append to the existing `*_water_budget_*.nc` files.** Those are audited
and checksummed in place and are the input to accepted, committed results; a new
family is cleaner and avoids invalidating them.

### Structure — must match the existing files exactly

| property | required value |
|---|---|
| dimensions | `time` = 288, `tile` = 112573, in that order |
| variable dims | `(time, tile)` |
| coordinates | `time` (datetime64), `lat` (tile, float32), `lon` (tile, float32) |
| time axis | monthly, `2000-06-01` through `2024-05-01`, first of month |
| tile ordering | identical to `LS_OLv8_M36.ldas_tilecoord.bin`; do not re-sort |
| dtype | float32 |
| compression | zlib, complevel 4, chunksizes `(12, 100000)` |
| fill value | NaN |

### Variable

Short name `PEATCLSM_FSWCHANGE`, with attributes in the same style as the
existing water-budget variables:

    long_name    : free_surface_water_on_peat_flux
    units        : kg m-2 s-1          # confirm against the daily files
    cell_methods : time: mean

Keep the native flux units (`kg m-2 s-1`) rather than pre-converting to monthly
totals — the local loader performs that conversion itself and requires the rate
units to trigger it. For reference, `WCHANGELAND` in the existing file carries
exactly `units = 'kg m-2 s-1'`, `cell_methods = 'time: mean'`.

Non-peat tiles should be **0.0**, not NaN, since the model sets `FSW_CHANGE = 0`
there. If they are written as NaN instead, say so explicitly — it changes the
masking on our side.

### Global attributes

Mirror the existing convention so provenance is traceable:

    source_root : <full Discover path to the cat/ens_avg directory used>
    file_prefix : LS_OLv8_M36   (or LS_DAv8_M36)
    note        : Monthly means from tavg24_1d_lnd_Nt.*; lat/lon attached on tile.
    Date        : <creation date>

### Transfer and verification

Provide the files plus a **SHA-256 checksum for each**. They will be placed in

    /Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2/

alongside the existing compressed products.

Please also report, so we can validate on arrival without guesswork:

1. total number of tiles with any non-zero `PEATCLSM_FSWCHANGE` over the record
   (we expect ~4,189, matching POROS >= 0.90);
2. the area-weighted 24-year mean, in kg m-2 yr-1, over all valid land and over
   peat tiles only;
3. the residual of the peat identity check from Task 1 step 4.

### What we will do with it

Add a `peat_fsw` dataset entry to `config/trend_breakpoint_inputs.json` and a
`PEATCLSM_FSWCHANGE` row to the variable-selection registry, then re-run the
water-year budget with a fifth closing term. If the identity holds, Figure 14's
residual should collapse from −2.67 to the storage-endpoint term alone,
with peatlands retained in the domain.

---

## 3. Task 2: restart-derived `TWLAND` at water-year boundaries

`TWLAND` appears only in time-averaged collections, so a true endpoint value is
not available from HISTORY output. It can, however, be reconstructed from
restarts, and that work is understood to be in progress.

**Deliver per tile, not as a domain mean.** Figure 14 uses the Northern
Hemisphere seasonal-snow domain (48,067 tiles: 20°N southern limit,
snow-cover-possible thresholds of 0.05 SCF and 5 kg m⁻² snow mass, and mean JJA
SCF < 0.20 to exclude permanent ice) — not global valid land. Per-tile delivery
lets us apply that mask, and any future mask, locally.

### Specification

| property | required value |
|---|---|
| variable | `TWLAND`, **instantaneous**, ensemble mean over all 24 members |
| dimensions | `time` = 7, `tile` = 112573, variable on `(time, tile)` |
| times | 0000 UTC on 1 October, 2000 through 2006 |
| tile ordering | identical to `LS_OLv8_M36.ldas_tilecoord.bin` |
| units | `kg m-2`, float32 |
| cell_methods | omit, or `time: point` — **not** `time: mean` |
| coordinates | `time`, `lat` (tile), `lon` (tile) |
| filenames | `OLv8_twland_wy_boundaries_compressed.nc`, `DAv8_twland_wy_boundaries_compressed.nc` |

Same global-attribute convention as the other products, plus a SHA-256 per file
and a note recording which restart files were used.

### Why per-year accuracy matters

Figure 14 panel (a) is a six-bar, per-water-year figure. Monthly-mean
approximation errors are largely random in sign and cancel over six years, so an
aggregate can look acceptable while individual years are poor — errors of ~18%
on a single year, and sign flips where the true value is small. Exact endpoints
are therefore required for the figure as drawn, not merely a refinement.

The monthly-mean fallback remains relevant only for periods where restarts do
not exist, should a future figure extend beyond WY2006.

---

## 4. Task 3: confirm the run configuration

Short checks, to close off assumptions we made from the source tree rather than
from the actual experiment:

1. Confirm **PEATCLSM is active** in `LS_OLv8_M36_v2` and `LS_DAv8_M36_v2`, and
   report the GEOSldas tag/version used.
2. Confirm `PEATCLSM_POROS_THRESHOLD = 0.90` in the version that was run.
3. Confirm the `soil_param.dat` used by the experiment matches our local copy at
   `GEOSldas_diagnostics/test_data/clsm/soil_param.dat` — compare the porosity
   column, or report a checksum. Our peat count is 4,189 of 112,573 tiles.
4. Confirm OL and DA share identical `soil_param.dat` and boundary conditions.

---

## 5. What to report back

1. Whether `PEATCLSM_FSWCHANGE` is in the daily archive; if so its short name,
   long name, units and sign convention.
2. Paths to any new monthly products created, and the code used.
3. The residual of the peat identity check in Task 1 step 4.
4. Which instantaneous prognostics are available, at what frequency (Task 2).
5. Answers to the four configuration questions in Task 3.
6. Anything that contradicts the diagnosis in Section 1 — particularly if the
   budget fails to close on non-peat tiles at source, which would mean the
   problem is in the monthly compression rather than in the model.

## 6. What is NOT being asked

- No model rerun.
- No change to any analysis, figure or manuscript text. This is a data and
  verification request only.
- No reprocessing of the existing compressed products beyond adding the new
  field; the current ones are audited and in use.
