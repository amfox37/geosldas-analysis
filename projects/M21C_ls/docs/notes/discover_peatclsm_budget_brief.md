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

**Term B is a storage-endpoint sampling artefact.** `dStorage` is taken as the
difference of *September monthly means* of `TWLAND`, while fluxes are integrated
1 Oct – 30 Sep. A monthly mean approximates mid-month, so the estimate is offset
by ~2 weeks at each end. Using the mean of the September and October monthly
means instead reduces the bias by 78% (−0.186 → −0.041 kg m⁻² yr⁻¹) and the RMS
error from 1.010 to 0.741.

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

3. **Deliver** in the same form and naming convention as the existing
   `*_water_budget_2000_2024_compressed.nc` files, ideally appended to those
   files as an extra variable, with tile ordering and time axis identical.

4. **Verify before delivering**: on peat tiles (POROS >= 0.90), confirm

       PRECTOTCORRLAND - EVLAND - RUNSURFLAND - BASEFLOWLAND - WCHANGELAND
           ≈ PEATCLSM_FSWCHANGE

   to within rounding, and that it is ~0 on non-peat tiles. Report the residual
   of that identity. This is the decisive test of the whole diagnosis: if it
   holds, the budget can be closed exactly with a fifth term.

---

## 3. Task 2: instantaneous `TWLAND` for the endpoint term

`TWLAND` appears **only** in time-averaged collections (`tavg24_1d_lnd_Nt`,
`tavg24_2d_lnd_Nx`). There is no instantaneous collection carrying it, so a
true endpoint value is not directly available.

Two options, in order of preference:

**(a) Reconstruct it** from instantaneous prognostics using the model's own
formula:

    TWLAND = cdcr2/(1-wpwet) - catdef + rzexc + srfexc + capac + sum(wesnn)

`cdcr2` and `wpwet` are static (in `clsm/` parameter files); the prognostics are
in the restart files and possibly in `inst3_1d_lndfcstana_Nt` or
`catch_progn_incr`. Please check which of `CATDEF`, `RZEXC`, `SRFEXC`, `CAPAC`,
`WESNN` are archived instantaneously and at what frequency, for both OL and DA.
A single value per 1 October and 30 September would be sufficient.

**(b) If reconstruction is not possible**, report that, and we will adopt the
September/October two-month mean approximation, which we have already shown
removes 78% of the bias.

Do **not** trigger a model rerun for this. It is worth ~7% of the residual once
peat is handled, and is not worth the cost on its own.

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
