# H121/H139 ASCAT reader QC updates (sensitivity, subsurface scattering, backscatter quality)

Status update, 2026-06-26: these changes are implemented in
`GEOSldas_GridComp` on branch `feature/amfox/ascat-hsaf-v8`. This note is
kept as the implementation checklist and rationale; older branches may still
lack this reader or carry the pre-update QC.

## File / subroutine

`src/Components/@GEOSldas_GridComp/GEOSlandassim_GridComp/clsm_ensupd_read_obs.F90`

`read_obs_sm_ASCAT_HSAF` — covers both H121 CDR and H139 ICDR (same file-format reader).

## Why

Hahn et al. (2026, ESSD discussion paper essd-2025-746, Sect. 4.1.1) state that
`surface_soil_moisture_sensitivity` values below 1 dB indicate unreliable
retrievals (typically dense vegetation with low backscatter signal
variation). Separately, the same paper's own validation methodology
(Sect. 3.5) masks `subsurface_scattering_probability > 5%`. The original
H121/H139 reader draft did not screen sensitivity and used `thr_subsfc = 10.`,
twice as permissive as what the dataset's authors validated against.

Tested locally against a 10-day global sample (all three Metop satellites,
2020-06-01 to 06-10): sensitivity <=1 dB removes ~5% of currently-passing
obs (single-day check); subsfc tightened to 5% removes another ~10.5%
(single-day check). Cross-checked against legacy BUFR as a (non-ground-truth)
baseline over the full 10-day window — that comparison was inconclusive
on its own (see `legacy_vs_h121_qc_flags.md` Sect. 7, 9), which is expected
given legacy isn't ground truth; none of the three changes below should be
read as validated-or-rejected by that comparison alone. See
`projects/ascat_da/report/legacy_vs_h121_qc_flags.md` for the full writeup
and rationale, including a third change (Sect. 9 below) decided on
assimilation-design grounds rather than that cross-check.

## Implemented Changes

### 1. Add a new threshold parameter

Add near the existing `thr_wetland`, `thr_topo`, `thr_subsfc` parameters:

```fortran
real, parameter :: thr_sens = 1.0    ! reject if SSM sensitivity <= 1 dB
```

### 2. Tighten the existing subsurface-scattering threshold

```fortran
real, parameter :: thr_subsfc = 5.   ! was 10.
```

### 3. Add a scale-factor parameter for the new variable

`surface_soil_moisture_sensitivity` is `NC_INT`, `scale_factor=1e-7`, units
dB, no offset (per PUM Table 4.12):

```fortran
real*8, parameter :: sens_scale = 1.0d-7
```

### 4. Declare a varid and a raw array

Alongside the existing `subsfc_varid` / `subsfc_raw`:

```fortran
integer :: sens_varid
integer, allocatable :: sens_raw(:)   ! surface_soil_moisture_sensitivity (NC_INT)
```

### 5. Variable-ID lookup

Add to the block that does `nf90_inq_varid` calls:

```fortran
ierr = nf90_inq_varid(ncid, 'surface_soil_moisture_sensitivity', sens_varid)
```

### 6. Allocate / read / deallocate

Alongside the other per-file arrays:

```fortran
allocate(sens_raw(N_obs_file))
...
ierr = nf90_get_var(ncid, sens_varid, sens_raw)
...
deallocate(sens_raw)   ! add to the existing
                       ! deallocate(wetland_raw, topo_raw, subsfc_raw) line
```

### 7. QC loop: add a new skip check

The fill value for this field is `-2^31`, so `sens_raw(ii) * sens_scale` for
a fill value is hugely negative and will naturally fail the `<= thr_sens`
test — no separate fill-value branch needed:

```fortran
! skip if sensitivity too low (or fill)
if (real(sens_raw(ii)) * real(sens_scale) <= thr_sens) cycle
```

### 8. Update the subroutine's header comment block

Near the existing "QC applied" list (around line ~2172), add:

```
!   surface_soil_moisture_sensitivity <= thr_sens (1 dB)      -> reject (new)
```

and update the existing `subsurface_scattering_probability >= thr_subsfc`
line to note the threshold is now 5%, not 10%. (Also add the
`backscatter40_flag` line described below.)

---

## 9. Third change: reject noisy backscatter (`backscatter40_flag` bit 4)

### Why

`backscatter40_flag` (PUM Table 4.15) carries three independent bits:
bit 1 = `sigma0_usable`, bit 2 = `sigma0_slightly_degraded`, bit 4 =
`sigma0_noise_out_of_limits`. None of these are currently read or screened.

This matters because **`read_obs_sm_ASCAT_HSAF` assigns every accepted
observation the same fixed `this_obs_param%errstd`** (line ~2618:
`ASCAT_sm_std(ii) = this_obs_param%errstd / 100.`), regardless of flag
status — there is no per-observation noise down-weighting anywhere in this
reader (or in any other obs-type reader in this file; `errstd` is a static
per-species namelist constant throughout). So a `noise_out_of_limits`
observation currently gets exactly the same weight in the Kalman update as
a clean one — hard rejection is the only available lever against it.

We tested this locally (10-day global sample, `projects/ascat_da/`) and
the legacy-cross-check was inconclusive on its own (see
`legacy_vs_h121_qc_flags.md` Sect. 10) — but given there is no fallback
down-weighting mechanism, we are recommending the conservative choice:
reject on `noise_out_of_limits` (bit 4) only. We are **not** recommending
also rejecting `slightly_degraded` (bit 2) — testing showed it discards
~3 more percentage points of global data for no measurable additional
benefit over bit 4 alone.

### Changes

Add a threshold/mask parameter near `thr_sens`/`thr_subsfc`:

```fortran
integer(1), parameter :: bsflag_bad_bits = 4_1   ! noise_out_of_limits
```

Declare a varid and raw array — `backscatter40_flag` is `NC_UBYTE` (same
type pattern as `sflag_raw`/`pflag_raw`, i.e. `integer(1)` in Fortran):

```fortran
integer :: bsflag_varid
integer(1), allocatable :: bsflag_raw(:)   ! backscatter40_flag (NC_UBYTE)
```

Variable-ID lookup, alongside the other `nf90_inq_varid` calls:

```fortran
ierr = nf90_inq_varid(ncid, 'backscatter40_flag', bsflag_varid)
```

Allocate / read / deallocate, alongside the other per-file arrays:

```fortran
allocate(bsflag_raw(N_obs_file))
...
ierr = nf90_get_var(ncid, bsflag_varid, bsflag_raw)
...
deallocate(bsflag_raw)   ! add to the existing deallocate(sflag_raw, pflag_raw) line
```

QC loop: add a new skip check, in the same `iand`-bitmask style already
used for `surface_flag`/`processing_flag`:

```fortran
! skip if backscatter at 40 deg is noisy
if (iand(bsflag_raw(ii), bsflag_bad_bits) /= 0_1) cycle
```

Update the header "QC applied" comment block to add:

```
!   backscatter40_flag bit 4 (noise_out_of_limits)            -> reject (new)
```

## Before promoting/merging

These three changes came from a mix of the paper's stated recommendations
(sensitivity, subsurface scattering) and an assimilation-design argument
specific to this codebase (backscatter noise has no down-weighting
fallback, so hard QC is the only lever) — not from tuning against GEOSldas
innovations directly. Please validate with the usual OFA-bias/match-fraction
check before promoting beyond the feature branch (the Fortran/production-side equivalent of
`check_global_superobs_vs_ofa.py`) in case any of the three thresholds need
adjusting for this use case. The Python-mirror side of all three changes is
implemented and validated in `projects/ascat_da/` (see
`legacy_vs_h121_qc_flags.md`); cross-check against those results before
either implementation ships.
