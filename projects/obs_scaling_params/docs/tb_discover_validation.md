# SMAP Tb Tile-Space Port — Discover Validation (2026-07-31)

Pre-production validation of `scripts/run_tb_scaling_params.py` (the SMAP Tb
tile-space scaling-parameter port, commit `cd3eb01`) on NASA Discover,
against real GEOS-LDAS experiment archives and the unmodified MATLAB
reference implementation, culminating in a real GEOS-LDAS run consuming the
generated scaling files.

**Result: parity confirmed (byte-identical MATLAB/Python pentad files, both
orbits), GEOS-LDAS consumption confirmed (correct QC gating and exact
scaling-formula application). No divergence found at any step.**

Supporting scripts referenced below are preserved under
[`tb_discover_validation_scripts/`](tb_discover_validation_scripts/). They
were run from a scratch working directory, not committed as part of the
package (see "Reproducing this validation" below) — this doc + that folder
is what makes the run reproducible, not repo state alone.

---

## 1. Checkout and unit tests

```bash
cd /discover/nobackup/projects/land_da/geosldas-analysis
git pull --ff-only origin main
git rev-parse --short HEAD   # cd3eb01

cd projects/obs_scaling_params
module purge && module load python/GEOSpyD/25.3.1-0/3.13
python -m pytest -q tests
# 10 passed
```

## 2. Experiment selection

No single live Discover archive satisfied every "preferred" criterion (M09
native domain + ≥3 live contiguous months + current obsparam schema +
`.bin` ObsFcstAna, all at once). Two real experiments were used, confirmed
to share the **identical M36 tile space** (`read_tilecoord` byte-for-byte
match on `tile_id`, `com_lon`, `com_lat`; 112,573 tiles in every case):

| Purpose | EXP_PATH / EXP_RUN | Domain | ObsFcstAna format | Obsparam schema |
|---|---|---|---|---|
| §3 CLI smoke test (live, current production format) | `hsaf_cdr_test/OLv7_M36_MULTI_type_13_H121` | `SMAP_EASEv2_M36_GLOBAL` | `.nc4` | current (`fcstvarname`/`fcstunits` present) |
| §4-6, §8 MATLAB parity | `Experiment_archive/M21C_land_sweeper_OLv8_M36/LS_OLv8_M36` (`.bin` ObsFcstAna extracted from `output/.../SMAP_EASEv2_M36_GLOBAL_ana_tars/Y2016.tar`, months `M01`-`M03`, into a private tree at `/discover/nobackup/projects/land_da/obs_scaling_test/`) | `SMAP_EASEv2_M36_GLOBAL` | `.bin` | oldest (no `fcstvarname`/`fcstunits`) |
| §7 GEOS-LDAS consumption | `hsaf_cdr_test/hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE` (new, cloned from `hsaf_cdr_test_DAv8_M36`), scaling stats sourced from `OLv7_M36_MULTI_type_13_H121` | `SMAP_EASEv2_M36_GLOBAL` | `.nc4` | current |

All live experiments found are **M36-native** (`gridtype='EASEv2-M36'` /
`'EASEv2_M36'`), not M09. `convert_grid='EASEv2_M36'` degenerates to a 1:1
identity mapping (`n_administering_tile == n_model_tile == 112573`) in this
case — the administering-tile *code path* is fully exercised, but not its
M09→M36 many-to-one behavior. No live M09 archive was found on Discover in
the time available. This is the one open item against the "preferred" spec.

**Species confirmed from each experiment's own obsparam** (never assumed):

| descr | species (OLv7 / M21C) | orbit | pol | angle |
|---|---|---|---|---|
| `SMAP_L1C_Tbh_A` | 1 / 5 | 1 (A) | 1 (H) | 40.0° |
| `SMAP_L1C_Tbh_D` | 2 / 6 | 2 (D) | 1 (H) | 40.0° |
| `SMAP_L1C_Tbv_A` | 3 / 7 | 1 (A) | 2 (V) | 40.0° |
| `SMAP_L1C_Tbv_D` | 4 / 8 | 2 (D) | 2 (V) | 40.0° |

## 3. Descending-orbit Python CLI test (live current-format archive)

```bash
cd projects/obs_scaling_params
module purge && module load python/GEOSpyD/25.3.1-0/3.13
python -u scripts/run_tb_scaling_params.py \
  --exp-path /discover/nobackup/projects/land_da/hsaf_cdr_test \
  --exp-run OLv7_M36_MULTI_type_13_H121 \
  --domain SMAP_EASEv2_M36_GLOBAL \
  --start 2016-01 --end 2016-03 --run-months 1,2,3 \
  --orbit D --angles 40 --description-contains SMAP_L1C \
  --window-days 75 --ndata-min 20 \
  --obsfcstana-format nc4 \
  --prefix PYTHON_TEST_SMAP_Tb_
```

(`--obsfcstana-format` switched from the example's `bin` to `nc4` — this
live archive stores `.nc4`, not `.bin`; confirmed by listing the archive
before running, not assumed.)

**Result:** 8 files written (4 DOY + 4 pentad), 112,573/112,573
administering tiles selected, no non-administering-tile error.
**Runtime: 39.4 s wall / 1.40 GB peak RSS** (`/usr/bin/time -v`).

## 4. MATLAB reference run (unmodified)

`matlab_reference/get_model_and_obs_clim_stats.m` and
`Run_get_L4_Tb_scale_SMAP.m` were run **byte-for-byte unmodified**, invoked
from a temporary wrapper
([`tb_discover_validation_scripts/matlab/run_matlab_test_D.m`](tb_discover_validation_scripts/matlab/run_matlab_test_D.m)),
against the M21C `.bin` archive, matching the Python config exactly (same 3
months/year, D orbit, H+V, 40°, 3-hr/`t0_assim=0`, 75-day window,
`Ndata_min=20`, `hscale=0`, `convert_grid='EASEv2_M36'`), output prefix
`MATLAB_TEST_SMAP_Tb_`.

```bash
module load matlab/R2023b
matlab -batch "run_matlab_test_D"   # from tb_discover_validation_scripts/matlab/
```

**Result:** 8 files written. **Runtime: 73.6 s wall (internal `toc`: 62.1 s)
/ 1.73 GB peak RSS.**

### Why M21C (`.bin`, oldest schema) instead of OLv7 for the MATLAB side

`matlab_reference/get_model_and_obs_clim_stats.m` hardcodes a `.bin` suffix
when building the ObsFcstAna filename (`clsm_ensupd_read_obs.F90`-style
Fortran-sequential reader) — it cannot read `.nc4`. Every live experiment
recent enough to carry the *current* obsparam schema (`fcstvarname`/
`fcstunits`, which `obs_scaling/io.py:read_obs_param` expects) also turned
out to write `.nc4` ObsFcstAna. M21C is the one real experiment on Discover
with `.bin` ObsFcstAna, but its obsparam predates *both* `fcstvarname`/
`fcstunits` **and** the `flistpath`/`flistname` fields the checked-in
`matlab2python/shared/matlab/read_obsparam.m` itself expects.

Empirically (see
[`check_m21c_schema.m`](tb_discover_validation_scripts/matlab/check_m21c_schema.m)):
the **unmodified** `read_obsparam.m` parses M21C's `SMAP_L1C_Tb*` entries
(species 5-8, i.e. exactly what this test needs) correctly — the schema gap
only manifests starting at the first ASCAT entry (species 9+), which this
test doesn't touch. So the MATLAB side of the comparison used the
production reader unmodified.

Python's `read_obs_param` cannot parse M21C's obsparam at all (crashes
immediately on the missing `fcstvarname`/`fcstunits` fields), and per the
validation ground rules the implementation was not to be modified. Since
`generate_tb_scaling_params()` accepts a pre-built `obs_params` list and
does not call `read_obs_param` internally, the Python side of the M21C
comparison used
[`run_python_m21c_test.py`](tb_discover_validation_scripts/python/run_python_m21c_test.py),
which constructs the four `SMAP_L1C_Tb*` `ObsParam` records directly from
values independently cross-checked against the raw obsparam text and
against MATLAB's own parse of the same file. Every other function under
test (`generate_tb_scaling_params`, `build_administering_tiles`,
`select_tb_species`, `read_obs_fcst_ana`, `write_tb_scaling_file`) ran
completely unmodified.

```bash
python -u tb_discover_validation_scripts/python/run_python_m21c_test.py
```
**Runtime: 22.3 s wall / 1.37 GB peak RSS.**

## 5. Pentad-file comparison (descending)

Compared with **two independently-implemented readers** — deliberately not
the same decoder on both sides:

- Python: `obs_scaling.seqbin.read_tb_scaling_file()` (production code)
- MATLAB output: a from-scratch big-endian Fortran-sequential binary reader,
  [`independent_matlab_reader.py`](tb_discover_validation_scripts/python/independent_matlab_reader.py)

```bash
python -u tb_discover_validation_scripts/python/compare_matlab_python.py
```

| Check | Result |
|---|---|
| `asc_flag` | 0 for both, all 4 pentads |
| 2nd header int (legacy version slot) | 0 for both |
| start/end times, angles | exact match |
| tile count | 112,573 = 112,573 |
| tile-ID sets | identical, 0 only-in-either |
| tile ordering | identical (no join needed) |
| duplicate tile IDs | 0 in both |
| lon/lat | max/median/p99 abs diff = **0.000e+00** |
| H/V no-data masks | identical |
| H/V `N_data` | exact integer match |
| `mean_obs`, `std_obs`, `mean_mod`, `std_mod` (×H,V) | max/median/p99 abs diff = **0.000e+00** (tolerance was 5e-5/5e-4) |
| debug fields 11-14 | max/median/p99 abs diff = **0.000e+00** |

All 4 pentad files (`p08`-`p11`) are **byte-identical** (`cmp -s`).

## 6. Administering-tile coverage audit

```bash
python -u tb_discover_validation_scripts/python/audit_tile_coverage.py
```

Across all 728 real `.bin` ObsFcstAna files (Jan-Mar 2016, D orbit, species
6+8):

| Metric | Value |
|---|---:|
| Total model tiles | 112,573 |
| Total M36 administering tiles | 112,573 |
| Unique observed tile numbers | 63,051 |
| Observed tiles absent from scaling file | **0** |
| Duplicate scaling tile IDs | **0** |

## 7. GEOS-LDAS consumption test

Three SLURM iterations on `hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE`
(`ldas_setup` config, `LDASsa_SPECIAL_inputs_ensupd.nml` diff, and exeinp
preserved under
[`tb_discover_validation_scripts/geosldas_smoke_test/`](tb_discover_validation_scripts/geosldas_smoke_test/)):

1. **Job 57618171 — failed**, unrelated to this port:
   `LDAS ERROR (3000) from scale_obs_Tb_zscore: could not open file` on
   species 1 (`SMOS_fit_Tbh_A`), processed before the target species. Cause:
   `hsaf_cdr_test/LDASsa_SPECIAL_inputs_ensupd.nml` (shared production
   namelist) has `scale=.true.` for SMOS (21-24) and ascending-Tb Tb (31,
   33) species pointing `scalepath` at
   `/discover/nobackup/amfox/Experiments/M21C_land_sweeper_OLv8_M36/...`,
   which does not exist. Pre-existing, unrelated to this port — worked
   around in a private nml copy by setting `scale=.false.` for just those
   species (see `nml_diff_vs_production.diff`).
2. **Job 57618394 — completed cleanly** (structural pass: file opened,
   header parsed, `asc_flag`/tile-count as expected) but the scaling file
   read had **zero valid tiles anywhere**. Cause: the window (15 days) was
   sized against `hsaf_cdr_test_DAv8_M36`'s archive, which turned out to
   have only **one day** (2020-01-01, 8 files) of live ObsFcstAna
   data — `Ndata_min=20` was mathematically unreachable at any window size
   from that source. Caught by explicitly checking valid-tile counts before
   declaring success (`check_raw_obs_counts.py` confirmed max observed
   samples per tile = 1).
3. Regenerated a genuinely valid scaling file (75-day window, real
   multi-year `OLv7_M36_MULTI_type_13_H121` archive, `GEOSLDAS_SMOKE_VALID_
   SMAP_Tb_...` prefix, 11 valid tiles for pentads 1-2) and reran.
   **Job 57619200 — completed cleanly**, 22:19 walltime, 48 ObsFcstAna
   files produced.

```bash
source /discover/nobackup/projects/land_da/GEOSldas_develop/GEOSldas/@env/g5_modules.sh
/discover/nobackup/projects/land_da/GEOSldas_develop/GEOSldas/install/bin/ldas_setup setup \
  /discover/nobackup/projects/land_da/hsaf_cdr_test/ \
  /discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE.txt \
  /discover/nobackup/projects/land_da/hsaf_cdr_test/bat_inp_debug.txt
cd /discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE/run
sbatch lenkf.j
```

**Verification against real output** (species 6 `Tbh_D`, 8 `Tbv_D`), via
[`verify_geosldas_consumption.py`](tb_discover_validation_scripts/python/verify_geosldas_consumption.py):

| Check | Result |
|---|---|
| Abort on orbit/angle/tile-ID/coordinate mismatch | none |
| Valid-stat obs → assimilable | 98/98 got `assim_flag=1` |
| Invalid-stat obs → do-not-assimilate | 258,911/258,911 got `assim_flag=0`; 0 invalid-stat obs ever assimilated |
| `obsvar == errstd² × (std_mod/std_obs)²` | 98/98 matched, **0 relative error** |
| Scaling formula (`mean_mod + (raw − mean_obs)·std_mod/std_obs`) | inverse-reconstructed raw Tb landed in [240, 263] K for all 98 — physically plausible |

## 8. Ascending orbit repeat

Same M21C archive, species 5/7, `MATLAB_TEST_SMAP_Tb_ASC_` /
`PYTHON_TEST_M21C_SMAP_Tb_ASC_` prefixes
([`run_python_m21c_test_A.py`](tb_discover_validation_scripts/python/run_python_m21c_test_A.py),
[`run_matlab_test_A.m`](tb_discover_validation_scripts/matlab/run_matlab_test_A.m),
[`compare_matlab_python_asc.py`](tb_discover_validation_scripts/python/compare_matlab_python_asc.py)).

**Result:** `asc_flag=1` confirmed for both tools, all 4 pentads, tile
sets/ordering identical, **byte-identical files**, 0 differences on every
field. Python: 38.5 s / 1.40 GB. MATLAB: 82.9 s wall / 1.74 GB.

---

## First point of divergence

**None found.** Every byte-level, field-level, and downstream-consumption
check passed exactly.

## Findings worth fixing separately (not blockers, not bugs in this port)

1. **`matlab2python/shared/matlab/read_obsparam.m` is stale.** It predates
   the `fcstvarname`/`fcstunits` fields current obsparam files carry (hard
   crash, `fscanf: Invalid size`, on the 2nd entry of any current-format
   file — see
   [`schema_bug_diagnosis.m`](tb_discover_validation_scripts/matlab/schema_bug_diagnosis.m)).
   Independently, it is not quote-aware like Python's `shlex`-based reader,
   so any quoted field with an embedded space (e.g. ASCAT's
   `fcstunits='m3 m-3'`) desyncs `fscanf('%s')` and corrupts every
   subsequent field in the file. A schema-aware, quote-aware replacement
   was prototyped for this validation
   ([`read_obsparam_newschema.m`](tb_discover_validation_scripts/matlab/read_obsparam_newschema.m)
   + [`tokenize_obsparam_file.m`](tb_discover_validation_scripts/matlab/tokenize_obsparam_file.m))
   but not applied to §4-6/§8 (the unmodified reader happened to work there
   — see §4's "why M21C" note) or committed to the repo.
2. **`hsaf_cdr_test/LDASsa_SPECIAL_inputs_ensupd.nml`** (shared production
   namelist) has dangling `scalepath` entries for SMOS/ascending-Tb species
   pointing at a path that doesn't exist — aborts any GEOSldas run using it
   unmodified, on the very first assimilation cycle.
3. **`hsaf_cdr_test_DAv8_M36`** (the experiment referenced in the top-level
   `CLAUDE.md`) only has 1 day of live ObsFcstAna archived
   (`ana/ens_avg/Y2020/M01/`, 8 files), despite the month directory
   structure suggesting more. Not usable as a scaling-climatology source at
   any `Ndata_min` above ~1.

## Artifacts left on Discover (not yet cleaned up)

- `/discover/nobackup/projects/land_da/obs_scaling_test/` — M21C `.bin`
  extraction + Python/MATLAB test output (§4-6, §8). ~2 GB.
- `/discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE/`
  — the GEOS-LDAS smoke-test experiment directory (§7), including restart
  links, build, and 48 output ObsFcstAna files.
- `/discover/nobackup/projects/land_da/hsaf_cdr_test/nml_GEOSLDAS_SMOKE_TEST/`
  — the patched namelist used to work around Finding 2.
- `/discover/nobackup/projects/land_da/hsaf_cdr_test/hsaf_cdr_test_DAv8_M36_GEOSLDAS_SMOKE.txt`
  — the exeinp file for the above.
- `GEOSLDAS_SMOKE_VALID_SMAP_Tb_*`, `GEOSLDAS_SMOKE_SMAP_Tb_*`,
  `GEOSLDAS_SMOKE2_SMAP_Tb_*`, `PYTHON_TEST_*`, `MATLAB_TEST_*` scaling
  files under the various experiments' `stats/z_score_clim/` directories.

None of these were required for or affect the production run this
validation gates. Deletion left to a separate, explicit cleanup pass.

## Reproducing this validation

The scripts under `tb_discover_validation_scripts/` hardcode the specific
Discover paths used above (module versions, experiment paths, scratch
output locations). They are preserved as an exact record of what was run,
not as a maintained/parameterized test suite — expect to edit paths before
re-running against a different experiment or date range.
