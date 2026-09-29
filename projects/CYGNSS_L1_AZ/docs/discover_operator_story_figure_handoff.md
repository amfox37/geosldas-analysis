# Discover task: stage inputs for a fixed-operator "operator story" figure

## Goal

The paper's methods schematic is a 3×3 single-observation figure made by
`CYGNSS_operator/scripts/plot_cygnss_good_obs_operator_story.py`:

- L1 DDM with the 3×5 crop
- DEM
- footprint
- vegetation opacity
- forecast SFMC
- C_t, R_t
- the SRFEXC increment and SFMC response for that assimilation cycle

The only existing version is from the 2019-11-01 two-day test, which was run
with the **pre-fix** GEOSldas L1 operator, so its increment panels cannot be
used. The fix (GEOSldas_GridComp `cdfaaf4` + `11dfdb1`) went in on 2026-09-25.

The script will be **run on the laptop**, which already has IGOT, GDAL and the
Copernicus DEM. Your job on Discover is to:

1. choose a good observation and cycle from the fixed-operator runs
2. stage the few files the script needs into one small bundle to copy back

Background: `projects/CYGNSS_L1_AZ/docs/cygl1_operator_test_project_README.md`
(experiment table, data layout) and `docs/paper_source_manifest.md`.

## Constraints

- **Read-only on existing experiment and obs directories.** Write only into
  the bundle directory below.
- **No `sbatch` / `scancel`.** If a rerun is needed (see step 1), stop and
  report what would be needed.
- **No git commits or pushes.**
- The bundle should be tens of MB. If it grows past about 500 MB, stop and
  report.

Bundle directory: `/discover/nobackup/amfox/cygl1_story_figure_bundle/`

Experiment root: `/discover/nobackup/projects/land_da/cygl1_operator_test/`

## Arms

The primary L1 arm is not decided yet, so stage **the same cycle from both
arms**:

- `DA_L1_full_xc015_coh040216_err39_fixedop` (best by O-F, errstd 3.9)
- `DA_L1_full_xc015_coh040216_fixedop` (best in situ, errstd 2.75)

Both arms assimilate the screened obs in `CYGNSS_L1_coh040_216/`. Both run
2020–2022.

## Steps

### 1. Check that the needed output collections exist

The script needs, per cycle:

- `*.inst3_1d_lndfcstana_Nt.<YYYYMMDD_HHMM>z.nc4` (forecast/analysis SFMC)
- `*.catch_progn_incr.<YYYYMMDD_HHMM>z.nc4` (SRFEXC increment)
- `*.ldas_mwRTMparam.*.nc4`, or whatever holds the mwRTM parameters for the domain

Look under `<EXP_ID>/output/<domain>/cat/ens_avg/Y<YYYY>/M<MM>/` and `rc_out/`.
If either arm did not write `inst3_1d_lndfcstana_Nt` or `catch_progn_incr`, stop here. Report:

- which collections are missing
- the HISTORY / exeinp settings that would turn them on
- the smallest rerun that would produce them (one day, restarted from that arm's restart)

### 2. Pick the observation

CYGNSS L1 is species **56** in the ObsFcstAna files (`ana/ens_avg/`). Use
tilecoord for tile lat/lon. OFA's own lat/lon is the superob centre and jitters
from cycle to cycle.

Avoid the spring degradation windows:

- May–Jun 2020
- Apr–Jun 2021
- Feb–May 2022

Good candidates are Jul–Oct 2020 (monsoon-season wetting) or Oct–Dec 2020.

Pick one cycle where the err39 arm has a species-56 obs meeting all of these:

- `assim_flag == 1`
- a clearly nonzero SFMC increment around the obs (top few percent of |incr| for that month)
- |innov| not an outlier (within about 2 × stdv of that month's O-F)
- the specular point is inland in Arizona, not at the domain edge
- the support tiles all fall inside the domain

List 3–5 candidates and choose one. For each, give:

- date and cycle
- tile ID, lat and lon
- obs, fcst, innov, incr
- incr for the same tile in the errstd 2.75 arm

### 3. Identify the obs in the coefficient product

Open the day's product in `CYGNSS_L1_coh040_216/`, the `cygnss_l1_..._m36_<date>_cygNN.nc4` files. Find the obs that matches the chosen OFA obs by:

- time within the 3-h window
- specular point closest to the OFA tile
- observed y

Record for that obs:

- spacecraft (`cygNN`)
- `sample_id`
- `ch_id`
- specular lon/lat
- incidence angle
- observed y (dB)
- support tile count

Also record these global attributes of the product:

- `vegopacity_file`
- `tile_index0_role`
- any attribute naming the tilecoord or experiment that `tile_index0` refers to

### 4. Find the raw L1 DDM file

The DDM panel is built by IGOT from the raw CYGNSS L1 v3.2 day file for that
spacecraft. Find it on Discover (year/DOY/`cygNN`) and stage it. The laptop has
raw L1 only for 2019.

### 5. Stage the bundle

Copy into the bundle directory. Do not symlink; the bundle will be copied off.

- For both arms, only for the chosen cycle: `inst3_1d_lndfcstana_Nt` and `catch_progn_incr`
- the mwRTM param file and tilecoord for the domain (one copy is enough if both arms share them)
- the day's coefficient product file for the spacecraft, from `CYGNSS_L1_coh040_216/`
- the raw L1 v3.2 day file
- the vegopacity file named by the product attribute, **only** if it is not the global M36 `vegopacity.bin`
- the product's `tile_index0` reference file, **only** if the attributes point to something other than the arm's own tilecoord

Then write `MANIFEST.txt` in the bundle. It should list:

- each file: source path, size, md5
- all the step 2–3 values
- the exact laptop command, filled in with bundle-relative paths:

```
python scripts/plot_cygnss_good_obs_operator_story.py \
  --sample-id <S> --ch-id <C> \
  --config <run_config.ini, see note> \
  --product <bundle>/<product>.nc4 \
  --inst3 <bundle>/<arm>/<inst3>.nc4 \
  --increment-file <bundle>/<arm>/<incr>.nc4 \
  --mwrtm-param <bundle>/<mwrtm>.nc4 \
  [--tile-index-reference <file>] \
  --bounds <W> <E> <S> <N>   # ~2° box centred on the specular point
```

Note on `--config`: this is an IGOT run-config `.ini` with year, day,
`sc_num`, `ch_id`, `sample_id` and the DEM options. If the Discover
preprocessor wrote one for this day and spacecraft, stage it. If not, say so;
one will be made on the laptop from the 2019 template.

## Report back (short)

- whether the collections existed (step 1), and if not, the rerun that would be needed
- the chosen obs and its runner-up candidates (the step 2 table)
- the bundle path, total size and `MANIFEST.txt` contents
- anything surprising, e.g.:
  - the product's `tile_index0` reference differs from the arm's tilecoord
  - the two arms disagree in sign on the increment
