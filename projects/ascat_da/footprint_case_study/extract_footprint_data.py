#!/usr/bin/env python3

"""
Stage 1 (Discover-only): build the entire self-contained data bundle for the
peat-QC footprint case study from raw sources -- no external/undocumented
CSV inputs. Everything downstream (stage 2, plot_footprint_case_study.py)
reads only the CSVs this script writes into the bundle's data directory.

This folds together, into one auditable pipeline, what were previously three
separate scripts in the parent project directory:
  - compare_peatqc_boundary_tiles_latlon.py  (raw OFA extraction + kept/lost)
  - check_ndata_closure_boundary_tiles.py    (N_data closure validation)
  - the original extract_footprint_data.py   (tile geometry + BCS porosity)

Experiments compared:
  baseline: DAv7_M36_ASCAT_type_13_H121_FOV12p5
  peat-QC:  DAv7_M36_ASCAT_type_13_H121_FOV12p5_peatlandqc
Period: 2015-04-01 through 2015-09-30 (6 months, 3-hourly OFA files).
Species: 8, 9 = ASCAT_HSAF_META_SM, ASCAT_HSAF_METB_SM (Metop-A/B).
Target tiles: 3823, 4338, 6548, 9449, 4813, 5230 (peat-boundary examples).

--------------------------------------------------------------------------
WHAT "kept" MEANS (the exact predicate, and why)
--------------------------------------------------------------------------
An OFA row for a given (timestamp, species, tilenum) key is called VALID if
it exists in the file AND its `obs` value is finite (not the netCDF fill
value). A baseline candidate is "kept" if a VALID row exists for the same
key in the peat-QC output; otherwise it is "lost".

This predicate was chosen empirically, not assumed: four candidate
predicates were tested --
  (a) row exists at all
  (b) row exists and obs is finite               <- the one used
  (c) row exists and obs AND fcst are finite
  (d) row exists and assim_flag == 1
-- against the production N_data counts (from the cached monthly
ObsFcstAna_sums, the same postproc pipeline used elsewhere in this project).
For these two ASCAT HSAF species, on these 6 tiles, **all four predicates
are numerically identical to N_data and to each other** (validated below,
hard-asserted, not just noted): an OFA row for an assimilated species is
apparently never written except when finite and assim_flag==1, so there is
no ambiguity to resolve here. If this assertion ever fails for a different
species/tile set, that would mean the reader can emit a "written but
invalid" row and the predicate would need to be revisited.

--------------------------------------------------------------------------
FOOTPRINT PEAT-FRACTION ALGORITHM (reconstructed for stage 2)
--------------------------------------------------------------------------
Reproduces clsm_ensupd_upd_routines.F90 on branch
feature/amfox/ascat-peatland-qc @ 95f0317 (the commit that built the
GEOSldas.x used for DAv7_M36_ASCAT_type_13_H121_FOV12p5_peatlandqc,
job 57695424):
  1. get_tile_num_in_ellipse() [LDAS_TileCoordRoutines.F90 line 1053]: for
     ASCAT HSAF species (FOV=12.5 km, FOV_units='km'), search radius =
     fac_search_FOV_km(=2.0) * FOV, converted km->deg via dist_km2deg()
     [clsm_ensupd_upd_routines.F90 line 6206] using MAPL_RADIUS=6371.0 km:
       dist_y_deg = dist_km * (180/pi) / 6371.0
       dist_x_deg = dist_y_deg / cos(lat_rad)
     Tiles (by center-of-mass com_lon/com_lat) within the ellipse
     ((dlon/r_x)^2 + (dlat/r_y)^2 <= 1) are candidates.
  2. Back in the caller [line ~1779]: the raw ellipse-normalized distance
     (relative to the *search* radius r_x,r_y) is rescaled by
     fac_search_FOV_km^2, which algebraically reduces to normalizing
     distance by the FOV itself (not the 2x-expanded search radius):
     rescaled_ndst2 = (dlon/FOV_x_deg)^2 + (dlat/FOV_y_deg)^2.
     Gaussian weight = exp(-0.5 * rescaled_ndst2); combined weight =
     Gaussian weight * tile area.
  3. Peat indicator per tile: 1 if POROS >= PEATCLSM_POROS_THRESHOLD(=0.90),
     else 0 [line 1557-1563].
  4. peatfrac = weighted-mean(peat indicator) over tiles in the ellipse
     [line 1789-1817]. Reject (tmpPeat=True) if peatfrac >= 0.10, or if no
     tile in the ellipse has a valid soil class.

POROS comes from the BCS catchment parameter file for this experiment's
BCS_PATH (from the exeinp: NL5, EASEv2_M36), *not* a generic peat map:
  /discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/NL5/land/EASEv2_M36/clsm/catch_params.nc4
Verified: N_tile there (112573) == tile_coord.bin N_tile, and tile_id is
sequential 1..N, so POROS[tilenum-1] is the tile's porosity directly.
"""

import glob
import os
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset

WORKTREE = '/gpfsm/dnb06/projects/p284/hsaf_cdr_test/.worktrees/obsfcstana-nc4-postproc/GEOSldas_App/util'
sys.path.insert(0, WORKTREE + '/shared/python')
from read_GEOSldas import read_tilecoord, read_obs_param

EXPDIR = '/discover/nobackup/projects/land_da/hsaf_cdr_test/'
DOMAIN = 'SMAP_EASEv2_M36_GLOBAL'
BASELINE_ID = 'DAv7_M36_ASCAT_type_13_H121_FOV12p5'
PEATQC_ID = 'DAv7_M36_ASCAT_type_13_H121_FOV12p5_peatlandqc'
EXPTAG = {BASELINE_ID: 'DA_baseline_FOV12p5', PEATQC_ID: 'DA_peatlandqc'}
SUMS_CACHE = '/gpfsm/dnb06/projects/p284/hsaf_cdr_test/omf_compare_sums'
BCS_POROS_FILE = '/discover/nobackup/projects/gmao/bcs_shared/fvInput/ExtData/esm/tiles/NL5/land/EASEv2_M36/clsm/catch_params.nc4'

SPECIES = (8, 9)
TARGET_TILES = [3823, 4338, 6548, 9449, 4813, 5230]
MONTHS = [f'2015{m:02d}' for m in range(4, 10)]  # 201504..201509
OBSPARAM_TIME = '20150401_0000'

BUNDLE_DIR = '/gpfsm/dnb06/projects/p284/hsaf_cdr_test/footprint_case_study'
DATA_DIR = f'{BUNDLE_DIR}/data'
os.makedirs(DATA_DIR, exist_ok=True)

# generous box around each target tile's com_lon/com_lat for the tile-geometry
# extract: comfortably covers the ~25 km (~0.2-0.5 deg, lat-dependent) search
# radius plus any offset between the obs footprint centre and tile centre.
BOX_LON_DEG = 1.5
BOX_LAT_DEG = 0.8

REPORT_LINES = []


def log(msg):
    print(msg, flush=True)
    REPORT_LINES.append(msg)


# ----------------------------------------------------------------------------------
# 1. raw OFA extraction, both experiments, with an explicit duplicate-key check
# ----------------------------------------------------------------------------------

def ofa_files(expid, yyyymm):
    yyyy, mm = yyyymm[:4], yyyymm[4:6]
    pattern = f'{EXPDIR}{expid}/output/{DOMAIN}/ana/ens_avg/Y{yyyy}/M{mm}/{expid}.ens_avg.ldas_ObsFcstAna.*z.nc4'
    return sorted(glob.glob(pattern))


def extract_experiment_rows(expid, target_tiles, species, months):
    """Return dict[(timestamp, species, tilenum)] -> row dict, for every row
    matching target_tiles/species across all files. Raises if the same key
    is ever seen twice (duplicate-key assertion)."""
    target_tiles_set = set(target_tiles)
    rows = {}
    n_files = 0
    n_seen = 0
    duplicates = []

    for yyyymm in months:
        for fpath in ofa_files(expid, yyyymm):
            n_files += 1
            fname = os.path.basename(fpath)
            ts = fname.split('.ldas_ObsFcstAna.')[1].split('z.nc4')[0]
            with Dataset(fpath) as d:
                tilenum = d.variables['tilenum'][:]
                sp = d.variables['species'][:]
                mask = np.isin(sp, species) & np.isin(tilenum, list(target_tiles_set))
                if not np.any(mask):
                    continue
                idx = np.nonzero(mask)[0]
                lon = d.variables['lon'][idx]
                lat = d.variables['lat'][idx]
                obs = d.variables['obs'][idx]
                fcst = d.variables['fcst'][idx]
                assim_flag = d.variables['assim_flag'][idx]
                tn = tilenum[idx]
                sp = sp[idx]

            for j in range(len(idx)):
                key = (ts, int(sp[j]), int(tn[j]))
                n_seen += 1
                if key in rows:
                    duplicates.append(key)
                    continue
                obs_v = obs[j]
                fcst_v = fcst[j]
                obs_finite = (not np.ma.is_masked(obs_v)) and np.isfinite(float(obs_v))
                fcst_finite = (not np.ma.is_masked(fcst_v)) and np.isfinite(float(fcst_v))
                rows[key] = {
                    'timestamp': ts, 'species': int(sp[j]), 'tilenum': int(tn[j]),
                    'lon': float(lon[j]), 'lat': float(lat[j]),
                    'obs': float(obs_v) if not np.ma.is_masked(obs_v) else np.nan,
                    'fcst': float(fcst_v) if not np.ma.is_masked(fcst_v) else np.nan,
                    'assim_flag': int(assim_flag[j]),
                    'a_exists': True, 'b_finite_obs': obs_finite,
                    'c_finite_obs_fcst': obs_finite and fcst_finite,
                    'd_assim_flag': int(assim_flag[j]) == 1,
                }

    assert not duplicates, (
        f'{expid}: {len(duplicates)} duplicate (timestamp,species,tilenum) keys found '
        f'-- comparison key is not unique, cannot proceed. First few: {duplicates[:5]}')

    log(f'{expid}: opened {n_files} files, {n_seen} rows matched target tiles/species, '
        f'0 duplicate keys (key uniqueness assertion PASSED)')
    return rows


log('=== Stage 1a: raw OFA extraction + key-uniqueness assertion ===')
base_rows = extract_experiment_rows(BASELINE_ID, TARGET_TILES, SPECIES, MONTHS)
peat_rows = extract_experiment_rows(PEATQC_ID, TARGET_TILES, SPECIES, MONTHS)

# ----------------------------------------------------------------------------------
# 1b. N_data closure assertion: do our four candidate predicates reproduce the
#     production N_data counts exactly, for every (tile, species) pair, both
#     experiments? This is what pins down the "kept" predicate above.
# ----------------------------------------------------------------------------------


def species_index(obsparam, species_id):
    for i, p in enumerate(obsparam):
        if int(p['species']) == species_id:
            return i
    raise ValueError(f'species {species_id} not found in obsparam')


def ndata_counts(expid, target_tiles, species_list, months):
    exptag = EXPTAG[expid]
    fop = f"{EXPDIR}{expid}/output/{DOMAIN}/rc_out/Y2015/M04/{expid}.ldas_obsparam.{OBSPARAM_TIME}z.txt"
    obsparam = read_obs_param(fop)
    idxs = {s: species_index(obsparam, s) for s in species_list}

    tc_local = read_tilecoord(f"{EXPDIR}{expid}/output/{DOMAIN}/rc_out/{expid}.ldas_tilecoord.bin")
    tile_pos = {int(t): i for i, t in enumerate(tc_local['tile_id'])}

    totals = {(t, s): 0 for t in target_tiles for s in species_list}
    for yyyymm in months:
        yyyy, mm = yyyymm[:4], yyyymm[4:6]
        fout = f"{SUMS_CACHE}/{exptag}/Y{yyyy}/M{mm}/{exptag}.ens_avg.ldas_ObsFcstAna_sums.{yyyymm}.nc4"
        with Dataset(fout) as nc:
            N_data = nc.variables['N_data'][:]
        for t in target_tiles:
            pos = tile_pos.get(t)
            if pos is None:
                continue
            for s in species_list:
                totals[(t, s)] += int(N_data[pos, idxs[s]])
    return totals


log('\n=== Stage 1b: N_data closure assertion ===')
predicates = ['a_exists', 'b_finite_obs', 'c_finite_obs_fcst', 'd_assim_flag']

closure_failures = []
for label, rows, expid in [('BASELINE', base_rows, BASELINE_ID), ('PEAT-QC', peat_rows, PEATQC_ID)]:
    nd = ndata_counts(expid, TARGET_TILES, SPECIES, MONTHS)
    raw_counts = {(t, s): {p: 0 for p in predicates} for t in TARGET_TILES for s in SPECIES}
    for (ts, sp, tn), row in rows.items():
        c = raw_counts[(tn, sp)]
        for p in predicates:
            if row[p]:
                c[p] += 1

    for t in TARGET_TILES:
        for s in SPECIES:
            n = nd[(t, s)]
            for p in predicates:
                v = raw_counts[(t, s)][p]
                if v != n:
                    closure_failures.append((label, t, s, p, v, n))

    match_str = ', '.join(f'{p}: {sum(1 for t in TARGET_TILES for s in SPECIES if raw_counts[(t,s)][p]==nd[(t,s)])}/{len(TARGET_TILES)*len(SPECIES)}' for p in predicates)
    log(f'{label}: N_data exact-match counts by predicate -- {match_str}')

assert not closure_failures, (
    f'N_data closure FAILED for {len(closure_failures)} (label,tile,species,predicate) combos: '
    f'{closure_failures[:10]}')
log('N_data closure assertion PASSED: all 4 predicates reproduce N_data exactly '
    f'for all {len(TARGET_TILES)} tiles x {len(SPECIES)} species x 2 experiments '
    '-- confirms "row exists with finite obs" == "finite obs+fcst" == "assim_flag==1" '
    'for these ASCAT HSAF species (no ambiguity in the kept/lost predicate).')

# ----------------------------------------------------------------------------------
# 1c. classify every valid baseline candidate as kept/lost using predicate (b)
# ----------------------------------------------------------------------------------

records = []
for key, brow in base_rows.items():
    if not brow['b_finite_obs']:
        continue
    prow = peat_rows.get(key)
    kept = (prow is not None) and prow['b_finite_obs']
    records.append({
        'timestamp': brow['timestamp'], 'species': brow['species'], 'tilenum': brow['tilenum'],
        'lon': brow['lon'], 'lat': brow['lat'], 'obs': brow['obs'], 'fcst': brow['fcst'],
        'assim_flag': brow['assim_flag'], 'kept': kept,
    })

df_obs = pd.DataFrame(records)
n_kept = int(df_obs['kept'].sum())
log(f'\nbaseline valid candidates: {len(df_obs)}  kept in peat-QC: {n_kept} '
    f'({100*n_kept/len(df_obs):.1f}%)  lost: {len(df_obs)-n_kept}')

obs_out = f'{DATA_DIR}/obs_candidates.csv'
df_obs.to_csv(obs_out, index=False)
log(f'wrote {obs_out}  ({len(df_obs)} rows)')

# ----------------------------------------------------------------------------------
# 2. tile geometry (com_lon, com_lat, area) + BCS porosity, neighborhoods around
#    the 6 target tiles
# ----------------------------------------------------------------------------------

log('\n=== Stage 2: tile geometry + BCS porosity extraction ===')

tc = read_tilecoord(
    f'{EXPDIR}{PEATQC_ID}/output/{DOMAIN}/rc_out/{PEATQC_ID}.ldas_tilecoord.bin')

tile_id = tc['tile_id']       # 1-based, sequential 1..N_tile
com_lon = tc['com_lon']
com_lat = tc['com_lat']
area = tc['area']

with Dataset(BCS_POROS_FILE) as d:
    poros = d.variables['POROS'][:]

assert len(poros) == tc['N_tile'], \
    f"BCS POROS length {len(poros)} != tile_coord N_tile {tc['N_tile']}"
assert np.array_equal(tile_id, np.arange(1, tc['N_tile'] + 1)), \
    'tile_id is not sequential 1..N -- POROS[tilenum-1] indexing assumption broken'
log('BCS POROS / tile_coord alignment assertions PASSED '
    f'(N_tile={tc["N_tile"]}, tile_id sequential 1..N)')

target_pos = {t: t - 1 for t in TARGET_TILES}

rows = []
for t in TARGET_TILES:
    clon = com_lon[target_pos[t]]
    clat = com_lat[target_pos[t]]
    in_box = (np.abs(com_lon - clon) <= BOX_LON_DEG) & (np.abs(com_lat - clat) <= BOX_LAT_DEG)
    idx = np.nonzero(in_box)[0]
    for i in idx:
        rows.append({
            'target_tile': t,
            'tile_id': int(tile_id[i]),
            'com_lon': float(com_lon[i]),
            'com_lat': float(com_lat[i]),
            'area': float(area[i]),
            'poros': float(poros[i]),
        })

df_tiles = pd.DataFrame(rows)
tiles_out = f'{DATA_DIR}/local_tile_geometry_porosity.csv'
df_tiles.to_csv(tiles_out, index=False)
log(f'wrote {tiles_out}  ({len(df_tiles)} rows, neighborhoods for {len(TARGET_TILES)} target tiles)')

# ----------------------------------------------------------------------------------
# 3. constants used by stage 2
# ----------------------------------------------------------------------------------

const_out = f'{DATA_DIR}/constants.csv'
pd.DataFrame([
    {'name': 'FOV_km', 'value': 12.5,
     'source': 'obsparam species 8/9 (ASCAT_HSAF_META/METB_SM), FOV field, FOV_units=km'},
    {'name': 'fac_search_FOV_km', 'value': 2.0,
     'source': 'clsm_ensupd_upd_routines.F90 line 1067'},
    {'name': 'PEATCLSM_POROS_THRESHOLD', 'value': 0.90,
     'source': 'catch_constants.f90 line 152'},
    {'name': 'ASCAT_max_peat_frac', 'value': 0.10,
     'source': 'clsm_ensupd_upd_routines.F90 line 1070'},
    {'name': 'MAPL_RADIUS_km', 'value': 6371.0,
     'source': 'PhysicalConstants.F90 line 29 (MAPL_RADIUS=6371.0E3 m)'},
]).to_csv(const_out, index=False)
log(f'wrote {const_out}')

# ----------------------------------------------------------------------------------
# 4. validation report
# ----------------------------------------------------------------------------------

report_out = f'{BUNDLE_DIR}/validation_report.txt'
with open(report_out, 'w') as f:
    f.write('\n'.join(REPORT_LINES) + '\n')
log(f'\nwrote {report_out}')

log('\nDONE -- footprint_case_study/ is a self-contained, portable data bundle, '
    'built from raw OFA files with no external/undocumented CSV inputs.')
log('scp this whole directory to run stage 2 (plot_footprint_case_study.py) '
    'anywhere (only needs numpy/pandas/matplotlib).')
