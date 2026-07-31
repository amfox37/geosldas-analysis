import sys
sys.path.insert(0, "/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params")
sys.path.insert(0, "/gpfsm/dnb34/tdirs/login/discover34.252777.amfox/claude-236601614/-gpfsm-dnb06-projects-p284-geosldas-analysis/bad106a1-1ee2-4c2d-96cc-d1c5e084a6fd/scratchpad")
import numpy as np
from obs_scaling.seqbin import read_tb_scaling_file
from independent_matlab_reader import read_matlab_scaling_file

OUT = "/discover/nobackup/projects/land_da/obs_scaling_test/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/stats/z_score_clim"
pentads = ["p08", "p09", "p10", "p11"]
py_base = f"{OUT}/PYTHON_TEST_M21C_SMAP_Tb_2016_p1_2016_p18_hscale_0.00_W_15p_Nmin_20_D_"
ml_base = f"{OUT}/MATLAB_TEST_SMAP_Tb_2016_p1_2016_p18_hscale_0.00_W_15p_Nmin_20_D_"

field_names = [
    "mean_obs_H", "std_obs_H", "mean_mod_H", "std_mod_H", "N_data_H",
    "mean_obs_V", "std_obs_V", "mean_mod_V", "std_mod_V", "N_data_V",
    "H_obs(dbg11)", "H_mod(dbg12)", "V_obs(dbg13)", "V_mod(dbg14)",
]
NODATA = -9999.0

overall_ok = True

for p in pentads:
    py_path = f"{py_base}{p}.bin"
    ml_path = f"{ml_base}{p}.bin"
    print(f"\n=== Pentad {p} ===")

    py = read_tb_scaling_file(py_path)
    ml = read_matlab_scaling_file(ml_path)

    print(f"Python asc_flag={py.asc_flag}  MATLAB asc_flag={ml.asc_flag} (expect 0, descending)")
    print(f"MATLAB header int2 (version)={ml.version} (expect 0)")
    print(f"Python angles={py.angles}  MATLAB angles={ml.angles}")
    print(f"Python start={py.start_time}  end={py.end_time}")
    print(f"MATLAB  start={ml.start_time}  end={ml.end_time}")
    print(f"Python n_tiles={py.tile_id.size}  MATLAB n_grid={ml.n_grid}")

    if py.asc_flag != 0 or ml.asc_flag != 0:
        print("  !! asc_flag mismatch from expected 0")
        overall_ok = False
    if ml.version != 0:
        print("  !! MATLAB header int2 != 0")
        overall_ok = False

    py_ids = py.tile_id.astype(np.int64)
    ml_ids = ml.tile_id.astype(np.int64)
    py_dup = py_ids.size != np.unique(py_ids).size
    ml_dup = ml_ids.size != np.unique(ml_ids).size
    print(f"Python duplicate tile IDs: {py_dup}   MATLAB duplicate tile IDs: {ml_dup}")

    set_py = set(py_ids.tolist())
    set_ml = set(ml_ids.tolist())
    only_py = set_py - set_ml
    only_ml = set_ml - set_py
    print(f"tile ID sets identical: {len(only_py) == 0 and len(only_ml) == 0}  "
          f"(only_py={len(only_py)}, only_ml={len(only_ml)})")
    if only_py or only_ml:
        overall_ok = False

    order_matches = py_ids.size == ml_ids.size and np.array_equal(py_ids, ml_ids)
    print(f"tile ordering identical: {order_matches}")

    # join by tile_id
    py_pos = {int(t): i for i, t in enumerate(py_ids)}
    ml_pos = {int(t): i for i, t in enumerate(ml_ids)}
    common = sorted(set_py & set_ml)
    py_idx = np.array([py_pos[t] for t in common])
    ml_idx = np.array([ml_pos[t] for t in common])

    py_lon = py.lon[py_idx]; ml_lon = ml.lon[ml_idx]
    py_lat = py.lat[py_idx]; ml_lat = ml.lat[ml_idx]
    dlon = np.abs(py_lon - ml_lon)
    dlat = np.abs(py_lat - ml_lat)
    print(f"lon diff: max={dlon.max():.3e} median={np.median(dlon):.3e} p99={np.percentile(dlon,99):.3e}")
    print(f"lat diff: max={dlat.max():.3e} median={np.median(dlat):.3e} p99={np.percentile(dlat,99):.3e}")

    # data shape: py.data (14, n_tiles, n_angles); ml.data (14, n_grid, n_angle)
    py_data = py.data[:, py_idx, 0]
    ml_data = ml.data[:, ml_idx, 0]

    for fi, fname in enumerate(field_names):
        pv = py_data[fi]
        mv = ml_data[fi]
        py_nodata = np.isclose(pv, NODATA, atol=1e-2)
        ml_nodata = np.isclose(mv, NODATA, atol=1e-2)
        mask_match = np.array_equal(py_nodata, ml_nodata)
        both_valid = (~py_nodata) & (~ml_nodata)
        if "N_data" in fname:
            same = np.array_equal(pv[both_valid], mv[both_valid])
            print(f"  {fname:16s} nodata-mask match={mask_match!s:5} exact-match(valid)={same!s:5} n_valid={both_valid.sum()}")
            if not same or not mask_match:
                overall_ok = False
        else:
            diff = np.abs(pv[both_valid] - mv[both_valid])
            if diff.size:
                dmax, dmed, dp99 = diff.max(), np.median(diff), np.percentile(diff, 99)
            else:
                dmax = dmed = dp99 = float("nan")
            tol = 5e-5 if "mean" in fname or "obs" in fname or "mod" in fname else 5e-4
            print(f"  {fname:16s} nodata-mask match={mask_match!s:5} n_valid={both_valid.sum():6d} "
                  f"max={dmax:.3e} median={dmed:.3e} p99={dp99:.3e}")
            if not mask_match:
                overall_ok = False

print("\n=== SUMMARY ===")
print("ALL CHECKS PASSED" if overall_ok else "SOME CHECKS FAILED (see !! markers above)")
