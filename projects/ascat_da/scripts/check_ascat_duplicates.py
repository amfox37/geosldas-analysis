#!/usr/bin/env python3
"""
Scan consecutive ObsFcstAna binary files and report identical ASCAT_METB_SM observations
(as defined by matching tile number + lon/lat (rounded) + obs value within tolerance).

Run from the repository root or from the Jupyter folder.
"""
import os
from datetime import datetime, timedelta
import statistics
import numpy as np

# Import readers from the local read_GEOSldas.py
from read_GEOSldas import read_ObsFcstAna, read_obs_param

# Configuration (match plot_obs_maps.py defaults)
EXPDIR = "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/"
EXPID = "LS_OLv8_M36"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
START_DATE = datetime(2020,5,1,0,0)
END_DATE = datetime(2020,5,3,0,0)
DA_DT = 3 * 3600
OBSPARAM_TIME = "20200101_0000"

# tolerance for comparing float observation values (absolute)
OBS_TOL = 1e-6
FCST_LARGE_TOL = 1e-3
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
PAIR_REPORT = os.path.join(SCRIPT_DIR, "ascat_duplicate_report.csv")
TRIP_REPORT = os.path.join(SCRIPT_DIR, "ascat_triplicate_report.csv")


def find_obsparam(rc_out_dir):
    for root, dirs, files in os.walk(rc_out_dir):
        for file in files:
            if file.endswith(f".ldas_obsparam.{OBSPARAM_TIME}z.txt"):
                return os.path.join(root, file)
    return None


def get_obs_file_path(dt):
    year_str = dt.strftime('%Y')
    month_str = dt.strftime('%m')
    datetime_str = dt.strftime('%Y%m%d_%H%M')
    return f"{EXPDIR}/{EXPID}/output/{DOMAIN}/ana/ens_avg/Y{year_str}/M{month_str}/{EXPID}.ens_avg.ldas_ObsFcstAna.{datetime_str}z.bin"


def load_species_for_time(fname, target_species_numbers):
    """Return structured arrays for observations matching any of the target species numbers.
    Returns None if file doesn't exist or no data for species."""
    if not os.path.exists(fname):
        return None
    ofa = read_ObsFcstAna(fname)
    if len(ofa['obs_lon']) == 0:
        return None
    mask = np.isin(ofa['obs_species'], target_species_numbers)
    if not np.any(mask):
        return None
    # Build an array with tilenum, lon, lat, obs, fcst, obsvar, assim
    arr = np.vstack((ofa['obs_tilenum'][mask].astype(int),
                     ofa['obs_lon'][mask].astype(float),
                     ofa['obs_lat'][mask].astype(float),
                     ofa['obs_obs'][mask].astype(float),
                     ofa['obs_fcst'][mask].astype(float),
                     ofa['obs_obsvar'][mask].astype(float),
                     ofa['obs_assim'][mask].astype(int))).T
    # Use structured dtype
    dtype = [('tilenum','i4'),('lon','f4'),('lat','f4'),('obs','f4'),('fcst','f4'),('obsvar','f4'),('assim','i4')]
    return np.array([tuple(row) for row in arr], dtype=dtype)


def make_keys(arr, round_to=6):
    lonr = np.round(arr["lon"], round_to)
    latr = np.round(arr["lat"], round_to)
    keys = list(zip(arr["tilenum"], lonr, latr))
    return keys


def analyze_duplicates(arr_a, arr_b, time_a, time_b, *, write_csv=True, verbose=True):
    """Compare two arrays and optionally write/report duplicates."""
    result = {
        "dup_count": 0,
        "frac_a": 0.0,
        "frac_b": 0.0,
        "fcst_stats": None,
        "dup_info": {},
    }
    if arr_a is None or arr_b is None:
        if verbose:
            print("  One of the files is missing or empty; skipping duplicate analysis.")
        return result

    keys_a = make_keys(arr_a)
    keys_b = make_keys(arr_b)
    obs_a = arr_a["obs"]
    obs_b = arr_b["obs"]

    prev_dict = {}
    for k, obs_val, fcst, obsvar, assim in zip(
        keys_a, obs_a, arr_a["fcst"], arr_a["obsvar"], arr_a["assim"]
    ):
        prev_dict.setdefault(k, []).append(
            (float(obs_val), float(fcst), float(obsvar), int(assim))
        )

    dup_count = 0
    examples = []
    fcst_diffs = []
    large_diff_count = 0
    csv_rows = []
    dup_info = {}

    for k, obs_val, fcst_curr, obsvar_curr, assim_curr in zip(
        keys_b, obs_b, arr_b["fcst"], arr_b["obsvar"], arr_b["assim"]
    ):
        if k not in prev_dict:
            continue
        for (prev_obs, fcst_prev, obsvar_prev, assim_prev) in prev_dict[k]:
            if np.isfinite(prev_obs) and np.isfinite(obs_val) and abs(prev_obs - obs_val) <= OBS_TOL:
                dup_count += 1
                if np.isfinite(fcst_prev) and np.isfinite(fcst_curr):
                    diff = float(fcst_curr) - float(fcst_prev)
                    fcst_diffs.append(diff)
                    if abs(diff) > FCST_LARGE_TOL:
                        large_diff_count += 1
                if len(examples) < 5:
                    examples.append(
                        (
                            k,
                            float(obs_val),
                            float(fcst_prev) if np.isfinite(fcst_prev) else None,
                            float(fcst_curr) if np.isfinite(fcst_curr) else None,
                        )
                    )
                if write_csv:
                    csv_rows.append(
                        (
                            time_a.strftime("%Y%m%d_%H%M"),
                            time_b.strftime("%Y%m%d_%H%M"),
                            k[0],
                            k[1],
                            k[2],
                            float(obs_val),
                            float(fcst_prev) if np.isfinite(fcst_prev) else "",
                            float(fcst_curr) if np.isfinite(fcst_curr) else "",
                            float(obsvar_prev) if np.isfinite(obsvar_prev) else "",
                            float(obsvar_curr) if np.isfinite(obsvar_curr) else "",
                            int(assim_prev),
                            int(assim_curr),
                        )
                    )
                sig = (
                    k[0],
                    round(k[1], 6),
                    round(k[2], 6),
                    round(float(obs_val), 6),
                )
                dup_info[sig] = {
                    "time_a": time_a,
                    "time_b": time_b,
                    "obs": float(obs_val),
                    "fcst_a": float(fcst_prev) if np.isfinite(fcst_prev) else None,
                    "fcst_b": float(fcst_curr) if np.isfinite(fcst_curr) else None,
                    "assim_a": int(assim_prev),
                    "assim_b": int(assim_curr),
                }
                break

    result["dup_count"] = dup_count
    result["frac_a"] = dup_count / len(arr_a) if len(arr_a) > 0 else 0.0
    result["frac_b"] = dup_count / len(arr_b) if len(arr_b) > 0 else 0.0
    result["dup_info"] = dup_info

    if verbose:
        print(f"  Duplicates between {time_a} and {time_b}: {dup_count}")
        print(
            f"    Fraction of file {time_b.strftime('%Y%m%d_%H%M')}: {result['frac_b']:.4f} ({result['frac_b']*100:.2f}%)"
        )
        print(
            f"    Fraction of file {time_a.strftime('%Y%m%d_%H%M')}: {result['frac_a']:.4f} ({result['frac_a']*100:.2f}%)"
        )

    if fcst_diffs:
        mean_diff = statistics.mean(fcst_diffs)
        med_diff = statistics.median(fcst_diffs)
        std_diff = statistics.pstdev(fcst_diffs)
        p95 = np.percentile(np.abs(fcst_diffs), 95)
        if verbose:
            print(
                f"    forecast diff (file2 - file1): mean={mean_diff:.6g}, median={med_diff:.6g}, std={std_diff:.6g}, 95pct_abs={p95:.6g}"
            )
            print(
                f"    large fcst diffs (>|{FCST_LARGE_TOL}|): {large_diff_count} ({large_diff_count/len(fcst_diffs):.3f} of compared)"
            )
        result["fcst_stats"] = {
            "mean": mean_diff,
            "median": med_diff,
            "std": std_diff,
            "p95_abs": p95,
            "large_fraction": large_diff_count / len(fcst_diffs),
        }
    elif verbose:
        print("    No finite forecast pairs found for duplicates")

    if verbose and examples:
        print("  Examples (tilenum, lon, lat, prev_fcst, curr_fcst, obs):")
        for k, obsval, pf, cf in examples:
            print(
                f"    tilenum={k[0]}, lon={k[1]}, lat={k[2]}, prev_fcst={pf}, curr_fcst={cf}, obs={obsval}"
            )

    if write_csv and csv_rows:
        write_header = not os.path.exists(PAIR_REPORT)
        MAX_ROWS_PER_PAIR = 10000
        with open(PAIR_REPORT, "a") as cf:
            if write_header:
                cf.write(
                    "prev_time,curr_time,tilenum,lon,lat,obs,prev_fcst,curr_fcst,prev_obsvar,curr_obsvar,prev_assim,curr_assim\n"
                )
            for row in csv_rows[:MAX_ROWS_PER_PAIR]:
                cf.write(",".join(map(str, row)) + "\n")

    return result


def write_triplicate_report(rows):
    if not rows:
        return
    write_header = not os.path.exists(TRIP_REPORT)
    with open(TRIP_REPORT, "a") as cf:
        if write_header:
            cf.write(
                "prev_time,curr_time,next_time,tilenum,lon,lat,obs,prev_fcst,curr_fcst,next_fcst,prev_assim,curr_assim,next_assim\n"
            )
        for row in rows:
            cf.write(",".join(map(str, row)) + "\n")


def main():
    rc_out_dir = f"{EXPDIR}/{EXPID}/output/{DOMAIN}/rc_out/"
    obsparam_file = find_obsparam(rc_out_dir)
    if not obsparam_file:
        print(f"Could not find obsparam file in {rc_out_dir}")
        return
    obs_param = read_obs_param(obsparam_file)
    # Find species numbers for ASCAT_METB_SM (including variants like _A/_D)
    target_species_numbers = []
    for sp in obs_param:
        descr = sp['descr']
        if 'ASCAT_METB_SM' in descr:
            target_species_numbers.append(int(sp['species']))
    if not target_species_numbers:
        print("No ASCAT_METB_SM species found in obs_param. Aborting.")
        return
    print(f"Target ASCAT_METB_SM species ids: {target_species_numbers}")

    # iterate times and load data
    current = START_DATE
    prev_dt = None
    prev_arr = None

    # Prefetch first/current and next arrays
    curr_arr = load_species_for_time(get_obs_file_path(current), target_species_numbers)
    print(f"{current} -> file: {'FOUND' if curr_arr is not None else 'MISSING/NO_DATA'}")

    next_dt = current + timedelta(seconds=DA_DT)
    next_arr = None
    if next_dt < END_DATE or next_dt == END_DATE:
        next_arr = load_species_for_time(get_obs_file_path(next_dt), target_species_numbers)
        print(f"{next_dt} -> file: {'FOUND' if next_arr is not None else 'MISSING/NO_DATA'}")

    while current < END_DATE:
        prev_result = None
        if prev_arr is not None and curr_arr is not None:
            prev_result = analyze_duplicates(prev_arr, curr_arr, prev_dt, current)

        if (
            prev_result
            and prev_result["dup_count"] > 0
            and curr_arr is not None
            and next_arr is not None
        ):
            next_result = analyze_duplicates(
                curr_arr,
                next_arr,
                current,
                next_dt,
                write_csv=False,
                verbose=False,
            )
            trip_sigs = set(prev_result["dup_info"]).intersection(next_result["dup_info"])
            if trip_sigs:
                print(
                    f"  Triplicate duplicates spanning {prev_dt}, {current}, {next_dt}: {len(trip_sigs)}"
                )
                rows = []
                for sig in list(trip_sigs)[:5]:
                    info_prev = prev_result["dup_info"][sig]
                    info_next = next_result["dup_info"][sig]
                    print(
                        f"    tilenum={sig[0]}, lon={sig[1]}, lat={sig[2]}, obs={info_prev['obs']}, "
                        f"prev_fcst={info_prev['fcst_a']}, curr_fcst={info_prev['fcst_b']}, next_fcst={info_next['fcst_b']}"
                    )
                for sig in trip_sigs:
                    info_prev = prev_result["dup_info"][sig]
                    info_next = next_result["dup_info"][sig]
                    rows.append(
                        (
                            info_prev["time_a"].strftime("%Y%m%d_%H%M"),
                            info_prev["time_b"].strftime("%Y%m%d_%H%M"),
                            info_next["time_b"].strftime("%Y%m%d_%H%M"),
                            sig[0],
                            sig[1],
                            sig[2],
                            info_prev["obs"],
                            info_prev["fcst_a"] if info_prev["fcst_a"] is not None else "",
                            info_prev["fcst_b"] if info_prev["fcst_b"] is not None else "",
                            info_next["fcst_b"] if info_next["fcst_b"] is not None else "",
                            info_prev["assim_a"],
                            info_prev["assim_b"],
                            info_next["assim_b"],
                        )
                    )
                write_triplicate_report(rows)
            else:
                print(
                    f"  No triplicate duplicates spanning {prev_dt}, {current}, {next_dt}"
                )

        prev_arr = curr_arr
        prev_dt = current
        current = next_dt
        curr_arr = next_arr
        next_dt = current + timedelta(seconds=DA_DT)
        if current < END_DATE:
            fname = get_obs_file_path(next_dt)
            next_arr = load_species_for_time(fname, target_species_numbers)
            print(f"{next_dt} -> file: {'FOUND' if next_arr is not None else 'MISSING/NO_DATA'}")
        else:
            next_arr = None

    print("Done.")


if __name__ == '__main__':
    main()
