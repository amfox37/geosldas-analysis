#!/usr/bin/env python3
"""
Scan consecutive ObsFcstAna binary files and report identical ASCAT_METB_SM observations
(as defined by matching tile number + lon/lat (rounded) + obs value within tolerance).

Run from the repository root or from the Jupyter folder.
"""
import os
from datetime import datetime, timedelta
import numpy as np

# Import readers from the local read_GEOSldas.py
from read_GEOSldas import read_ObsFcstAna, read_obs_param

# Configuration (match plot_obs_maps.py defaults)
EXPDIR = "/discover/nobackup/projects/land_da/CYGNSS_Experiments/DAv8_M36_cd_all"
EXPID = "DAv8_M36_cd_all"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
START_DATE = datetime(2018,8,1,0,0)
END_DATE = datetime(2018,8,5,0,0)
DA_DT = 3 * 3600
OBSPARAM_TIME = "20180801_0000"

# tolerance for comparing float observation values (absolute)
OBS_TOL = 1e-6


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

    # Prepare CSV report rows
    report_rows = []

    while current < END_DATE:
        fname = get_obs_file_path(current)
        arr = load_species_for_time(fname, target_species_numbers)
        print(f"{current} -> file: {'FOUND' if arr is not None else 'MISSING/NO_DATA'}")

        if prev_arr is not None and arr is not None:
            # Compare prev_arr and arr for duplicates
            # We'll match on tilenum and rounded lon/lat (6 decimal places) and obs equal within tol
            # Create keys
            def make_keys(a):
                lonr = np.round(a['lon'],6)
                latr = np.round(a['lat'],6)
                obsr = a['obs']
                keys = list(zip(a['tilenum'], lonr, latr))
                return keys, obsr

            keys_prev, obs_prev = make_keys(prev_arr)
            keys_curr, obs_curr = make_keys(arr)

            # Build dict for prev including fcst/obsvar/assim
            prev_dict = {}
            for k,v,fcst,obsvar,assim in zip(keys_prev, obs_prev, prev_arr['fcst'], prev_arr['obsvar'], prev_arr['assim']):
                prev_dict.setdefault(k, []).append((float(v), float(fcst), float(obsvar), int(assim)))

            # Check matches and compare forecasts
            dup_count = 0
            examples = []
            fcst_diffs = []
            large_diff_count = 0
            FCST_LARGE_TOL = 1e-3
            csv_rows = []

            for k, v, fcst_cur, obsvar_cur, assim_cur in zip(keys_curr, obs_curr, arr['fcst'], arr['obsvar'], arr['assim']):
                if k in prev_dict:
                    prev_entries = prev_dict[k]
                    for (pv, fcst_prev, obsvar_prev, assim_prev) in prev_entries:
                        if np.isfinite(pv) and np.isfinite(v) and abs(float(pv)-float(v)) <= OBS_TOL:
                            dup_count += 1
                            # compute fcst diff if both finite
                            if np.isfinite(fcst_prev) and np.isfinite(fcst_cur):
                                diff = float(fcst_cur) - float(fcst_prev)
                                fcst_diffs.append(diff)
                                if abs(diff) > FCST_LARGE_TOL:
                                    large_diff_count += 1
                            # collect small number of examples
                            if len(examples) < 5:
                                examples.append((k, float(v), float(fcst_prev) if np.isfinite(fcst_prev) else None, float(fcst_cur) if np.isfinite(fcst_cur) else None))
                            # save csv row (cap will be applied later)
                            csv_rows.append((prev_dt.strftime('%Y%m%d_%H%M'), current.strftime('%Y%m%d_%H%M'), k[0], k[1], k[2], float(v), float(fcst_prev) if np.isfinite(fcst_prev) else '', float(fcst_cur) if np.isfinite(fcst_cur) else '', float(obsvar_prev) if np.isfinite(obsvar_prev) else '', float(obsvar_cur) if np.isfinite(obsvar_cur) else '', int(assim_prev), int(assim_cur)))
                            break

            n_prev = len(prev_arr)
            n_curr = len(arr)
            frac_curr = dup_count / n_curr if n_curr > 0 else 0.0
            frac_prev = dup_count / n_prev if n_prev > 0 else 0.0


            print(f"  Number of obs in {current}: {n_curr}") 
            print(f"  Duplicates between {prev_dt} and {current}: {dup_count}")
            print(f"    Fraction of current file duplicated: {frac_curr:.4f} ({frac_curr*100:.2f}%)")
            print(f"    Fraction of previous file duplicated: {frac_prev:.4f} ({frac_prev*100:.2f}%)")

            if fcst_diffs:
                import statistics
                mean_diff = statistics.mean(fcst_diffs)
                med_diff = statistics.median(fcst_diffs)
                std_diff = statistics.pstdev(fcst_diffs)
                p95 = np.percentile(np.abs(fcst_diffs),95)
                print(f"    forecast diff (curr - prev): mean={mean_diff:.6g}, median={med_diff:.6g}, std={std_diff:.6g}, 95pct_abs={p95:.6g}")
                print(f"    large fcst diffs (>|{FCST_LARGE_TOL}|): {large_diff_count} ({large_diff_count/len(fcst_diffs):.3f} of compared)")
            else:
                print("    No finite forecast pairs found for duplicates")

            if examples:
                print("  Examples (tilenum, lon, lat, prev_fcst, curr_fcst, obs):")
                for ex in examples:
                    k, obsval, pf, cf = ex[0], ex[1], ex[2], ex[3]
                    print(f"    tilenum={k[0]}, lon={k[1]}, lat={k[2]}, prev_fcst={pf}, curr_fcst={cf}, obs={obsval}")

            # write CSV rows for this pair (append to global file)
            if csv_rows:
                out_csv = os.path.join(os.path.dirname(__file__), 'ascat_duplicate_report.csv')
                write_header = not os.path.exists(out_csv)
                # cap total appended rows per pair to avoid huge files
                MAX_ROWS_PER_PAIR = 10000
                with open(out_csv, 'a') as cf:
                    if write_header:
                        cf.write('prev_time,curr_time,tilenum,lon,lat,obs,prev_fcst,curr_fcst,prev_obsvar,curr_obsvar,prev_assim,curr_assim\n')
                    for row in csv_rows[:MAX_ROWS_PER_PAIR]:
                        cf.write(','.join(map(str,row)) + '\n')

        prev_arr = arr
        prev_dt = current
        current += timedelta(seconds=DA_DT)

    print("Done.")

    # Write CSV report if we collected any rows
    if report_rows:
        import csv
        out_csv = os.path.join(os.path.dirname(__file__), 'ascat_duplicate_report.csv')
        with open(out_csv, 'w', newline='') as csvfile:
            fieldnames = ['prev_time','curr_time','n_prev','n_curr','n_duplicates','frac_curr','frac_prev']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for r in report_rows:
                writer.writerow(r)
        print(f"Wrote duplicate report: {out_csv}")


if __name__ == '__main__':
    main()
