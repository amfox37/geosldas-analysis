"""Validate and merge per-cell IMS contingency counts: replace the 2024 contribution of the original 2000-2024 counts file
with a recomputed (flip-corrected) 2024-only counts file. Pooled scopes are additive in A,B,C,D,N:
  ALL_PERIOD            = old - old[YEAR 2024]            + new[YEAR 2024]
  SEASON_ALL_YEARS(s)   = old - old[YEAR_SEASON 2024, s]  + new[YEAR_SEASON 2024, s]
  YEAR 2024 / YEAR_SEASON 2024 = new
Usage: python3 -I merge_counts.py check   (validation only)
       python3 -I merge_counts.py merge <out.nc4>"""
import sys, shutil, numpy as np, netCDF4 as nc
W = "/discover/nobackup/projects/land_da/VIIRS_testing/ims_m21c_fix/"
OLD = "/discover/nobackup/projects/land_da/geosldas-analysis/projects/IMS/output/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2000_2024_thr0p50.nc4"
NEWF = W + "out_new2024/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2024_2024_thr0p50.nc4"
OLDF = W + "out_old2024/ims_ol_da_cell_counts_metrics_SMAP_EASEv2_M36_GLOBAL_2024_2024_thr0p50.nc4"
CNT = ["A", "B", "C", "D", "N"]

def meta(d):
    return list(zip(d["scope_type_code"][:].tolist(), np.ma.filled(d["scope_year"][:], -1).tolist(), np.ma.filled(d["scope_season_code"][:], -1).tolist()))

def idx(m, t, y, s):
    return m.index((t, y, s))

def load(d, k):
    return np.ma.filled(d[k][:], 0).astype(np.int64)

old, new, olx = nc.Dataset(OLD), nc.Dataset(NEWF), nc.Dataset(OLDF)
mo, mn, mx = meta(old), meta(new), meta(olx)
print("old scopes", len(mo), "new scopes", len(mn), mn)
assert (old["cell_i"][:] == new["cell_i"][:]).all() and (old["cell_j"][:] == new["cell_j"][:]).all()
assert mn == mx
# map new scope -> old scope index
pairs = {i: idx(mo, *m) for i, m in enumerate(mn)}
ok = True
for k in CNT:
    o, x, n = load(old, k), load(olx, k), load(new, k)
    for i, j in pairs.items():
        t = mn[i][0]
        if t in (2, 3):  # YEAR / YEAR_SEASON 2024: old-IMS rerun must equal the original
            same = np.array_equal(o[:, j], x[:, i])
            if not same: ok = False; print("MISMATCH rerun-vs-original", k, mn[i], int(np.abs(o[:, j] - x[:, i]).sum()))
    # additivity of the original file: ALL_PERIOD == sum of YEAR scopes
    yrs = [a for a, m in enumerate(mo) if m[0] == 2]
    for e in range(2):
        s = o[e, yrs].sum(0) if o.ndim == 3 else None
    if o.ndim == 3:
        tot = o[:, yrs, :].sum(1); ap = o[:, idx(mo, 0, -1, -1), :]
        if not np.array_equal(tot, ap): ok = False; print("NOT ADDITIVE ALL_PERIOD", k, int(np.abs(tot - ap).sum()))
        for sc in range(4):
            ys = [a for a, m in enumerate(mo) if m[0] == 3 and m[2] == sc]
            if not np.array_equal(o[:, ys, :].sum(1), o[:, idx(mo, 1, -1, sc), :]): ok = False; print("NOT ADDITIVE SEASON", k, sc)
print("VALIDATION", "PASSED" if ok else "FAILED")
if sys.argv[1] == "check" or not ok:
    sys.exit(0 if ok else 1)

out = sys.argv[2]
shutil.copyfile(OLD, out)
d = nc.Dataset(out, "a")
for k in CNT:
    o, n = load(old, k), load(new, k)
    m = o.copy()
    ap, y24 = idx(mo, 0, -1, -1), idx(mo, 2, 2024, -1)
    ny24 = idx(mn, 2, 2024, -1)
    m[:, ap, :] = o[:, ap, :] - o[:, y24, :] + n[:, ny24, :]
    m[:, y24, :] = n[:, ny24, :]
    for sc in range(4):
        ys = idx(mo, 3, 2024, sc); ns = idx(mn, 3, 2024, sc); sa = idx(mo, 1, -1, sc)
        m[:, sa, :] = o[:, sa, :] - o[:, ys, :] + n[:, ns, :]
        m[:, ys, :] = n[:, ns, :]
    # the new file's own ALL_PERIOD / SEASON_ALL_YEARS are 2024-only and are not used
    d[k][:] = m
    print(k, "total change in ALL_PERIOD:", int((m[:, ap, :] - o[:, ap, :]).sum()))
d.history = "2024 contribution replaced with flip-corrected IMS (merge_counts.py); metric variables (accuracy etc.) NOT yet recomputed -- regenerate with --reuse-cell-counts-nc"
d.close()
print("wrote", out)
