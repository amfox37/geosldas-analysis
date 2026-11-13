#!/usr/bin/env python3
"""Per-tile TC / IV maps and summaries for the CYGNSS L1 AZ arms.

Reads the sparse step-4 files written by ``projects/iv_tc`` (keyed by the M36 cell
index ``idx0``), joins them to the 909 AZ domain tiles via
``output/iv_tc_az/tile_map.csv`` (see ``build_iv_tc_tile_map.py``), and draws maps
in the style of the CYGNSS L3 DA paper (figs. 6, S5; IV delta maps):

* TC:  R^2_mod and model error variance sigma^2_mod per run; arm-minus-OL
  dR^2_mod and % d sigma^2 (drawn so that positive = better).
* IV:  anomaly correlation R(model, obs) per reference sensor; arm-minus-OL dR for
  each reference, a sign-agreement map, and the mean of the two.

Nothing is collapsed before the maps: every tile keeps its own value, tiles are
shown on one fixed common mask (the L1 domain, ``in_l1_mask``), and tiles that
fail the sample-count / validity rule are drawn grey, not dropped.  Domain means
are computed afterwards, always with the number of tiles that contributed.

Used from the notebook ``notebooks/iv_tc_az_maps.ipynb`` or as a script:
  python -I scripts/iv_tc_az_maps.py --root <iv_tc_output> --out <figdir>
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parent.parent
DEFAULT_TILE_MAP = PROJ / "output/iv_tc_az/tile_map.csv"
DEFAULT_ROOT = Path("/discover/nobackup/projects/land_da/cygl1_operator_test/iv_tc_output")
DOMAIN = (-118.0, -106.0, 29.0, 40.0)  # lon0, lon1, lat0, lat1
IV_UNDEFINED = -9999.0
GREY = "0.82"

# run label -> (title, short label), in plotting order.  OL is the reference.
RUNS = {"OL": "OL", "L1_filt275": "L1 (filter, errstd 2.75)", "L3": "L3-only"}
REFERENCE = "OL"


# --------------------------------------------------------------------------- paths
def tc_file(root, run, start="20200101", end="20221231", primary="smap_l3", secondary="ascat_h121",
            season="annual", nmin=20):
    trip = f"{primary}__model__{secondary}"
    return Path(root) / "step4_tc" / trip / run / f"{start}_{end}_{season}_n{nmin}.npz"


def iv_file(root, sensor, run, start="20200101", end="20221231", lag=2, nmin=100):
    return Path(root) / "step4_iv" / sensor / run / f"{start}_{end}_lag{lag}_n{nmin}.npz"


# --------------------------------------------------------------------------- loaders
def load_tile_map(path=DEFAULT_TILE_MAP) -> pd.DataFrame:
    tm = pd.read_csv(path)
    return tm.sort_values("tile_axis_index").reset_index(drop=True)


def load_tc(path, tm: pd.DataFrame) -> pd.DataFrame:
    """One row per domain tile.  Columns: N, r2_{pri,mod,sec}, sig2_{pri,mod,sec}.

    The three TC columns are ordered (primary, model, secondary), as written by
    ``compute_tc.py``.  Tiles with no TC estimate are NaN (N = 0)."""
    with np.load(path) as z:
        d = pd.DataFrame({"idx0": z["idx0"].astype(np.int64), "N": z["N_sm"].astype(float)})
        for k, name in enumerate(("pri", "mod", "sec")):
            d[f"r2_{name}"] = z["R2_TC"][:, k]
            d[f"sig2_{name}"] = z["sigma2"][:, k]
        d.attrs["Nmin"] = int(z["Nmin"])
        d.attrs["names"] = tuple(str(s) for s in z["source_names"])
    out = tm[["tile_axis_index", "idx0"]].merge(d, on="idx0", how="left")
    out.attrs.update(d.attrs)
    return out.sort_values("tile_axis_index").reset_index(drop=True)


def load_iv(path, tm: pd.DataFrame) -> pd.DataFrame:
    """One row per domain tile.  R = anomaly correlation of model vs obs."""
    with np.load(path) as z:
        r = z["R_mod_obs"].astype(float)
        r[r <= IV_UNDEFINED + 1] = np.nan
        d = pd.DataFrame({"idx0": z["idx0"].astype(np.int64), "N": z["N_sm"].astype(float), "R": r})
        d.attrs["Nmin"] = int(z["Nmin"])
    out = tm[["tile_axis_index", "idx0"]].merge(d, on="idx0", how="left")
    out.attrs.update(d.attrs)
    return out.sort_values("tile_axis_index").reset_index(drop=True)


# --------------------------------------------------------------------------- masks / metrics
def clamp01(a):
    a = np.asarray(a, float).copy()
    a[(a <= 0) | (a > 1)] = np.nan
    return a


@dataclass
class TCComparison:
    """arm-minus-reference TC fields on one common tile set."""
    tm: pd.DataFrame
    valid: np.ndarray                 # common-valid tiles (inside L1 mask, enough N, finite, both runs)
    d_r2: np.ndarray                  # R2_mod(arm) - R2_mod(ref); positive = better
    pct_sig2_inv: np.ndarray          # -100*(sig2_arm/sig2_ref - 1); positive = better
    r2_ref: np.ndarray
    r2_arm: np.ndarray
    sig2_ref: np.ndarray
    sig2_arm: np.ndarray
    floor: float
    counts: dict = field(default_factory=dict)


def compare_tc(ref: pd.DataFrame, arm: pd.DataFrame, tm: pd.DataFrame, nmin=None, floor_pct=10,
               use_mask=True) -> TCComparison:
    nmin = ref.attrs.get("Nmin", 20) if nmin is None else nmin
    inside = tm.in_l1_mask.to_numpy() if use_mask else np.ones(len(tm), bool)
    r0, r1 = clamp01(ref.r2_mod), clamp01(arm.r2_mod)
    s0, s1 = ref.sig2_mod.to_numpy(float), arm.sig2_mod.to_numpy(float)
    enough = (ref.N.to_numpy() >= nmin) & (arm.N.to_numpy() >= nmin)
    ok = inside & enough & np.isfinite(r0) & np.isfinite(r1) & np.isfinite(s0) & np.isfinite(s1) & (s0 > 0) & (s1 > 0)
    # Paper: mask tiny baseline sigma^2 before a percentage change.
    floor = float(np.nanpercentile(s0[ok], floor_pct)) if ok.any() else np.nan
    pct_ok = ok & (s0 > floor)
    d_r2 = np.where(ok, r1 - r0, np.nan)
    pct = np.where(pct_ok, -100.0 * (s1 / s0 - 1.0), np.nan)
    counts = {
        "tiles_in_mask": int(inside.sum()),
        "tiles_with_both_TC": int((inside & np.isfinite(ref.r2_mod) & np.isfinite(arm.r2_mod)).sum()),
        "tiles_enough_N": int((inside & enough).sum()),
        "tiles_valid": int(ok.sum()),
        "tiles_pct_valid": int(pct_ok.sum()),
    }
    return TCComparison(tm, ok, d_r2, pct, np.where(ok, r0, np.nan), np.where(ok, r1, np.nan),
                        np.where(ok, s0, np.nan), np.where(ok, s1, np.nan), floor, counts)


def compare_iv(ref: pd.DataFrame, arm: pd.DataFrame, tm: pd.DataFrame, nmin=None, use_mask=True):
    """Return (valid, dR) with dR = R(arm) - R(ref); positive = better."""
    nmin = ref.attrs.get("Nmin", 100) if nmin is None else nmin
    inside = tm.in_l1_mask.to_numpy() if use_mask else np.ones(len(tm), bool)
    r0, r1 = ref.R.to_numpy(), arm.R.to_numpy()
    ok = inside & (ref.N.to_numpy() >= nmin) & (arm.N.to_numpy() >= nmin) & np.isfinite(r0) & np.isfinite(r1)
    return ok, np.where(ok, r1 - r0, np.nan)


def summarize(values, valid=None, label=""):
    """Domain statistics of a per-tile field, always with how many tiles contributed."""
    v = np.asarray(values, float)
    m = np.isfinite(v) if valid is None else (np.asarray(valid) & np.isfinite(v))
    x = v[m]
    if x.size == 0:
        return {"label": label, "n_tiles": 0}
    q1, med, q3 = np.nanpercentile(x, [25, 50, 75])
    return {"label": label, "n_tiles": int(x.size), "mean": float(x.mean()), "std": float(x.std()),
            "median": float(med), "IQR": float(q3 - q1), "frac_positive": float((x > 0).mean())}


# --------------------------------------------------------------------------- drawing
def _edges(centers):
    c = np.unique(np.round(np.asarray(centers, dtype=float), 4))
    mid = 0.5 * (c[:-1] + c[1:])
    return c, np.concatenate([[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]])


def map_axes(fig, pos, extent=DOMAIN):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    args = (pos,) if isinstance(pos, int) else tuple(pos)
    ax = fig.add_subplot(*args, projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.5, edgecolor="0.35")
    ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.6, edgecolor="0.25")
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5, edgecolor="0.25")
    return ax


def tile_squares(ax, tm, values, mask_background=True, **kw):
    """Draw every domain tile as its EASE cell.  Tiles with a finite value are coloured; tiles that
    are inside the L1 mask but have no valid value are grey; tiles outside the mask are very light."""
    import cartopy.crs as ccrs
    from matplotlib.collections import PolyCollection
    xc, xe = _edges(tm.com_lon)
    yc, ye = _edges(tm.com_lat)
    ix = np.searchsorted(xc, np.round(tm.com_lon.to_numpy(), 4))
    iy = np.searchsorted(yc, np.round(tm.com_lat.to_numpy(), 4))
    verts = np.array([[(xe[i], ye[j]), (xe[i + 1], ye[j]), (xe[i + 1], ye[j + 1]), (xe[i], ye[j + 1])]
                      for i, j in zip(ix, iy)])
    vals = np.asarray(values, float)
    fin = np.isfinite(vals)
    inside = tm.in_l1_mask.to_numpy()
    if mask_background:
        for sel, col in ((~inside, "0.95"), (inside & ~fin, GREY)):
            if sel.any():
                pc = PolyCollection(verts[sel], facecolor=col, linewidths=0, transform=ccrs.PlateCarree(), zorder=1)
                ax.add_collection(pc)
    vmin, vmax = kw.pop("vmin", None), kw.pop("vmax", None)
    pc = PolyCollection(verts[fin], linewidths=0, transform=ccrs.PlateCarree(), zorder=2, **kw)
    pc.set_array(vals[fin])
    pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    return pc


def draw_field(fig, pos, tm, values, title, cmap, vmin, vmax, cbar_label="", levels=20, stats=True, extend="both"):
    """One map panel with a discrete colour scale and mean/std/n in the title."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm
    ax = map_axes(fig, pos)
    edges = np.linspace(vmin, vmax, levels + 1)
    n_ext = {"both": 2, "min": 1, "max": 1}.get(extend, 0)
    cm = plt.get_cmap(cmap, levels + n_ext)
    pc = tile_squares(ax, tm, values, cmap=cm, norm=BoundaryNorm(edges, cm.N, extend=extend))
    if stats:
        s = summarize(values)
        sub = f"\nn={s['n_tiles']}" + (f", mean {s['mean']:.3g} ± {s['std']:.3g}" if s["n_tiles"] else "")
    else:
        sub = ""
    ax.set_title(title + sub, fontsize=9)
    cb = fig.colorbar(pc, ax=ax, orientation="horizontal", pad=0.05, fraction=0.05, extend=extend)
    cb.set_label(cbar_label, fontsize=9)
    cb.ax.tick_params(labelsize=8)
    return ax


# --------------------------------------------------------------------------- figures
def fig_tc_baseline(tc, tm, run="OL", sig_max=2e-3, savepath=None):
    """Reference run: sigma^2_mod, R^2_mod and sample count per tile."""
    import matplotlib.pyplot as plt
    r = tc[run]
    ok = (r.N >= r.attrs.get("Nmin", 20)).to_numpy()
    fig = plt.figure(figsize=(15, 5.6))
    draw_field(fig, 131, tm, np.where(ok, r.sig2_mod, np.nan), f"{RUNS.get(run, run)}: model error variance σ²_ε",
               "YlGnBu", 0, sig_max, "m⁶ m⁻⁶", extend="max")
    draw_field(fig, 132, tm, np.where(ok, clamp01(r.r2_mod), np.nan), f"{RUNS.get(run, run)}: R²_mod",
               "viridis", 0, 1, "R²", extend="neither")
    draw_field(fig, 133, tm, r.N.where(r.N > 0), "TC sample count N", "magma_r", 0, float(np.nanmax(r.N)),
               "days", extend="neither", stats=False)
    _finish(fig, savepath)
    return fig


def fig_tc_delta(cmp: TCComparison, arm_label, ref_label="OL", d_r2_lim=0.2, pct_lim=60, savepath=None):
    """arm - reference.  Both panels: positive (red) = arm better."""
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11, 5.6))
    draw_field(fig, 121, cmp.tm, cmp.d_r2, f"{arm_label} − {ref_label}: ΔR²_mod", "RdBu_r", -d_r2_lim, d_r2_lim, "ΔR²")
    draw_field(fig, 122, cmp.tm, cmp.pct_sig2_inv, f"{ref_label} − {arm_label}: % Δσ²_ε", "RdBu_r", -pct_lim, pct_lim, "% (positive = lower error)")
    fig.text(0.5, 0.005, f"valid {cmp.counts['tiles_valid']} of {cmp.counts['tiles_in_mask']} mask tiles; "
             f"σ² floor {cmp.floor:.2e}; improved R²: {np.nanmean(cmp.d_r2[cmp.valid] > 0):.0%}", ha="center", fontsize=9)
    _finish(fig, savepath)
    return fig


def fig_iv_delta(iv, tm, arm, ref="OL", arm_label="", sensors=("smap_l3", "ascat_h121"), lim=0.2, savepath=None):
    """ΔanomR vs each reference sensor, the sign-agreement class map, and their mean."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    dr = {}
    for s in sensors:
        _, dr[s] = compare_iv(iv[(s, ref)], iv[(s, arm)], tm)
    a, b = dr[sensors[0]], dr[sensors[1]]
    fig = plt.figure(figsize=(19, 5.6))
    for k, s in enumerate(sensors):
        draw_field(fig, (1, 4, k + 1), tm, dr[s], f"{arm_label} − {ref}: ΔanomR vs {s}", "RdBu_r", -lim, lim, "ΔanomR")
    both = np.isfinite(a) & np.isfinite(b)
    cls = np.full(len(tm), np.nan)
    cls[both & (a > 0) & (b > 0)] = 3
    cls[both & (a > 0) & (b <= 0)] = 2
    cls[both & (a <= 0) & (b > 0)] = 1
    cls[both & (a <= 0) & (b <= 0)] = 0
    ax = map_axes(fig, (1, 4, 3))
    cm = ListedColormap(["#4575b4", "#abd9e9", "#fdae61", "#d73027"])
    pc = tile_squares(ax, tm, cls, cmap=cm, norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cm.N))
    cb = fig.colorbar(pc, ax=ax, orientation="horizontal", pad=0.05, fraction=0.05, ticks=[0, 1, 2, 3])
    cb.ax.set_xticklabels([f"both ↓", f"{sensors[1]} ↑", f"{sensors[0]} ↑", "both ↑"], fontsize=7)
    ax.set_title(f"sign agreement\nn={int(both.sum())}", fontsize=9)
    mean = np.where(both, 0.5 * (a + b), np.nan)
    draw_field(fig, (1, 4, 4), tm, mean, "mean of the two ΔanomR", "RdBu_r", -lim, lim, "ΔanomR")
    _finish(fig, savepath)
    return fig, dr


def _finish(fig, savepath):
    import matplotlib.pyplot as plt
    fig.subplots_adjust(wspace=0.08, bottom=0.12, top=0.88)
    if savepath:
        Path(savepath).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(savepath, dpi=170, bbox_inches="tight", facecolor="white")


# --------------------------------------------------------------------------- driver
def load_all(root, tm, runs=tuple(RUNS), start="20200101", end="20221231", tc_nmin=20, iv_nmin=100,
             primary="smap_l3", secondary="ascat_h121", sensors=("smap_l3", "ascat_h121")):
    tc = {r: load_tc(tc_file(root, r, start, end, primary, secondary, nmin=tc_nmin), tm) for r in runs}
    iv = {(s, r): load_iv(iv_file(root, s, r, start, end, nmin=iv_nmin), tm) for s in sensors for r in runs}
    return tc, iv


def summary_table(tc, iv, tm, runs, ref=REFERENCE, tc_nmin=None, iv_nmin=None):
    """One row per arm: per-tile deltas summarized over the common mask, with tile counts."""
    rows = []
    for r in runs:
        if r == ref:
            continue
        cmp = compare_tc(tc[ref], tc[r], tm, nmin=tc_nmin)
        row = {"arm": r, "tc_tiles_valid": cmp.counts["tiles_valid"], "tc_tiles_in_mask": cmp.counts["tiles_in_mask"],
               "R2_mod_ref_mean": np.nanmean(cmp.r2_ref), "R2_mod_arm_mean": np.nanmean(cmp.r2_arm),
               "dR2_mean": np.nanmean(cmp.d_r2), "dR2_median": np.nanmedian(cmp.d_r2),
               "dR2_frac_improved": float(np.nanmean(cmp.d_r2[cmp.valid] > 0)) if cmp.valid.any() else np.nan,
               "pct_sig2_inv_median": np.nanmedian(cmp.pct_sig2_inv)}
        for (s, rr), _ in list(iv.items()):
            if rr != r:
                continue
            ok, d = compare_iv(iv[(s, ref)], iv[(s, r)], tm, nmin=iv_nmin)
            row[f"dAnomR_{s}_mean"] = np.nanmean(d)
            row[f"dAnomR_{s}_frac_improved"] = float(np.nanmean(d[ok] > 0)) if ok.any() else np.nan
            row[f"dAnomR_{s}_n"] = int(ok.sum())
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", type=Path, default=DEFAULT_ROOT, help="iv_tc output root (step4_tc / step4_iv under it)")
    p.add_argument("--out", type=Path, required=True, help="directory for figures and the summary CSV")
    p.add_argument("--tile-map", type=Path, default=DEFAULT_TILE_MAP)
    p.add_argument("--start", default="20200101")
    p.add_argument("--end", default="20221231")
    p.add_argument("--tc-nmin", type=int, default=20)
    p.add_argument("--iv-nmin", type=int, default=100)
    a = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    tm = load_tile_map(a.tile_map)
    runs = tuple(RUNS)
    tc, iv = load_all(a.root, tm, runs, a.start, a.end, a.tc_nmin, a.iv_nmin)
    a.out.mkdir(parents=True, exist_ok=True)
    fig_tc_baseline(tc, tm, REFERENCE, savepath=a.out / "tc_baseline_OL.png")
    for r in runs:
        if r == REFERENCE:
            continue
        fig_tc_delta(compare_tc(tc[REFERENCE], tc[r], tm, nmin=a.tc_nmin), RUNS[r], savepath=a.out / f"tc_delta_{r}.png")
        fig_iv_delta(iv, tm, r, REFERENCE, RUNS[r], savepath=a.out / f"iv_delta_{r}.png")
    tab = summary_table(tc, iv, tm, runs, tc_nmin=a.tc_nmin, iv_nmin=a.iv_nmin)
    tab.to_csv(a.out / "iv_tc_summary.csv", index=False)
    print(tab.to_string(index=False))
    print(f"figures and summary -> {a.out}")


if __name__ == "__main__":
    main()
