#!/usr/bin/env python3
"""Per-tile IGBP land-cover fractions for the 909 AZ-domain EASEv2 M36 tiles.

Source: MODIS MCD12C1 v061 (0.05 deg CMG, IGBP = Land_Cover_Type_1_Percent), raw HDF4
files in <repo>/data/landcover_igbp/ (download from LP DAAC / Earthdata, not done here).
HDF4 is read by shelling out to `hdp dumpsds` (cygnss env) -> flat uint8 cache in
data/landcover_igbp/cache/ (git-ignored).  Layer 0 = water (IGBP 17), layer k = IGBP k.

Footprint: the EASEv2 M36 cell is exactly a rectangle in EPSG:6933, and a 0.05 deg CMG
pixel is also an exact rectangle there (lon -> x linear; lat -> y monotonic), so the
overlap area of every CMG pixel with the cell is computed exactly (no sampling).
Per-tile class fraction = sum(overlap * pct_k) / sum(overlap * sum_k pct_k) over non-fill
pixels; valid_frac = overlap-weighted share of the cell covered by non-fill pixels.

  python -I projects/CYGNSS_L1_AZ/scripts/benefit_drivers/build_igbp_landcover.py [--years 2020 2021 2022]
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[1]
REPO = PROJ.parents[1]
RAW = REPO / "data/landcover_igbp"
OUT = PROJ / "output/benefit_drivers"
HDP = "/Users/amfox/mamba/envs/cygnss/bin/hdp"

NX, NY = 964, 406
CS = 36032.220840584  # EASEv2 M36 cell size (m) = 4 * 9008.055210146
X0, Y0 = -NX * CS / 2, NY * CS / 2  # upper-left corner (x, y)
R_X = 6378137.0 * np.cos(np.pi / 6) / np.sqrt(1 - 0.0066943799901413165 * 0.25)  # EPSG:6933 x = R_X * lon_rad (lat_ts 30)
NLAT, NLON, NC = 3600, 7200, 17
RES = 0.05

NAMES = {1: "ENF", 2: "EBF", 3: "DNF", 4: "DBF", 5: "MF", 6: "CSH", 7: "OSH", 8: "WSV", 9: "SAV",
         10: "GRS", 11: "WET", 12: "CRP", 13: "URB", 14: "CVM", 15: "SNO", 16: "BAR", 17: "WAT"}
GROUPS = {"forest": [1, 2, 3, 4, 5], "shrub": [6, 7], "savanna": [8, 9], "grass": [10],
          "crop_mosaic": [12, 14], "barren": [16], "water": [17], "wetland": [11],
          "urban": [13], "snowice": [15]}


def cell_corners(i, j):
    x0 = X0 + i * CS
    y1 = Y0 - j * CS
    return x0, x0 + CS, y1 - CS, y1  # xmin, xmax, ymin, ymax


def load_percent(year):
    """(3600, 7200, 17) uint8 view; layer order as stored (0 = water)."""
    hdfs = sorted(RAW.glob(f"MCD12C1.A{year}001.061.*.hdf"))
    if not hdfs:
        raise FileNotFoundError(f"no MCD12C1 {year} HDF in {RAW}")
    cache = RAW / "cache" / f"p{year}.bin"
    if not cache.exists() or cache.stat().st_size != NLAT * NLON * NC:
        cache.parent.mkdir(exist_ok=True)
        subprocess.run([HDP, "dumpsds", "-b", "-d", "-n", "Land_Cover_Type_1_Percent",
                        "-o", str(cache), str(hdfs[0])], check=True, stdout=subprocess.DEVNULL)
    return np.memmap(cache, dtype=np.uint8, mode="r", shape=(NLAT, NLON, NC)), hdfs[0].name


def aggregate(pct, i, j, fwd):
    xmin, xmax, ymin, ymax = cell_corners(i, j)
    # cell -> lon/lat box (exact: x <-> lon linear, y <-> lat monotonic)
    lon_l, lat_lo = fwd(xmin, ymin)
    lon_r, lat_hi = fwd(xmax, ymax)
    c0 = int(np.floor((lon_l + 180) / RES)) - 1
    c1 = int(np.ceil((lon_r + 180) / RES)) + 1
    r0 = int(np.floor((90 - lat_hi) / RES)) - 1
    r1 = int(np.ceil((90 - lat_lo) / RES)) + 1
    r0, c0, r1, c1 = max(r0, 0), max(c0, 0), min(r1, NLAT), min(c1, NLON)
    rows = np.arange(r0, r1)
    cols = np.arange(c0, c1)
    # CMG pixel edges in projected space
    lon_e = -180 + RES * np.arange(c0, c1 + 1)
    lat_e = 90 - RES * np.arange(r0, r1 + 1)
    xe = np.deg2rad(lon_e) * R_X
    ye = np.array([inv(0.0, la)[1] for la in lat_e])  # y of each lat edge
    ox = np.clip(np.minimum(xe[1:], xmax) - np.maximum(xe[:-1], xmin), 0, None)
    oy = np.clip(np.minimum(ye[:-1], ymax) - np.maximum(ye[1:], ymin), 0, None)
    w = oy[:, None] * ox[None, :] / (CS * CS)  # fraction of the cell covered by each pixel
    blk = np.asarray(pct[r0:r1, c0:c1, :]).astype(np.float64)
    fill = blk[..., 0] == 255
    blk[blk == 255] = 0
    wv = np.where(fill, 0.0, w)
    cls = np.tensordot(wv, blk, axes=([0, 1], [0, 1]))  # (17,) layer order 0..16
    tot = cls.sum()
    frac = np.zeros(17)
    if tot > 0:
        # IGBP 1..16 = layers 1..16 ; IGBP 17 (water) = layer 0
        frac = np.concatenate([cls[1:], cls[:1]]) / tot
    return frac, wv.sum(), w.sum(), tot / max(wv.sum(), 1e-12) / 100.0


_t = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)
_ti = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)
fwd = _t.transform
inv = _ti.transform


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--years", nargs="+", type=int, default=[2020, 2021, 2022])
    ap.add_argument("--main-year", type=int, default=2021)
    a = ap.parse_args()

    cov = pd.read_csv(OUT / "covariate_table.csv",
                      usecols=["tile_axis_index", "tile_id", "idx0", "com_lon", "com_lat", "in_l1_mask",
                               "stat_modis_landcover_mode", "stat_mwrtm_vegcls"])
    assert len(cov) == 909
    cov["i"] = cov.idx0 % NX
    cov["j"] = cov.idx0 // NX

    # step 1: verify geometry against com_lon/com_lat
    cx = X0 + (cov.i + 0.5) * CS
    cy = Y0 - (cov.j + 0.5) * CS
    clon, clat = fwd(cx.values, cy.values)
    dlon = clon - cov.com_lon
    dlat = clat - cov.com_lat
    km = np.hypot(dlon * 111.32 * np.cos(np.deg2rad(cov.com_lat)), dlat * 111.2)
    geo = (f"cell-centre vs com: max {km.max():.1f} km, median {np.median(km):.1f} km "
           f"(cell = 36 km); max |dlon|={np.abs(dlon).max():.3f} deg, max |dlat|={np.abs(dlat).max():.3f} deg; "
           f"com inside cell for {int(((np.abs(dlon*np.cos(np.deg2rad(cov.com_lat)))*111.32<18)&(np.abs(dlat)*111.2<18)).sum())}/909")
    print(geo)

    res = {}
    for y in a.years:
        pct, name = load_percent(y)
        F = np.zeros((909, 17))
        V = np.zeros(909)
        S = np.zeros(909)
        for n, (i, j) in enumerate(zip(cov.i, cov.j)):
            F[n], V[n], _, S[n] = aggregate(pct, int(i), int(j), fwd)
        res[y] = (F, V, S, name)
        print(y, name, "sum frac min/max", F.sum(1).min(), F.sum(1).max(), "valid min", V.min())

    F, V, S, name = res[a.main_year]
    out = cov[["tile_axis_index", "tile_id", "idx0"]].copy()
    for k in range(17):
        out[f"igbp_{k+1:02d}"] = F[:, k]
    out["valid_frac"] = V
    out["pixel_pct_sum_mean"] = S
    out["dominant_class"] = F.argmax(1) + 1
    out["dominant_name"] = [NAMES[c] for c in out.dominant_class]
    out["dominant_frac"] = F.max(1)
    for g, cl in GROUPS.items():
        out[f"frac_{g}"] = F[:, [c - 1 for c in cl]].sum(1)
    out["igbp_year"] = a.main_year
    out.to_csv(OUT / "igbp_landcover_per_tile.csv", index=False, float_format="%.5f")

    # stability across years
    lines = []
    for y in a.years:
        if y == a.main_year:
            continue
        Fy = res[y][0]
        same = (Fy.argmax(1) == F.argmax(1)).mean()
        lines.append(f"- {y} vs {a.main_year}: dominant class identical in {same*100:.1f}% of tiles; "
                     f"max |dfrac| per tile-class = {np.abs(Fy-F).max():.3f}; mean abs dfrac = {np.abs(Fy-F).mean():.4f}")

    # validation vs existing columns
    chk = out.copy()
    chk["modis_mode"] = cov.stat_modis_landcover_mode.values
    chk["mwrtm"] = cov.stat_mwrtm_vegcls.values
    chk["in_l1"] = cov.in_l1_mask.values
    vl = []
    for col in ["modis_mode", "mwrtm"]:
        s = chk.dropna(subset=[col])
        agree = (s.dominant_class == s[col]).mean()
        fr = [s[f"igbp_{int(c):02d}"].iloc[n] if 1 <= c <= 17 else np.nan for n, c in enumerate(s[col])]
        vl.append(f"- vs `{col}`: {len(s)} tiles with a value; dominant IGBP class agrees in {agree*100:.1f}% "
                  f"(median IGBP fraction of the existing class = {np.nanmedian(fr):.2f})")
        ct = pd.crosstab(s[col].astype(int).rename("existing"), s.dominant_class.rename("igbp_dom"))
        vl.append("\n```\n" + ct.to_string() + "\n```")
    # class counts over L1 mask
    m = chk[chk.in_l1]
    cc = m.groupby(["dominant_class", "dominant_name"]).size().rename("n_tiles").reset_index()
    cc["mean_dom_frac"] = m.groupby("dominant_class").dominant_frac.mean().reindex(cc.dominant_class).values
    cc.to_csv(OUT / "igbp_class_counts_l1mask.csv", index=False, float_format="%.3f")

    cc_txt = cc.to_string(index=False, float_format=lambda x: f'{x:.2f}')
    sums_ok = np.allclose(F.sum(1), 1, atol=1e-9)
    make_map(out, chk.in_l1.values)

    readme = f"""# IGBP land cover per AZ M36 tile

Output: `igbp_landcover_per_tile.csv` (909 tiles; key tile_axis_index, tile_id, idx0), built by
`scripts/benefit_drivers/build_igbp_landcover.py`.

## Source
MODIS/Terra+Aqua MCD12C1 **v061**, CMG 0.05 deg, IGBP scheme (`Land_Cover_Type_1_Percent`, uint8 percent,
17 layers; layer 0 = water -> IGBP 17, layer k = IGBP k, 255 = fill). Main year **{a.main_year}**
(`{name}`), also processed {', '.join(str(y) for y in a.years if y != a.main_year)}.
Downloaded from LP DAAC (Earthdata login via ~/.netrc worked) into `data/landcover_igbp/` (git-ignored).
HDF4 read with `hdp dumpsds` (no Python HDF4 reader installed).

## Aggregation
- Cell footprint: EASEv2 M36 cell (964 x 406, cell = {CS:.3f} m, EPSG:6933, upper-left corner
  x={X0:.2f}, y={Y0:.2f}) from idx0 = j*964+i. {geo}.
  The tile centre of mass is NOT assumed to be the cell centre; the footprint comes from the (i,j) index.
- Because x is linear in longitude and y monotonic in latitude, each 0.05 deg CMG pixel is an exact rectangle
  in EPSG:6933; its overlap area with the cell is the exact rectangle intersection (about 7x7 pixels plus edges;
  partial pixels weighted by overlap, so no sampling error).
- Fraction of class k = sum(overlap * pct_k) / sum(overlap * sum of pct over classes) over non-fill pixels, so
  fractions sum to 1 (max |sum-1| = {np.abs(F.sum(1)-1).max():.1e}). `valid_frac` = share of the cell area covered by
  non-fill pixels (min {V.min():.3f}, mean {V.mean():.3f}); `pixel_pct_sum_mean` = mean per-pixel percent total/100
  (min {S.min():.3f}; below 1 means unclassified remainder was renormalised away).
- Dominant class = argmax of the 17 fractions. Groups: forest 1-5, shrub 6-7, savanna 8-9 (WSV+SAV), grass 10,
  crop_mosaic 12+14, barren 16, water 17, wetland 11, urban 13, snowice 15 (`frac_*` columns).

## Validation
- 909 tiles, fractions sum to 1: {sums_ok}. Valid fraction min {V.min():.3f}.
- Year stability:
{chr(10).join(lines)}
- Agreement with existing land-cover columns (those exist only for tiles with L1 observations):
{chr(10).join(vl)}
  Note `stat_modis_landcover_mode` is the mode of per-observation (specular point) land cover, so it
  reflects where CYGNSS samples within a cell, not the full cell area; `stat_mwrtm_vegcls` is the model's
  MWRTM class. Disagreements are expected to be mostly shrub/grass/savanna neighbours.
- Dominant-class counts over the 552 L1-mask tiles: `igbp_class_counts_l1mask.csv`:

```
{cc_txt}
```
- Map: `igbp_dominant_class_map.png`.
"""
    (OUT / "igbp_landcover_README.md").write_text(readme)
    print(cc.to_string(index=False))
    print("\n".join(lines + vl))


def make_map(out, in_l1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cov = pd.read_csv(OUT / "covariate_table.csv", usecols=["com_lon", "com_lat"])
    cols = ["#1b5e20", "#2e7d32", "#558b2f", "#7cb342", "#9ccc65", "#a1887f", "#d7b377", "#c0ca33",
            "#ffd54f", "#aed581", "#4db6ac", "#ffb74d", "#e53935", "#ff8a65", "#e0f7fa", "#bdbdbd", "#1e88e5"]
    cm = ListedColormap(cols)
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(cov.com_lon, cov.com_lat, c=out.dominant_class, cmap=cm, norm=BoundaryNorm(np.arange(.5, 18.5), 17),
                    marker="s", s=70, edgecolor="none")
    ax.scatter(cov.com_lon[in_l1], cov.com_lat[in_l1], facecolor="none", edgecolor="k", marker="s", s=70, lw=0.3)
    cb = fig.colorbar(sc, ticks=range(1, 18), shrink=0.8)
    cb.ax.set_yticklabels([f"{k} {NAMES[k]}" for k in range(1, 18)], fontsize=7)
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(f"Dominant IGBP class per M36 tile (MCD12C1 v061, {out.igbp_year.iloc[0]})\nblack outline = L1 mask")
    fig.tight_layout()
    fig.savefig(OUT / "igbp_dominant_class_map.png", dpi=150)


if __name__ == "__main__":
    main()
