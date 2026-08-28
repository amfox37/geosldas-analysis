#!/usr/bin/env python3
"""Assimilated MODIS snow-cover observation counts are banded in longitude.

read_obs_MODIS_SCF reads MODIS SCF for each assimilation window over the
longitude band where local solar time matches the satellite overpass. With
3-hourly windows those band edges fall on multiples of 45 degrees, and the
band is deliberately widened by 3x the maximum tile longitude extent so that
tiles straddling an edge get all their CMG cells.

The assimilated counts show that widening is not compensated: at several band
edges a single grid column carries almost exactly twice its neighbours' count,
i.e. those tiles are assimilated from two consecutive windows. Other edges
instead show a sustained level shift between adjacent bands. Over the Tibetan
Plateau the 90E edge produces a threefold step in the snow analysis increment
while the model background runs smoothly through it.
"""
from __future__ import annotations
from pathlib import Path
import sys
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT.parent.parent / "common/python/io"))
from read_GEOSldas import read_tilecoord

DATA = Path("/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2")
STATS = DATA / "temporal_stats_DA_20000601_20070531.nc4"    # MODIS-only period
NCOL, SPECIES = 964, {"Terra (11)": [11], "Aqua (12)": [12], "MODIS total": [11, 12]}
EDGES = np.arange(-180, 181, 45)

tc = read_tilecoord(str(DATA / "LS_OLv8_M36.ldas_tilecoord.bin"))
icol = np.asarray(tc["i_indg"]).astype(int)
tlat = np.asarray(tc["com_lat"])
counts = np.asarray(Dataset(STATS)["N_data"][:])
counts[~np.isfinite(counts)] = np.nan
band = (tlat >= 25) & (tlat <= 70)
lon_centre = -180.0 + np.arange(NCOL) * (360.0 / NCOL) + 0.5 * (360.0 / NCOL)


def column_profile(values: np.ndarray, minimum: int = 10) -> np.ndarray:
    out = np.full(NCOL, np.nan)
    for c in range(NCOL):
        keep = band & (icol == c)
        if keep.sum() >= minimum:
            out[c] = np.nanmean(values[keep])
    return out


profiles = {name: column_profile(np.nansum(counts[:, sp], axis=1)) for name, sp in SPECIES.items()}
total = profiles["MODIS total"]

width = 3.0 * float(np.nanmax(np.asarray(tc["max_lon"]) - np.asarray(tc["min_lon"])))
print(f"band widening in read_obs_MODIS_SCF = 3 x max tile lon extent = {width:.3f} deg "
      f"({width / (360 / NCOL):.1f} grid columns)\n")
print(f"  {'edge':>6s} {'col':>5s} {'edge count':>11s} {'neighbours':>11s} {'excess':>8s} {'level shift':>12s}")
for lon in EDGES[:-1]:
    c0 = int(round((lon + 180) / (360.0 / NCOL)))
    window = [c for c in (c0 - 1, c0, c0 + 1) if np.isfinite(total[c])]
    if not window:
        print(f"  {lon:6.0f} {'':>5s} {'insufficient land':>34s}")
        continue
    c = max(window, key=lambda k: total[k])                 # edge column, allowing +-1
    left, right = total[c - 4:c - 1], total[c + 2:c + 5]
    neighbours = np.nanmean(np.concatenate([left, right]))
    print(f"  {lon:6.0f} {c:5d} {total[c]:11.0f} {neighbours:11.0f} "
          f"{100 * (total[c] / neighbours - 1):+7.0f}% {np.nanmean(right) - np.nanmean(left):+12.0f}")

# Two distinct artifacts live at these edges. A doubled column is confined to one
# grid column. A step changes the level on either side with no doubled column, so
# every tile east of the edge is sampled differently from every tile west of it.
# The 90E edge is a pure step, and it reproduces across unrelated terrain.
print("\n90.0E edge, by latitude band: step or doubled column?")
edge = int(round((90.0 + 180) / (360.0 / NCOL)))
print(f"  {'band':>22s} {'west':>7s} {'east':>7s} {'step':>7s} {'vs local grad':>14s}  {'edge column':>14s}")
for lo, hi, label in ((28, 38, "Tibetan Plateau"), (40, 50, "Kazakh / Altai"),
                      (50, 60, "W and C Siberia"), (60, 70, "Northern Siberia")):
    keep = (tlat >= lo) & (tlat <= hi)
    prof = np.array([np.nanmean(np.nansum(counts[:, [11, 12]], axis=1)[keep & (icol == c)])
                     if (keep & (icol == c)).sum() >= 8 else np.nan for c in range(NCOL)])
    if not np.isfinite(prof[edge - 8:edge + 8]).all():
        continue
    west, east = prof[edge - 6:edge], prof[edge:edge + 6]
    grad = abs(np.polyfit(lon_centre[edge - 6:edge], west, 1)[0]) * (360.0 / NCOL)
    step = east[0] - west[-1]
    doubled = prof[edge] > 1.5 * np.nanmean(prof[edge + 1:edge + 4])
    print(f"  {label:>22s} {west.mean():7.0f} {east.mean():7.0f} {step:+7.0f} "
          f"{abs(step) / max(grad, 1e-9):11.0f}x  {'doubled' if doubled else 'pure step':>14s}")

fig, ax = plt.subplots(figsize=(13, 4.8))
for name, colour in (("Terra (11)", "#e07b39"), ("Aqua (12)", "#2b6cb0"), ("MODIS total", "#333333")):
    ax.plot(lon_centre, profiles[name], lw=1.1, color=colour, label=name)
for e in EDGES:
    ax.axvline(e, color="#c0392b", lw=1.0, ls="--", alpha=0.8)
ax.set_xlim(-180, 180); ax.set_xticks(EDGES)
ax.set_xlabel("longitude (°E)")
ax.set_ylabel("assimilated obs per tile\n(2000-06 to 2007-05)")
ax.set_title("Assimilated MODIS snow-cover observation count is banded in longitude, "
             "with discontinuities on 45° multiples (dashed)", loc="left", fontsize=11.5)
ax.legend(frameon=False, ncols=3, fontsize=9.5)
ax.grid(color="0.92"); ax.set_axisbelow(True)
out = ROOT / "output/monthly_synthesis_diagnostics/figures/modis_obs_count_longitude_bands.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=200, bbox_inches="tight")
print("\nwrote", out.relative_to(ROOT))
