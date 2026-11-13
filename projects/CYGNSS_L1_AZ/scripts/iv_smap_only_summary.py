#!/usr/bin/env python3
"""SMAP-L3-only IV summary: per-tile dR = R(arm) - R(OL) vs SMAP L3, with 2-degree block-bootstrap CIs.

Tiles are spatially correlated, so the bootstrap resamples 2x2-degree blocks, not tiles.
Reports the full SMAP-valid tile set, latitude bands, and the split between tiles that are
also valid for ASCAT and tiles that are SMAP-only (to show the tile set is not driving the result).

  python -I scripts/iv_smap_only_summary.py --root <iv_tc_output>
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import iv_tc_az_maps as m  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=m.DEFAULT_ROOT)
    p.add_argument("--nboot", type=int, default=2000)
    a = p.parse_args()

    tm = m.load_tile_map()
    iv = {(s, r): m.load_iv(m.iv_file(a.root, s, r), tm) for s in ("smap_l3", "ascat_h121") for r in m.RUNS}
    rng = np.random.default_rng(0)
    lat = tm.com_lat.to_numpy()
    blk = (np.floor(tm.com_lon / 2).astype(int).astype(str) + "_" + np.floor(tm.com_lat / 2).astype(int).astype(str)).to_numpy()
    ub = np.unique(blk)
    by_block = {k: np.flatnonzero(blk == k) for k in ub}

    def ci(d, sel):
        v = np.where(sel, d, np.nan)
        boots = []
        for _ in range(a.nboot):
            idx = np.concatenate([by_block[k] for k in rng.choice(ub, len(ub))])
            x = v[idx]
            x = x[np.isfinite(x)]
            if x.size:
                boots.append(x.mean())
        lo, hi = np.percentile(boots, [2.5, 97.5])
        return np.nanmean(v), lo, hi, int(np.isfinite(v).sum()), float(np.nanmean(v[np.isfinite(v)] > 0))

    def line(label, d, sel):
        mu, lo, hi, n, fp = ci(d, sel)
        print(f"   {label:<28s} n={n:3d}  mean dR={mu:+.4f}  95% CI [{lo:+.4f}, {hi:+.4f}]  frac>0={fp:.2f}")

    okA, _ = m.compare_iv(iv[("ascat_h121", "OL")], iv[("ascat_h121", "L1_filt275")], tm)
    for arm in [r for r in m.RUNS if r != m.REFERENCE]:
        ok, d = m.compare_iv(iv[("smap_l3", "OL")], iv[("smap_l3", arm)], tm)
        print(f"{arm}  (OL mean R vs SMAP = {np.nanmean(np.where(ok, iv[('smap_l3', 'OL')].R, np.nan)):.3f})")
        line("all SMAP-valid tiles", d, ok)
        for lo_, hi_, lab in ((29, 32, "29-32N"), (32, 35, "32-35N"), (35, 40, "35-40N")):
            s = ok & (lat >= lo_) & (lat < hi_)
            if s.sum() > 5:
                line(lab, d, s)
        line("also ASCAT-valid", d, ok & okA)
        line("SMAP-only tiles", d, ok & ~okA)


if __name__ == "__main__":
    main()
