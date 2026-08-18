#!/usr/bin/env python3
"""Simple regional RZMC DA-OL figure: monthly series + period-mean steps."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.dates as mdates

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from changepoint_detection import seasonal_adjustment  # noqa: E402
from m21c_periods import load_period_frames  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "output" / "regional_rzmc_transitions"
FIGDIR = ROOT / "docs" / "regional_rzmc_diagnostic_figures"
REG = json.loads((ROOT/"config"/"regional_rzmc_regions.json").read_text())["regions"]
PCOL = ["#ece7d8","#e2eaf3","#e3efe3","#f8e8da","#eae2f0","#dff0ec","#f9e4ec","#e3e9f3","#f0ece0"]


def bh(p):
    p = np.asarray(p, float); n = p.size; o = np.argsort(p)
    q = np.empty(n); run = 1.0
    for k in range(n-1, -1, -1):
        run = min(run, p[o[k]]*n/(k+1)); q[o[k]] = run
    return q


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)
    ds = xr.open_dataset(OUT/"regional_rzmc_monthly.nc"); t = pd.DatetimeIndex(ds.time.values)
    _, fine, _, _ = load_period_frames()
    label_of = {r["region_id"]: r["label"] for r in REG}
    order = [r["region_id"] for r in REG]
    ids = [r.period_id for r in fine.itertuples()]

    stats = {}
    for rid in order:
        v = ds["delta"].sel(region=rid).values.astype("float64")
        adj, _, _ = seasonal_adjustment(v, t)
        s = pd.Series(adj, index=t)
        resid = np.concatenate([(s[(s.index>=r.start)&(s.index<=r.end)].values -
                                 s[(s.index>=r.start)&(s.index<=r.end)].mean())
                                for r in fine.itertuples()])
        x = resid - resid.mean()
        rho = float(np.clip((x[1:]*x[:-1]).sum()/max((x[:-1]**2).sum(), 1e-30), 0.0, 0.98))
        sd = resid.std(ddof=1)
        pm = {}
        for r in fine.itertuples():
            seg = s[(s.index>=r.start)&(s.index<=r.end)]
            neff = max(len(seg)*(1-rho)/(1+rho), 1.0)
            pm[r.period_id] = (float(seg.mean()), float(sd/np.sqrt(neff)), len(seg))
        stats[rid] = {"pm": pm, "rho": rho, "raw": pd.Series(v, index=t)}

    # adjacent differences, FDR one family per transition across regions
    rows = []
    from scipy.stats import norm
    for a, b in zip(ids[:-1], ids[1:]):
        fam = []
        for rid in order:
            m1, e1, _ = stats[rid]["pm"][a]; m2, e2, _ = stats[rid]["pm"][b]
            d = m2 - m1; se = float(np.hypot(e1, e2))
            fam.append({"region": label_of[rid], "transition": f"{b}−{a}",
                        "diff": d, "se": se, "p": 2*(1-norm.cdf(abs(d)/se))})
        q = bh([f["p"] for f in fam])
        for f, qq in zip(fam, q):
            f["q"] = float(qq); f["sig_fdr"] = bool(qq < 0.05); rows.append(f)
    diffs = pd.DataFrame(rows)
    diffs.to_csv(OUT/"regional_period_diffs_fdr.csv", index=False)

    t0, t1 = fine.start.iloc[0], fine.end.iloc[-1]
    sig_lookup = {(r.region, r.transition): bool(r.sig_fdr) for _, r in diffs.iterrows()}

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 7.2), sharex=True)
    for k, rid in enumerate(order):
        ax = axes.flat[k]
        for j, r in enumerate(fine.itertuples()):
            ax.axvspan(r.start, r.end, color=PCOL[j % len(PCOL)], zorder=0, lw=0)
            if j > 0:
                tr = f"{ids[j]}−{ids[j-1]}"
                if sig_lookup.get((label_of[rid], tr), False):
                    ax.axvline(r.start, color="black", lw=1.6, ls="--", zorder=6)
                else:
                    ax.axvline(r.start, color="0.45", lw=0.8, ls="--", zorder=3)
            ax.text(r.start + (r.end-r.start)/2, 0.972, r.period_id,
                    transform=ax.get_xaxis_transform(), ha="center", va="top",
                    fontsize=6.4, fontweight="bold", color="0.25", zorder=7)
        ax.axhline(0, color="0.4", lw=0.7, zorder=1)
        ax.plot(t, stats[rid]["raw"].values, color="#c07a72", lw=0.75, alpha=0.6, zorder=2)
        for r in fine.itertuples():
            m = stats[rid]["pm"][r.period_id][0]
            ax.plot([r.start, r.end], [m, m], color="#8c1d16", lw=1.6,
                    solid_capstyle="butt", zorder=5)
        ax.set_xlim(t0, t1)                      # precise fit to the record
        ax.margins(y=0.16)                       # headroom for the P labels
        ax.set_title(f"({chr(97+k)}) {label_of[rid]}", fontsize=10.5)
        ax.grid(axis="y", color="0.9", lw=0.6); ax.set_axisbelow(True)
        if k % 3 == 0: ax.set_ylabel(r"RZMC DA $-$ OL (m$^3$ m$^{-3}$)")
        ax.xaxis.set_major_locator(mdates.YearLocator(4))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    h = [plt.Line2D([],[],color="#c07a72",lw=1.0,alpha=0.6),
         plt.Line2D([],[],color="#8c1d16",lw=1.6),
         plt.Line2D([],[],color="black",lw=1.6,ls="--"),
         plt.Line2D([],[],color="0.45",lw=0.8,ls="--")]
    fig.legend(h, ["monthly DA $-$ OL", "observing-system period mean",
                   "transition resolved (FDR $q<0.05$)", "period boundary"],
               loc="lower center", ncols=4, frameon=False, fontsize=9.5,
               bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Regional area-weighted RZMC DA $-$ OL by observing-system period, "
                 "June 2000 – May 2024", fontsize=12.5, y=0.98)
    fig.tight_layout(rect=[0, 0.04, 1, 0.955])
    fig.savefig(FIGDIR/"regional_rzmc_period_means.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    pd.set_option("display.width", 220)
    print("=== adjacent differences x10^-3 m3 m-3 (bold = FDR significant per transition family) ===")
    for reg in diffs.region.unique():
        sub = diffs[diffs.region == reg]
        cells = [f"{r.transition}:{r['diff']*1000:+.2f}{'*' if r.sig_fdr else ' '}"
                 for _, r in sub.iterrows()]
        print(f"  {reg:28s} " + "  ".join(cells))
    print(f"\nFDR-significant: {int(diffs.sig_fdr.sum())} of {len(diffs)} tests")
    print("wrote regional_rzmc_period_means.png")


if __name__ == "__main__":
    main()
