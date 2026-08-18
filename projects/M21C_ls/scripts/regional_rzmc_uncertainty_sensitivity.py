#!/usr/bin/env python3
"""Compare three serial-correlation-aware intervals for adjacent-period differences.

A  effective-n AR(1) analytic
B  moving-block bootstrap of residuals about period means
C  fitted-AR(1) parametric bootstrap of residuals about period means
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, xarray as xr
from scipy.stats import norm

HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
from changepoint_detection import seasonal_adjustment  # noqa: E402
from m21c_periods import load_period_frames  # noqa: E402

ROOT = HERE.parent
OUT = ROOT / "output" / "regional_rzmc_transitions"
REG = json.loads((ROOT/"config"/"regional_rzmc_regions.json").read_text())["regions"]
BLOCK, NREP, SEED = 24, 2000, 20260818


def bh(p):
    p = np.asarray(p, float); n = p.size; o = np.argsort(p); q = np.empty(n); run = 1.0
    for k in range(n-1, -1, -1):
        run = min(run, p[o[k]]*n/(k+1)); q[o[k]] = run
    return q


def ar1_of(x):
    x = x - x.mean()
    den = (x[:-1]**2).sum()
    return float(np.clip((x[1:]*x[:-1]).sum()/den, 0.0, 0.98)) if den > 0 else 0.0


def main():
    ds = xr.open_dataset(OUT/"regional_rzmc_monthly.nc"); t = pd.DatetimeIndex(ds.time.values)
    _, fine, _, _ = load_period_frames()
    label_of = {r["region_id"]: r["label"] for r in REG}
    order = [r["region_id"] for r in REG]
    ids = [r.period_id for r in fine.itertuples()]
    pidx = {r.period_id: np.where((t >= r.start) & (t <= r.end))[0] for r in fine.itertuples()}
    rng = np.random.default_rng(SEED)
    rows = []

    for rid in order:
        v = ds["delta"].sel(region=rid).values.astype("float64")
        adj, _, _ = seasonal_adjustment(v, t)
        pm = {p: adj[ix].mean() for p, ix in pidx.items()}
        fitted = np.empty_like(adj)
        for p, ix in pidx.items():
            fitted[ix] = pm[p]
        resid = adj - fitted
        rho = ar1_of(resid); sd = resid.std(ddof=1); n = resid.size

        # --- B: moving-block bootstrap ---
        starts = np.arange(0, n - BLOCK + 1)
        nblk = int(np.ceil(n / BLOCK))
        bootB = np.empty((NREP, len(ids)))
        for b in range(NREP):
            pick = rng.choice(starts, size=nblk, replace=True)
            e = np.concatenate([resid[s:s+BLOCK] for s in pick])[:n]
            y = fitted + e
            bootB[b] = [y[pidx[p]].mean() for p in ids]

        # --- C: fitted-AR(1) parametric bootstrap ---
        burn = 120
        bootC = np.empty((NREP, len(ids)))
        innov_sd = sd*np.sqrt(max(1.0 - rho**2, 1e-6))
        for b in range(NREP):
            z = rng.normal(0.0, innov_sd, size=n + burn)
            e = np.empty(n + burn); e[0] = z[0]
            for i in range(1, n + burn):
                e[i] = rho*e[i-1] + z[i]
            y = fitted + e[burn:]
            bootC[b] = [y[pidx[p]].mean() for p in ids]

        for a, bb in zip(ids[:-1], ids[1:]):
            ia, ib = ids.index(a), ids.index(bb)
            d = pm[bb] - pm[a]
            # A: effective-n
            se = 0.0
            for p in (a, bb):
                neff = max(len(pidx[p])*(1-rho)/(1+rho), 1.0)
                se += (sd/np.sqrt(neff))**2
            seA = np.sqrt(se); pA = 2*(1-norm.cdf(abs(d)/seA))
            # B / C: centered bootstrap distributions of the difference
            dB = bootB[:, ib] - bootB[:, ia]; sB = dB.std(ddof=1)
            dC = bootC[:, ib] - bootC[:, ia]; sC = dC.std(ddof=1)
            pB = 2*(1-norm.cdf(abs(d)/sB)); pC = 2*(1-norm.cdf(abs(d)/sC))
            rows.append({"region": label_of[rid], "transition": f"{bb}-{a}", "diff": d,
                         "se_effn": seA, "se_block": sB, "se_ar1boot": sC,
                         "p_effn": pA, "p_block": pB, "p_ar1boot": pC})
        print(f"  {label_of[rid]:28s} rho={rho:.2f}", flush=True)

    df = pd.DataFrame(rows)
    for m in ("effn", "block", "ar1boot"):
        q = np.empty(len(df))
        for tr in df.transition.unique():
            s = df.transition == tr
            q[s.values] = bh(df.loc[s, f"p_{m}"].values)
        df[f"q_{m}"] = q; df[f"sig_{m}"] = q < 0.05
    df.to_csv(OUT/"regional_uncertainty_sensitivity.csv", index=False)

    print("\n=== median SE ratio vs effective-n ===")
    print(f"  block  / effn : {(df.se_block/df.se_effn).median():.2f}")
    print(f"  ar1boot/ effn : {(df.se_ar1boot/df.se_effn).median():.2f}")
    print("\n=== FDR-significant counts (of 48) ===")
    for m in ("effn", "block", "ar1boot"):
        print(f"  {m:8s}: {int(df[f'sig_{m}'].sum())}")
    dis = df[(df.sig_effn != df.sig_block) | (df.sig_effn != df.sig_ar1boot)]
    print(f"\n=== disagreements ({len(dis)}) ===")
    if len(dis):
        print(dis[["region","transition","diff","sig_effn","sig_block","sig_ar1boot"]]
              .to_string(index=False, float_format=lambda x: f"{x:+.5f}"))
    else:
        print("  none - all three methods give the identical significance pattern")


if __name__ == "__main__":
    main()
