#!/usr/bin/env python3
"""First-look tables and figures for analysis A (obs availability), from the per-arm NetCDFs written by
obs_availability.py. Reads config/arms.json for arms, snow years, species and the switch date.

Writes
  output/obs_availability/summary_obs_per_day.csv       arm, snow year, species: mean obs/day, tiles/day, assimilated tiles/day
  output/obs_availability/summary_overlap.csv           arm, snow year: tiles/day with both sensors, unions, only-X
  docs/viirs_transition_report_figures/obs_availability_summary.md    the same tables, as markdown
  docs/viirs_transition_report_figures/figA_*.png       three figures

Counts are domain totals over the 112,573 M36 tiles (southern hemisphere included); n_obs equals n_tiles because a
tile gets at most one obs of a species per day. Snow years are 1 Aug - 31 Jul (arms.json). Run:
  /usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3 obs_availability_summary.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
CFG = json.load(open(PROJECT / "config" / "arms.json"))
OUT = PROJECT / "output" / "obs_availability"
FIG = PROJECT / "docs" / "viirs_transition_report_figures"
FIG.mkdir(parents=True, exist_ok=True)

ARMS = list(CFG["arms"])
SPECIES = [CFG["obsfcstana_species_ids"][k]["name"] for k in "1234"]            # Aqua Terra J1 J2
# fixed categorical order (reference palette slots 1-4); entity keeps its colour in every figure.
COLOR = {"Terra": "#2a78d6", "Aqua": "#eb6834", "J1": "#1baf7a", "J2": "#6250d6"}
STYLE = {"Terra": "-", "Aqua": "-", "J1": "--", "J2": "--"}                     # MODIS solid, VIIRS dashed
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SY = CFG["snow_years"]
SWITCH = np.datetime64(CFG["switch"])


def load(arm):
    d = netCDF4.Dataset(OUT / f"obs_availability_{arm}.nc")
    d.set_auto_mask(False)
    date = np.array([np.datetime64(f"{str(x)[:4]}-{str(x)[4:6]}-{str(x)[6:]}") for x in d["date"][:]])
    r = dict(date=date, steps=d["n_steps"][:].astype(int), bands=list(d["band_name"][:]),
             n_obs=d["n_obs"][:], n_tiles=d["n_tiles"][:], n_assim=d["n_assim"][:])
    for k in ("tiles_both_terra_aqua", "tiles_both_j1_j2", "tiles_any_modis", "tiles_any_viirs", "tiles_terra_only",
              "tiles_aqua_only", "tiles_j1_only", "tiles_j2_only"):
        r[k] = d[k][:]
    d.close()
    return r


def full_days(r):
    """days with all 8 steps (the first and last day of a run have fewer)"""
    return r["steps"] == 8


def sy_mask(r, sy):
    return (r["date"] >= np.datetime64(sy["beg"])) & (r["date"] <= np.datetime64(sy["end"])) & full_days(r)


def monthly_mean(date, x, mask):
    m = date.astype("datetime64[M]")
    keys = np.unique(m[mask])
    return keys.astype("datetime64[D]"), np.array([x[mask & (m == k)].mean() for k in keys])


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def main():
    R = {a: load(a) for a in ARMS if (OUT / f"obs_availability_{a}.nc").exists()}
    missing = [a for a in ARMS if a not in R]
    if missing:
        print("missing arms (not summarised):", missing)

    # ---- tables
    rows_obs, rows_ov = [], []
    for a, r in R.items():
        for sy in SY:
            m = sy_mask(r, sy)
            if m.sum() < 100:
                continue
            for si, sp in enumerate(SPECIES):
                no = r["n_obs"][m, si].sum(axis=1).mean()
                nt = r["n_tiles"][m, si].sum(axis=1).mean()
                na = r["n_assim"][m, si].sum(axis=1).mean()
                rows_obs.append((a, sy["id"], sp, int(m.sum()), no, nt, na))
            rows_ov.append((a, sy["id"], int(m.sum()), *[r[k][m].mean() for k in (
                "tiles_both_terra_aqua", "tiles_any_modis", "tiles_terra_only", "tiles_aqua_only",
                "tiles_both_j1_j2", "tiles_any_viirs", "tiles_j1_only", "tiles_j2_only")]))
    with open(OUT / "summary_obs_per_day.csv", "w") as f:
        f.write("arm,snow_year,species,n_days,mean_obs_per_day,mean_tiles_per_day,mean_assimilated_obs_per_day\n")
        for x in rows_obs:
            f.write(f"{x[0]},{x[1]},{x[2]},{x[3]},{x[4]:.1f},{x[5]:.1f},{x[6]:.1f}\n")
    with open(OUT / "summary_overlap.csv", "w") as f:
        f.write("arm,snow_year,n_days,both_terra_aqua,any_modis,terra_only,aqua_only,both_j1_j2,any_viirs,j1_only,j2_only\n")
        for x in rows_ov:
            f.write(f"{x[0]},{x[1]},{x[2]}," + ",".join(f"{v:.0f}" for v in x[3:]) + "\n")

    md = ["# Analysis A first look: SCF obs availability", "",
          "Domain totals per day over the 112,573 M36 tiles, averaged over the days of each snow year (1 Aug - 31 Jul). "
          "A tile has at most one obs per species per day, so obs/day = tiles/day. 'assim' counts obs with assim_flag = 1. "
          "Zero MODIS in 2026 for the OL and VIIRS arms is the known monitoring gap (their MODIS path ends 2025-12-31).", "",
          "## Obs (= tiles) per day by species", "",
          "| arm | snow year | days | Aqua | Terra | J1 | J2 | assimilated per day |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    by = {}
    for x in rows_obs:
        by.setdefault((x[0], x[1], x[3]), {})[x[2]] = (x[4], x[6])
    for (a, s, n), v in by.items():
        asm = sum(v[sp][1] for sp in SPECIES)
        md.append(f"| {a} | {s} | {n} | " + " | ".join(f"{v[sp][0]:,.0f}" for sp in SPECIES) + f" | {asm:,.0f} |")
    md += ["", "## Overlap per day (tiles)", "",
           "| arm | snow year | Terra and Aqua | any MODIS | Terra only | Aqua only | J1 and J2 | any VIIRS | J1 only | J2 only |",
           "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for x in rows_ov:
        md.append(f"| {x[0]} | {x[1]} | " + " | ".join(f"{v:,.0f}" for v in x[3:]) + " |")
    (FIG / "obs_availability_summary.md").write_text("\n".join(md) + "\n")

    # ---- checks against arms.json, and zero-obs (gap) days
    lines = ["# Analysis A checks", ""]
    names_idx = {sp: i for i, sp in enumerate(SPECIES)}
    gaps = []
    for a, r in R.items():
        date_i = np.array([int(str(x).replace("-", "")) for x in r["date"]])
        tot = r["n_obs"].sum(axis=2)                                   # day x species
        asm = r["n_assim"].sum(axis=2)
        full = full_days(r)
        exp = np.zeros((len(date_i), len(SPECIES)), bool)
        for seg in CFG["arms"][a]["segments"]:
            b, e = int(seg["beg"].replace("-", "")), int(seg["end"].replace("-", ""))
            for spid in seg["assim_species"]:
                nm = CFG["species"][str(spid)]["name"]
                exp[(date_i >= b) & (date_i < e), names_idx[nm]] = True
        # (1) a tile has at most one obs of a species per day
        ok1 = bool((tot == r["n_tiles"].sum(axis=2)).all())
        # (2) assimilated species: every obs is assimilated; others: none is (the segment-end 00z file of the
        #     switch day, 2023-08-01, still belongs to the old segment, so that day is skipped)
        sw = int(CFG["switch"].replace("-", ""))
        chk = full & (date_i != sw)
        bad_a = int((((asm < tot) & exp & (tot > 0)) & chk[:, None]).sum())
        bad_b = int((((asm > 0) & ~exp) & chk[:, None]).sum())
        lines.append(f"- {a}: n_obs == n_tiles on all days: {ok1}; days x species where an assimilated species has "
                     f"unassimilated obs: {bad_a}; where a monitor-only species has assimilated obs: {bad_b}")
        # (3) zero-obs ranges per species, for full days of the arm (data gaps / the 2026 MODIS monitoring gap)
        for si, sp in enumerate(SPECIES):
            z = (tot[:, si] == 0) & full
            if not z.any():
                continue
            idx = np.flatnonzero(z)
            p0 = 0                                                   # position in idx where the current run starts
            for k in range(1, len(idx) + 1):
                if k == len(idx) or idx[k] != idx[k - 1] + 1:
                    run = idx[p0:k]
                    gaps.append((a, sp, int(date_i[run[0]]), int(date_i[run[-1]]), len(run),
                                 bool(exp[run, si].any())))
                    p0 = k
    with open(OUT / "zero_obs_days.csv", "w") as f:
        f.write("arm,species,first_day,last_day,n_days,species_assimilated_in_arm_on_these_days\n")
        for g in gaps:
            f.write(",".join(str(x) for x in g) + "\n")
    lines += ["", "## Zero-obs periods (full days with no obs of the species)", "",
              "| arm | species | first | last | days | assimilated there |", "|---|---|---|---|---:|---|"]
    for g in gaps:
        if g[4] >= 1:
            lines.append(f"| {g[0]} | {g[1]} | {g[2]} | {g[3]} | {g[4]} | {'yes' if g[5] else 'no'} |")
    (FIG / "obs_availability_checks.md").write_text("\n".join(lines) + "\n")

    # ---- Figure 1: monthly mean tiles/day per species, one panel per arm
    fig, axs = plt.subplots(len(R), 1, figsize=(9, 1.9 * len(R) + 0.8), sharex=True, sharey=True)
    axs = np.atleast_1d(axs)
    for ax, (a, r) in zip(axs, R.items()):
        for si, sp in enumerate(SPECIES):
            x = r["n_tiles"][:, si, :].sum(axis=1).astype(float)
            ok = full_days(r)
            if not ok.any():
                continue
            t, y = monthly_mean(r["date"], x, ok)
            ax.plot(t, y / 1000, color=COLOR[sp], ls=STYLE[sp], lw=1.6, label=sp)
        ax.axvline(SWITCH, color=INK2, lw=0.8, ls=":")
        ax.set_ylabel("tiles/day (k)", fontsize=8, color=INK2)
        ax.set_title(f"{a}  ({CFG['arms'][a]['role']})", fontsize=8.5, loc="left", color=INK)
        style(ax)
    axs[0].legend(ncol=4, fontsize=8, frameon=False, loc="upper left", bbox_to_anchor=(0, 1.55))
    axs[-1].annotate("dotted line: 2023-08-01 switch", (SWITCH, 0), xytext=(4, 4), textcoords="offset points",
                     fontsize=7, color=INK2)
    fig.suptitle("SCF obs per day (tiles with an obs), monthly mean; MODIS solid, VIIRS dashed", fontsize=9,
                 color=INK, x=0.01, ha="left", y=0.995)
    fig.tight_layout()
    fig.savefig(FIG / "figA_1_tiles_per_day_by_species.png", dpi=150)
    plt.close(fig)

    # ---- Figure 2: overlap fractions. MODIS from the continuous T+A chain, VIIRS from TA_to_J1J2
    def chain_series(key, num, den):
        out_t, out_y = [], []
        for arm, beg, end in CFG["continuous_series"]["MODIS_TA"]["parts"]:
            if arm not in R:
                continue
            r = R[arm]
            m = full_days(r) & (r["date"] >= np.datetime64(beg)) & (r["date"] <= np.datetime64(end))
            x = r[num].astype(float) / np.maximum(r[den].astype(float), 1)
            t, y = monthly_mean(r["date"], x, m)
            out_t.append(t); out_y.append(y)
        return np.concatenate(out_t), np.concatenate(out_y)

    fig, ax = plt.subplots(figsize=(9, 3.4))
    if "TA_to_J1J2" in R and "TA_control" in R:
        t, y = chain_series("modis", "tiles_both_terra_aqua", "tiles_any_modis")
        ax.plot(t, y, color=COLOR["Terra"], lw=1.8, label="Terra and Aqua / any MODIS (T+A chain)")
        r = R["TA_to_J1J2"]
        m = full_days(r)
        x = r["tiles_both_j1_j2"].astype(float) / np.maximum(r["tiles_any_viirs"], 1)
        t2, y2 = monthly_mean(r["date"], x, m & (r["tiles_any_viirs"] > 0) & (r["date"] >= np.datetime64("2023-03-01")))
        ax.plot(t2, y2, color=COLOR["J1"], lw=1.8, ls="--", label="J1 and J2 / any VIIRS (TA_to_J1J2; J2 from 2023-02)")
    ax.axvline(SWITCH, color=INK2, lw=0.8, ls=":")
    ax.set_ylim(0, 1)
    ax.set_ylabel("share of tiles with an obs", fontsize=8, color=INK2)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    ax.set_title("Same-day overlap of the two sensors of a constellation, monthly mean", fontsize=9, loc="left", color=INK)
    style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "figA_2_constellation_overlap.png", dpi=150)
    plt.close(fig)

    # ---- Figure 3: tiles/day by latitude band, SY4-SY5 (full post-switch years; SY6 has the 2026 MODIS gap in 3 arms)
    fig, ax = plt.subplots(figsize=(9, 3.6))
    sy_use = [s for s in SY if s["id"] in ("SY4", "SY5")]
    bands = R[next(iter(R))]["bands"]
    src = {"Terra": "TA_control", "Aqua": "TA_control", "J1": "TA_to_J1J2", "J2": "TA_to_J1J2"}
    w = 0.2
    for si, sp in enumerate(SPECIES):
        arm = src[sp]
        if arm not in R:
            continue
        r = R[arm]
        m = np.zeros(len(r["date"]), bool)
        for s in sy_use:
            m |= sy_mask(r, s)
        vals = [r["n_tiles"][m, si, bi].mean() / 1000 for bi in range(len(bands))]
        k = ["Terra", "Aqua", "J1", "J2"].index(sp)
        ax.bar(np.arange(len(bands)) + (k - 1.5) * w, vals, width=w - 0.02, color=COLOR[sp], label=f"{sp} ({arm})")
    ax.set_xticks(range(len(bands)))
    ax.set_xticklabels([b.replace("_", " ") for b in bands], fontsize=8)
    ax.set_ylabel("tiles/day (k)", fontsize=8, color=INK2)
    ax.legend(fontsize=8, frameon=False, ncol=2)
    ax.set_title("Tiles with an obs per day by latitude band, mean of SY4-SY5", fontsize=9, loc="left", color=INK)
    style(ax)
    fig.tight_layout()
    fig.savefig(FIG / "figA_3_tiles_by_latitude_band.png", dpi=150)
    plt.close(fig)
    print("wrote tables to", OUT, "and figures + markdown to", FIG)


if __name__ == "__main__":
    main()
