#!/usr/bin/env python3
"""Build one tile-level covariate table (key tile_axis_index, 909 AZ M36 tiles) from LOCAL files only.

Run:  /Users/amfox/mamba/envs/regrid/bin/python -I projects/CYGNSS_L1_AZ/scripts/benefit_drivers/build_covariate_table.py
Writes output/benefit_drivers/covariate_table.csv.  See covariate_table_dictionary.md for column meanings.

Column prefixes: geo_/stat_ static, sup_ obs supply, op_ operator skill proxy (O vs H(x) in the OL),
sn_ signal/noise, y_ outcomes (DA-vs-OL change; NOT covariates, kept for the later regression).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset

PROJ = Path(__file__).resolve().parents[2]
OUT = PROJ / "output/benefit_drivers/covariate_table.csv"
STATS = PROJ / "output/cygl1_fixedop_2020_2022_stats"
sys.path.insert(0, str(PROJ.parent / "iv_tc"))
from iv_tc.readers import M36_NX, read_tilecoord  # noqa: E402

SEAS = {12: "djf", 1: "djf", 2: "djf", 3: "mam", 4: "mam", 5: "mam", 6: "jja", 7: "jja", 8: "jja",
        9: "son", 10: "son", 11: "son"}
COH_LO, COH_HI = 0.40, 2.16
MIN_N = 30  # min obs for a rank correlation


def base():
    t = pd.read_csv(STATS / "tile_coords.csv")
    tc = read_tilecoord(PROJ / "OLv8_M36_all_sensors_AZ_describe/OLv8_M36_all_sensors_AZ.ldas_tilecoord.bin")
    row = {int(x): k for k, x in enumerate(tc["tile_id"])}
    r = np.array([row[int(x)] for x in t.tile_id])
    t["idx0"] = tc["j_indg"][r].astype(int) * M36_NX + tc["i_indg"][r].astype(int)
    t["geo_elev_tilecoord_m"] = tc["elev"][r]
    t["geo_area_km2"] = tc["area"][r]
    t["geo_frac_cell"] = tc["frac_cell"][r]
    c = pd.read_csv(PROJ / "output/ismn_fixedop_20200101_20221231/tile_l1_obs_count_2020_2022.csv")
    t = t.merge(c.rename(columns={"tile_index": "tile_axis_index", "n_l1_obs": "sup_n_l1_obs_ismncount"}),
                on="tile_axis_index", how="left")
    t["in_l1_mask"] = t.sup_n_l1_obs_ismncount >= 100
    return t


def srho(x, y):
    if len(x) < MIN_N:
        return np.nan
    return pd.Series(x).rank().corr(pd.Series(y).rank())


def per_obs(t):
    """OL per-obs L1 table (Jan 2020-Dec 2022): supply, QC, operator-skill proxy, static fields."""
    q = pd.read_csv(PROJ / "output/cygnss_da_quality_bundle/per_obs/cygnss_ol_innovations_with_quality_all.csv")
    q["tile_axis_index"] = q.tilenum.astype(int) - 1  # verified: nearest-tile match 99.99%
    d = pd.to_datetime(q.ddm_time_utc, format="ISO8601", utc=True)
    q["year"] = d.dt.year
    q["seas"] = d.dt.month.map(SEAS)
    q["coh_pass"] = q.coherency_ratio.between(COH_LO, COH_HI)
    q["inn"] = q.obs - q.fcst
    g = q.groupby("tile_axis_index")
    o = pd.DataFrame(index=g.size().index)
    o["sup_n_obs_ol"] = g.size()
    for y in (2020, 2021, 2022):
        o[f"sup_n_obs_{y}"] = q[q.year == y].groupby("tile_axis_index").size()
    for s in ("djf", "mam", "jja", "son"):
        o[f"sup_n_obs_{s}"] = q[q.seas == s].groupby("tile_axis_index").size()
    o["sup_frac_coh_pass"] = g.coh_pass.mean()
    o["sup_frac_quality_flags_clean"] = g.quality_flags_2_clean.mean()
    o["sup_coh_ratio_mean"] = g.coherency_ratio.mean()
    o["sup_inc_angle_mean_deg"] = g.sp_inc_angle.mean()
    o["sup_inc_angle_std_deg"] = g.sp_inc_angle.std()
    o["sup_ddm_snr_mean_db"] = g.ddm_snr.mean()
    o["sup_rx_gain_mean_db"] = g.sp_rx_gain.mean()
    o["sup_quality_score_mean"] = g.quality_score.mean()
    for sc in range(1, 9):
        o[f"sup_frac_cyg0{sc}"] = (q.sc_num == sc).groupby(q.tile_axis_index).mean()
    o["sup_sc_entropy"] = o[[f"sup_frac_cyg0{s}" for s in range(1, 9)]].apply(
        lambda r: -(r[r > 0] * np.log(r[r > 0])).sum(), axis=1)
    # operator skill proxy: Spearman(O, H(x)) in the OL (raw dB, unscaled), whole period and per season
    o["op_spearman_OH"] = pd.Series({k: srho(v.obs.values, v.fcst.values) for k, v in g})
    for s in ("djf", "mam", "jja", "son"):
        o[f"op_spearman_OH_{s}"] = pd.Series({k: srho(v.obs.values, v.fcst.values)
                                              for k, v in q[q.seas == s].groupby("tile_axis_index")})
    o["op_spearman_OH_cohpass"] = pd.Series({k: srho(v.obs.values, v.fcst.values)
                                             for k, v in q[q.coh_pass].groupby("tile_axis_index")})
    # signal/noise (OL, raw dB)
    o["sn_H_mean_db"] = g.fcst.mean()
    o["sn_H_std_db"] = g.fcst.std()
    o["sn_H_range_p5_p95_db"] = g.fcst.quantile(0.95) - g.fcst.quantile(0.05)
    o["sn_O_std_db"] = g.obs.std()
    o["sn_O_over_H_std"] = o.sn_O_std_db / o.sn_H_std_db
    o["sn_innov_mean_db"] = g.inn.mean()
    o["sn_innov_std_db"] = g.inn.std()
    o["sn_spread_H_mean_db"] = np.sqrt(g.fcstvar.mean())
    o["sn_innov_over_spread"] = o.sn_innov_std_db / np.sqrt(g.fcstvar.mean() + g.obsvar.mean())
    o["sn_innov_over_Hspread"] = o.sn_innov_std_db / o.sn_spread_H_mean_db
    o["sn_K_mean"] = g.K.mean()
    o["sn_K_p90"] = g.K.quantile(0.9)
    o["sn_obsvar_mean_db2"] = g.obsvar.mean()
    o["sn_implied_incr_rms_db"] = np.sqrt((q.K ** 2 * q.inn ** 2).groupby(q.tile_axis_index).mean())
    for s in ("djf", "mam", "jja", "son"):
        gs = q[q.seas == s].groupby("tile_axis_index")
        o[f"sn_H_std_db_{s}"] = gs.fcst.std()
        o[f"sn_spread_H_db_{s}"] = np.sqrt(gs.fcstvar.mean())
        o[f"sn_K_mean_{s}"] = gs.K.mean()
    # static fields carried on the obs (tile constant or averaged)
    o["stat_dem_obs_mean_m"] = g.srtm_dem_alt.mean()
    o["stat_dem_obs_std_m"] = g.srtm_dem_alt.std()
    o["stat_slope_obs_mean_deg"] = g.srtm_slope.mean()
    o["stat_slope_obs_p90_deg"] = g.srtm_slope.quantile(0.9)
    o["stat_modis_landcover_mode"] = g.modis_land_cover.agg(lambda x: x.mode().iloc[0] if x.notna().any() else np.nan)
    o["stat_mwrtm_vegcls"] = g.MWRTM_VEGCLS.first()
    o["stat_mwrtm_soilcls"] = g.MWRTM_SOILCLS.first()
    o["stat_sand_frac"] = g.MWRTM_SAND.first()
    o["stat_clay_frac"] = g.MWRTM_CLAY.first()
    o["stat_porosity"] = g.MWRTM_POROS.first()
    o["stat_opacity_mean"] = g.coefficient_weighted_opacity.mean()
    o["stat_veg_geometry_factor_mean"] = g.vegetation_geometry_factor.mean()
    o["stat_attenuation_mean"] = g.coefficient_weighted_attenuation.mean()
    o["stat_brcs_crop_sum_median"] = g.brcs_crop_sum.median()
    return o.reset_index()


def paired_scaled(t):
    """2020 OL scaled L1 and L3 obs vs F (tile-resolved paired OFA csv): scaled-space skill."""
    d = pd.read_csv(PROJ / "output/cygnss_l1_l3_paired_ofa_2020.csv.gz", comment="#")
    ax = dict(zip(t.tile_id, t.tile_axis_index))
    d["tile_axis_index"] = d.tile_id.map(ax)
    d["inn"] = d.obs - d.fcst
    o = pd.DataFrame(index=sorted(d.tile_axis_index.dropna().astype(int).unique()))
    for sp, nm in ((13, "l1"), (12, "l3")):
        s = d[d.species == sp]
        g = s.groupby("tile_axis_index")
        o[f"sup_n_{nm}_scaled_2020"] = g.size()
        o[f"sn_{nm}_scaled_O_std_2020"] = g.obs.std()
        o[f"sn_{nm}_scaled_F_std_2020"] = g.fcst.std()
        o[f"sn_{nm}_scaled_innov_std_2020"] = g.inn.std()
        o[f"op_{nm}_scaled_pearson_OF_2020"] = g.apply(
            lambda v: np.corrcoef(v.obs, v.fcst)[0, 1] if len(v) >= MIN_N else np.nan)
        o[f"sn_{nm}_scaled_F_mean_2020"] = g.fcst.mean()
    o.index.name = "tile_axis_index"
    return o.reset_index()


def nc(path):
    with Dataset(path) as ds:
        out = {}
        for k, v in ds.variables.items():
            a = np.ma.filled(v[:], np.nan).astype(float)
            a[a <= -9998] = np.nan
            out[k] = a
    return out


GROUPS = {"smos": [0, 1, 2, 3], "smap": [4, 5, 6, 7], "ascat": [8, 9, 10], "l3": [11], "l1": [12]}


def stats_cols():
    sd = STATS / "stats"
    ts = "_20200101_20221231.nc4"
    prim = "full_xc015_coh040216_err39_fixedop"
    ol = nc(sd / f"temporal_stats_OL_paired_monitor_xmask_{prim}{ts}")
    da = nc(sd / f"temporal_stats_DA_paired_{prim}{ts}")
    o = {}
    # OL signal/noise in scaled obs space (xmask: same obs population as the primary arm)
    for nm, sp in (("l1", 12), ("l3", 11)):
        o[f"sn_ol_{nm}_F_mean_scaled"] = ol["F_mean"][:, sp]
        o[f"sn_ol_{nm}_F_std_scaled"] = ol["F_stdv"][:, sp]
        o[f"sn_ol_{nm}_O_std_scaled"] = ol["O_stdv"][:, sp]
        o[f"sn_ol_{nm}_OmF_std_scaled"] = ol["OmF_stdv"][:, sp]
        o[f"sn_ol_{nm}_OmF_norm_std"] = ol["OmF_norm_stdv"][:, sp]
        o[f"sup_ol_{nm}_N_data_{prim[:8]}"] = ol["N_data"][:, sp]
    o["sn_ol_smap_Tb_F_std_K"] = np.nanmean(ol["F_stdv"][:, 4:8], axis=1)
    o["sn_ol_smap_Tb_OmF_std_K"] = np.nanmean(ol["OmF_stdv"][:, 4:8], axis=1)
    o["sn_ol_smap_Tb_OmF_norm_std"] = np.nanmean(ol["OmF_norm_stdv"][:, 4:8], axis=1)
    o["sn_ol_smap_Tb_N"] = np.nansum(ol["N_data"][:, 4:8], axis=1)
    o["sn_ol_smos_Tb_OmF_std_K"] = np.nanmean(ol["OmF_stdv"][:, 0:4], axis=1)
    # L1 own-obs response in the primary DA arm
    o["sn_da_l1_OmA_over_OmF_std"] = da["OmA_stdv"][:, 12] / da["OmF_stdv"][:, 12]
    o["sn_da_l1_F_std_change_pct"] = 100 * (da["F_stdv"][:, 12] / ol["F_stdv"][:, 12] - 1)
    # outcomes: % change of O-F stdv, DA vs same-tag OL cross-mask, per species group
    def pct(daf, olf, sps):
        v = 100 * (daf["OmF_stdv"][:, sps] / olf["OmF_stdv"][:, sps] - 1)
        v[(daf["N_data"][:, sps] < 20) | (olf["N_data"][:, sps] < 20)] = np.nan
        return np.nanmean(v, axis=1) if v.shape[1] > 1 else v[:, 0]
    import warnings; warnings.simplefilter("ignore")
    for arm, lab in ((prim, "l1err39"), ("l3_fixedop", "l3only"), ("smap_fixedop", "smaponly")):
        a_da = nc(sd / f"temporal_stats_DA_paired_{arm}{ts}")
        a_ol = nc(sd / f"temporal_stats_OL_paired_monitor_xmask_{arm}{ts}")
        for gname, sps in GROUPS.items():
            o[f"y_{lab}_pct_omf_std_{gname}"] = pct(a_da, a_ol, sps)
    for yr in (2020, 2021, 2022):
        pdir = sd / "per_year"
        f = f"_{yr}0101_{yr}1231.nc4"
        a_da = nc(pdir / f"temporal_stats_DA_paired_{prim}{f}")
        a_ol = nc(pdir / f"temporal_stats_OL_paired_monitor_xmask_{prim}{f}")
        for gname in ("smap", "smos", "l3", "l1"):
            o[f"y_l1err39_{yr}_pct_omf_std_{gname}"] = pct(a_da, a_ol, GROUPS[gname])
    return pd.DataFrame(o).rename_axis("tile_axis_index").reset_index()


def main():
    t = base()
    for part in (per_obs(t), paired_scaled(t), stats_cols()):
        t = t.merge(part, on="tile_axis_index", how="left")
    # tile-level obs count used for mask definition cross-check
    t["y_l1err39_pct_omf_std_smap_minus_l3only"] = t.y_l1err39_pct_omf_std_smap - t.y_l3only_pct_omf_std_smap
    t.to_csv(OUT, index=False)
    print(t.shape, "->", OUT)
    print("in_l1_mask:", int(t.in_l1_mask.sum()), "; cols with all-NaN:", [c for c in t if t[c].isna().all()])


if __name__ == "__main__":
    main()
