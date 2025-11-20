#!/usr/bin/env python3
import sys, argparse
import numpy as np
import xarray as xr

CLIM_START = "2001-01-01"
CLIM_END   = "2020-12-31"
CHUNKS = {"time": 36, "tile": 100_000}

# thresholds
TEMP_THRESH_K = 275.15   # 2°C
SNOW_EPS      = 1e-2     # 1% snow cover

def deseasonalize(da: xr.DataArray, clim_start=CLIM_START, clim_end=CLIM_END) -> xr.DataArray:
    clim = da.sel(time=slice(clim_start, clim_end)).groupby("time.month").mean("time")
    return da.groupby("time.month") - clim

def trend_per_decade(da: xr.DataArray) -> xr.DataArray:
    coeffs = da.polyfit(dim="time", deg=1, skipna=True)
    slope_per_ns = coeffs.polyfit_coefficients.sel(degree=1)
    ns_per_decade = np.float64(10 * 365.25 * 24 * 3600 * 1e9)
    out = slope_per_ns * ns_per_decade
    out.name = (da.name or "var") + "_trend_per_decade"
    return out

def area_weights(lat: xr.DataArray) -> xr.DataArray:
    return np.cos(np.deg2rad(lat))

def apply_frozen_snow_mask(sm_da: xr.DataArray, tsoil: xr.DataArray, frsnow: xr.DataArray) -> xr.DataArray:
    """Mask SM where ground is frozen OR snow-covered, month-by-month."""
    mask = (tsoil < TEMP_THRESH_K) | (frsnow > SNOW_EPS)
    sm_masked = sm_da.where(~mask)
    return sm_masked

def align_common_time(*das):
    return xr.align(*das, join="inner")

def run_analysis(
    ol_states_nc: str,
    da_states_nc: str,
    incr_nc: str,
    out_nc: str,
    use_anomalies: bool = True,
    mask_increments: bool = True,   # apply same mask to SM increments
):
    ds_ol = xr.open_dataset(ol_states_nc, chunks=CHUNKS)
    ds_da = xr.open_dataset(da_states_nc, chunks=CHUNKS)
    ds_inc = xr.open_dataset(incr_nc, chunks=CHUNKS)

    # required vars
    for v in ["SFMC","RZMC","PRECTOTCORRLAND","SNOMASLAND","TSOIL1","FRLANDSNO"]:
        if v not in ds_ol or v not in ds_da:
            raise KeyError(f"Missing {v} in state files.")

    for v in ["SFMC_INC","RZMC_INC","SNOWMASS_INCR"]:
        if v not in ds_inc:
            raise KeyError(f"Missing {v} in increment file.")

    # Align states first (time,tile)
    SFMC_ol, SFMC_da = align_common_time(ds_ol["SFMC"], ds_da["SFMC"])
    RZMC_ol,  RZMC_da  = align_common_time(ds_ol["RZMC"], ds_da["RZMC"])
    P_ol, P_da         = align_common_time(ds_ol["PRECTOTCORRLAND"], ds_da["PRECTOTCORRLAND"])
    SNO_ol, SNO_da     = align_common_time(ds_ol["SNOMASLAND"],     ds_da["SNOMASLAND"])
    T_ol,  T_da        = align_common_time(ds_ol["TSOIL1"],         ds_da["TSOIL1"])
    FRS_ol, FRS_da     = align_common_time(ds_ol["FRLANDSNO"],      ds_da["FRLANDSNO"])

    # Make a common time/tile grid for states
    SFMC_ol, SFMC_da, RZMC_ol, RZMC_da, P_ol, P_da, SNO_ol, SNO_da, T_ol, T_da, FRS_ol, FRS_da = xr.align(
        SFMC_ol, SFMC_da, RZMC_ol, RZMC_da, P_ol, P_da, SNO_ol, SNO_da, T_ol, T_da, FRS_ol, FRS_da, join="inner"
    )

    # Apply frozen/snow mask to soil moisture (CNTL and DA)
    SFMC_ol_m = apply_frozen_snow_mask(SFMC_ol, T_ol, FRS_ol)
    SFMC_da_m = apply_frozen_snow_mask(SFMC_da, T_da, FRS_da)
    RZMC_ol_m = apply_frozen_snow_mask(RZMC_ol, T_ol, FRS_ol)
    RZMC_da_m = apply_frozen_snow_mask(RZMC_da, T_da, FRS_da)

    # Increments aligned and (optionally) masked using DA mask (or combined mask)
    SFMC_INC, RZMC_INC, SNO_INC = align_common_time(ds_inc["SFMC_INC"], ds_inc["RZMC_INC"], ds_inc["SNOWMASS_INCR"])
    # Bring increments onto same final (time,tile) intersection
    SFMC_ol_m, SFMC_da_m, RZMC_ol_m, RZMC_da_m, P_ol, P_da, SNO_ol, SNO_da, SFMC_INC, RZMC_INC, SNO_INC = xr.align(
        SFMC_ol_m, SFMC_da_m, RZMC_ol_m, RZMC_da_m, P_ol, P_da, SNO_ol, SNO_da, SFMC_INC, RZMC_INC, SNO_INC, join="inner"
    )

    if mask_increments:
        # mask increments where DA soil moisture is masked (frozen/snow)
        da_mask_sf  = SFMC_da_m.isnull()
        da_mask_rz  = RZMC_da_m.isnull()
        SFMC_INC = SFMC_INC.where(~da_mask_sf)
        RZMC_INC = RZMC_INC.where(~da_mask_rz)

    # Deseasonalize (recommended). If raw trends desired, pass use_anomalies=False
    if use_anomalies:
        SFMC_ol_a = deseasonalize(SFMC_ol_m); SFMC_da_a = deseasonalize(SFMC_da_m)
        RZMC_ol_a = deseasonalize(RZMC_ol_m); RZMC_da_a = deseasonalize(RZMC_da_m)
        P_ol_a    = deseasonalize(P_ol);      P_da_a    = deseasonalize(P_da)
        SNO_ol_a  = deseasonalize(SNO_ol);    SNO_da_a  = deseasonalize(SNO_da)
    else:
        SFMC_ol_a, SFMC_da_a = SFMC_ol_m, SFMC_da_m
        RZMC_ol_a, RZMC_da_a = RZMC_ol_m, RZMC_da_m
        P_ol_a,    P_da_a    = P_ol,    P_da
        SNO_ol_a,  SNO_da_a  = SNO_ol,  SNO_da

    # Trends (per decade)
    t_SFMC_ol = trend_per_decade(SFMC_ol_a).rename("SFMC_trend_CNTL")
    t_SFMC_da = trend_per_decade(SFMC_da_a).rename("SFMC_trend_DA")
    t_RZMC_ol = trend_per_decade(RZMC_ol_a).rename("RZMC_trend_CNTL")
    t_RZMC_da = trend_per_decade(RZMC_da_a).rename("RZMC_trend_DA")

    t_P_ol    = trend_per_decade(P_ol_a).rename("PREC_trend_CNTL")
    t_P_da    = trend_per_decade(P_da_a).rename("PREC_trend_DA")
    t_SNO_ol  = trend_per_decade(SNO_ol_a).rename("SNOW_trend_CNTL")
    t_SNO_da  = trend_per_decade(SNO_da_a).rename("SNOW_trend_DA")

    # ΔTrend and increments diagnostics
    d_SFMC = (t_SFMC_da - t_SFMC_ol).rename("SFMC_trend_delta")
    d_RZMC = (t_RZMC_da - t_RZMC_ol).rename("RZMC_trend_delta")
    d_SNO  = (t_SNO_da  - t_SNO_ol ).rename("SNOW_trend_delta")

    Cum_SFMC = SFMC_INC.cumsum("time")
    Cum_RZMC = RZMC_INC.cumsum("time")
    Cum_SNO  = SNO_INC.cumsum("time")

    t_Cum_SFMC = trend_per_decade(Cum_SFMC).rename("CumInc_SFMC_trend_per_decade")
    t_Cum_RZMC = trend_per_decade(Cum_RZMC).rename("CumInc_RZMC_trend_per_decade")
    t_Cum_SNO  = trend_per_decade(Cum_SNO ).rename("CumInc_SNOW_trend_per_decade")

    mean_SFMC_inc = SFMC_INC.mean("time").rename("SFMC_increment_mean")
    mean_RZMC_inc = RZMC_INC.mean("time").rename("RZMC_increment_mean")
    mean_SNO_inc  = SNO_INC.mean("time").rename("SNOW_increment_mean")

    lat = ds_da["lat"]; lon = ds_da["lon"]
    out = xr.Dataset(
        data_vars=dict(
            SFMC_trend_CNTL=t_SFMC_ol, SFMC_trend_DA=t_SFMC_da, SFMC_trend_delta=d_SFMC,
            RZMC_trend_CNTL=t_RZMC_ol, RZMC_trend_DA=t_RZMC_da, RZMC_trend_delta=d_RZMC,
            PREC_trend_CNTL=t_P_ol, PREC_trend_DA=t_P_da,
            SNOW_trend_CNTL=t_SNO_ol, SNOW_trend_DA=t_SNO_da, SNOW_trend_delta=d_SNO,
            CumInc_SFMC_trend=t_Cum_SFMC, CumInc_RZMC_trend=t_Cum_RZMC, CumInc_SNOW_trend=t_Cum_SNO,
            SFMC_increment_mean=mean_SFMC_inc, RZMC_increment_mean=mean_RZMC_inc, SNOW_increment_mean=mean_SNO_inc,
        ),
        coords=dict(tile=lat.tile, lat=(("tile",), lat.data), lon=(("tile",), lon.data)),
        attrs=dict(note=f"Trends with monthly SM masked where TSOIL1<{TEMP_THRESH_K}K or FRLANDSNO>{SNOW_EPS}")
    )

    comp = dict(zlib=True, complevel=4)
    enc = {v: {**comp, "chunksizes": (CHUNKS["tile"],)} for v in out.data_vars}
    out.to_netcdf(out_nc, encoding=enc)
    return out

def parse_args(argv=None):
    import argparse
    p = argparse.ArgumentParser(description="Trend analysis with frozen/snow masking.")
    p.add_argument("--ol-states", default="/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/OLv8_land_variables_2000_2024_compressed.nc")
    p.add_argument("--da-states", default="/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/DAv8_land_variables_2000_2024_compressed.nc")
    p.add_argument("--increments", default="/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/LS_monthly_increments_2000_2024.nc")
    p.add_argument("--out-nc", default="LS_trends_masked.nc")
    p.add_argument("--raw", action="store_true", help="Use raw monthly means (no deseasonalization).")
    p.add_argument("--no-mask-increments", action="store_true", help="Do not mask SM increments.")
    return p.parse_args(argv)

def main(argv=None):
    args = parse_args(argv)
    out = run_analysis(
        ol_states_nc=args.ol_states,
        da_states_nc=args.da_states,
        incr_nc=args.increments,
        out_nc=args.out_nc,
        use_anomalies=not args.raw,
        mask_increments=not args.no_mask_increments,
    )
    print(f"Saved: {args.out_nc}")
    print("Vars:", list(out.data_vars))

if __name__ == "__main__":
    sys.exit(main())