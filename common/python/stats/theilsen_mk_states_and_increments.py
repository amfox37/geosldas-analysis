#!/usr/bin/env python3
import sys, argparse, time, logging
import numpy as np
import dask, xarray as xr

# ---------- Dask / Xarray options ----------
dask.config.set(scheduler="threads")
xr.set_options(file_cache_maxsize=1)

try:
    from dask.diagnostics import ProgressBar
    print("ProgressBar OK")
    HAVE_PBAR = True
except Exception as e:
    HAVE_PBAR = False
    print("No ProgressBar:", e)

# -----------------------------
# Config
# -----------------------------
CLIM_START = "2001-01-01"
CLIM_END   = "2020-12-31"

TEMP_THRESH_K     = 275.15   # 2 °C
SNOW_EPS          = 1e-2     # 1% snow cover
MIN_VALID_FRAC    = 0.25     # require >=25% valid months per gridcell (after masking)
SNOW_CLIM_THRESH  = 0.5      # kg m^-2; month is "snowy" if mean SWE > this
MIN_N_DEFAULT     = 10       # min valid months for general trend calc
MIN_N_SNOW        = 36       # stricter: min valid months for snow trend calc (~3 winters)

# For speed: single chunk along time + big tile chunks (warnings OK)
CHUNKS = {"time": -1, "tile": 50_000}

# Optional robust methods
try:
    from scipy.stats import theilslopes, kendalltau
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

try:
    import pymannkendall as mk
    HAVE_PMK = True
except Exception:
    HAVE_PMK = False

# -----------------------------
# Logging helpers
# -----------------------------
def stamp(msg: str):
    logging.info(f"[{time.strftime('%H:%M:%S')}] {msg}")

def maybe_compute(*delayed, progress=False):
    if progress and HAVE_PBAR:
        with ProgressBar():
            return dask.compute(*delayed)
    return dask.compute(*delayed)

def compute_with_progress(objs, progress=False, label="compute"):
    stamp(f"Starting dask {label} (progress={progress})")
    if progress and HAVE_PBAR:
        with ProgressBar():
            out = dask.compute(*objs)
    else:
        out = dask.compute(*objs)
    stamp(f"Finished dask {label}")
    return out

# -----------------------------
# Helpers
# -----------------------------
def deseasonalize(da: xr.DataArray, clim_start=CLIM_START, clim_end=CLIM_END):
    clim = da.sel(time=slice(clim_start, clim_end)).groupby("time.month").mean("time", skipna=True)
    return da.groupby("time.month") - clim

def apply_frozen_snow_mask(sm_da: xr.DataArray, tsoil: xr.DataArray, frsnow: xr.DataArray) -> xr.DataArray:
    # time-varying mask for soil moisture only
    mask = (tsoil < TEMP_THRESH_K) | (frsnow > SNOW_EPS)
    return sm_da.where(~mask)

def mask_sparse(da: xr.DataArray, min_frac=MIN_VALID_FRAC) -> xr.DataArray:
    n = da.sizes["time"]
    valid = da.notnull().sum("time") >= int(np.ceil(min_frac * n))
    return da.where(valid)

# --- Snow: data-driven selection of locally snowy months (both hemispheres, mountains) ---
def select_snowy_months(da_snow: xr.DataArray,
                        clim_start=CLIM_START,
                        clim_end=CLIM_END,
                        thresh=SNOW_CLIM_THRESH) -> xr.DataArray:
    """
    Keep only months that are climatologically snowy for each tile.
    A month is "snowy" at a tile if its SWE climatology > thresh (kg m^-2) for that month.
    """
    clim = (da_snow
            .sel(time=slice(clim_start, clim_end))
            .groupby("time.month")
            .mean("time", skipna=True))     # dims: (month, tile)
    snowy = clim > thresh                   # (month, tile) boolean
    keep = snowy.sel(month=da_snow.time.dt.month)  # (time, tile)
    return da_snow.where(keep)

# --- Theil–Sen + MK with configurable min_n ---
from functools import partial

def _theil_sen_mk_1d_with_min(y: np.ndarray, t: np.ndarray, min_n: int):
    m = np.isfinite(y)
    if m.sum() < min_n:
        return np.nan, np.nan

    t_ns0 = t[m][0].astype("datetime64[ns]").astype("int64")
    x_years = (t[m].astype("datetime64[ns]").astype("int64") - t_ns0) / 1e9 / (365.25 * 24 * 3600)
    yv = y[m].astype(float)

    if HAVE_SCIPY:
        slope_per_year, _, _, _ = theilslopes(yv, x=x_years)
    else:
        x = x_years
        xmean = x.mean(); ymean = yv.mean()
        denom = np.sum((x - xmean) ** 2)
        slope_per_year = np.nan if denom == 0 else np.sum((x - xmean) * (yv - ymean)) / denom

    slope_per_decade = slope_per_year * 10.0

    if HAVE_PMK:
        try:
            res = mk.original_test(yv)
            pval = res.p
        except Exception:
            pval = np.nan
    elif HAVE_SCIPY:
        _, pval = kendalltau(x_years, yv)
    else:
        pval = np.nan

    return slope_per_decade, pval

def theil_sen_mk(da: xr.DataArray, min_n: int = MIN_N_DEFAULT) -> tuple[xr.DataArray, xr.DataArray]:
    """
    Vectorized Theil–Sen slope (per decade) + MK p-value along 'time'.
    min_n sets the minimum number of valid months required per tile.
    """
    if hasattr(da, "chunks") and da.chunks is not None and "time" in da.dims:
        da = da.chunk({"time": -1})
    tvals = np.asarray(da["time"].values)

    f = partial(_theil_sen_mk_1d_with_min, min_n=min_n)
    slope, p = xr.apply_ufunc(
        f,
        da,
        xr.DataArray(tvals, dims=("time",)),
        input_core_dims=[["time"], ["time"]],
        output_core_dims=[[], []],
        dask="parallelized",
        vectorize=True,
        output_dtypes=[np.float64, np.float64],
        dask_gufunc_kwargs={"allow_rechunk": False},
    )
    return slope, p

def align_all(join="inner", *args):
    return xr.align(*args, join=join)

# -----------------------------
# Staged write helper
# -----------------------------
def write_batch(out_nc: str, mode: str, data_vars: dict, coords: dict, attrs: dict,
                encoding_template: dict, progress=False):
    ds = xr.Dataset(data_vars=data_vars, coords=coords, attrs=attrs)

    # Make enc per var
    enc = {name: dict(encoding_template) for name in data_vars}

    # Persist this batch so writes are fast and one-time
    if progress and HAVE_PBAR:
        with ProgressBar():
            ds = ds.persist()
    else:
        ds = ds.persist()

    delayed = ds.to_netcdf(out_nc, mode=mode, encoding=enc, compute=False)

    if progress and HAVE_PBAR:
        with ProgressBar():
            dask.compute(delayed)
    else:
        dask.compute(delayed)

    try:
        ds.close()
    except Exception:
        pass

# -----------------------------
# Pipeline
# -----------------------------
def run(
    ol_states_nc: str,
    da_states_nc: str,
    incr_nc: str,
    out_nc: str,
    use_anomalies: bool = True,
    mask_sm_increments: bool = True,   # SM increments masked like DA SM
    progress: bool = False,
    min_valid_frac: float = MIN_VALID_FRAC,
    snow_clim_thresh: float = SNOW_CLIM_THRESH,
    min_n_snow: int = MIN_N_SNOW,
):
    stamp("Opening datasets")
    ds_ol  = xr.open_dataset(ol_states_nc, chunks=CHUNKS)
    ds_da  = xr.open_dataset(da_states_nc, chunks=CHUNKS)
    ds_inc = xr.open_dataset(incr_nc,      chunks=CHUNKS)

    for v in ["SFMC","RZMC","TSOIL1","FRLANDSNO","PRECTOTCORRLAND","SNOMASLAND"]:
        if v not in ds_ol or v not in ds_da:
            raise KeyError(f"Missing {v} in state files.")
    for v in ["SFMC_INC","RZMC_INC","SNOWMASS_INCR"]:
        if v not in ds_inc:
            raise KeyError(f"Missing {v} in increments file.")

    stamp("Aligning state variables")
    SFMC_ol, SFMC_da = align_all("inner", ds_ol["SFMC"], ds_da["SFMC"])
    RZMC_ol, RZMC_da = align_all("inner", ds_ol["RZMC"], ds_da["RZMC"])
    T_ol, T_da       = align_all("inner", ds_ol["TSOIL1"], ds_da["TSOIL1"])
    F_ol, F_da       = align_all("inner", ds_ol["FRLANDSNO"], ds_da["FRLANDSNO"])
    P_ol, P_da       = align_all("inner", ds_ol["PRECTOTCORRLAND"], ds_da["PRECTOTCORRLAND"])
    SNO_ol, SNO_da   = align_all("inner", ds_ol["SNOMASLAND"],     ds_da["SNOMASLAND"])

    SFMC_ol, SFMC_da, RZMC_ol, RZMC_da, T_ol, T_da, F_ol, F_da, P_ol, P_da, SNO_ol, SNO_da = xr.align(
        SFMC_ol, SFMC_da, RZMC_ol, RZMC_da, T_ol, T_da, F_ol, F_da, P_ol, P_da, SNO_ol, SNO_da, join="inner"
    )

    stamp("Applying frozen/snow mask to SM")
    SFMC_ol = apply_frozen_snow_mask(SFMC_ol, T_ol, F_ol)
    SFMC_da = apply_frozen_snow_mask(SFMC_da, T_da, F_da)
    RZMC_ol = apply_frozen_snow_mask(RZMC_ol,  T_ol, F_ol)
    RZMC_da = apply_frozen_snow_mask(RZMC_da,  T_da, F_da)

    stamp(f"Dropping sparse SM tiles (min_valid_frac={min_valid_frac})")
    SFMC_ol = mask_sparse(SFMC_ol, min_frac=min_valid_frac); SFMC_da = mask_sparse(SFMC_da, min_frac=min_valid_frac)
    RZMC_ol = mask_sparse(RZMC_ol, min_frac=min_valid_frac); RZMC_da = mask_sparse(RZMC_da, min_frac=min_valid_frac)

    stamp("Aligning increments")
    SFMC_INC, RZMC_INC, SNO_INC = align_all("inner", ds_inc["SFMC_INC"], ds_inc["RZMC_INC"], ds_inc["SNOWMASS_INCR"])
    SFMC_ol, SFMC_da, RZMC_ol, RZMC_da, P_ol, P_da, SNO_ol, SNO_da, SFMC_INC, RZMC_INC, SNO_INC = xr.align(
        SFMC_ol, SFMC_da, RZMC_ol, RZMC_da, P_ol, P_da, SNO_ol, SNO_da, SFMC_INC, RZMC_INC, SNO_INC, join="inner"
    )
    if mask_sm_increments:
        stamp("Masking SM increments where DA SM is masked")
        da_mask_sf  = SFMC_da.isnull()
        da_mask_rz  = RZMC_da.isnull()
        SFMC_INC = SFMC_INC.where(~da_mask_sf)
        RZMC_INC = RZMC_INC.where(~da_mask_rz)

    # --- Snow: retain locally snowy months (data-driven, both hemispheres) ---
    stamp(f"Selecting snowy months per tile (SWE climatology > {snow_clim_thresh} kg m^-2)")
    SNO_ol = select_snowy_months(SNO_ol, thresh=snow_clim_thresh)
    SNO_da = select_snowy_months(SNO_da, thresh=snow_clim_thresh)

    stamp(f"Deseasonalizing (anomalies={use_anomalies})")
    if use_anomalies:
        SFMC_ol_a = deseasonalize(SFMC_ol); SFMC_da_a = deseasonalize(SFMC_da)
        RZMC_ol_a = deseasonalize(RZMC_ol); RZMC_da_a = deseasonalize(RZMC_da)
        P_ol_a    = deseasonalize(P_ol);    P_da_a    = deseasonalize(P_da)
        SNO_ol_a  = deseasonalize(SNO_ol);  SNO_da_a  = deseasonalize(SNO_da)
    else:
        SFMC_ol_a, SFMC_da_a = SFMC_ol, SFMC_da
        RZMC_ol_a, RZMC_da_a = RZMC_ol, RZMC_da
        P_ol_a,    P_da_a    = P_ol,    P_da
        SNO_ol_a,  SNO_da_a  = SNO_ol,  SNO_da

    # ---------- Trends: SM ----------
    stamp("Computing Theil–Sen + MK for SM (CNTL)")
    sfmc_slope_cntl, sfmc_p_cntl = theil_sen_mk(SFMC_ol_a, min_n=MIN_N_DEFAULT)
    stamp("Computing Theil–Sen + MK for SM (DA)")
    sfmc_slope_da,   sfmc_p_da   = theil_sen_mk(SFMC_da_a, min_n=MIN_N_DEFAULT)
    stamp("Computing Theil–Sen + MK for RZ (CNTL)")
    rzmc_slope_cntl, rzmc_p_cntl = theil_sen_mk(RZMC_ol_a, min_n=MIN_N_DEFAULT)
    stamp("Computing Theil–Sen + MK for RZ (DA)")
    rzmc_slope_da,   rzmc_p_da   = theil_sen_mk(RZMC_da_a, min_n=MIN_N_DEFAULT)

    # ---------- Trends: precip & snow ----------
    stamp("Computing Theil–Sen + MK for PREC (CNTL/DA)")
    p_slope_cntl, p_p_cntl = theil_sen_mk(P_ol_a, min_n=MIN_N_DEFAULT)
    p_slope_da,   p_p_da   = theil_sen_mk(P_da_a, min_n=MIN_N_DEFAULT)

    stamp(f"Computing Theil–Sen + MK for SNOW (CNTL/DA), min_n={min_n_snow}")
    sno_slope_cntl, sno_p_cntl = theil_sen_mk(SNO_ol_a, min_n=min_n_snow)
    sno_slope_da,   sno_p_da   = theil_sen_mk(SNO_da_a, min_n=min_n_snow)

    # ---------- Cum increments ----------
    stamp("Building cumulative increments and trends")
    Cum_SFMC = SFMC_INC.fillna(0).cumsum("time")
    Cum_RZMC = RZMC_INC.fillna(0).cumsum("time")
    Cum_SNO  = SNO_INC.fillna(0).cumsum("time")

    cum_sfmc_slope, cum_sfmc_p = theil_sen_mk(Cum_SFMC, min_n=MIN_N_DEFAULT)
    cum_rzmc_slope, cum_rzmc_p = theil_sen_mk(Cum_RZMC, min_n=MIN_N_DEFAULT)
    cum_sno_slope,  cum_sno_p  = theil_sen_mk(Cum_SNO,  min_n=MIN_N_DEFAULT)

    mean_sfmc_inc = SFMC_INC.mean("time")
    mean_rzmc_inc = RZMC_INC.mean("time")
    mean_sno_inc  = SNO_INC.mean("time")

    # ---------- ONE BIG COMPUTE ----------
    (sfmc_slope_cntl, sfmc_p_cntl, sfmc_slope_da, sfmc_p_da,
     rzmc_slope_cntl, rzmc_p_cntl, rzmc_slope_da, rzmc_p_da,
     p_slope_cntl,    p_p_cntl,    p_slope_da,    p_p_da,
     sno_slope_cntl,  sno_p_cntl,  sno_slope_da,  sno_p_da,
     cum_sfmc_slope,  cum_sfmc_p,  cum_rzmc_slope, cum_rzmc_p,
     cum_sno_slope,   cum_sno_p,
     mean_sfmc_inc,   mean_rzmc_inc, mean_sno_inc) = compute_with_progress(
        [
            sfmc_slope_cntl, sfmc_p_cntl, sfmc_slope_da, sfmc_p_da,
            rzmc_slope_cntl, rzmc_p_cntl, rzmc_slope_da, rzmc_p_da,
            p_slope_cntl,    p_p_cntl,    p_slope_da,    p_p_da,
            sno_slope_cntl,  sno_p_cntl,  sno_slope_da,  sno_p_da,
            cum_sfmc_slope,  cum_sfmc_p,  cum_rzmc_slope, cum_rzmc_p,
            cum_sno_slope,   cum_sno_p,
            mean_sfmc_inc,   mean_rzmc_inc, mean_sno_inc,
        ],
        progress=progress,
        label="trend & increment metrics"
    )

    # ---------- Prepare metadata ----------
    lat = ds_da["lat"]; lon = ds_da["lon"]
    coords = dict(tile=np.arange(lat.size, dtype=np.int64),
                  lat=(("tile",), lat.data),
                  lon=(("tile",), lon.data))
    attrs = dict(
        note=("Theil–Sen slope (per decade) & MK p-value. SM masked where TSOIL1<275.15K or FRLANDSNO>1e-2; "
              "PREC unmasked; SNOW is restricted to locally snowy months (data-driven) and then deseasonalized. "
              "CumInc slopes/p for SFMC,RZMC (masked like DA SM) and SNOW (unmasked cumulative).")
    )
    comp = dict(zlib=True, complevel=4, chunksizes=(CHUNKS["tile"],))

    # ---------- STAGED WRITES ----------
    stamp(f"Writing batch 1/4 → {out_nc} (SM trends)")
    write_batch(
        out_nc, mode="w",
        data_vars=dict(
            SFMC_slope_CNTL = sfmc_slope_cntl,
            SFMC_p_CNTL     = sfmc_p_cntl,
            SFMC_slope_DA   = sfmc_slope_da,
            SFMC_p_DA       = sfmc_p_da,
            RZMC_slope_CNTL = rzmc_slope_cntl,
            RZMC_p_CNTL     = rzmc_p_cntl,
            RZMC_slope_DA   = rzmc_slope_da,
            RZMC_p_DA       = rzmc_p_da,
        ),
        coords=coords, attrs=attrs, encoding_template=comp, progress=progress
    )

    stamp("Writing batch 2/4 (PREC & SNOW trends)")
    write_batch(
        out_nc, mode="a",
        data_vars=dict(
            PREC_slope_CNTL = p_slope_cntl,    PREC_p_CNTL = p_p_cntl,
            PREC_slope_DA   = p_slope_da,      PREC_p_DA   = p_p_da,
            SNOW_slope_CNTL = sno_slope_cntl,  SNOW_p_CNTL = sno_p_cntl,
            SNOW_slope_DA   = sno_slope_da,    SNOW_p_DA   = sno_p_da,
        ),
        coords=coords, attrs=attrs, encoding_template=comp, progress=progress
    )

    stamp("Writing batch 3/4 (CumInc slopes & p)")
    write_batch(
        out_nc, mode="a",
        data_vars=dict(
            CumInc_SFMC_slope = cum_sfmc_slope,  CumInc_SFMC_p = cum_sfmc_p,
            CumInc_RZMC_slope = cum_rzmc_slope,  CumInc_RZMC_p = cum_rzmc_p,
            CumInc_SNOW_slope = cum_sno_slope,   CumInc_SNOW_p = cum_sno_p,
        ),
        coords=coords, attrs=attrs, encoding_template=comp, progress=progress
    )

    stamp("Writing batch 4/4 (mean increments)")
    write_batch(
        out_nc, mode="a",
        data_vars=dict(
            SFMC_increment_mean = mean_sfmc_inc,
            RZMC_increment_mean = mean_rzmc_inc,
            SNOW_increment_mean = mean_sno_inc,
        ),
        coords=coords, attrs=attrs, encoding_template=comp, progress=progress
    )

    # Close inputs to release file handles
    for _ds in (ds_ol, ds_da, ds_inc):
        try: _ds.close()
        except Exception: pass

    stamp("All done")
    return xr.open_dataset(out_nc)

# -----------------------------
# CLI
# -----------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Theil–Sen + MK trends for SM (masked), precip (unmasked) & snow (snowy-months filtered), "
                    "plus cumulative increments, with progress."
    )
    p.add_argument("--ol-states",
        default="/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/OLv8_land_variables_2000_2024_compressed.nc")
    p.add_argument("--da-states",
        default="/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/DAv8_land_variables_2000_2024_compressed.nc")
    p.add_argument("--increments",
        default="/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper/LS_DAv8_M36_v2/LS_DAv8_M36/output/SMAP_EASEv2_M36_GLOBAL/cat/LS_monthly_increments_2000_2024.nc")
    p.add_argument("--out-nc", default="LS_theilsen_mk_states_increments.nc")
    p.add_argument("--raw", action="store_true", help="Use raw monthly means instead of anomalies.")
    p.add_argument("--no-mask-sm-increments", action="store_true", help="Do not mask SM increments.")
    p.add_argument("--min-valid-frac", type=float, default=MIN_VALID_FRAC,
                   help="Minimum fraction of non-missing months required per tile for SM/RZ (default: 0.25).")
    p.add_argument("--snow-clim-thresh", type=float, default=SNOW_CLIM_THRESH,
                   help="SWE climatology threshold (kg m^-2) to define snowy months per tile (default: 0.5).")
    p.add_argument("--snow-min-months", type=int, default=MIN_N_SNOW,
                   help="Minimum valid months required for snow trend (default: 36).")
    p.add_argument("--progress", action="store_true", help="Show live progress bars for Dask computes.")
    p.add_argument("--log", default="INFO", help="Log level (DEBUG, INFO, WARNING, ERROR).")
    return p.parse_args(argv)

def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log.upper(), logging.INFO),
                        format="%(message)s")
    out = run(
        ol_states_nc=args.ol_states,
        da_states_nc=args.da_states,
        incr_nc=args.increments,
        out_nc=args.out_nc,
        use_anomalies=not args.raw,
        mask_sm_increments=not args.no_mask_sm_increments,
        progress=args.progress,
        min_valid_frac=args.min_valid_frac,
        snow_clim_thresh=args.snow_clim_thresh,
        min_n_snow=args.snow_min_months,
    )
    print(f"Saved: {args.out_nc}")
    print("Vars:", list(out.data_vars))
    print("Dims:", dict(out.dims))

if __name__ == "__main__":
    sys.exit(main())
