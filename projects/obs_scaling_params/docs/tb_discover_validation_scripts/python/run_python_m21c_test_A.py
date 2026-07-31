"""
Same as run_python_m21c_test.py but ascending orbit (orbit=1), species 5,7.
"""
from __future__ import annotations
import sys
from pathlib import Path

REPO = Path("/discover/nobackup/projects/land_da/geosldas-analysis/projects/obs_scaling_params")
sys.path.insert(0, str(REPO))

import numpy as np
from obs_scaling.io import ObsParam
from obs_scaling.tb_tile_stats import generate_tb_scaling_params


def make_param(descr, species, orbit, pol):
    return ObsParam(
        descr=descr, species=species, orbit=orbit, pol=pol, n_ang=1,
        ang=np.array([40.0]), freq=1.41e9, fov=20.0, fov_units="km",
        assim="F", scale="F", getinnov="T", rtm_id=4,
        bias_npar=0.0, bias_trel=864000.0, bias_tcut=432000.0, nodata=-9999.0,
        varname="Tb", units="K", fcstvarname="Tb", fcstunits="K",
        path="", name="", maskpath="", maskname="",
        scalepath="", scalename="", flistpath="", flistname="",
        errstd=4.0, std_normal_max=2.5, zeromean="T", coarsen_pert="T",
        xcorr=0.25, ycorr=0.25, adapt=0.0,
    )


obs_params = [
    make_param("SMAP_L1C_Tbh_A", 5, 1, 1),
    make_param("SMAP_L1C_Tbh_D", 6, 2, 1),
    make_param("SMAP_L1C_Tbv_A", 7, 1, 2),
    make_param("SMAP_L1C_Tbv_D", 8, 2, 2),
]

written = generate_tb_scaling_params(
    run_months=[1, 2, 3],
    exp_path="/discover/nobackup/projects/land_da/obs_scaling_test",
    exp_run="LS_OLv8_M36",
    domain="SMAP_EASEv2_M36_GLOBAL",
    start_year=[2016, 2016, 2016],
    end_year=[2016, 2016, 2016],
    dt_assim=3 * 60 * 60,
    t0_assim=0,
    obs_params=obs_params,
    description_contains="SMAP_L1C",
    orbit=1,
    angles=[40.0],
    window_days=75,
    ndata_min=20,
    prefix="PYTHON_TEST_M21C_SMAP_Tb_ASC_",
    convert_grid="EASEv2_M36",
    obsfcstana_format="bin",
)
print(f"Wrote {len(written)} files")
for w in written:
    print(" ", w)
