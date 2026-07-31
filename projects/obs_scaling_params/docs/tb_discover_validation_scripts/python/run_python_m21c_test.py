"""
Temporary parity-test driver: runs the UNMODIFIED
obs_scaling.tb_tile_stats.generate_tb_scaling_params against the
M21C_land_sweeper_OLv8_M36 (.bin) archive.

obs_scaling.io.read_obs_param cannot parse this archive's obsparam file: it
predates the fcstvarname/fcstunits fields the shipped reader (cd3eb01)
expects (see check_obsparam_schema.m / read_obsparam_newschema.m for the
MATLAB-side finding of the same schema evolution). generate_tb_scaling_params
itself does not call read_obs_param internally -- it only consumes an
obs_params list -- so this script builds that list directly from values
independently verified against the raw obsparam text and against MATLAB's
own (unmodified) read_obsparam.m output for this exact file:

  i=5 SMAP_L1C_Tbh_A species=5 orbit=1 pol=1 ang=40 varname=Tb units=K
  i=6 SMAP_L1C_Tbh_D species=6 orbit=2 pol=1 ang=40 varname=Tb units=K
  i=7 SMAP_L1C_Tbv_A species=7 orbit=1 pol=2 ang=40 varname=Tb units=K
  i=8 SMAP_L1C_Tbv_D species=8 orbit=2 pol=2 ang=40 varname=Tb units=K

This does not modify obs_scaling/io.py or any other implementation file.
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
    orbit=2,
    angles=[40.0],
    window_days=75,
    ndata_min=20,
    prefix="PYTHON_TEST_M21C_SMAP_Tb_",
    convert_grid="EASEv2_M36",
    obsfcstana_format="bin",
)
print(f"Wrote {len(written)} files")
for w in written:
    print(" ", w)
