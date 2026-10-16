from pathlib import Path
import shutil

import conftest
import netCDF4 as nc
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]      # test data; scripts: see conftest.py

from obs_scaling_utils import (  # noqa: E402
    get_doy_pentad_range,
    get_ObsFcstAna_fname,
    get_year_bounds,
    read_obs_param,
    read_ObsFcstAna_nc4,
    read_tilecoord,
)


INPUTS = PROJECT_ROOT / "test_data" / "inputs"
OFA_NC4 = INPUTS / "LS_DAv8_M36_as_test2.ens_avg.ldas_ObsFcstAna.20200502_0900z.nc4"


def test_read_obs_param_handles_empty_quoted_fields_and_fcst_metadata():
    ls_params = read_obs_param(INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt")
    ol_params = read_obs_param(INPUTS / "OLv7_M36_MULTI_type_13_H121.ldas_obsparam.20150401_0000z.txt")

    assert len(ls_params) == 14
    assert [p["descr"] for p in ls_params[:3]] == ["SMOS_fit_Tbh_A", "SMOS_fit_Tbh_D", "SMOS_fit_Tbv_A"]
    assert ls_params[0]["species"] == 1
    assert ls_params[0]["orbit"] == 1
    np.testing.assert_array_equal(ls_params[0]["ang"], [40.0])
    assert ls_params[0]["FOV_units"] == "km"
    assert ls_params[-1]["descr"] == "CYGNSS_SM_6hr"
    assert ls_params[-1]["varname"] == "sfmc"
    assert ls_params[-1]["units"] == "m3 m-3"
    assert ls_params[-1]["fcstvarname"] == "sfmc"
    assert ls_params[-1]["fcstunits"] == "m3 m-3"
    assert ls_params[0]["name"] == ""
    assert ls_params[0]["maskname"] == ""
    assert ls_params[0]["flistname"] == ""
    assert ls_params[11]["name"] == "MYD10C1.Ayyyyddd.061.hdf"
    assert ls_params[8]["maskname"] == "ascat_combined_mask_p1.nc"

    assert len(ol_params) == 10
    assert ol_params[-1]["descr"] == "ASCAT_HSAF_METC_SM"
    assert ol_params[-1]["varname"] == "sfds"
    assert ol_params[-1]["units"] == "%"
    assert ol_params[-1]["fcstvarname"] == "sfmc"
    assert ol_params[0]["flistname"] == "SMAP_L1C_TB_A_list.txt"


def test_read_ObsFcstAna_nc4_matches_file():
    ofa = read_ObsFcstAna_nc4(str(OFA_NC4))
    with nc.Dataset(OFA_NC4) as ds:
        np.testing.assert_array_equal(ofa["obs_species"], ds["species"][:])
        np.testing.assert_array_equal(ofa["obs_tilenum"], ds["tilenum"][:])
        for name in ("lon", "lat", "obs", "obsvar", "fcst", "fcstvar", "ana", "anavar"):
            expected = np.ma.filled(ds[name][:], np.nan).astype(np.float32)
            if name not in ("lon", "lat"):
                expected[expected == -9999.0] = np.nan
            np.testing.assert_array_equal(ofa["obs_" + name], expected)
    assert ofa["obs_lon"].dtype == np.float32


@pytest.mark.skipif(conftest.GRIDCOMP_READERS, reason="the GridComp reader returns empty arrays; the scripts check os.path.exists")
def test_read_ObsFcstAna_nc4_missing_file_returns_None():
    assert read_ObsFcstAna_nc4(str(OFA_NC4) + ".missing") is None


def test_read_ObsFcstAna_nc4_obs_nodata_becomes_nan(tmp_path):
    target = tmp_path / OFA_NC4.name
    shutil.copy2(OFA_NC4, target)
    with nc.Dataset(target, "r+") as ds:
        ds["obs"][0] = -9999.0
    assert np.isnan(read_ObsFcstAna_nc4(str(target))["obs_obs"][0])


def test_read_tilecoord():
    tile_coord = read_tilecoord(INPUTS / "LS_DAv8_M36_as_test2.ldas_tilecoord.bin")
    assert tile_coord["N_tile"] == 112573
    assert tile_coord["tile_id"].shape == (112573,)
    assert np.all(np.isfinite(tile_coord["com_lat"]))


def test_time_helpers():
    assert get_ObsFcstAna_fname("/x/out", "EXP", 2020, 5, 2, 9, 0) == (
        "/x/out/ana/ens_avg/Y2020/M05/EXP.ens_avg.ldas_ObsFcstAna.20200502_0900z.nc4"
    )
    # Feb 2019 - Aug 2020: Jan is only available in 2020, Sep-Dec only in 2019
    run_months = [1, 2, 8, 9]
    first, last = get_year_bounds(run_months, 2019, 2, 2020, 8)
    assert first == [2020, 2019, 2019, 2019]
    assert last == [2020, 2020, 2020, 2019]
    assert get_doy_pentad_range([5, 6, 7], [2019] * 3, [2020] * 3) == (121, 212, 25, 42)
