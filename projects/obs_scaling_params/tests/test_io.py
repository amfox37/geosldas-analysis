from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from obs_scaling.io import read_obs_fcst_ana, read_obs_param  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"


def test_read_obs_param_handles_empty_quoted_fields_and_fcst_metadata():
    ls_params = read_obs_param(
        INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt"
    )
    ol_params = read_obs_param(
        INPUTS / "OLv7_M36_MULTI_type_13_H121.ldas_obsparam.20150401_0000z.txt"
    )

    assert len(ls_params) == 14
    assert [param.descr for param in ls_params[:3]] == [
        "SMOS_fit_Tbh_A",
        "SMOS_fit_Tbh_D",
        "SMOS_fit_Tbv_A",
    ]
    assert ls_params[-1].descr == "CYGNSS_SM_6hr"
    assert ls_params[-1].varname == "sfmc"
    assert ls_params[-1].units == "m3 m-3"
    assert ls_params[-1].fcstvarname == "sfmc"
    assert ls_params[-1].fcstunits == "m3 m-3"
    assert ls_params[0].name == ""
    assert ls_params[0].maskname == ""
    assert ls_params[0].flistname == ""
    assert ls_params[11].name == "MYD10C1.Ayyyyddd.061.hdf"
    assert ls_params[8].maskname == "ascat_combined_mask_p1.nc"

    assert len(ol_params) == 10
    assert ol_params[-1].descr == "ASCAT_HSAF_METC_SM"
    assert ol_params[-1].varname == "sfds"
    assert ol_params[-1].units == "%"
    assert ol_params[-1].fcstvarname == "sfmc"
    assert ol_params[-1].fcstunits == "m3 m-3"
    assert ol_params[0].flistname == "SMAP_L1C_TB_A_list.txt"


def test_obsfcstana_binary_matches_netcdf_sample():
    binary_path = (
        INPUTS / "LS_DAv8_M36_as_test2.ens_avg.ldas_ObsFcstAna.20200502_0900z.bin"
    )
    netcdf_path = (
        INPUTS / "LS_DAv8_M36_as_test2.ens_avg.ldas_ObsFcstAna.20200502_0900z.nc4"
    )

    record = read_obs_fcst_ana(binary_path)
    netcdf_record = read_obs_fcst_ana(netcdf_path)
    assert record is not None
    assert netcdf_record is not None
    assert record.date_time.year == 2020
    assert record.date_time.month == 5
    assert record.date_time.day == 2
    assert record.date_time.hour == 9
    assert record.date_time.dofyr == 123
    assert record.date_time.pentad == 25
    assert netcdf_record.date_time == record.date_time
    np.testing.assert_array_equal(netcdf_record.obs_assim, record.obs_assim)
    np.testing.assert_array_equal(netcdf_record.obs_species, record.obs_species)
    np.testing.assert_array_equal(netcdf_record.obs_tilenum, record.obs_tilenum)

    with nc.Dataset(netcdf_path) as dataset:
        assert record.obs_obs.size == len(dataset.dimensions["n_obs"])

        np.testing.assert_array_equal(
            record.obs_assim.astype(np.int32), dataset.variables["assim_flag"][:]
        )
        np.testing.assert_array_equal(record.obs_species, dataset.variables["species"][:])
        np.testing.assert_array_equal(record.obs_tilenum, dataset.variables["tilenum"][:])

        float_pairs = {
            "lon": record.obs_lon,
            "lat": record.obs_lat,
            "obs": record.obs_obs,
            "obsvar": record.obs_obsvar,
            "fcst": record.obs_fcst,
            "fcstvar": record.obs_fcstvar,
            "ana": record.obs_ana,
            "anavar": record.obs_anavar,
        }
        for name, values in float_pairs.items():
            expected = np.ma.filled(dataset.variables[name][:], np.nan)
            if name in {"obsvar", "fcst", "fcstvar", "ana", "anavar"}:
                expected = np.where(expected == -9999.0, np.nan, expected)
            np.testing.assert_allclose(
                values, expected, rtol=0.0, atol=0.0, equal_nan=True
            )

        record_pairs = {
            "lon": (netcdf_record.obs_lon, record.obs_lon),
            "lat": (netcdf_record.obs_lat, record.obs_lat),
            "obs": (netcdf_record.obs_obs, record.obs_obs),
            "obsvar": (netcdf_record.obs_obsvar, record.obs_obsvar),
            "fcst": (netcdf_record.obs_fcst, record.obs_fcst),
            "fcstvar": (netcdf_record.obs_fcstvar, record.obs_fcstvar),
            "ana": (netcdf_record.obs_ana, record.obs_ana),
            "anavar": (netcdf_record.obs_anavar, record.obs_anavar),
        }
        for _, (values, expected) in record_pairs.items():
            np.testing.assert_allclose(
                values, expected, rtol=0.0, atol=0.0, equal_nan=True
            )
