from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np
import pytest


CYGNSS_L1_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = CYGNSS_L1_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(CYGNSS_L1_DIR))

from get_obs_scaling_owner_tile import (  # noqa: E402
    cygnss_l1_descr,
    get_obs_scaling_owner_tile,
    get_owner_grid_ind,
    select_cygnss_l1_species,
    write_owner_grid_nc4,
)
from obs_scaling_utils import read_obs_param  # noqa: E402
from run_cygnss_l1_scaling_params import main, parse_args  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"
EXP_RUN = "TINY_CYGL1"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"


def _grid():
    return dict(
        gridtype="EASEv2_M36",
        ind_base=1,
        i_dir=1,
        j_dir=-1,
        N_lon=4,
        N_lat=3,
        i_offg=10,
        j_offg=20,
        ll_lon=-180.0,
        ll_lat=-90.0,
        ur_lon=180.0,
        ur_lat=90.0,
        dlon=0.0,
        dlat=0.0,
    )


def _tile_coord():
    integer = np.asarray([1, 2], dtype=np.int32)
    zeros = np.zeros(2, dtype=np.float32)
    ones = np.ones(2, dtype=np.float32)
    return dict(
        N_tile=2,
        tile_id=np.asarray([101, 202], dtype=np.int32),
        typ=integer,
        pfaf=integer,
        com_lon=np.asarray([-110.0, -109.0], dtype=np.float32),
        com_lat=np.asarray([32.0, 32.0], dtype=np.float32),
        min_lon=zeros,
        max_lon=zeros,
        min_lat=zeros,
        max_lat=zeros,
        i_indg=np.asarray([11, 14], dtype=np.int32),
        j_indg=np.asarray([21, 23], dtype=np.int32),
        frac_cell=ones,
        frac_pfaf=ones,
        area=ones,
        elev=zeros,
    )


def _cygnss_param(scale="F"):
    param = read_obs_param(INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt")[0]
    param.update(
        descr=cygnss_l1_descr,
        species=13.0,
        scale=scale,
        varname="cygl1scal",
        units="dB",
        fcstvarname="cygl1scal",
        fcstunits="dB",
    )
    return param


def test_owner_grid_ind_removes_offset_and_index_base():
    i, j = get_owner_grid_ind(_tile_coord(), _grid())
    np.testing.assert_array_equal(i, [0, 3])
    np.testing.assert_array_equal(j, [0, 2])


def test_owner_grid_ind_requires_easev2_m36():
    grid = _grid()
    grid["gridtype"] = "EASEv2_M09"
    with pytest.raises(ValueError, match="requires EASEv2_M36"):
        get_owner_grid_ind(_tile_coord(), grid)


def test_cygnss_species_is_selected_by_descr_and_must_be_unscaled():
    assert select_cygnss_l1_species([_cygnss_param()])["species"] == 13

    with pytest.raises(ValueError, match="unscaled"):
        select_cygnss_l1_species([_cygnss_param(scale="T")])


def test_owner_grid_writer_uses_dense_zero_based_cells(tmp_path):
    from datetime import datetime

    data = np.full((7, 2, 2), np.nan)
    data[:, 0, 0] = [-20.0, 2.0, -19.0, 1.0, 5.0, -22.0, -17.0]
    data[:, 1, 1] = [-30.0, 3.0, -28.0, 2.0, 6.0, -33.0, -24.0]
    times = [datetime(2014, 1, 1), datetime(2014, 1, 6)]
    fname = str(tmp_path / "owner.nc4")
    write_owner_grid_nc4(fname, data, np.asarray([0, 3]), np.asarray([0, 2]), _grid(), times, times,
                         5, 1, 1e-6, 13, cygnss_l1_descr)

    with nc.Dataset(fname) as ds:
        assert {name: len(dim) for name, dim in ds.dimensions.items()} == {"pentad": 2, "y": 3, "x": 4}
        np.testing.assert_array_equal(ds["x"][:], [0, 1, 2, 3])
        np.testing.assert_array_equal(ds["y"][:], [0, 1, 2])
        assert ds.source_ind_base == 1
        assert ds.source_i_offg == 10
        assert ds["o_mean"][0, 0, 0] == -20.0
        assert ds["o_mean"][1, 2, 3] == -30.0
        assert np.ma.is_masked(ds["o_mean"][0, 1, 1])
        assert ds["m_min"][0, 0] == -22.0
        assert ds["m_max"][2, 3] == -24.0
        assert np.ma.is_masked(ds["m_min"][1, 1])


def _write_record(fp, values, dtype):
    payload = np.asarray(values, dtype=dtype).tobytes()
    tag = np.asarray([len(payload)], dtype="<i4").tobytes()
    fp.write(tag + payload + tag)


def _write_tile_files(rc_dir):
    rc_dir.mkdir(parents=True, exist_ok=True)
    tile_coord = _tile_coord()
    fields = ("tile_id", "typ", "pfaf", "com_lon", "com_lat", "min_lon", "max_lon", "min_lat",
              "max_lat", "i_indg", "j_indg", "frac_cell", "frac_pfaf", "area", "elev")
    with (rc_dir / f"{EXP_RUN}.ldas_tilecoord.bin").open("wb") as fp:
        _write_record(fp, [tile_coord["N_tile"]], "<i4")
        for field in fields:
            values = tile_coord[field]
            _write_record(fp, values, "<i4" if values.dtype.kind in "iu" else "<f4")

    grid = _grid()
    with (rc_dir / f"{EXP_RUN}.ldas_tilegrids.bin").open("wb") as fp:
        for _ in range(2):
            payload = (
                grid["gridtype"].encode("ascii").ljust(40, b"\x00")
                + np.asarray([grid[k] for k in ("ind_base", "i_dir", "j_dir", "N_lon", "N_lat",
                                                "i_offg", "j_offg")], dtype="<i4").tobytes()
                + np.asarray([grid[k] for k in ("ll_lon", "ll_lat", "ur_lon", "ur_lat", "dlon",
                                                "dlat")], dtype="<f4").tobytes()
            )
            tag = np.asarray([len(payload)], dtype="<i4").tobytes()
            fp.write(tag + payload + tag)


def _write_ofa(fname, day, scale=0, duplicate=False):
    fname.parent.mkdir(parents=True, exist_ok=True)
    if duplicate:
        tile = [1, 1, 2]
        obs = [-20.0 + day, -19.5 + day, -30.0 + day]
        fcst = [-10.0 + 2 * day, -9.5 + 2 * day, -15.0 + 3 * day]
    else:
        tile = [1, 2]
        obs = [-20.0 + day, -30.0 + day]
        fcst = [-10.0 + 2 * day, -9999.0 if day == 3 else -15.0 + 3 * day]
    n_obs = len(tile)
    with nc.Dataset(fname, "w") as ds:
        ds.createDimension("n_obs", n_obs)
        ds.createDimension("n_species", 1)
        arrays = {
            "assim_flag": ("i4", np.zeros(n_obs, dtype=np.int32)),
            "species": ("i4", np.full(n_obs, 13, dtype=np.int32)),
            "tilenum": ("i4", tile),
            "lon": ("f4", np.zeros(n_obs)),
            "lat": ("f4", np.zeros(n_obs)),
            "obs": ("f4", obs),
            "obsvar": ("f4", np.ones(n_obs)),
            "fcst": ("f4", fcst),
            "fcstvar": ("f4", np.ones(n_obs)),
            "ana": ("f4", np.zeros(n_obs)),
            "anavar": ("f4", np.ones(n_obs)),
        }
        for name, (dtype, values) in arrays.items():
            ds.createVariable(name, dtype, ("n_obs",))[:] = values
        metadata = {
            "obsparam_species_id": ("i4", 13),
            "obsparam_assim": ("i4", 0),
            "obsparam_scale": ("i4", scale),
            "obsparam_errstd": ("f4", 3.0),
            "obsparam_varname": (str, "cygl1scal"),
            "obsparam_units": (str, "dB"),
            "obsparam_fcstvarname": (str, "cygl1scal"),
            "obsparam_fcstunits": (str, "dB"),
            "obsparam_descr": (str, cygnss_l1_descr),
        }
        for name, (dtype, value) in metadata.items():
            ds.createVariable(name, dtype, ("n_species",))[0] = value


def _ofa_fname(inpath, day):
    return inpath / "ana" / "ens_avg" / "Y2020" / "M01" / f"{EXP_RUN}.ens_avg.ldas_ObsFcstAna.202001{day:02d}_0000z.nc4"


def _generate(tmp_path, **kwargs):
    return get_obs_scaling_owner_tile(
        str(tmp_path), EXP_RUN, DOMAIN, [_cygnss_param()], cygnss_l1_descr,
        [1], [2020], [2020], 86400, 0, 5, 1, "tiny_", "cygnss_l1_z_score_clim",
        std_epsilon=0.0, **kwargs,
    )


def test_tiny_generation_uses_only_paired_samples(tmp_path):
    inpath = tmp_path / EXP_RUN / "output" / DOMAIN
    _write_tile_files(inpath / "rc_out")
    for day in range(1, 7):
        _write_ofa(_ofa_fname(inpath, day), day)

    fname = _generate(tmp_path)
    assert fname.endswith("tiny_2020_doy1_2020_doy31_W_5d_Nmin_1_sp_CYGNSS_L1_DDM3X5_CROP_SCALAR_all_pentads.nc4")
    with nc.Dataset(fname) as ds:
        np.testing.assert_allclose(ds["o_mean"][0, 0, 0], -16.0)
        np.testing.assert_allclose(ds["o_std"][0, 0, 0], np.sqrt(2.0))
        np.testing.assert_allclose(ds["m_mean"][0, 0, 0], -2.0)
        np.testing.assert_allclose(ds["m_std"][0, 0, 0], np.sqrt(8.0))
        assert ds["n_data"][0, 0, 0] == 5.0
        np.testing.assert_allclose(ds["o_mean"][0, 2, 3], -25.75)
        np.testing.assert_allclose(ds["m_mean"][0, 2, 3], -2.25)
        assert ds["n_data"][0, 2, 3] == 4.0
        assert ds["m_min"][2, 3] == -9.0
        assert ds["m_max"][2, 3] == 3.0


def test_obs_keep_drops_obs(tmp_path):
    inpath = tmp_path / EXP_RUN / "output" / DOMAIN
    _write_tile_files(inpath / "rc_out")
    for day in range(1, 7):
        _write_ofa(_ofa_fname(inpath, day), day)

    def obs_keep(time, lon, lat):
        return np.full(len(lon), time.day != 2)  # pentad 1 window is Jan 2-6

    with nc.Dataset(_generate(tmp_path, obs_keep=obs_keep)) as ds:
        assert ds["n_data"][0, 0, 0] == 4.0
        np.testing.assert_allclose(ds["o_mean"][0, 0, 0], -15.5)


@pytest.mark.parametrize(
    ("scale", "duplicate", "message"),
    [(1, False, "already scaled"), (0, True, "more than one CYGNSS L1 obs per tile")],
)
def test_generation_rejects_scaled_or_duplicate_obs(tmp_path, scale, duplicate, message):
    inpath = tmp_path / EXP_RUN / "output" / DOMAIN
    _write_tile_files(inpath / "rc_out")
    _write_ofa(_ofa_fname(inpath, 1), 1, scale=scale, duplicate=duplicate)
    with pytest.raises(ValueError, match=message):
        _generate(tmp_path)


def test_cli_defaults():
    args = parse_args(["--exp-path", "/tmp/exps", "--exp-run", "OL", "--start", "2020-01", "--end", "2022-12"])
    assert args.descr == cygnss_l1_descr
    assert args.domain == DOMAIN
    assert args.window_days == 75
    assert args.ndata_min == 20
    assert args.std_epsilon == 1e-6
    assert args.out_dir == "cygnss_l1_z_score_clim"
    assert args.prefix == "CYGNSS_L1_owner_zscore_"
    assert args.overwrite


def test_cli_tiny_run(tmp_path):
    inpath = tmp_path / EXP_RUN / "output" / DOMAIN
    _write_tile_files(inpath / "rc_out")
    for day in range(1, 7):
        _write_ofa(_ofa_fname(inpath, day), day)
    obsparam = tmp_path / "obsparam.txt"
    obsparam.write_text(
        (INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt").read_text()
    )
    # obsparam with the CYGNSS L1 species (first entry rewritten)
    lines = obsparam.read_text().splitlines()
    lines[0] = "1"
    param_lines = lines[1:37]
    param_lines[0] = f"'{cygnss_l1_descr}'"
    param_lines[1] = "13"
    param_lines[10] = "'F'"
    param_lines[17] = "'cygl1scal'"
    param_lines[18] = "'dB'"
    param_lines[19] = "'cygl1scal'"
    param_lines[20] = "'dB'"
    obsparam.write_text("\n".join([lines[0]] + param_lines) + "\n")

    main(["--exp-path", str(tmp_path), "--exp-run", EXP_RUN, "--start", "2020-01", "--end", "2020-01",
          "--run-months", "1", "--window-days", "5", "--ndata-min", "1", "--std-epsilon", "0",
          "--dt-assim", "86400", "--obsparam", str(obsparam), "--prefix", "cli_"])

    fname = inpath / "stats" / "cygnss_l1_z_score_clim" / (
        "cli_2020_doy1_2020_doy31_W_5d_Nmin_1_sp_CYGNSS_L1_DDM3X5_CROP_SCALAR_all_pentads.nc4")
    with nc.Dataset(fname) as ds:
        assert ds["n_data"][0, 0, 0] == 5.0
