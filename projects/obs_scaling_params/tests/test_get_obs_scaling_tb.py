from pathlib import Path

import netCDF4 as nc
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]      # test data; scripts: see conftest.py

from get_obs_scaling_tb import (  # noqa: E402
    EASEv2_ind2latlon,
    EASEv2_latlon2ind,
    get_administering_tiles,
    get_obs_scaling_tb,
    select_tb_species,
)
from obs_scaling_utils import read_obs_param, read_tilecoord  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"
OBSPARAM = INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt"


def read_tb_scaling_bin(fname):
    """Minimal reader for the big-endian Tb scaling file (as GEOSldas read_Tb_scale_param.py)."""

    def record(fid, dtype):
        n = int(np.frombuffer(fid.read(4), ">i4")[0])
        values = np.frombuffer(fid.read(n), dtype)
        fid.read(4)
        return values

    with open(fname, "rb") as fid:
        header = record(fid, ">i4")
        start = record(fid, ">i4")
        end = record(fid, ">i4")
        n_out, n_ang, _ = record(fid, ">i4")
        angles = record(fid, ">f4")
        lon = record(fid, ">f4")
        lat = record(fid, ">f4")
        tile_id = record(fid, ">i4")
        data = np.array([[record(fid, ">f4") for _ in range(n_ang)] for _ in range(14)]).transpose(0, 2, 1)
        assert fid.read(1) == b""
    return dict(header=header, start=start, end=end, angles=angles, lon=lon, lat=lat, tile_id=tile_id, data=data)


def test_easev2_m36_forward_inverse_at_cell_centers():
    rows = np.array([0, 100, 202, 405])
    cols = np.array([0, 200, 481, 963])
    lat, lon = EASEv2_ind2latlon(rows, cols, "M36")
    actual_rows, actual_cols = EASEv2_latlon2ind(lat, lon, "M36")
    np.testing.assert_array_equal(actual_rows, rows)
    np.testing.assert_array_equal(actual_cols, cols)


def test_administering_tiles_on_m36_model_grid_are_all_tiles():
    tile_coord = read_tilecoord(INPUTS / "LS_DAv8_M36_as_test2.ldas_tilecoord.bin")
    grid = dict(gridtype="EASEv2_M36", ind_base=0, i_offg=0, j_offg=0, N_lon=964, N_lat=406)
    tiles = get_administering_tiles(tile_coord, grid)
    np.testing.assert_array_equal(np.sort(tiles), np.arange(tile_coord["N_tile"]))   # in M36 cell order


def test_smap_species_selection_is_by_orbit_pol_and_angle():
    tb_species = select_tb_species(read_obs_param(OBSPARAM), "SMAP_L1C", 2, [40.0])
    params = {int(p["species"]): p["descr"] for p in read_obs_param(OBSPARAM)}
    assert {params[s]: v for s, v in tb_species.items()} == {
        "SMAP_L1C_Tbh_D": (0, 0),
        "SMAP_L1C_Tbv_D": (1, 0),
    }


def _write_le_record(fp, values, dtype):
    payload = np.asarray(values, dtype=dtype).tobytes()
    tag = np.asarray([len(payload)], dtype="<i4").tobytes()
    fp.write(tag + payload + tag)


def _write_tiny_tile_files(rc_dir: Path, exp_run: str):
    int_values = {
        "tile_id": [101, 202], "typ": [100, 100], "pfaf": [1, 2],
        "i_indg": [200, 201], "j_indg": [100, 100],
    }
    float_values = {
        "com_lon": [-110.0, -109.0], "com_lat": [32.0, 32.0],
        "min_lon": [-110.5, -109.5], "max_lon": [-109.5, -108.5],
        "min_lat": [31.5, 31.5], "max_lat": [32.5, 32.5],
        "frac_cell": [1.0, 1.0], "frac_pfaf": [1.0, 1.0],
        "area": [1.0, 1.0], "elev": [0.0, 0.0],
    }
    fields = [
        "tile_id", "typ", "pfaf", "com_lon", "com_lat", "min_lon", "max_lon",
        "min_lat", "max_lat", "i_indg", "j_indg", "frac_cell", "frac_pfaf",
        "area", "elev",
    ]
    with (rc_dir / f"{exp_run}.ldas_tilecoord.bin").open("wb") as fp:
        _write_le_record(fp, [2], "<i4")
        for field in fields:
            values = int_values.get(field, float_values.get(field))
            _write_le_record(fp, values, "<i4" if field in int_values else "<f4")

    with (rc_dir / f"{exp_run}.ldas_tilegrids.bin").open("wb") as fp:
        for _ in range(2):
            payload = (
                b"EASEv2_M36".ljust(40, b"\x00")
                + np.asarray([0, 1, 1, 964, 406, 0, 0], dtype="<i4").tobytes()
                + np.asarray([-180, -90, 180, 90, 0, 0], dtype="<f4").tobytes()
            )
            tag = np.asarray([len(payload)], dtype="<i4").tobytes()
            fp.write(tag + payload + tag)


def _write_tiny_ofa(path: Path, h_species: int, v_species: int, day: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("n_obs", 4)
        values = {
            "assim_flag": ("i4", [1, 1, 1, 1]),
            "species": ("i4", [h_species, h_species, v_species, v_species]),
            "tilenum": ("i4", [1, 2, 1, 2]),
            "lon": ("f4", [-110.0, -109.0, -110.0, -109.0]),
            "lat": ("f4", [32.0, 32.0, 32.0, 32.0]),
            "obs": ("f4", [200 + day, 220 + day, 250 + day, 270 + day]),
            "obsvar": ("f4", [1, 1, 1, 1]),
            "fcst": ("f4", [210 + day, 230 + day, 260 + day, 280 + day]),
            "fcstvar": ("f4", [1, 1, 1, 1]),
            "ana": ("f4", [0, 0, 0, 0]),
            "anavar": ("f4", [1, 1, 1, 1]),
        }
        for name, (dtype, data) in values.items():
            ds.createVariable(name, dtype, ("n_obs",))[:] = data


def test_tiny_tb_uses_only_obs_fcst_pairs(tmp_path):
    # tile 2 H-pol obs is no-data on day 5 -> its fcst (235) is not used either
    exp_run = "TINY_OL"
    domain = "SMAP_EASEv2_M36_GLOBAL"
    output = tmp_path / exp_run / "output" / domain
    rc_dir = output / "rc_out"
    rc_dir.mkdir(parents=True)
    _write_tiny_tile_files(rc_dir, exp_run)
    obs_param = read_obs_param(OBSPARAM)
    h_species, v_species = select_tb_species(obs_param, "SMAP_L1C", 2, [40.0])
    for day in range(1, 6):
        fname = (output / "ana" / "ens_avg" / "Y2020" / "M01"
                 / f"{exp_run}.ens_avg.ldas_ObsFcstAna.202001{day:02d}_0000z.nc4")
        _write_tiny_ofa(fname, h_species, v_species, day)
        if day == 5:
            with nc.Dataset(fname, "a") as ds:
                ds["obs"][1] = -9999.0

    get_obs_scaling_tb(
        str(tmp_path), exp_run, domain, obs_param, "SMAP_L1C", 2, [40.0],
        [1], [2020], [2020], 3 * 60 * 60, 0, 5, 1, "tiny_", "z_score_clim", convert_grid=False,
    )

    data = read_tb_scaling_bin(
        output / "stats" / "z_score_clim"
        / "tiny_2020_p1_2020_p6_hscale_0.00_W_1p_Nmin_1_D_p01.bin"
    )["data"]
    np.testing.assert_allclose(data[4, :, 0], [5.0, 4.0])          # N, H-pol
    np.testing.assert_allclose(data[0, :, 0], [203.0, 222.5])      # o_mean
    np.testing.assert_allclose(data[2, :, 0], [213.0, 232.5])      # m_mean: days 1-4 only for tile 2


def test_tiny_tb_generation_produces_expected_pentad_statistics(tmp_path):
    exp_run = "TINY_OL"
    domain = "SMAP_EASEv2_M36_GLOBAL"
    output = tmp_path / exp_run / "output" / domain
    rc_dir = output / "rc_out"
    rc_dir.mkdir(parents=True)
    _write_tiny_tile_files(rc_dir, exp_run)
    obs_param = read_obs_param(OBSPARAM)
    h_species, v_species = select_tb_species(obs_param, "SMAP_L1C", 2, [40.0])
    for day in range(1, 6):
        _write_tiny_ofa(
            output / "ana" / "ens_avg" / "Y2020" / "M01"
            / f"{exp_run}.ens_avg.ldas_ObsFcstAna.202001{day:02d}_0000z.nc4",
            h_species, v_species, day,
        )

    get_obs_scaling_tb(
        str(tmp_path), exp_run, domain, obs_param, "SMAP_L1C", 2, [40.0],
        [1], [2020], [2020], 3 * 60 * 60, 0, 5, 1, "tiny_", "z_score_clim", convert_grid=False,
    )

    result = read_tb_scaling_bin(
        output / "stats" / "z_score_clim"
        / "tiny_2020_p1_2020_p6_hscale_0.00_W_1p_Nmin_1_D_p01.bin"
    )
    np.testing.assert_array_equal(result["header"], [0, 0])
    np.testing.assert_array_equal(result["tile_id"], [101, 202])
    data = result["data"]
    np.testing.assert_allclose(data[0, :, 0], [203.0, 223.0])
    np.testing.assert_allclose(data[1, :, 0], np.sqrt(2.0), atol=1e-6)
    np.testing.assert_allclose(data[2, :, 0], [213.0, 233.0])
    np.testing.assert_allclose(data[4, :, 0], 5.0)
    np.testing.assert_allclose(data[5, :, 0], [253.0, 273.0])
    np.testing.assert_allclose(data[7, :, 0], [263.0, 283.0])
    np.testing.assert_allclose(data[10, :, 0], [203.0, 223.0])
    np.testing.assert_allclose(data[13, :, 0], [263.0, 283.0])
