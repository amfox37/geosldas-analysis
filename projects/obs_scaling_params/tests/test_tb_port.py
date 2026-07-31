from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from obs_scaling.easev2 import ind_to_latlon, latlon_to_ind  # noqa: E402
from obs_scaling.io import DateTime, read_obs_param  # noqa: E402
from obs_scaling.seqbin import read_tb_scaling_file, write_tb_scaling_file  # noqa: E402
from obs_scaling.tb_tile_stats import (  # noqa: E402
    build_administering_tiles,
    generate_tb_scaling_params,
    select_tb_species,
)
from obs_scaling.tile_io import TileGrid, read_tilecoord  # noqa: E402
from scripts.run_tb_scaling_params import parse_args  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"


def test_easev2_m36_forward_inverse_at_cell_centers():
    rows = np.array([0, 100, 202, 405])
    cols = np.array([0, 200, 481, 963])
    lat, lon = ind_to_latlon(rows, cols, "M36")
    actual_rows, actual_cols = latlon_to_ind(lat, lon, "M36", rounded=True)
    np.testing.assert_array_equal(actual_rows, rows)
    np.testing.assert_array_equal(actual_cols, cols)


def test_tilecoord_reader_parses_real_fixture_and_identity_admin_mapping():
    tile_coord = read_tilecoord(INPUTS / "LS_DAv8_M36_as_test2.ldas_tilecoord.bin")
    assert tile_coord.n_tile == 112573
    assert tile_coord.tile_id.shape == (tile_coord.n_tile,)
    assert np.all(np.isfinite(tile_coord.com_lon))
    assert np.all(np.isfinite(tile_coord.com_lat))

    grid = TileGrid(
        gridtype="EASEv2_M36",
        ind_base=0,
        i_dir=1,
        j_dir=1,
        n_lon=964,
        n_lat=406,
        i_offg=0,
        j_offg=0,
        ll_lon=-180.0,
        ll_lat=-90.0,
        ur_lon=180.0,
        ur_lat=90.0,
        dlon=0.0,
        dlat=0.0,
    )
    administering = build_administering_tiles(tile_coord, grid, convert_grid=None)
    np.testing.assert_array_equal(
        administering.model_indices, np.arange(tile_coord.n_tile)
    )
    np.testing.assert_array_equal(administering.tile_id, tile_coord.tile_id)


def test_smap_species_selection_is_by_orbit_pol_and_angle():
    params = read_obs_param(
        INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt"
    )
    selected = select_tb_species(
        params, description_contains="SMAP_L1C", orbit=2, angles=[40.0]
    )
    assert [(item.description, item.pol_index, item.angle_index) for item in selected] == [
        ("SMAP_L1C_Tbh_D", 0, 0),
        ("SMAP_L1C_Tbv_D", 1, 0),
    ]


def test_tb_scaling_seqbin_round_trip(tmp_path):
    path = tmp_path / "scaling_D_p25.bin"
    rng = np.random.default_rng(4)
    data = rng.normal(size=(14, 3, 2)).astype(np.float32)
    data[2, 1, 0] = -9999.0
    write_tb_scaling_file(
        path,
        asc_flag=0,
        ndata_min=20,
        start_time=DateTime(2014, 2, 1, 21, 0, 0),
        end_time=DateTime(2014, 4, 17, 21, 0, 0),
        angles=np.array([35.0, 40.0]),
        lon=np.array([-110.0, -109.5, -109.0]),
        lat=np.array([32.0, 32.5, 33.0]),
        tile_id=np.array([101, 202, 303]),
        data=data,
    )
    actual = read_tb_scaling_file(path)
    assert actual.asc_flag == 0
    assert actual.ndata_min == 20
    assert actual.start_time == DateTime(2014, 2, 1, 21, 0, 0)
    assert actual.end_time == DateTime(2014, 4, 17, 21, 0, 0)
    np.testing.assert_array_equal(actual.angles, [35.0, 40.0])
    np.testing.assert_array_equal(actual.tile_id, [101, 202, 303])
    np.testing.assert_array_equal(actual.data, data)


def test_tb_cli_has_production_smap_defaults():
    args = parse_args(
        [
            "--exp-path", "/tmp/archive", "--exp-run", "OL", "--start", "2015-04",
            "--end", "2021-03", "--prefix", "SMAP_stats_",
        ]
    )
    assert args.domain == "SMAP_EASEv2_M09_GLOBAL"
    assert args.orbit == "D"
    assert args.angles == "40"
    assert args.description_contains == "SMAP_L1C"
    assert args.window_days == 75
    assert args.ndata_min == 20
    assert args.obsfcstana_format == "bin"


def _write_le_record(fp, values, dtype):
    payload = np.asarray(values, dtype=dtype).tobytes()
    tag = np.asarray([len(payload)], dtype="<i4").tobytes()
    fp.write(tag)
    fp.write(payload)
    fp.write(tag)


def _write_tiny_tile_files(rc_dir: Path, exp_run: str):
    tilecoord = rc_dir / f"{exp_run}.ldas_tilecoord.bin"
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
    with tilecoord.open("wb") as fp:
        _write_le_record(fp, [2], "<i4")
        for field in fields:
            values = int_values.get(field, float_values.get(field))
            _write_le_record(fp, values, "<i4" if field in int_values else "<f4")

    tilegrids = rc_dir / f"{exp_run}.ldas_tilegrids.bin"
    with tilegrids.open("wb") as fp:
        for _ in range(2):
            gridtype = b"EASEv2_M36".ljust(40, b"\x00")
            integers = np.asarray([0, 1, 1, 964, 406, 0, 0], dtype="<i4").tobytes()
            floats = np.asarray([-180, -90, 180, 90, 0, 0], dtype="<f4").tobytes()
            payload = gridtype + integers + floats
            tag = np.asarray([len(payload)], dtype="<i4").tobytes()
            fp.write(tag + payload + tag)


def _write_tiny_ofa(path: Path, h_species: int, v_species: int, day: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(path, "w") as dataset:
        dataset.createDimension("n_obs", 4)
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
            dataset.createVariable(name, dtype, ("n_obs",))[:] = data


def test_tiny_tb_generation_produces_expected_pentad_statistics(tmp_path):
    exp_run = "TINY_OL"
    domain = "SMAP_EASEv2_M36_GLOBAL"
    output = tmp_path / exp_run / "output" / domain
    rc_dir = output / "rc_out"
    rc_dir.mkdir(parents=True)
    _write_tiny_tile_files(rc_dir, exp_run)
    params = read_obs_param(
        INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt"
    )
    species = select_tb_species(
        params, description_contains="SMAP_L1C", orbit=2, angles=[40.0]
    )
    h_species, v_species = (item.species for item in species)
    for day in range(1, 6):
        path = (
            output / "ana" / "ens_avg" / "Y2020" / "M01"
            / f"{exp_run}.ens_avg.ldas_ObsFcstAna.202001{day:02d}_0000z.nc4"
        )
        _write_tiny_ofa(path, h_species, v_species, day)

    written = generate_tb_scaling_params(
        run_months=[1],
        exp_path=tmp_path,
        exp_run=exp_run,
        domain=domain,
        start_year=[2020],
        end_year=[2020],
        dt_assim=3 * 60 * 60,
        t0_assim=0,
        obs_params=params,
        description_contains="SMAP_L1C",
        orbit=2,
        angles=[40.0],
        window_days=5,
        ndata_min=1,
        prefix="tiny_",
        convert_grid=None,
        obsfcstana_format="nc4",
    )
    pentad_path = next(path for path in written if path.name.endswith("_p01.bin"))
    result = read_tb_scaling_file(pentad_path)
    assert result.ndata_min == 0
    np.testing.assert_allclose(result.data[0, :, 0], [203.0, 223.0])
    np.testing.assert_allclose(result.data[1, :, 0], np.sqrt(2.0), atol=1e-6)
    np.testing.assert_allclose(result.data[2, :, 0], [213.0, 233.0])
    np.testing.assert_allclose(result.data[4, :, 0], 5.0)
    np.testing.assert_allclose(result.data[5, :, 0], [253.0, 273.0])
    np.testing.assert_allclose(result.data[7, :, 0], [263.0, 283.0])
    np.testing.assert_allclose(result.data[10, :, 0], [203.0, 223.0])
    np.testing.assert_allclose(result.data[13, :, 0], [263.0, 283.0])
