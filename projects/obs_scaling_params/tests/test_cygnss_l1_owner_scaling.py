from dataclasses import replace
from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from obs_scaling.io import DateTime, read_obs_param  # noqa: E402
from obs_scaling.owner_grid_io import write_netcdf_owner_grid  # noqa: E402
from obs_scaling.owner_tile_stats import (  # noqa: E402
    CYGNSS_L1_DESCRIPTION,
    build_owner_grid_mapping,
    generate_cygnss_l1_scaling_params,
    select_cygnss_l1_species,
)
from obs_scaling.tile_io import TileCoordinates, TileGrid  # noqa: E402
from scripts.run_cygnss_l1_scaling_params import parse_args  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"


def _grid() -> TileGrid:
    return TileGrid(
        gridtype="EASEv2_M36",
        ind_base=1,
        i_dir=1,
        j_dir=-1,
        n_lon=4,
        n_lat=3,
        i_offg=10,
        j_offg=20,
        ll_lon=-180.0,
        ll_lat=-90.0,
        ur_lon=180.0,
        ur_lat=90.0,
        dlon=0.0,
        dlat=0.0,
    )


def _tile_coordinates() -> TileCoordinates:
    integer = np.asarray([1, 2], dtype=np.int32)
    zeros = np.zeros(2, dtype=np.float32)
    ones = np.ones(2, dtype=np.float32)
    return TileCoordinates(
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


def _cygnss_param():
    base = read_obs_param(
        INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt"
    )[0]
    return replace(
        base,
        descr=CYGNSS_L1_DESCRIPTION,
        species=13,
        scale="F",
        varname="cygl1scal",
        units="dB",
        fcstvarname="cygl1scal",
        fcstunits="dB",
    )


def test_owner_grid_mapping_round_trips_offset_and_index_base():
    mapping = build_owner_grid_mapping(_tile_coordinates(), _grid())
    np.testing.assert_array_equal(mapping.i, [0, 3])
    np.testing.assert_array_equal(mapping.j, [0, 2])
    np.testing.assert_array_equal(mapping.linear, [0, 11])


def test_cygnss_species_is_resolved_by_description_and_requires_unscaled_db():
    selected = select_cygnss_l1_species([_cygnss_param()])
    assert selected.species == 13

    with pytest.raises(ValueError, match="unscaled monitoring run"):
        select_cygnss_l1_species([replace(selected, scale="T")])


def test_owner_grid_writer_uses_dense_zero_based_cells(tmp_path):
    data = np.full((7, 2, 2), np.nan)
    data[:, 0, 0] = [-20.0, 2.0, -19.0, 1.0, 5.0, -22.0, -17.0]
    data[:, 1, 1] = [-30.0, 3.0, -28.0, 2.0, 6.0, -33.0, -24.0]
    times = [DateTime(2014, 1, 1, 0, 0, 0), DateTime(2014, 1, 6, 0, 0, 0)]
    path = tmp_path / "owner.nc4"
    write_netcdf_owner_grid(
        path,
        data=data,
        tile_i=np.asarray([0, 3]),
        tile_j=np.asarray([0, 2]),
        grid=_grid(),
        pentads=[1, 2],
        start_time=times,
        end_time=times,
        window_days=5,
        ndata_min=1,
        std_epsilon=1e-6,
        species_id=13,
        species_description=CYGNSS_L1_DESCRIPTION,
    )

    with nc.Dataset(path) as dataset:
        assert {name: len(dim) for name, dim in dataset.dimensions.items()} == {
            "pentad": 2,
            "y": 3,
            "x": 4,
        }
        np.testing.assert_array_equal(dataset["x"][:], [0, 1, 2, 3])
        np.testing.assert_array_equal(dataset["y"][:], [0, 1, 2])
        assert dataset.source_ind_base == 1
        assert dataset.source_i_offg == 10
        assert dataset["o_mean"][0, 0, 0] == -20.0
        assert dataset["o_mean"][1, 2, 3] == -30.0
        assert np.ma.is_masked(dataset["o_mean"][0, 1, 1])
        assert dataset["m_min"][0, 0] == -22.0
        assert dataset["m_max"][2, 3] == -24.0


def _write_record(fp, values, dtype):
    payload = np.asarray(values, dtype=dtype).tobytes()
    tag = np.asarray([len(payload)], dtype="<i4").tobytes()
    fp.write(tag)
    fp.write(payload)
    fp.write(tag)


def _write_tile_files(rc_dir: Path, exp_run: str) -> None:
    rc_dir.mkdir(parents=True, exist_ok=True)
    coords = _tile_coordinates()
    fields = (
        "tile_id",
        "typ",
        "pfaf",
        "com_lon",
        "com_lat",
        "min_lon",
        "max_lon",
        "min_lat",
        "max_lat",
        "i_indg",
        "j_indg",
        "frac_cell",
        "frac_pfaf",
        "area",
        "elev",
    )
    with (rc_dir / f"{exp_run}.ldas_tilecoord.bin").open("wb") as fp:
        _write_record(fp, [coords.n_tile], "<i4")
        for field in fields:
            values = getattr(coords, field)
            dtype = "<i4" if values.dtype.kind in "iu" else "<f4"
            _write_record(fp, values, dtype)

    with (rc_dir / f"{exp_run}.ldas_tilegrids.bin").open("wb") as fp:
        for grid in (_grid(), _grid()):
            payload = (
                grid.gridtype.encode("ascii").ljust(40, b"\x00")
                + np.asarray(
                    [
                        grid.ind_base,
                        grid.i_dir,
                        grid.j_dir,
                        grid.n_lon,
                        grid.n_lat,
                        grid.i_offg,
                        grid.j_offg,
                    ],
                    dtype="<i4",
                ).tobytes()
                + np.asarray(
                    [grid.ll_lon, grid.ll_lat, grid.ur_lon, grid.ur_lat, grid.dlon, grid.dlat],
                    dtype="<f4",
                ).tobytes()
            )
            tag = np.asarray([len(payload)], dtype="<i4").tobytes()
            fp.write(tag + payload + tag)


def _write_ofa(
    path: Path,
    day: int,
    *,
    scale: int = 0,
    duplicate: bool = False,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tile = [1, 1, 2] if duplicate else [1, 2]
    obs = [-20.0 + day, -19.5 + day, -30.0 + day] if duplicate else [-20.0 + day, -30.0 + day]
    model = [-10.0 + 2 * day, -9.5 + 2 * day, -15.0 + 3 * day] if duplicate else [
        -10.0 + 2 * day,
        -9999.0 if day == 3 else -15.0 + 3 * day,
    ]
    n_obs = len(tile)
    with nc.Dataset(path, "w") as dataset:
        dataset.createDimension("n_obs", n_obs)
        dataset.createDimension("n_species", 1)
        arrays = {
            "assim_flag": ("i4", np.zeros(n_obs, dtype=np.int32)),
            "species": ("i4", np.full(n_obs, 13, dtype=np.int32)),
            "tilenum": ("i4", tile),
            "lon": ("f4", np.zeros(n_obs)),
            "lat": ("f4", np.zeros(n_obs)),
            "obs": ("f4", obs),
            "obsvar": ("f4", np.ones(n_obs)),
            "fcst": ("f4", model),
            "fcstvar": ("f4", np.ones(n_obs)),
            "ana": ("f4", np.zeros(n_obs)),
            "anavar": ("f4", np.ones(n_obs)),
        }
        for name, (dtype, values) in arrays.items():
            dataset.createVariable(name, dtype, ("n_obs",))[:] = values
        metadata = {
            "obsparam_species_id": ("i4", 13),
            "obsparam_assim": ("i4", 0),
            "obsparam_scale": ("i4", scale),
            "obsparam_errstd": ("f4", 3.0),
            "obsparam_varname": (str, "cygl1scal"),
            "obsparam_units": (str, "dB"),
            "obsparam_fcstvarname": (str, "cygl1scal"),
            "obsparam_fcstunits": (str, "dB"),
            "obsparam_descr": (str, CYGNSS_L1_DESCRIPTION),
        }
        for name, (dtype, value) in metadata.items():
            variable = dataset.createVariable(name, dtype, ("n_species",))
            variable[0] = value


def _generation_args(tmp_path: Path) -> dict:
    return {
        "run_months": [1],
        "exp_path": tmp_path,
        "exp_run": "TINY_CYGL1",
        "domain": "SMAP_EASEv2_M36_GLOBAL",
        "start_year": [2020],
        "end_year": [2020],
        "dt_assim": 86400,
        "t0_assim": 0,
        "obs_params": [_cygnss_param()],
        "window_days": 5,
        "ndata_min": 1,
        "prefix": "tiny_",
        "std_epsilon": 0.0,
    }


def _tiny_output(tmp_path: Path) -> Path:
    args = _generation_args(tmp_path)
    return (
        tmp_path
        / args["exp_run"]
        / "output"
        / args["domain"]
    )


def test_tiny_generation_uses_strictly_paired_samples(tmp_path):
    output = _tiny_output(tmp_path)
    _write_tile_files(output / "rc_out", "TINY_CYGL1")
    for day in range(1, 7):
        _write_ofa(
            output
            / "ana"
            / "ens_avg"
            / "Y2020"
            / "M01"
            / f"TINY_CYGL1.ens_avg.ldas_ObsFcstAna.202001{day:02d}_0000z.nc4",
            day,
        )

    path = generate_cygnss_l1_scaling_params(**_generation_args(tmp_path))
    with nc.Dataset(path) as dataset:
        np.testing.assert_allclose(dataset["o_mean"][0, 0, 0], -16.0)
        np.testing.assert_allclose(dataset["o_std"][0, 0, 0], np.sqrt(2.0))
        np.testing.assert_allclose(dataset["m_mean"][0, 0, 0], -2.0)
        np.testing.assert_allclose(dataset["m_std"][0, 0, 0], np.sqrt(8.0))
        assert dataset["n_data"][0, 0, 0] == 5.0
        np.testing.assert_allclose(dataset["o_mean"][0, 2, 3], -25.75)
        np.testing.assert_allclose(dataset["m_mean"][0, 2, 3], -2.25)
        assert dataset["n_data"][0, 2, 3] == 4.0
        assert dataset["m_min"][2, 3] == -9.0
        assert dataset["m_max"][2, 3] == 3.0


@pytest.mark.parametrize(
    ("scale", "duplicate", "message"),
    [(1, False, "already scaled"), (0, True, "Duplicate CYGNSS L1 owner tiles")],
)
def test_generation_rejects_invalid_ofa_contract(tmp_path, scale, duplicate, message):
    output = _tiny_output(tmp_path)
    _write_tile_files(output / "rc_out", "TINY_CYGL1")
    _write_ofa(
        output
        / "ana"
        / "ens_avg"
        / "Y2020"
        / "M01"
        / "TINY_CYGL1.ens_avg.ldas_ObsFcstAna.20200101_0000z.nc4",
        1,
        scale=scale,
        duplicate=duplicate,
    )
    with pytest.raises(ValueError, match=message):
        generate_cygnss_l1_scaling_params(**_generation_args(tmp_path))


def test_cygnss_cli_defaults_to_owner_grid_contract():
    args = parse_args(
        [
            "--exp-path",
            "/tmp/archive",
            "--exp-run",
            "OL",
            "--start",
            "2020-01",
            "--end",
            "2022-12",
        ]
    )
    assert args.description == CYGNSS_L1_DESCRIPTION
    assert args.window_days == 75
    assert args.ndata_min == 20
    assert args.std_epsilon == 1e-6
    assert args.out_dir == "cygnss_l1_z_score_clim"
