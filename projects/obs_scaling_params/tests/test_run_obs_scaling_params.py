from pathlib import Path

import netCDF4 as nc
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]      # test data; scripts: see conftest.py

from run_obs_scaling_params import main, parse_args  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"
REQUIRED = ["--exp-path", "/x", "--exp-run", "OL", "--domain", "D", "--start", "2019-02",
            "--end", "2020-08", "--prefix", "p_"]


def test_cli_defaults():
    args = parse_args(["--method", "latlon", "--species", "A,B"] + REQUIRED)
    assert args.out_dir == "z_score_clim"
    assert args.window_days == 75
    assert args.ndata_min == 20
    assert args.grid_resolution == 0.25
    assert args.combine_species and args.all_pentads
    assert not args.each_doy and not args.each_pentad

    args = parse_args(["--method", "tb"] + REQUIRED)
    assert args.orbit == "D"
    assert args.angles == [40.0]
    assert args.descr_contains == "SMAP_L1C"
    assert args.convert_grid


def test_cli_latlon_requires_species():
    with pytest.raises(SystemExit):
        parse_args(["--method", "latlon"] + REQUIRED)


def _write_ofa(path, species, lon, lat, obs, fcst, tilenum):
    path.parent.mkdir(parents=True, exist_ok=True)
    with nc.Dataset(path, "w") as ds:
        ds.createDimension("n_obs", len(obs))
        for name, dtype, data in (
            ("assim_flag", "i4", [1] * len(obs)), ("species", "i4", species), ("tilenum", "i4", tilenum),
            ("lon", "f4", lon), ("lat", "f4", lat), ("obs", "f4", obs), ("obsvar", "f4", [1] * len(obs)),
            ("fcst", "f4", fcst), ("fcstvar", "f4", [1] * len(obs)), ("ana", "f4", fcst),
            ("anavar", "f4", [1] * len(obs)),
        ):
            ds.createVariable(name, dtype, ("n_obs",))[:] = data


def test_tiny_latlon_run(tmp_path):
    # 10-deg grid, 5 daily files, 5-day window -> pentad 1.
    # Obs 1-2: same grid cell (both must count).  Obs 3: on the cell edge lon=-170, lat=-80
    # (ceil -> cell i=0, j=0, the same cell as obs 1-2).  Obs 4: other species.
    # Obs 5: excluded tile.  Obs 6: no-data obs (obs and fcst not counted).
    exp_run, domain = "TINY", "TEST_DOMAIN"
    out = tmp_path / exp_run / "output" / domain
    rc = out / "rc_out" / "Y2020" / "M01"
    rc.mkdir(parents=True)
    (rc / f"{exp_run}.ldas_obsparam.20200101_0000z.txt").write_text(
        (INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt").read_text()
    )
    for day in range(1, 6):
        _write_ofa(
            out / "ana" / "ens_avg" / "Y2020" / "M01" / f"{exp_run}.ens_avg.ldas_ObsFcstAna.202001{day:02d}_0000z.nc4",
            species=[14, 14, 14, 13, 14, 14],
            lon=[-175, -172, -170, -175, -175, 5],
            lat=[-85, -85, -80, -85, -85, 5],
            obs=[0.1 * day, 0.2 * day, 0.3 * day, 9, 9, -9999],
            fcst=[1, 2, 3, 9, 9, 4],
            tilenum=[1, 2, 3, 4, 5, 6],
        )
    csv = tmp_path / "excl.csv"
    csv.write_text("tilenum,reject\n1,False\n5,True\n6,False\n")

    main(["--method", "latlon", "--exp-path", str(tmp_path), "--exp-run", exp_run, "--domain", domain,
          "--start", "2020-01", "--end", "2020-01", "--run-months", "1", "--species", "CYGNSS_SM_6hr",
          "--grid-resolution", "10", "--window-days", "5", "--ndata-min", "1", "--prefix", "t_",
          "--exclude-tiles-csv", str(csv)])

    fname = out / "stats" / "z_score_clim" / "t_2020_doy1_2020_doy31_W_5d_Nmin_1_sp_ALL_all_pentads.nc4"
    with nc.Dataset(fname) as ds:
        assert ds["o_mean"].shape == (73, 36, 18)
        n = ds["n_data"][0]
        assert n[0, 0] == 15                             # 3 obs x 5 days, same cell
        assert np.isnan(n[18, 9])                        # no-data obs only: N=0 -> NaN
        obs = np.array([[0.1 * d, 0.2 * d, 0.3 * d] for d in range(1, 6)], dtype=np.float32).astype(float)
        np.testing.assert_allclose(ds["o_mean"][0][0, 0], obs.mean())
        np.testing.assert_allclose(ds["o_std"][0][0, 0], obs.std())
        np.testing.assert_allclose(ds["m_mean"][0][0, 0], 2.0)
        assert ds["m_min"][0, 0] == 1 and ds["m_max"][0, 0] == 3
        assert np.isnan(ds["o_mean"][1][0, 0])           # pentad 2 not computed
