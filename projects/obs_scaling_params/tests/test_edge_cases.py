from pathlib import Path
import sys

import netCDF4 as nc
import numpy as np
import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from get_obs_scaling_tb import get_obs_scaling_tb  # noqa: E402
from obs_scaling_utils import read_obs_param  # noqa: E402
from run_obs_scaling_params import main, parse_args  # noqa: E402


INPUTS = PROJECT_ROOT / "test_data" / "inputs"
OBSPARAM = INPUTS / "LS_DAv8_M36_as_test2.ldas_obsparam.20200501_0000z.txt"
EXP_RUN, DOMAIN = "TINY", "TEST_DOMAIN"
CYGNSS_SM, CYGNSS_SM_ID = "CYGNSS_SM_6hr", 14
REQUIRED = ["--exp-path", "/x", "--exp-run", "OL", "--domain", "D", "--prefix", "p_"]


def _write_ofa(fname, species, lon, lat, obs, fcst):
    fname.parent.mkdir(parents=True, exist_ok=True)
    n = len(obs)
    with nc.Dataset(fname, "w") as ds:
        ds.createDimension("n_obs", n)
        for name, dtype, data in (
            ("assim_flag", "i4", [1] * n), ("species", "i4", species), ("tilenum", "i4", range(1, n + 1)),
            ("lon", "f4", lon), ("lat", "f4", lat), ("obs", "f4", obs), ("obsvar", "f4", [1] * n),
            ("fcst", "f4", fcst), ("fcstvar", "f4", [1] * n), ("ana", "f4", fcst), ("anavar", "f4", [1] * n),
        ):
            ds.createVariable(name, dtype, ("n_obs",))[:] = data


def _ofa_fname(out, year, month, day):
    return (out / "ana" / "ens_avg" / f"Y{year}" / f"M{month:02d}" /
            f"{EXP_RUN}.ens_avg.ldas_ObsFcstAna.{year}{month:02d}{day:02d}_0000z.nc4")


def _setup(tmp_path, year=2020, month=1):
    out = tmp_path / EXP_RUN / "output" / DOMAIN
    rc = out / "rc_out" / f"Y{year}" / f"M{month:02d}"
    rc.mkdir(parents=True)
    (rc / f"{EXP_RUN}.ldas_obsparam.{year}{month:02d}01_0000z.txt").write_text(OBSPARAM.read_text())
    return out


def _run_latlon(tmp_path, start, end, run_months, *extra):
    main(["--method", "latlon", "--exp-path", str(tmp_path), "--exp-run", EXP_RUN, "--domain", DOMAIN,
          "--start", start, "--end", end, "--run-months", run_months, "--species", CYGNSS_SM,
          "--grid-resolution", "10", "--window-days", "5", "--ndata-min", "1", "--dt-assim", "86400",
          "--prefix", "t_", *extra])


# ---- driver input checks ---------------------------------------------------------------

@pytest.mark.parametrize("bad", [
    ["--start", "2020-13", "--end", "2020-14"],
    ["--start", "2020-1x", "--end", "2020-12"],
    ["--start", "2021-01", "--end", "2020-12"],
    ["--start", "2020-01", "--end", "2020-12", "--run-months", "0,13"],
    ["--start", "2020-01", "--end", "2020-12", "--run-months", "1,a"],
])
def test_cli_rejects_bad_dates_and_months(bad):
    with pytest.raises(SystemExit):
        parse_args(["--method", "latlon", "--species", "A"] + REQUIRED + bad)


def test_cli_strips_species_and_rejects_empty_or_duplicate_selections():
    dates = ["--start", "2020-01", "--end", "2020-12"]
    args = parse_args(["--method", "latlon", "--species", " A, B ,"] + REQUIRED + dates)
    assert args.species == ["A", "B"]
    with pytest.raises(SystemExit):
        parse_args(["--method", "latlon", "--species", " , "] + REQUIRED + dates)
    with pytest.raises(SystemExit):
        parse_args(["--method", "tb", "--angles", "40,40"] + REQUIRED + dates)


def test_tb_rejects_duplicate_angles(tmp_path):
    with pytest.raises(ValueError, match="duplicate angles"):
        get_obs_scaling_tb(str(tmp_path), EXP_RUN, DOMAIN, read_obs_param(OBSPARAM), "SMAP_L1C", 2, [40.0, 40.0],
                           [1], [2020], [2020], 86400, 0, 5, 1, "t_", "z_score_clim")


# ---- no input ---------------------------------------------------------------------------

def test_latlon_without_ObsFcstAna_files_fails_and_writes_nothing(tmp_path):
    out = _setup(tmp_path)
    with pytest.raises(ValueError, match="no ObsFcstAna files"):
        _run_latlon(tmp_path, "2020-01", "2020-01", "1")
    assert not list((out / "stats").rglob("*.nc4")) if (out / "stats").exists() else True


def test_latlon_without_selected_obs_fails_before_all_pentads_file(tmp_path):
    out = _setup(tmp_path)
    for day in range(1, 6):
        _write_ofa(_ofa_fname(out, 2020, 1, day), [13], [-175], [-85], [0.5], [0.3])     # other species
    with pytest.raises(ValueError, match="no obs used"):
        _run_latlon(tmp_path, "2020-01", "2020-01", "1")
    assert not list((out / "stats").rglob("*all_pentads.nc4"))


def test_tb_without_ObsFcstAna_files_fails(tmp_path):
    with pytest.raises(ValueError, match="no ObsFcstAna files"):
        get_obs_scaling_tb(str(tmp_path), EXP_RUN, DOMAIN, read_obs_param(OBSPARAM), "SMAP_L1C", 2, [40.0],
                           [1], [2020], [2020], 86400, 0, 5, 1, "t_", "z_score_clim")


# ---- calendar ---------------------------------------------------------------------------

def test_latlon_window_wraps_the_year(tmp_path):
    # run months 12,1 for Dec 2019 - Jan 2020; pentad 73 (center DOY 363) uses Dec 28 - Jan 1
    out = _setup(tmp_path, 2019, 12)
    days = [(2019, 12, d) for d in (27, 28, 29, 30, 31)] + [(2020, 1, d) for d in (1, 2)]
    for y, m, d in days:
        _write_ofa(_ofa_fname(out, y, m, d), [CYGNSS_SM_ID], [-175], [-85], [d / 100], [0.3])
    _run_latlon(tmp_path, "2019-12", "2020-01", "12,1")

    fname = out / "stats" / "z_score_clim" / "t_2019_doy335_2020_doy31_W_5d_Nmin_1_sp_ALL_all_pentads.nc4"
    with nc.Dataset(fname) as ds:
        assert ds["n_data"][72][0, 0] == 5
        expected = np.array([28, 29, 30, 31, 1], dtype=np.float32) / np.float32(100)
        np.testing.assert_allclose(ds["o_mean"][72][0, 0], expected.astype(float).mean())


def test_latlon_skips_feb29(tmp_path):
    # the climatological year has no Feb 29; pentad 12 (center DOY 59 = Feb 28) uses
    # Feb 26, 27, 28, Mar 1, 2 of 2020 and not the Feb 29 file
    out = _setup(tmp_path, 2020, 2)
    for m, d in [(2, 26), (2, 27), (2, 28), (2, 29), (3, 1), (3, 2)]:
        _write_ofa(_ofa_fname(out, 2020, m, d), [CYGNSS_SM_ID], [-175], [-85],
                   [9.0 if d == 29 else 0.5], [0.3])
    _run_latlon(tmp_path, "2020-02", "2020-03", "2,3")

    fname = out / "stats" / "z_score_clim" / "t_2020_doy32_2020_doy90_W_5d_Nmin_1_sp_ALL_all_pentads.nc4"
    with nc.Dataset(fname) as ds:
        assert ds["n_data"][11][0, 0] == 5
        assert ds["o_mean"][11][0, 0] == 0.5


# ---- obs/fcst pairs --------------------------------------------------------------------

def test_latlon_uses_only_obs_fcst_pairs(tmp_path):
    # same cell: valid pair (obs 0.1, fcst 0.2) and no-data obs with fcst 10 -> fcst 10 not used
    out = _setup(tmp_path)
    for day in range(2, 7):
        _write_ofa(_ofa_fname(out, 2020, 1, day), [CYGNSS_SM_ID] * 2, [-175, -172], [-85, -85],
                   [0.1, -9999.], [0.2, 10.])
    _run_latlon(tmp_path, "2020-01", "2020-01", "1")

    fname = out / "stats" / "z_score_clim" / "t_2020_doy1_2020_doy31_W_5d_Nmin_1_sp_ALL_all_pentads.nc4"
    with nc.Dataset(fname) as ds:
        assert ds["n_data"][0][0, 0] == 5
        np.testing.assert_allclose(ds["m_mean"][0][0, 0], np.float32(0.2))
        np.testing.assert_allclose(ds["m_std"][0][0, 0], 0., atol=1e-7)
        np.testing.assert_allclose(ds["m_max"][0, 0], np.float32(0.2))


# ---- species ----------------------------------------------------------------------------

def test_latlon_separate_species_files_have_their_own_statistics(tmp_path):
    # --species not in species-ID order: MOD10C1 (13) obs 0.1, CYGNSS_SM_6hr (14) obs 0.9
    out = _setup(tmp_path)
    for day in range(2, 7):
        _write_ofa(_ofa_fname(out, 2020, 1, day), [13, 14], [-175, -175], [-85, -85], [0.1, 0.9], [0.3, 0.3])
    main(["--method", "latlon", "--exp-path", str(tmp_path), "--exp-run", EXP_RUN, "--domain", DOMAIN,
          "--start", "2020-01", "--end", "2020-01", "--run-months", "1", "--species", "CYGNSS_SM_6hr,MOD10C1",
          "--no-combine-species", "--grid-resolution", "10", "--window-days", "5", "--ndata-min", "1",
          "--dt-assim", "86400", "--prefix", "t_"])

    base = out / "stats" / "z_score_clim" / "t_2020_doy1_2020_doy31_W_5d_Nmin_1_sp_"
    for name, value in (("CYGNSS_SM_6hr", 0.9), ("MOD10C1", 0.1)):
        with nc.Dataset(f"{base}{name}_all_pentads.nc4") as ds:
            np.testing.assert_allclose(ds["o_mean"][0][0, 0], np.float32(value))


def test_latlon_rejects_duplicate_species(tmp_path):
    out = _setup(tmp_path)
    _write_ofa(_ofa_fname(out, 2020, 1, 2), [14], [-175], [-85], [0.5], [0.3])
    with pytest.raises(ValueError, match="duplicate species"):
        main(["--method", "latlon", "--exp-path", str(tmp_path), "--exp-run", EXP_RUN, "--domain", DOMAIN,
              "--start", "2020-01", "--end", "2020-01", "--run-months", "1",
              "--species", "CYGNSS_SM_6hr,CYGNSS_SM_6hr", "--prefix", "t_"])


# ---- output options ---------------------------------------------------------------------

def test_latlon_each_pentad_without_all_pentads_and_chunking(tmp_path):
    out = _setup(tmp_path)
    for day in range(2, 7):                 # daily files at 0z: pentad 1 window is Jan 2-6
        _write_ofa(_ofa_fname(out, 2020, 1, day), [CYGNSS_SM_ID], [-175], [-85], [0.5], [0.3])
    _run_latlon(tmp_path, "2020-01", "2020-01", "1", "--each-pentad", "--no-all-pentads")

    files = sorted(p.name for p in (out / "stats" / "z_score_clim").glob("*.nc4"))
    assert files == [f"t_2020_p1_2020_p6_W_1p_Nmin_1_sp_ALL_p{p:02d}.nc4" for p in range(1, 7)]
    with nc.Dataset(out / "stats" / "z_score_clim" / files[0]) as ds:
        assert ds["n_data"][0][0, 0] == 5
        assert ds["o_mean"].chunking() == [1, 36, 18]
        assert ds["m_min"].chunking() == [36, 18]


def test_latlon_header_only_exclusion_csv_excludes_nothing(tmp_path):
    out = _setup(tmp_path)
    for day in range(2, 7):
        _write_ofa(_ofa_fname(out, 2020, 1, day), [CYGNSS_SM_ID], [-175], [-85], [0.5], [0.3])
    csv = tmp_path / "excl.csv"
    csv.write_text("tilenum,reject\n")
    _run_latlon(tmp_path, "2020-01", "2020-01", "1", "--exclude-tiles-csv", str(csv))

    fname = out / "stats" / "z_score_clim" / "t_2020_doy1_2020_doy31_W_5d_Nmin_1_sp_ALL_all_pentads.nc4"
    with nc.Dataset(fname) as ds:
        assert ds["n_data"][0][0, 0] == 5
