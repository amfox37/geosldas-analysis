"""
Python driver that mirrors `Run_get_model_and_obs_clim_stats_latlon_grid.m`.

It prepares the same experiment configuration and invokes the Python port of
`get_model_and_obs_clim_stats_latlon_grid`.
"""

from __future__ import annotations

from datetime import datetime
from math import ceil
from pathlib import Path
from typing import Iterable

try:
    from .get_model_and_obs_clim_stats_latlon_grid import (
        get_model_and_obs_clim_stats_latlon_grid,
    )
    from .obs_scaling_utils import read_obs_param
except ImportError:
    # Allow running as a standalone script
    import sys
    from pathlib import Path
    THIS_DIR = Path(__file__).resolve().parent
    if str(THIS_DIR) not in sys.path:
        sys.path.append(str(THIS_DIR))
    from get_model_and_obs_clim_stats_latlon_grid import (
        get_model_and_obs_clim_stats_latlon_grid,
    )
    from obs_scaling_utils import read_obs_param

def main() -> None:
    exp_path = Path("/discover/nobackup/projects/land_da/Experiment_archive/M21C_land_sweeper_OLv8_M36")
    exp_run = "LS_OLv8_M36"
    domain = "SMAP_EASEv2_M36_GLOBAL"
    prefix_out = "M36_python_dedup_zscore_stats_"

    start_month = 6
    start_year = 2007
    end_month = 5
    end_year = 2024

    species_names = ["ASCAT_META_SM", "ASCAT_METB_SM", "ASCAT_METC_SM"]
    combine_species_stats = True

    grid_resolution = 0.25
    w_days = 75
    ndata_min = 20
    dt_assim = 3 * 60 * 60
    t0_assim = 0

    print_each_doy = True
    print_each_pentad = False
    print_all_pentads = True

    out_dir = "python_z_score_dedup_clim_quarter_degree"

    run_months = list(range(1, 13)) + list(range(1, ceil(w_days / 30) + 1))
    earliest_year, latest_year = _compute_year_bounds(
        run_months, start_year, start_month, end_year, end_month
    )

    obs_param_fname = (
        exp_path
        / exp_run
        / "output"
        / domain
        / f"rc_out/Y{start_year:04d}/M{start_month:02d}"
        / f"{exp_run}.ldas_obsparam.{start_year:04d}{start_month:02d}01_0000z.txt"
    )

    obs_params = read_obs_param(obs_param_fname)
    species = _collect_species_ids(obs_params, species_names)

    get_model_and_obs_clim_stats_latlon_grid(
        species_names=species_names,
        run_months=run_months,
        exp_path=str(exp_path),
        exp_run=exp_run,
        domain=domain,
        start_year=earliest_year,
        end_year=latest_year,
        dt_assim=dt_assim,
        t0_assim=t0_assim,
        species=species,
        combine_species_stats=combine_species_stats,
        resol=grid_resolution,
        w_days=w_days,
        ndata_min=ndata_min,
        prefix=prefix_out,
        print_each_DOY=print_each_doy,
        print_each_pentad=print_each_pentad,
        print_all_pentads=print_all_pentads,
        out_dir=out_dir,
    )


def _compute_year_bounds(
    run_months: Iterable[int],
    start_year: int,
    start_month: int,
    end_year: int,
    end_month: int,
) -> tuple[list[int], list[int]]:
    earliest = []
    latest = []
    start_ref = datetime(start_year, start_month, 1)
    end_ref = datetime(end_year, end_month, 1)
    for month in run_months:
        current_start = datetime(start_year, month, 1)
        earliest.append(start_year + 1 if current_start < start_ref else start_year)
        current_end = datetime(end_year, month, 1)
        latest.append(end_year - 1 if current_end > end_ref else end_year)
    return earliest, latest


def _collect_species_ids(obs_params, species_names: list[str]) -> list[int]:
    ids: list[int] = []
    for name in species_names:
        matches = [param.species for param in obs_params if param.descr == name]
        if not matches:
            raise ValueError(f"Species '{name}' not found in ldas_obsparam file")
        ids.extend(matches)
    return sorted(set(ids))


if __name__ == "__main__":
    main()
