from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat, savemat

# Make the shared EASE grid utilities available
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(REPO_ROOT / "shared" / "python"))
from EASEv2 import EASEv2_ind2latlon  # noqa: E402


def _as_vector(mat_array: np.ndarray) -> np.ndarray:
    """Flatten MATLAB-loaded arrays using column-major order."""
    return np.asarray(mat_array, dtype=float).ravel(order="F")


def load_skill_file(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load R2 fields from a MATLAB stats file."""
    data = loadmat(path, simplify_cells=True)
    return _as_vector(data["R2_ivs_obs"]), _as_vector(data["R2_ivs_mod"])


def build_ease_grid() -> tuple[np.ndarray, np.ndarray]:
    """Build latitude/longitude grids for the M36 EASEv2 layout."""
    rows = np.arange(406)
    cols = np.arange(964)
    lat, lon = EASEv2_ind2latlon(rows, cols, "M36")
    lon_grid = np.tile(lon[:, None], (1, lat.size))
    lat_grid = np.tile(lat, (lon.size, 1))
    return lat_grid, lon_grid


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute R differences between two IVS skill files and save a .mat output.",
    )
    parser.add_argument(
        "--data-path",
        type=Path,
        default=Path.home() / "for_Andy/IVs",
        help="Directory containing the IVS stats .mat files.",
    )
    parser.add_argument(
        "--d1-version",
        default="OLv7_M36_MULTI_type_13_comb_fp_scaled",
        help="Label for the first dataset (D1).",
    )
    parser.add_argument(
        "--d2-version",
        default="DAv7_M36_ASCAT_type_13_comb_fp_scaled",
        help="Label for the second dataset (D2).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        help="Directory to write the output .mat file.",
    )
    parser.add_argument(
        "--plot",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Whether to display a pcolor plot of the R difference.",
    )
    args = parser.parse_args()

    data_path = args.data_path.expanduser()
    f_D1 = data_path / f"SMPL3_{args.d1_version}_IVD_IVS_stats_lag2day_201504_202103.mat"
    f_D2 = data_path / f"SMPL3_{args.d2_version}_IVD_IVS_stats_lag2day_201504_202103.mat"

    R2_obs_D1, R2_mod_D1 = load_skill_file(f_D1)
    R2_obs_D2, R2_mod_D2 = load_skill_file(f_D2)

    R_D1 = np.sqrt(R2_mod_D1)
    R_D2 = np.sqrt(R2_mod_D2)
    R_OBS = np.sqrt(R2_obs_D2)

    R_D1[R_D1 < 0.1] = np.nan
    R_D2[R_D2 < 0.1] = np.nan
    R_OBS[R_OBS < 0.1] = np.nan

    R_D1[np.isnan(R_OBS)] = np.nan
    R_D2[np.isnan(R_OBS)] = np.nan

    Rdiff_vector = R_D2 - R_D1

    lat_grid, lon_grid = build_ease_grid()
    Rdiff_grid = np.reshape(Rdiff_vector, lon_grid.shape, order="F")

    if args.plot:
        fig, ax = plt.subplots(figsize=(12, 5))
        mesh = ax.pcolormesh(lon_grid, lat_grid, Rdiff_grid, shading="flat")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title(f"R diff: {args.d2_version} minus {args.d1_version}")
        fig.colorbar(mesh, ax=ax, label="R difference (D2 - D1)")
        plt.tight_layout()
        plt.show()

    outname = f"Rdiff_{args.d2_version}_minus_{args.d1_version}.mat"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    savemat(
        args.output_dir / outname,
        {
            "Rdiff_vector": Rdiff_vector,
            "lats": lat_grid.reshape(-1, order="F"),
            "lons": lon_grid.reshape(-1, order="F"),
        },
    )


if __name__ == "__main__":
    main()
