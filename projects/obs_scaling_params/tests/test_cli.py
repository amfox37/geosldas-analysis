from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.run_scaling_params import parse_args  # noqa: E402


def test_dedup_is_opt_in():
    args = parse_args([])
    assert args.enable_dedup is False
    assert args.prefix == "M36_python_zscore_stats_"
    assert args.out_dir == "python_z_score_clim_quarter_degree"

    args = parse_args(["--dedup"])
    assert args.enable_dedup is True

