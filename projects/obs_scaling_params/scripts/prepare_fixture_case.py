"""Create a tiny GEOS-LDAS-like fixture tree for scaling parity runs."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
INPUTS = PROJECT_ROOT / "test_data" / "inputs"

EXP_RUN = "LS_DAv8_M36_as_test2"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
TIMESTAMP = "20200502_0900z"
OBS_PARAM_TIMESTAMP = "20200501_0000z"


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    fixture_root = args.fixture_root.resolve()

    input_obsfcstana_bin = INPUTS / f"{EXP_RUN}.ens_avg.ldas_ObsFcstAna.{TIMESTAMP}.bin"
    input_obsfcstana_nc4 = INPUTS / f"{EXP_RUN}.ens_avg.ldas_ObsFcstAna.{TIMESTAMP}.nc4"
    input_obsparam = INPUTS / f"{EXP_RUN}.ldas_obsparam.{OBS_PARAM_TIMESTAMP}.txt"
    input_tilecoord = INPUTS / f"{EXP_RUN}.ldas_tilecoord.bin"
    obsfcstana_inputs = []
    if args.obsfcstana_format in {"both", "bin"}:
        obsfcstana_inputs.append(input_obsfcstana_bin)
    if args.obsfcstana_format in {"both", "nc4"}:
        obsfcstana_inputs.append(input_obsfcstana_nc4)

    for path in (*obsfcstana_inputs, input_obsparam, input_tilecoord):
        if not path.exists():
            raise FileNotFoundError(f"Required fixture input is missing: {path}")

    output_domain = fixture_root / EXP_RUN / "output" / DOMAIN
    ana_dir = output_domain / "ana" / "ens_avg" / "Y2020" / "M05"
    rc_dir = output_domain / "rc_out" / "Y2020" / "M05"
    stats_dir = output_domain / "stats"
    ana_dir.mkdir(parents=True, exist_ok=True)
    rc_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    for path in obsfcstana_inputs:
        shutil.copy2(path, ana_dir / path.name)
    shutil.copy2(input_obsparam, rc_dir / input_obsparam.name)
    shutil.copy2(input_tilecoord, rc_dir / input_tilecoord.name)

    print("Created obs-scaling fixture")
    print(f"  fixture_root: {fixture_root}")
    print(f"  exp_run:      {EXP_RUN}")
    print(f"  domain:       {DOMAIN}")
    print(f"  timestamp:    {TIMESTAMP}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fixture-root",
        type=Path,
        default=Path("/tmp/obs_scaling_fixture"),
        help="Directory that will contain EXP_RUN/output/DOMAIN.",
    )
    parser.add_argument(
        "--obsfcstana-format",
        choices=("both", "bin", "nc4"),
        default="both",
        help="ObsFcstAna fixture input formats to require and copy.",
    )
    return parser.parse_args(argv)


if __name__ == "__main__":
    main()
