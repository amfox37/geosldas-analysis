# Run the tests against the scripts in this project (default) or against another copy,
# e.g. the copy in GEOSldas_GridComp:
#
#   OBS_SCALING_DIR=<GEOSldas_GridComp>/GEOSldas_App/util/inputs/obs_scaling_params python -m pytest -q tests
#
# The GridComp copy reads the GEOSldas files with ../../shared/python/read_GEOSldas.py
# instead of the readers in obs_scaling_utils.py; the tests then use those readers.

import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(os.environ.get("OBS_SCALING_DIR", PROJECT_ROOT)).resolve()

sys.path.insert(0, str(SCRIPT_DIR))

import obs_scaling_utils  # noqa: E402

GRIDCOMP_READERS = not hasattr(obs_scaling_utils, "read_obs_param")

if GRIDCOMP_READERS:
    sys.path.insert(1, str((SCRIPT_DIR / "../../shared/python").resolve()))
    import read_GEOSldas  # noqa: E402

    for name in ("read_obs_param", "read_tilecoord", "read_tilegrids", "read_ObsFcstAna_nc4"):
        setattr(obs_scaling_utils, name, getattr(read_GEOSldas, name))


def pytest_report_header(config):
    return f"obs scaling scripts: {SCRIPT_DIR}" + (" (shared GridComp readers)" if GRIDCOMP_READERS else "")
