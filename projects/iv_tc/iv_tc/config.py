"""Configuration objects for the IV/TC workflow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


DEFAULT_DOMAIN = "SMAP_EASEv2_M36_GLOBAL"
DEFAULT_COLLECTION = ".tavg24_1d_lnd_Nt."


@dataclass(frozen=True)
class ProductRoots:
    """Discover roots for observation products used by IV/TC fixtures."""

    ascat_h121_root: Path = Path("/discover/nobackup/projects/land_da/ASCAT_SSM_CDR")
    ascat_h119_h120_root: Path = Path("/discover/nobackup/qliu/merra_land/DATA/ASCAT_HSAF")
    smosic_root: Path = Path("/discover/nobackup/projects/land_da/SMOS_IC/preprocessed_m36_daily")
    smap_l3_root: Path = Path("/discover/nobackup/projects/land_da/Evaluation/IVs/data/SPL3SMP_v009")
    cygnss_l3_root: Path = Path("/discover/nobackup/projects/gmao/smap/SMAP_Nature/CYGNSS")
    domain: str = DEFAULT_DOMAIN
    collection: str = DEFAULT_COLLECTION


@dataclass(frozen=True)
class RunConfig:
    """A GEOSldas experiment/run root."""

    name: str
    root: Path


def parse_run_spec(spec: str) -> RunConfig:
    """Parse NAME=/path/to/run into a RunConfig."""

    if "=" not in spec:
        raise ValueError(f"Run must be NAME=PATH, got {spec!r}")
    name, root = spec.split("=", 1)
    name = name.strip()
    root = root.strip()
    if not name or not root:
        raise ValueError(f"Run must be NAME=PATH, got {spec!r}")
    return RunConfig(name=name, root=Path(root))
