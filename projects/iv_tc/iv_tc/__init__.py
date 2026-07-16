"""Small IV/TC helpers for GEOSldas validation workflows."""

from .config import ProductRoots, RunConfig
from .generation import (
    PairGenerationResult,
    PairInputPaths,
    date_range,
    generate_daily_pairs,
    load_daily_pair_npz,
    pair_npz_is_valid,
    pair_output_path,
    read_daily_pair,
    save_daily_pair_npz,
    write_daily_pair,
)
from .pairs import DailyPair, SparseObservation
from .readers import (
    read_ascat_h119_h120_model_pair,
    read_ascat_h119_h120_sparse,
    read_ascat_h121_model_pair,
    read_ascat_h121_sparse,
    read_smap_l3_model_pair,
    read_smap_l3_sparse,
    read_smosic_model_pair,
    read_smosic_sparse,
)

__all__ = [
    "DailyPair",
    "PairGenerationResult",
    "PairInputPaths",
    "ProductRoots",
    "RunConfig",
    "SparseObservation",
    "date_range",
    "generate_daily_pairs",
    "load_daily_pair_npz",
    "pair_npz_is_valid",
    "pair_output_path",
    "read_ascat_h119_h120_model_pair",
    "read_ascat_h119_h120_sparse",
    "read_ascat_h121_model_pair",
    "read_ascat_h121_sparse",
    "read_smap_l3_model_pair",
    "read_smap_l3_sparse",
    "read_daily_pair",
    "read_smosic_model_pair",
    "read_smosic_sparse",
    "save_daily_pair_npz",
    "write_daily_pair",
]
