#!/usr/bin/env python3
"""Load audited M21C monthly fields for trend and breakpoint analysis."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import xarray as xr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT.parents[1]
DEFAULT_DATA_DIR = Path(
    "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/M21C_land_sweeper_v2"
)
DEFAULT_INPUT_CONTRACT = PROJECT_ROOT / "config" / "trend_breakpoint_inputs.json"
DEFAULT_VARIABLE_SELECTION = (
    PROJECT_ROOT / "config" / "trend_breakpoint_variable_selection.json"
)

SNOW_POSSIBLE_SCF_THRESHOLD = 0.05
SNOW_POSSIBLE_SWE_THRESHOLD = 5.0
PERMANENT_SNOW_JJA_MEAN_MAX = 0.20
NH_MIN_LAT = 20.0
WARM_SNOWFREE_SCF_MAX = 0.05
WARM_SNOWFREE_TSOIL_MIN_K = 277.15
WATER_RATE_UNITS = "kg m-2 s-1"
MONTHLY_WATER_UNITS = "kg m-2 month-1"


def _load_json(path: str | Path) -> dict[str, Any]:
    return json.loads(Path(path).read_text())


def find_legacy_increment_variables(values: Any) -> set[str]:
    """Return legacy monthly increment names, without rejecting raw diagnostics."""

    names = {str(value) for value in values}
    return {name for name in names if name.upper().endswith(("_INC", "_INCR"))}


def load_variable_selection(
    path: str | Path = DEFAULT_VARIABLE_SELECTION,
) -> pd.DataFrame:
    """Load the audited variable selection table."""

    payload = _load_json(path)
    frame = pd.DataFrame(payload["variables"], columns=payload["columns"])
    if frame["analysis_variable"].duplicated().any():
        raise ValueError("Duplicate analysis_variable entries in variable selection")
    legacy_variables = find_legacy_increment_variables(frame["source_variable"])
    if legacy_variables:
        raise ValueError(f"Legacy monthly increment variables selected: {legacy_variables}")
    return frame.set_index("analysis_variable", drop=False)


def _normalise_units(value: Any) -> str:
    return " ".join(str(value or "").split())


def read_tile_area(path: Path) -> np.ndarray:
    io_dir = REPO_ROOT / "common" / "python" / "io"
    if str(io_dir) not in sys.path:
        sys.path.insert(0, str(io_dir))
    from read_GEOSldas import read_tilecoord

    return np.asarray(read_tilecoord(str(path))["area"], dtype="float64")


class MonthlySeriesLoader:
    """Provide paired tile and area-weighted domain monthly series.

    The loader reads the monthly products already defined by the M21C synthesis
    workflow. It does not aggregate raw GEOS LDAS output or reconstruct any
    increment variable.
    """

    def __init__(
        self,
        data_dir: str | Path = DEFAULT_DATA_DIR,
        *,
        input_contract: str | Path = DEFAULT_INPUT_CONTRACT,
        variable_selection: str | Path = DEFAULT_VARIABLE_SELECTION,
        tile_area: xr.DataArray | np.ndarray | None = None,
    ) -> None:
        self.data_dir = Path(data_dir)
        self.input_contract_path = Path(input_contract)
        self.variable_selection_path = Path(variable_selection)
        self.contract = _load_json(self.input_contract_path)
        self.dataset_contracts = self.contract["datasets"]
        self.selection = load_variable_selection(self.variable_selection_path)
        self._datasets: dict[str, xr.Dataset] = {}
        self._mask_dataset: xr.Dataset | None = None

        self.expected_time = pd.date_range(
            self.contract["record_start"], self.contract["record_end"], freq="MS"
        )
        if len(self.expected_time) != int(self.contract["n_months"]):
            raise ValueError("Input contract month count does not match its date range")

        ol_land = self._open_dataset("ol_land")
        da_land = self._open_dataset("da_land")
        if not np.array_equal(ol_land.time.values, da_land.time.values):
            raise ValueError("OL and DA land time coordinates differ")

        self.n_tile = int(ol_land.sizes["tile"])
        self.tile = xr.DataArray(np.arange(self.n_tile), dims="tile", name="tile")
        self.lat = self._tile_coordinate(ol_land, "lat")
        self.lon = self._tile_coordinate(ol_land, "lon")

        if tile_area is None:
            tilecoord_path = self.data_dir / self.contract["tilecoord_file"]
            if not tilecoord_path.exists():
                raise FileNotFoundError(f"Missing tile-coordinate file: {tilecoord_path}")
            area_values = read_tile_area(tilecoord_path)
        else:
            area_values = np.asarray(tile_area, dtype="float64")
        if area_values.ndim != 1 or area_values.size != self.n_tile:
            raise ValueError(
                f"Tile area has shape {area_values.shape}; expected ({self.n_tile},)"
            )
        if not np.any(np.isfinite(area_values) & (area_values > 0)):
            raise ValueError("Tile area contains no finite positive weights")
        area_sum = float(np.nansum(np.where(area_values > 0, area_values, np.nan)))
        area_sum_min = self.contract.get("tile_area_sum_min")
        area_sum_max = self.contract.get("tile_area_sum_max")
        if area_sum_min is not None and area_sum < float(area_sum_min):
            raise ValueError(
                f"Tile area sum {area_sum:.6g} is below contracted minimum {area_sum_min}"
            )
        if area_sum_max is not None and area_sum > float(area_sum_max):
            raise ValueError(
                f"Tile area sum {area_sum:.6g} is above contracted maximum {area_sum_max}"
            )
        self.tile_area = xr.DataArray(
            area_values,
            dims="tile",
            coords={"tile": self.tile.values},
            name="tile_area",
            attrs={
                "long_name": "GEOS LDAS land-tile area",
                "units": str(self.contract.get("tile_area_units", "km2")),
            },
        ).where(np.isfinite(area_values) & (area_values > 0))
        self.valid_land_mask = (
            np.isfinite(self.lat) & np.isfinite(self.lon) & np.isfinite(self.tile_area)
        ).rename("valid_land")

    def __enter__(self) -> "MonthlySeriesLoader":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()

    def close(self) -> None:
        """Close all cached xarray datasets."""

        for dataset in self._datasets.values():
            dataset.close()
        self._datasets.clear()

    def _tile_coordinate(self, dataset: xr.Dataset, name: str) -> xr.DataArray:
        if name not in dataset:
            raise KeyError(f"{name} is absent from ol_land")
        values = dataset[name]
        if values.dims != ("tile",):
            raise ValueError(f"ol_land/{name} has dimensions {values.dims}, expected ('tile',)")
        return xr.DataArray(
            np.asarray(values.values, dtype="float64"),
            dims="tile",
            coords={"tile": self.tile.values},
            name=name,
            attrs=dict(values.attrs),
        )

    def _open_dataset(self, key: str) -> xr.Dataset:
        if key in self._datasets:
            return self._datasets[key]
        if key not in self.dataset_contracts:
            raise KeyError(f"Unknown dataset key: {key}")

        spec = self.dataset_contracts[key]
        path = self.data_dir / spec["filename"]
        if not path.exists():
            raise FileNotFoundError(f"Missing monthly dataset for {key}: {path}")
        dataset = xr.open_dataset(path, decode_times=True)
        try:
            self._validate_dataset(key, dataset)
        except Exception:
            dataset.close()
            raise
        self._datasets[key] = dataset
        return dataset

    def _validate_dataset(self, key: str, dataset: xr.Dataset) -> None:
        expected_tiles = int(self.contract["n_tiles"])
        if dataset.sizes.get("time") != len(self.expected_time):
            raise ValueError(
                f"{key} has {dataset.sizes.get('time')} months; expected {len(self.expected_time)}"
            )
        if dataset.sizes.get("tile") != expected_tiles:
            raise ValueError(
                f"{key} has {dataset.sizes.get('tile')} tiles; expected {expected_tiles}"
            )
        actual_time = pd.DatetimeIndex(dataset.time.values)
        if not actual_time.equals(self.expected_time):
            raise ValueError(f"{key} does not match the contracted monthly time coordinate")
        source_token = str(self.dataset_contracts[key].get("source_token", ""))
        source_root = str(dataset.attrs.get("source_root", ""))
        if source_token and source_token not in source_root:
            raise ValueError(
                f"{key} source_root does not contain {source_token!r}: {source_root!r}"
            )

    def _source_values(self, dataset_key: str, variable: str) -> xr.DataArray:
        dataset = self._open_dataset(dataset_key)
        if variable not in dataset:
            raise KeyError(f"{variable} is absent from {dataset_key}")
        values = dataset[variable]
        if values.dims != ("time", "tile"):
            raise ValueError(
                f"{dataset_key}/{variable} has dimensions {values.dims}; "
                "expected ('time', 'tile')"
            )
        # Dataset-specific lat/lon auxiliaries can differ at floating precision
        # even when the tile order is identical. Analysis geometry always comes
        # from the audited OL land/tile-coordinate contract.
        return values.reset_coords(drop=True).assign_coords(tile=self.tile.values)

    def _monthly_values(self, dataset_key: str, variable: str) -> xr.DataArray:
        source = self._source_values(dataset_key, variable)
        source_attrs = dict(source.attrs)
        source_units = _normalise_units(source_attrs.get("units", ""))
        if source_units == WATER_RATE_UNITS:
            output = source * (source.time.dt.days_in_month.astype("float64") * 86400.0)
            output.attrs = {
                **source_attrs,
                "units": MONTHLY_WATER_UNITS,
                "source_units": source_units,
                "monthly_transform": "monthly mean rate multiplied by seconds in month",
            }
            return output
        output = source.copy(deep=False)
        output.attrs = {
            **source_attrs,
            "units": source_units,
            "monthly_transform": "none; existing monthly product used directly",
        }
        return output

    def selected_variables(self, phase: int | None = None) -> pd.DataFrame:
        """Return all selected variables, optionally restricted to a phase."""

        selected = self.selection
        if phase is not None:
            selected = selected.loc[selected["phase"].astype(int) == int(phase)]
        return selected.copy()

    def masks(self) -> xr.Dataset:
        """Return the synthesis notebook's static and monthly analysis masks."""

        if self._mask_dataset is not None:
            return self._mask_dataset

        scf_ol = self._source_values("ol_land", "FRLANDSNO")
        scf_da = self._source_values("da_land", "FRLANDSNO")
        swe_ol = self._source_values("ol_land", "SNOMASLAND")
        swe_da = self._source_values("da_land", "SNOMASLAND")
        tsoil_ol = self._source_values("ol_land", "TSOIL1")
        tsoil_da = self._source_values("da_land", "TSOIL1")

        scf_pair_max = xr.concat([scf_ol, scf_da], dim="experiment").max(
            "experiment", skipna=True
        )
        swe_pair_max = xr.concat([swe_ol, swe_da], dim="experiment").max(
            "experiment", skipna=True
        )
        scf_any = scf_pair_max.max("time", skipna=True)
        swe_any = swe_pair_max.max("time", skipna=True)
        jja_scf = scf_pair_max.sel(time=scf_pair_max.time.dt.month.isin([6, 7, 8])).mean(
            "time", skipna=True
        )

        snow_possible = (
            (self.lat > NH_MIN_LAT)
            & (
                (scf_any > SNOW_POSSIBLE_SCF_THRESHOLD)
                | (swe_any > SNOW_POSSIBLE_SWE_THRESHOLD)
            )
            & self.valid_land_mask
        )
        seasonal_snow = snow_possible & (jja_scf < PERMANENT_SNOW_JJA_MEAN_MAX)
        warm_static = (
            (np.abs(self.lat) < 60.0)
            & (scf_any < WARM_SNOWFREE_SCF_MAX)
            & self.valid_land_mask
        )
        warm_snowfree_monthly = (
            (scf_ol < WARM_SNOWFREE_SCF_MAX)
            & (scf_da < WARM_SNOWFREE_SCF_MAX)
            & (tsoil_ol > WARM_SNOWFREE_TSOIL_MIN_K)
            & (tsoil_da > WARM_SNOWFREE_TSOIL_MIN_K)
            & self.valid_land_mask
        )
        locally_snowy_monthly = (
            seasonal_snow
            & (
                (scf_pair_max > SNOW_POSSIBLE_SCF_THRESHOLD)
                | (swe_pair_max > SNOW_POSSIBLE_SWE_THRESHOLD)
            )
        )

        self._mask_dataset = xr.Dataset(
            {
                "valid_land": self.valid_land_mask,
                "snow_possible": snow_possible.rename("snow_possible"),
                "seasonal_snow": seasonal_snow.rename("seasonal_snow"),
                "warm_static": warm_static.rename("warm_static"),
                "warm_snowfree_monthly": warm_snowfree_monthly.rename(
                    "warm_snowfree_monthly"
                ),
                "locally_snowy_monthly": locally_snowy_monthly.rename(
                    "locally_snowy_monthly"
                ),
            },
            attrs={
                "snow_possible_scf_threshold": SNOW_POSSIBLE_SCF_THRESHOLD,
                "snow_possible_swe_threshold_kg_m-2": SNOW_POSSIBLE_SWE_THRESHOLD,
                "permanent_snow_jja_mean_max": PERMANENT_SNOW_JJA_MEAN_MAX,
                "nh_min_lat_degrees": NH_MIN_LAT,
                "warm_snowfree_scf_max": WARM_SNOWFREE_SCF_MAX,
                "warm_snowfree_tsoil_min_K": WARM_SNOWFREE_TSOIL_MIN_K,
            },
        ).assign_coords(lat=self.lat, lon=self.lon, tile_area=self.tile_area)
        return self._mask_dataset

    def mask(self, name: str) -> xr.DataArray:
        """Return one named analysis mask."""

        if name == "valid_land":
            return self.valid_land_mask
        masks = self.masks()
        if name not in masks:
            raise KeyError(f"Unknown mask {name!r}; available masks: {list(masks.data_vars)}")
        return masks[name]

    def load_tile_series(
        self, analysis_variable: str, *, mask: str = "valid_land"
    ) -> xr.Dataset:
        """Return monthly tile fields with paired OL/DA sampling where applicable."""

        if analysis_variable not in self.selection.index:
            raise KeyError(f"Unselected analysis variable: {analysis_variable}")
        row = self.selection.loc[analysis_variable]
        source_variable = str(row["source_variable"])
        ol_key = str(row["ol_dataset"] or "")
        da_key = str(row["da_dataset"] or "")
        analysis_mask = self.mask(mask)
        eligible = analysis_mask.broadcast_like(
            self._source_values(ol_key or da_key, source_variable)
        ).astype(bool).rename("eligible")
        eligible.attrs = {
            "long_name": "months and tiles eligible under the selected analysis mask",
            "units": "1",
            "mask": mask,
        }

        if ol_key and da_key:
            ol = self._monthly_values(ol_key, source_variable)
            da = self._monthly_values(da_key, source_variable)
            if ol.dims != da.dims or not np.array_equal(ol.time.values, da.time.values):
                raise ValueError(f"OL and DA fields do not align for {analysis_variable}")
            ol_units = _normalise_units(ol.attrs.get("units", ""))
            da_units = _normalise_units(da.attrs.get("units", ""))
            if ol_units != da_units:
                raise ValueError(
                    f"OL and DA units differ for {analysis_variable}: {ol_units!r} vs {da_units!r}"
                )
            paired_valid = np.isfinite(ol) & np.isfinite(da) & analysis_mask
            ol = ol.where(paired_valid).rename("ol")
            da = da.where(paired_valid).rename("da")
            delta = (da - ol).rename("delta")
            common_attrs = {
                "units": ol_units,
                "analysis_variable": analysis_variable,
                "source_variable": source_variable,
                "sample_contract": "strict paired finite OL and DA values",
            }
            ol.attrs = {**ol.attrs, **common_attrs, "experiment": "OL"}
            da.attrs = {**da.attrs, **common_attrs, "experiment": "DA"}
            delta.attrs = {
                **common_attrs,
                "long_name": f"paired DA minus OL {source_variable}",
                "experiment": "DA - OL",
                "monthly_transform": ol.attrs["monthly_transform"],
            }
            output = xr.Dataset({"ol": ol, "da": da, "delta": delta, "eligible": eligible})
            paired = True
        elif da_key:
            value = self._monthly_values(da_key, source_variable)
            value = value.where(np.isfinite(value) & analysis_mask).rename("value")
            value.attrs = {
                **value.attrs,
                "analysis_variable": analysis_variable,
                "source_variable": source_variable,
                "experiment": "DA diagnostic",
                "sample_contract": "finite values from the existing monthly DA product",
            }
            output = xr.Dataset({"value": value, "eligible": eligible})
            paired = False
        else:
            raise ValueError(
                f"{analysis_variable} has neither a paired OL/DA source nor a DA source"
            )

        output = output.assign_coords(
            tile=self.tile.values,
            lat=self.lat,
            lon=self.lon,
            tile_area=self.tile_area,
        )
        output.attrs = {
            "analysis_variable": analysis_variable,
            "source_variable": source_variable,
            "kind": str(row["kind"]),
            "phase": int(row["phase"]),
            "primary_estimand": str(row["primary_estimand"]),
            "mask": mask,
            "paired": int(paired),
            "input_contract": str(self.input_contract_path),
            "variable_selection": str(self.variable_selection_path),
        }
        return output

    def domain_mean(
        self, analysis_variable: str, *, mask: str = "valid_land"
    ) -> xr.Dataset:
        """Return area-weighted monthly means and their changing sample support."""

        fields = self.load_tile_series(analysis_variable, mask=mask)
        weights = fields["tile_area"].fillna(0.0)
        series_names = [name for name in ("ol", "da", "delta", "value") if name in fields]
        means: dict[str, xr.DataArray] = {}
        for name in series_names:
            values = fields[name]
            valid_weights = weights.where(np.isfinite(values), 0.0)
            denominator = valid_weights.sum("tile")
            mean = ((values.fillna(0.0) * valid_weights).sum("tile") / denominator).where(
                denominator > 0
            )
            mean.name = name
            mean.attrs = {
                **values.attrs,
                "spatial_aggregation": "GEOS LDAS tile-area-weighted mean",
            }
            means[name] = mean

        reference_name = "delta" if "delta" in fields else "value"
        reference_valid = np.isfinite(fields[reference_name])
        n_tiles = reference_valid.sum("tile").rename("n_tiles")
        area_sum = weights.where(reference_valid, 0.0).sum("tile").rename("area_sum")
        n_tiles.attrs = {"long_name": "number of finite contributing land tiles", "units": "1"}
        area_sum.attrs = {
            "long_name": "area represented by contributing tiles",
            "units": str(fields["tile_area"].attrs.get("units", "km2")),
        }

        output = xr.Dataset({**means, "n_tiles": n_tiles, "area_sum": area_sum})
        output.attrs = {
            **fields.attrs,
            "spatial_aggregation": "GEOS LDAS tile-area-weighted mean",
            "support_reference": reference_name,
        }
        return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--input-contract", type=Path, default=DEFAULT_INPUT_CONTRACT)
    parser.add_argument("--variable-selection", type=Path, default=DEFAULT_VARIABLE_SELECTION)
    parser.add_argument("--variable", default="SFMC")
    parser.add_argument("--mask", default="valid_land")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with MonthlySeriesLoader(
        args.data_dir,
        input_contract=args.input_contract,
        variable_selection=args.variable_selection,
    ) as loader:
        result = loader.domain_mean(args.variable, mask=args.mask).load()
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
