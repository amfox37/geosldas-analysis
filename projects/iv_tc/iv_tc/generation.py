"""Daily pair generation for IV/TC step-2 inputs."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path

import numpy as np

from .config import ProductRoots, RunConfig
from .pairs import DailyPair
from .readers import (
    M36_NX,
    read_ascat_h119_h120_model_pair,
    read_ascat_h121_model_pair,
    read_cygnss_l3_model_pair,
    read_smap_l3_model_pair,
    read_smosic_model_pair,
)


DEFAULT_SENSORS = ("smosic", "smap_l3", "cygnss_l3", "ascat_h121", "ascat_h119_h120")

H121_MANIFESTS = {
    "metop_a": "H121_METOPA.txt",
    "metop_b": "H121_H139_METOPB.txt",
    "metop_c": "H121_H139_METOPC.txt",
}

SENSOR_ALIASES = {
    "smosic": "smosic",
    "smos-ic": "smosic",
    "smos_ic": "smosic",
    "smap_l3": "smap_l3",
    "smap-l3": "smap_l3",
    "smap": "smap_l3",
    "cygnss_l3": "cygnss_l3",
    "cygnss-l3": "cygnss_l3",
    "cygnss": "cygnss_l3",
    "cygl3": "cygnss_l3",
    "ascat_h121": "ascat_h121",
    "h121": "ascat_h121",
    "h121_h139": "ascat_h121",
    "ascat_h121_h139": "ascat_h121",
    "ascat_h119_h120": "ascat_h119_h120",
    "h119_h120": "ascat_h119_h120",
    "h119-h120": "ascat_h119_h120",
}


@dataclass(frozen=True)
class PairInputPaths:
    """Resolved inputs for one sensor/run/day pair."""

    observation_paths: tuple[Path, ...]
    model_path: Path
    tilecoord_path: Path
    auxiliary_path: Path | None = None


@dataclass(frozen=True)
class PairGenerationResult:
    """One output decision from daily pair generation."""

    date: date
    sensor: str
    run: str
    output_path: Path
    status: str
    n_points: int = 0
    message: str = ""


def date_range(start: date, end: date) -> list[date]:
    """Return inclusive daily dates from start through end."""

    if end < start:
        raise ValueError(f"end date {end} is before start date {start}")

    days: list[date] = []
    day = start
    while day <= end:
        days.append(day)
        day += timedelta(days=1)
    return days


def canonical_sensor(sensor: str) -> str:
    """Return the internal sensor token used in pair output paths."""

    key = sensor.strip().lower()
    if key == "ascat":
        raise ValueError("ASCAT is ambiguous here; use ascat_h121 or ascat_h119_h120")
    try:
        return SENSOR_ALIASES[key]
    except KeyError as exc:
        allowed = ", ".join(DEFAULT_SENSORS)
        raise ValueError(f"Unsupported sensor {sensor!r}; expected one of {allowed}") from exc


def pair_output_path(output_root: Path | str, sensor: str, run_name: str, day: date) -> Path:
    """Return the step-2 pair path for one sensor/run/day."""

    return Path(output_root) / "step2_pairs" / canonical_sensor(sensor) / run_name / f"{day:%Y%m%d}.npz"


def pair_npz_is_valid(path: Path | str) -> bool:
    """Return True when path has the compact daily pair arrays."""

    path = Path(path)
    if not path.exists():
        return False

    try:
        with np.load(path) as z:
            if not all(key in z for key in ("idx0", "sm_obs", "sm_mod")):
                return False
            idx0 = np.asarray(z["idx0"])
            sm_obs = np.asarray(z["sm_obs"])
            sm_mod = np.asarray(z["sm_mod"])
    except Exception:
        return False

    if idx0.ndim != 1 or sm_obs.ndim != 1 or sm_mod.ndim != 1:
        return False
    return bool(idx0.size == sm_obs.size == sm_mod.size)


def save_daily_pair_npz(pair: DailyPair, path: Path | str) -> None:
    """Write a DailyPair using the step-2 npz format consumed downstream."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        idx0=np.asarray(pair.idx, dtype=np.int32),
        sm_obs=np.asarray(pair.obs, dtype=np.float32),
        sm_mod=np.asarray(pair.model, dtype=np.float32),
        date=np.asarray(pair.date.strftime("%Y%m%d")),
        sensor=np.asarray(pair.sensor),
        run=np.asarray(pair.run),
        obs_units=np.asarray(pair.obs_units),
        model_units=np.asarray(pair.model_units),
    )


def load_daily_pair_npz(path: Path | str) -> DailyPair:
    """Load a step-2 pair npz as a DailyPair record."""

    path = Path(path)
    with np.load(path) as z:
        idx = np.asarray(z["idx0"], dtype=np.int64)
        obs = np.asarray(z["sm_obs"], dtype=np.float64)
        model = np.asarray(z["sm_mod"], dtype=np.float64)
        day = _parse_yyyymmdd(str(z["date"])) if "date" in z else _parse_yyyymmdd(path.stem)
        sensor = str(z["sensor"]) if "sensor" in z else ""
        run = str(z["run"]) if "run" in z else ""
        obs_units = str(z["obs_units"]) if "obs_units" in z else ""
        model_units = str(z["model_units"]) if "model_units" in z else ""

    return DailyPair(
        date=day,
        sensor=sensor,
        run=run,
        idx=idx,
        obs=obs,
        model=model,
        obs_units=obs_units,
        model_units=model_units,
    )


def resolve_daily_pair_inputs(
    sensor: str,
    day: date,
    run: RunConfig,
    roots: ProductRoots,
) -> PairInputPaths:
    """Resolve all input paths needed for one pair."""

    sensor = canonical_sensor(sensor)
    model_path = resolve_model_path(day, run, roots)
    tilecoord_path = resolve_tilecoord_path(run, roots)

    if sensor == "smosic":
        obs_path = _require_existing(roots.smosic_root / f"smos_ic_sm_m36_{day:%Y%m%d}.nc")
        return PairInputPaths((obs_path,), model_path, tilecoord_path)

    if sensor == "smap_l3":
        obs_path = _resolve_smap_l3_path(day, roots)
        return PairInputPaths((obs_path,), model_path, tilecoord_path)

    if sensor == "cygnss_l3":
        obs_path = _resolve_cygnss_l3_path(day, roots)
        return PairInputPaths((obs_path,), model_path, tilecoord_path)

    if sensor == "ascat_h121":
        obs_paths = _resolve_h121_files(day, roots)
        return PairInputPaths(tuple(obs_paths), model_path, tilecoord_path)

    if sensor == "ascat_h119_h120":
        mat_path = _require_existing(
            roots.ascat_h119_h120_root
            / "H119_H120_processed"
            / f"Y{day.year:04d}"
            / f"M{day.month:02d}"
            / f"ASCAT_HSAF_H119_SM_{day:%Y%m%d}_AD.mat"
        )
        auxiliary_path = _require_existing(
            roots.ascat_h119_h120_root / "Auxiliary" / "TUW_WARP5_grid_info_2_2.nc"
        )
        return PairInputPaths((mat_path,), model_path, tilecoord_path, auxiliary_path)

    raise AssertionError(f"Unhandled sensor {sensor!r}")


def resolve_model_path(day: date, run: RunConfig, roots: ProductRoots) -> Path:
    """Resolve a GEOSldas daily tile-space model file for one run/day."""

    exp_id = run.root.name
    base = run.root / "output" / roots.domain
    yyyy = f"{day.year:04d}"
    mm = f"{day.month:02d}"
    yyyymmdd = f"{day:%Y%m%d}"
    candidates: list[Path] = []
    for member in ("ens_avg", "ens0000"):
        model_dir = base / "cat" / member / f"Y{yyyy}" / f"M{mm}"
        candidates.extend(
            [
                model_dir / f"{exp_id}{roots.collection}{yyyymmdd}_1200z.nc4",
                model_dir / f"{exp_id}{roots.collection}{yyyymmdd}.nc4",
            ]
        )
    try:
        return _first_existing(candidates)
    except FileNotFoundError:
        pass

    glob_matches: list[Path] = []
    for member in ("ens_avg", "ens0000"):
        model_dir = base / "cat" / member / f"Y{yyyy}" / f"M{mm}"
        glob_matches.extend(sorted(model_dir.glob(f"*{roots.collection}{yyyymmdd}_1200z.nc4")))
        glob_matches.extend(sorted(model_dir.glob(f"*{roots.collection}{yyyymmdd}.nc4")))
    return _only_match(glob_matches, candidates[0])


def resolve_tilecoord_path(run: RunConfig, roots: ProductRoots) -> Path:
    """Resolve a GEOSldas tilecoord file for one run."""

    exp_id = run.root.name
    rc_out = run.root / "output" / roots.domain / "rc_out"
    expected = rc_out / f"{exp_id}.ldas_tilecoord.bin"
    if expected.exists():
        return expected
    return _only_match(sorted(rc_out.glob("*.ldas_tilecoord.bin")), expected)


def read_daily_pair(
    sensor: str,
    day: date,
    run: RunConfig,
    roots: ProductRoots,
    model_variable: str = "SFMC",
    h119_h120_method: str = "linear",
    h119_h120_fill_nearest: bool = False,
    nx: int = M36_NX,
) -> DailyPair:
    """Read one daily obs/model pair after resolving product paths."""

    sensor = canonical_sensor(sensor)
    paths = resolve_daily_pair_inputs(sensor, day, run, roots)

    if sensor == "smosic":
        return read_smosic_model_pair(
            paths.observation_paths[0],
            paths.model_path,
            paths.tilecoord_path,
            run=run.name,
            model_variable=model_variable,
            nx=nx,
        )

    if sensor == "smap_l3":
        return read_smap_l3_model_pair(
            paths.observation_paths[0],
            paths.model_path,
            paths.tilecoord_path,
            run=run.name,
            model_variable=model_variable,
            nx=nx,
        )

    if sensor == "cygnss_l3":
        return read_cygnss_l3_model_pair(
            paths.observation_paths[0],
            paths.model_path,
            paths.tilecoord_path,
            run=run.name,
            model_variable=model_variable,
            nx=nx,
        )

    if sensor == "ascat_h121":
        return read_ascat_h121_model_pair(
            list(paths.observation_paths),
            day,
            paths.model_path,
            paths.tilecoord_path,
            run=run.name,
            model_variable=model_variable,
            nx=nx,
        )

    if sensor == "ascat_h119_h120":
        if paths.auxiliary_path is None:
            raise FileNotFoundError("ASCAT H119/H120 auxiliary path was not resolved")
        return read_ascat_h119_h120_model_pair(
            paths.observation_paths[0],
            paths.auxiliary_path,
            paths.model_path,
            paths.tilecoord_path,
            run=run.name,
            method=h119_h120_method,
            fill_nearest=h119_h120_fill_nearest,
            model_variable=model_variable,
            nx=nx,
        )

    raise AssertionError(f"Unhandled sensor {sensor!r}")


def write_daily_pair(
    sensor: str,
    day: date,
    run: RunConfig,
    roots: ProductRoots,
    output_root: Path | str,
    skip_existing: bool = True,
    model_variable: str = "SFMC",
    h119_h120_method: str = "linear",
    h119_h120_fill_nearest: bool = False,
    nx: int = M36_NX,
) -> PairGenerationResult:
    """Generate one daily pair file, returning a status record."""

    sensor = canonical_sensor(sensor)
    output_path = pair_output_path(output_root, sensor, run.name, day)
    if skip_existing and pair_npz_is_valid(output_path):
        return PairGenerationResult(day, sensor, run.name, output_path, "exists")

    try:
        pair = read_daily_pair(
            sensor=sensor,
            day=day,
            run=run,
            roots=roots,
            model_variable=model_variable,
            h119_h120_method=h119_h120_method,
            h119_h120_fill_nearest=h119_h120_fill_nearest,
            nx=nx,
        )
    except FileNotFoundError as exc:
        return PairGenerationResult(day, sensor, run.name, output_path, "missing", message=str(exc))
    except (RuntimeError, OSError) as exc:
        return PairGenerationResult(day, sensor, run.name, output_path, "failed", message=str(exc))

    if pair.idx.size == 0:
        return PairGenerationResult(day, sensor, run.name, output_path, "empty")

    save_daily_pair_npz(pair, output_path)
    return PairGenerationResult(day, sensor, run.name, output_path, "written", n_points=int(pair.idx.size))


def generate_daily_pairs(
    dates: list[date],
    sensors: list[str],
    runs: list[RunConfig],
    roots: ProductRoots,
    output_root: Path | str,
    skip_existing: bool = True,
    model_variable: str = "SFMC",
    h119_h120_method: str = "linear",
    h119_h120_fill_nearest: bool = False,
    nx: int = M36_NX,
) -> list[PairGenerationResult]:
    """Loop over dates, sensors, and runs and write daily pair npz files."""

    results: list[PairGenerationResult] = []
    for run in runs:
        for sensor in sensors:
            sensor_token = canonical_sensor(sensor)
            for day in dates:
                results.append(
                    write_daily_pair(
                        sensor=sensor_token,
                        day=day,
                        run=run,
                        roots=roots,
                        output_root=output_root,
                        skip_existing=skip_existing,
                        model_variable=model_variable,
                        h119_h120_method=h119_h120_method,
                        h119_h120_fill_nearest=h119_h120_fill_nearest,
                        nx=nx,
                    )
                )
    return results


def _resolve_smap_l3_path(day: date, roots: ProductRoots) -> Path:
    flat_dir = roots.smap_l3_root / f"Y{day.year:04d}"
    pattern = f"SMAP_L3_SM_P_{day:%Y%m%d}_R19240_*.h5"
    matches = sorted(path for path in flat_dir.glob(pattern) if not path.name.startswith("._"))
    if not matches:
        raise FileNotFoundError(str(flat_dir / pattern))
    return matches[0]


def _resolve_cygnss_l3_path(day: date, roots: ProductRoots) -> Path:
    name = (
        f"cyg.ddmi.s{day:%Y%m%d}-030000-e{day:%Y%m%d}-210000."
        "l3.grid-soil-moisture-36km.a32.d33.nc"
    )
    return _require_existing(roots.cygnss_l3_root / f"Y{day.year:04d}" / f"M{day.month:02d}" / name)


def _resolve_h121_files(day: date, roots: ProductRoots) -> list[Path]:
    yyyy = f"{day.year:04d}"
    mm = f"{day.month:02d}"
    dd = f"{day.day:02d}"
    files: list[Path] = []
    missing_manifests: list[Path] = []
    missing_raw: list[Path] = []

    for platform, manifest_name in H121_MANIFESTS.items():
        manifest = roots.ascat_h121_root / "flists" / f"Y{yyyy}" / f"M{mm}" / f"D{dd}" / manifest_name
        if not manifest.exists():
            missing_manifests.append(manifest)
            continue

        raw_root = roots.ascat_h121_root / "H121" / platform / f"Y{yyyy}" / f"M{mm}"
        for raw_name in _read_manifest_names(manifest):
            raw = raw_root / raw_name
            if raw.exists():
                files.append(raw)
            else:
                missing_raw.append(raw)

    if missing_raw:
        raise FileNotFoundError(str(missing_raw[0]))
    if not files:
        if missing_manifests:
            raise FileNotFoundError(str(missing_manifests[0]))
        raise FileNotFoundError(str(roots.ascat_h121_root / "flists" / f"Y{yyyy}" / f"M{mm}" / f"D{dd}"))

    return files


def _read_manifest_names(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text().splitlines():
        name = line.strip()
        if name:
            names.append(name)
    return names


def _first_existing(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    raise FileNotFoundError(str(paths[0]))


def _only_match(paths: list[Path], fallback: Path) -> Path:
    if len(paths) == 1:
        return paths[0]
    if not paths:
        raise FileNotFoundError(str(fallback))
    names = ", ".join(str(path) for path in paths[:5])
    raise RuntimeError(f"Multiple candidate files found; use a real run root or a narrower fixture root: {names}")


def _require_existing(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(str(path))
    return path


def _parse_yyyymmdd(value: str) -> date:
    value = value.strip()
    return date(int(value[:4]), int(value[4:6]), int(value[6:8]))
