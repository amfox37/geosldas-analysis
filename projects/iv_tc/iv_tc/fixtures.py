"""Resolve and copy small Discover fixture inputs for IV/TC development."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
import shutil

from .config import ProductRoots, RunConfig


H121_MANIFESTS = {
    "metop_a": "H121_METOPA.txt",
    "metop_b": "H121_H139_METOPB.txt",
    "metop_c": "H121_H139_METOPC.txt",
}


@dataclass(frozen=True)
class FixtureFile:
    """One source file and its portable fixture destination."""

    label: str
    source: Path
    dest: Path

    @property
    def exists(self) -> bool:
        return self.source.exists()


def parse_date(value: str) -> date:
    """Parse YYYY-MM-DD or YYYYMMDD."""

    for fmt in ("%Y-%m-%d", "%Y%m%d"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            pass
    raise ValueError(f"Date must be YYYY-MM-DD or YYYYMMDD, got {value!r}")


def collect_fixture_files(
    dates: list[date],
    output_root: Path,
    roots: ProductRoots,
    runs: list[RunConfig],
) -> list[FixtureFile]:
    """Return all fixture files known for the requested dates and runs."""

    files: list[FixtureFile] = []
    for day in dates:
        files.extend(_h121_files(day, output_root, roots))
        files.extend(_legacy_ascat_files(day, output_root, roots))
        files.extend(_smosic_files(day, output_root, roots))
        files.extend(_smap_l3_files(day, output_root, roots))
        for run in runs:
            files.extend(_model_files(day, output_root, roots, run))
    return files


def copy_fixture_files(files: list[FixtureFile]) -> tuple[int, int]:
    """Copy existing files. Return (copied, missing)."""

    copied = 0
    missing = 0
    for item in files:
        if not item.exists:
            missing += 1
            continue
        item.dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(item.source, item.dest)
        copied += 1
    return copied, missing


def print_fixture_plan(files: list[FixtureFile]) -> None:
    """Print a dry-run friendly fixture plan."""

    for item in files:
        status = "FOUND" if item.exists else "MISS "
        print(f"[{status}] {item.label}")
        print(f"  src: {item.source}")
        print(f"  dst: {item.dest}")


def _date_parts(day: date) -> tuple[str, str, str, str]:
    yyyy = f"{day.year:04d}"
    mm = f"{day.month:02d}"
    dd = f"{day.day:02d}"
    yyyymmdd = f"{yyyy}{mm}{dd}"
    return yyyy, mm, dd, yyyymmdd


def _h121_files(day: date, output_root: Path, roots: ProductRoots) -> list[FixtureFile]:
    yyyy, mm, dd, _ = _date_parts(day)
    files: list[FixtureFile] = []

    for platform, manifest_name in H121_MANIFESTS.items():
        manifest = roots.ascat_cdr_root / "flists" / f"Y{yyyy}" / f"M{mm}" / f"D{dd}" / manifest_name
        manifest_dest = output_root / "ASCAT_SSM_CDR" / "flists" / f"Y{yyyy}" / f"M{mm}" / f"D{dd}" / manifest_name
        files.append(FixtureFile(f"H121 manifest {platform} {day}", manifest, manifest_dest))

        if not manifest.exists():
            continue

        raw_root = roots.ascat_cdr_root / "H121" / platform / f"Y{yyyy}" / f"M{mm}"
        raw_dest_root = output_root / "ASCAT_SSM_CDR" / "H121" / platform / f"Y{yyyy}" / f"M{mm}"
        for raw_name in _read_manifest_names(manifest):
            files.append(
                FixtureFile(
                    f"H121 raw {platform} {day}",
                    raw_root / raw_name,
                    raw_dest_root / raw_name,
                )
            )

    return files


def _read_manifest_names(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text().splitlines():
        name = line.strip()
        if name:
            names.append(name)
    return names


def _legacy_ascat_files(day: date, output_root: Path, roots: ProductRoots) -> list[FixtureFile]:
    yyyy, mm, _, yyyymmdd = _date_parts(day)
    name = f"ASCAT_HSAF_H119_SM_{yyyymmdd}_AD.mat"
    source = roots.legacy_ascat_root / "H119_H120_processed" / f"Y{yyyy}" / f"M{mm}" / name
    dest = output_root / "ASCAT_HSAF" / "H119_H120_processed" / f"Y{yyyy}" / f"M{mm}" / name
    return [FixtureFile(f"legacy ASCAT {day}", source, dest)]


def _smosic_files(day: date, output_root: Path, roots: ProductRoots) -> list[FixtureFile]:
    _, _, _, yyyymmdd = _date_parts(day)
    name = f"smos_ic_sm_m36_{yyyymmdd}.nc"
    source = roots.smosic_root / name
    dest = output_root / "SMOS_IC" / "preprocessed_m36_daily" / name
    return [FixtureFile(f"SMOS-IC {day}", source, dest)]


def _smap_l3_files(day: date, output_root: Path, roots: ProductRoots) -> list[FixtureFile]:
    yyyy, _, _, yyyymmdd = _date_parts(day)
    flat_dir = roots.smap_l3_root / f"Y{yyyy}"
    pattern = f"SMAP_L3_SM_P_{yyyymmdd}_R19240_*.h5"
    matches = sorted(path for path in flat_dir.glob(pattern) if not path.name.startswith("._"))

    if not matches:
        missing_source = flat_dir / pattern
        missing_dest = output_root / "SPL3SMP_v009" / f"Y{yyyy}" / pattern
        return [FixtureFile(f"SMAP L3 flat {day}", missing_source, missing_dest)]

    files: list[FixtureFile] = []
    for source in matches:
        dest = output_root / "SPL3SMP_v009" / f"Y{yyyy}" / source.name
        files.append(FixtureFile(f"SMAP L3 flat {day}", source, dest))
    return files


def _model_files(day: date, output_root: Path, roots: ProductRoots, run: RunConfig) -> list[FixtureFile]:
    # GEOSldas names output files after EXP_ID, which by convention equals the
    # run root's directory basename -- not the CLI label used to organize
    # fixture output. Using the label here silently produces "found" nothing
    # for any run whose label doesn't happen to match its EXP_ID exactly.
    exp_id = run.root.name
    yyyy, mm, _, yyyymmdd = _date_parts(day)
    base = run.root / "output" / roots.domain
    model_dir = base / "cat" / "ens_avg" / f"Y{yyyy}" / f"M{mm}"
    candidates = [
        model_dir / f"{exp_id}{roots.collection}{yyyymmdd}_1200z.nc4",
        model_dir / f"{exp_id}{roots.collection}{yyyymmdd}.nc4",
    ]
    model_source = _first_existing_or_first(candidates)
    model_dest = (
        output_root
        / "runs"
        / run.name
        / "output"
        / roots.domain
        / "cat"
        / "ens_avg"
        / f"Y{yyyy}"
        / f"M{mm}"
        / model_source.name
    )

    tilecoord = base / "rc_out" / f"{exp_id}.ldas_tilecoord.bin"
    tilecoord_dest = (
        output_root
        / "runs"
        / run.name
        / "output"
        / roots.domain
        / "rc_out"
        / tilecoord.name
    )

    return [
        FixtureFile(f"model daily {run.name} {day}", model_source, model_dest),
        FixtureFile(f"tilecoord {run.name}", tilecoord, tilecoord_dest),
    ]


def _first_existing_or_first(paths: list[Path]) -> Path:
    for path in paths:
        if path.exists():
            return path
    return paths[0]

