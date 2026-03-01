#!/usr/bin/env python3
"""
Download daily IMS ASCII snow-cover files (.asc.gz) and write one NetCDF file.

This script is intentionally simple but robust:
- configurable year/paths via CLI
- retry+backoff for HTTP downloads
- no fixed ASCII header length assumption
- grid shape validation before write
- basic NetCDF metadata for easier downstream use
"""

from __future__ import annotations

import argparse
import calendar
import gzip
import time
from datetime import date, datetime, timedelta
from pathlib import Path

import numpy as np
import requests
from netCDF4 import Dataset


# Root folder for IMS files. Per-year data are typically under:
#   <BASE_ROOT>/<grid>/<year>/
DEFAULT_BASE_ROOT = "https://noaadata.apps.nsidc.org/NOAA/G02156"
DEFAULT_GRID = "24km"
DEFAULT_VERSION = "v1.3"

# Grid size for 24 km IMS.
DEFAULT_NX = 1024
DEFAULT_NY = 1024


def parse_args() -> argparse.Namespace:
    """Parse command-line settings."""
    now_year = datetime.utcnow().year

    ap = argparse.ArgumentParser(
        description="Download IMS ASCII .asc.gz files and convert to NetCDF."
    )
    ap.add_argument("--year", type=int, default=now_year, help="Year to process (default: current UTC year)")
    ap.add_argument("--base-root", default=DEFAULT_BASE_ROOT, help="Base IMS root URL")
    ap.add_argument("--grid", default=DEFAULT_GRID, help="IMS grid tag in filenames/URL (default: 24km)")
    ap.add_argument("--version", default=DEFAULT_VERSION, help="IMS file version tag (default: v1.3)")
    ap.add_argument("--latlon-file", default="ims_24km_latlon.nc4", help="Path to lat/lon NetCDF file")
    ap.add_argument("--out-file", default=None, help="Output NetCDF path (default: ims_snowcover_<grid>_<year>.nc4)")
    ap.add_argument("--cache-dir", default="gz", help="Directory to store downloaded .gz files")
    ap.add_argument("--nx", type=int, default=DEFAULT_NX, help="Grid x size (default: 1024)")
    ap.add_argument("--ny", type=int, default=DEFAULT_NY, help="Grid y size (default: 1024)")

    # Optional day subset can be useful for testing or restart workflows.
    ap.add_argument("--start-doy", type=int, default=1, help="Start day-of-year (default: 1)")
    ap.add_argument("--end-doy", type=int, default=None, help="End day-of-year (default: last day of year)")

    # Download robustness controls.
    ap.add_argument("--timeout", type=int, default=60, help="HTTP timeout in seconds (default: 60)")
    ap.add_argument("--retries", type=int, default=3, help="HTTP retries per file (default: 3)")
    ap.add_argument("--backoff", type=float, default=1.0, help="Backoff seconds base (default: 1.0)")

    # If enabled, missing file days are left as NaN in output instead of aborting.
    ap.add_argument("--allow-missing-days", action="store_true", help="Continue when a day is missing or unreadable")
    ap.add_argument("--overwrite-download", action="store_true", help="Re-download .gz files even if already cached")
    ap.add_argument("--overwrite-output", action="store_true", help="Overwrite output NetCDF if it exists")

    return ap.parse_args()


def days_in_year(year: int) -> int:
    """Return 365 or 366 for the target year."""
    return 366 if calendar.isleap(year) else 365


def build_filename(year: int, doy: int, grid: str, version: str) -> str:
    """Build expected IMS ASCII gzip filename."""
    # IMS daily ASCII naming convention includes cycle (00UTC) in filename.
    # Example: ims2024001_00UTC_24km_v1.3.asc.gz
    return f"ims{year}{doy:03d}_00UTC_{grid}_{version}.asc.gz"


def day_from_year_doy(year: int, doy: int) -> date:
    """Convert year + day-of-year into calendar date."""
    return (datetime(year, 1, 1) + timedelta(days=doy - 1)).date()


def download_file_with_retries(
    session: requests.Session,
    url: str,
    local_path: Path,
    timeout: int,
    retries: int,
    backoff: float,
) -> None:
    """
    Download one file to local path with retry/backoff.

    We raise on persistent failures so the caller can decide whether to skip
    (allow-missing-days) or abort.
    """
    last_err: Exception | None = None

    for attempt in range(retries):
        try:
            r = session.get(url, timeout=timeout)
            if r.status_code == 200:
                local_path.write_bytes(r.content)
                return
            if r.status_code == 404:
                raise FileNotFoundError(f"Remote file not found: {url}")
            raise RuntimeError(f"HTTP {r.status_code} while downloading: {url}")
        except Exception as exc:
            last_err = exc
            if attempt == retries - 1:
                break
            sleep_s = backoff * (2 ** attempt)
            time.sleep(sleep_s)

    assert last_err is not None
    raise last_err


def parse_ims_ascii_gz(path: Path, nx: int, ny: int) -> np.ndarray:
    """
    Parse one IMS .asc.gz file into a (ny, nx) float32 array.

    Important: we do NOT assume a fixed header length.

    We support common IMS row layouts:
    1) whitespace-delimited rows with nx numeric tokens
    2) fixed-width rows with no whitespace (or compacted whitespace), where
       each row is represented by nx values packed as:
       - 1-char codes (len == nx), or
       - constant-width fields (len == nx * field_width)
    """
    def _try_parse_row(line: str) -> np.ndarray | None:
        s = line.strip()
        if not s:
            return None

        # Case 1: whitespace-delimited tokens.
        toks = s.split()
        if len(toks) == nx:
            try:
                return np.asarray(toks, dtype=np.float32)
            except ValueError:
                return None

        # Case 2: fixed-width rows (possibly with internal spaces removed).
        # Some IMS ASCII files store each row as contiguous characters without
        # separators, so a split() call returns only one token.
        compact = "".join(toks) if len(toks) > 1 else s
        if not compact:
            return None

        # 2a) One-character-per-cell format (most common for categorical grids).
        if len(compact) == nx:
            try:
                # Fast digit parse for simple code grids.
                return np.fromiter((float(ch) for ch in compact), dtype=np.float32, count=nx)
            except ValueError:
                return None

        # 2b) Generic fixed-width field parse.
        if len(compact) % nx == 0:
            width = len(compact) // nx
            if 1 < width <= 6:
                try:
                    return np.asarray(
                        [float(compact[i * width : (i + 1) * width]) for i in range(nx)],
                        dtype=np.float32,
                    )
                except ValueError:
                    return None

        return None

    rows: list[np.ndarray] = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            arr = _try_parse_row(line)
            if arr is None:
                # Header/metadata or non-numeric line.
                continue

            rows.append(arr)
            if len(rows) == ny:
                # We found a full grid; ignore any trailing content.
                break

    if len(rows) != ny:
        raise ValueError(f"Parsed {len(rows)} data rows from {path.name}, expected {ny}")

    grid = np.vstack(rows).astype(np.float32, copy=False)
    if grid.shape != (ny, nx):
        raise ValueError(f"Grid shape {grid.shape} does not match expected {(ny, nx)}")

    return grid


def read_latlon(latlon_file: Path, nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Read lat/lon arrays and validate shape.

    We expect both arrays to be 2D and match (ny, nx).
    """
    if not latlon_file.exists():
        raise FileNotFoundError(f"Lat/lon file not found: {latlon_file}")

    with Dataset(latlon_file, "r") as ds:
        if "lat" not in ds.variables or "lon" not in ds.variables:
            raise KeyError(f"Expected variables 'lat' and 'lon' in {latlon_file}")

        lat = np.asarray(ds["lat"][:], dtype=np.float64)
        lon = np.asarray(ds["lon"][:], dtype=np.float64)

    if lat.shape != (ny, nx):
        raise ValueError(f"lat shape {lat.shape} != expected {(ny, nx)}")
    if lon.shape != (ny, nx):
        raise ValueError(f"lon shape {lon.shape} != expected {(ny, nx)}")

    return lat, lon


def init_output_netcdf(
    out_file: Path,
    n_time: int,
    ny: int,
    nx: int,
    year: int,
    lat: np.ndarray,
    lon: np.ndarray,
) -> tuple[Dataset, np.ndarray]:
    """
    Create output NetCDF and return dataset handle + writable snow variable.

    We write each day directly to disk so memory usage stays modest.
    """
    ds = Dataset(out_file, "w", format="NETCDF4")

    ds.createDimension("time", n_time)
    ds.createDimension("y", ny)
    ds.createDimension("x", nx)

    time_var = ds.createVariable("time", "f8", ("time",))
    time_var.long_name = "time"
    time_var.standard_name = "time"
    time_var.units = f"days since {year}-01-01 00:00:00"
    time_var.calendar = "standard"

    lat_var = ds.createVariable("lat", "f8", ("y", "x"))
    lon_var = ds.createVariable("lon", "f8", ("y", "x"))
    lat_var.long_name = "latitude"
    lat_var.units = "degrees_north"
    lon_var.long_name = "longitude"
    lon_var.units = "degrees_east"
    lat_var[:, :] = lat
    lon_var[:, :] = lon

    # Store raw IMS values as float32.
    # We use NaN as fill for missing/unavailable days.
    snow_var = ds.createVariable(
        "ims_snowcover",
        "f4",
        ("time", "y", "x"),
        zlib=True,
        complevel=4,
        fill_value=np.float32(np.nan),
    )
    snow_var.long_name = "IMS daily snow-cover code (raw ASCII values)"
    snow_var.units = "1"
    snow_var.comment = (
        "Raw values from IMS ASCII grids. Interpretation of specific code values "
        "should follow IMS/NSIDC documentation for the selected product version."
    )

    ds.title = "IMS daily snow-cover from ASCII .asc.gz files"
    ds.source = "NOAA/NSIDC IMS G02156 (downloaded by script)"
    ds.history = f"Created {datetime.utcnow().isoformat()}Z"
    ds.Conventions = "CF-1.8"

    return ds, snow_var


def main() -> None:
    args = parse_args()

    year = int(args.year)
    nx = int(args.nx)
    ny = int(args.ny)

    n_year_days = days_in_year(year)
    start_doy = int(args.start_doy)
    end_doy = int(args.end_doy) if args.end_doy is not None else n_year_days

    if start_doy < 1 or start_doy > n_year_days:
        raise ValueError(f"--start-doy must be in [1, {n_year_days}]")
    if end_doy < start_doy or end_doy > n_year_days:
        raise ValueError(f"--end-doy must be in [{start_doy}, {n_year_days}]")

    out_file = Path(args.out_file) if args.out_file else Path(f"ims_snowcover_{args.grid}_{year}.nc4")
    latlon_file = Path(args.latlon_file)
    cache_dir = Path(args.cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    if out_file.exists() and not args.overwrite_output:
        raise FileExistsError(f"Output file exists (use --overwrite-output): {out_file}")

    base_url = f"{args.base_root.rstrip('/')}/{args.grid}/{year}"

    print(f"Year: {year}")
    print(f"DOY range: {start_doy}..{end_doy}")
    print(f"Base URL: {base_url}")
    print(f"Lat/Lon file: {latlon_file}")
    print(f"Output file: {out_file}")
    print(f"Cache dir: {cache_dir}")
    print(f"Grid shape: ny={ny}, nx={nx}")

    # Load static geolocation once.
    print("Reading lat/lon...")
    lat, lon = read_latlon(latlon_file, nx=nx, ny=ny)

    # Initialize output NetCDF and write slices day-by-day.
    n_time = end_doy - start_doy + 1
    ds_out, snow_var = init_output_netcdf(
        out_file=out_file,
        n_time=n_time,
        ny=ny,
        nx=nx,
        year=year,
        lat=lat,
        lon=lon,
    )

    # Use one session for connection reuse.
    session = requests.Session()

    n_ok = 0
    n_missing = 0
    n_failed = 0

    try:
        # Iterate over requested day-of-year window.
        for ti, doy in enumerate(range(start_doy, end_doy + 1)):
            this_date = day_from_year_doy(year, doy)
            filename = build_filename(year, doy, args.grid, args.version)
            url = f"{base_url}/{filename}"
            local_path = cache_dir / filename

            print(f"[{ti + 1:04d}/{n_time:04d}] {this_date.isoformat()}  {filename}")

            # Keep a continuous time axis for all requested days.
            # 0 corresponds to Jan-01 of target year.
            ds_out["time"][ti] = float((this_date - date(year, 1, 1)).days)

            # Download if needed.
            if args.overwrite_download or not local_path.exists():
                try:
                    download_file_with_retries(
                        session=session,
                        url=url,
                        local_path=local_path,
                        timeout=int(args.timeout),
                        retries=int(args.retries),
                        backoff=float(args.backoff),
                    )
                except FileNotFoundError as exc:
                    if args.allow_missing_days:
                        print(f"  missing (404), writing NaN slice: {exc}")
                        snow_var[ti, :, :] = np.nan
                        n_missing += 1
                        continue
                    raise

            # Parse the gz ASCII and validate shape.
            try:
                grid = parse_ims_ascii_gz(local_path, nx=nx, ny=ny)
            except Exception as exc:
                if args.allow_missing_days:
                    print(f"  parse/read failed, writing NaN slice: {exc}")
                    snow_var[ti, :, :] = np.nan
                    n_failed += 1
                    continue
                raise

            # Write this day directly to NetCDF.
            snow_var[ti, :, :] = grid
            n_ok += 1

    finally:
        ds_out.close()
        session.close()

    print("Finished.")
    print(f"Output: {out_file}")
    print(f"Days written: {n_ok}")
    print(f"Days missing (404, skipped with NaN): {n_missing}")
    print(f"Days failed parse/read (skipped with NaN): {n_failed}")


if __name__ == "__main__":
    main()
