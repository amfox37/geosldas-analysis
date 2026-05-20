#!/usr/bin/env python3
"""
Subset daily LDAS tavg24 NetCDF files to required variables and copy ObsFcstAna
binary files into a flat date-organized tree with exp_run_alias naming.
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
from pathlib import Path

from netCDF4 import Dataset


DEFAULT_VARS = [
    "SFMC",            # surface soil moisture
    "RZMC",            # root-zone soil moisture
    "TSURFLAND",       # surface temperature
    "TSOIL1",          # soil temperature level 1
    "PRECTOTCORRLAND", # precipitation
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Subset daily tavg24 NetCDF to required vars and copy ObsFcstAna bins."
    )
    p.add_argument("--exp-path", required=True, help="Experiment root path (e.g., .../CYGNSS_Experiments/DAv8_M36_cd)")
    p.add_argument("--exp-run", required=True, help="Experiment run name (e.g., DAv8_M36_cd)")
    p.add_argument("--exp-run-alias", required=True, help="Alias name for outputs (e.g., CYG_DA)")
    p.add_argument("--out-root", required=True, help="Output root (e.g., .../output)")
    p.add_argument(
        "--vars",
        default=",".join(DEFAULT_VARS),
        help="Comma-separated list of NetCDF variables to keep",
    )
    p.add_argument(
        "--hours",
        default="00,03,06,09,12,15,18,21",
        help="Comma-separated hours for ObsFcstAna (HH, e.g., 00,03,...,21)",
    )
    p.add_argument("--dry-run", action="store_true", help="Print actions without writing files")
    return p.parse_args()


def ensure_dir(path: Path, dry_run: bool) -> None:
    if dry_run:
        return
    path.mkdir(parents=True, exist_ok=True)


def copy_variable(src_var, dst, name: str) -> None:
    # Preserve dtype, dimensions, and compression/chunking where possible
    kwargs = {}
    try:
        filters = src_var.filters()
        if filters:
            kwargs.update(filters)
    except Exception:
        pass

    fill_value = getattr(src_var, "_FillValue", None)
    if fill_value is not None:
        kwargs["fill_value"] = fill_value

    if "chunksizes" not in kwargs and src_var.chunking() not in (None, "contiguous"):
        try:
            kwargs["chunksizes"] = src_var.chunking()
        except Exception:
            pass

    dst_var = dst.createVariable(name, src_var.datatype, src_var.dimensions, **kwargs)
    dst_var.setncatts({k: src_var.getncattr(k) for k in src_var.ncattrs()})
    dst_var[:] = src_var[:]


def subset_nc(src: Path, dst: Path, keep_vars: list[str], dry_run: bool) -> None:
    if dry_run:
        print(f"[dry-run] subset {src} -> {dst}")
        return

    with Dataset(src, "r") as ds_in:
        with Dataset(dst, "w", format=ds_in.file_format) as ds_out:
            # copy global attrs
            ds_out.setncatts({k: ds_in.getncattr(k) for k in ds_in.ncattrs()})

            # copy dimensions
            for name, dim in ds_in.dimensions.items():
                ds_out.createDimension(name, (len(dim) if not dim.isunlimited() else None))

            # copy coordinate variables used by kept vars
            coord_vars = set()
            for vname in keep_vars:
                if vname in ds_in.variables:
                    coord_vars.update(ds_in.variables[vname].dimensions)

            # Always keep time/lat/lon if present (safe for plotting)
            for extra in ["time", "lat", "lon", "longitude", "latitude"]:
                if extra in ds_in.variables:
                    coord_vars.add(extra)

            # copy coords first
            for vname in coord_vars:
                if vname in ds_in.variables:
                    copy_variable(ds_in.variables[vname], ds_out, vname)

            # copy requested variables
            for vname in keep_vars:
                if vname not in ds_in.variables:
                    raise KeyError(f"Variable not found in {src}: {vname}")
                copy_variable(ds_in.variables[vname], ds_out, vname)


def main() -> None:
    args = parse_args()
    exp_path = Path(args.exp_path)
    exp_run = args.exp_run
    exp_run_alias = args.exp_run_alias
    out_root = Path(args.out_root)
    keep_vars = [v.strip() for v in args.vars.split(",") if v.strip()]
    hours = [h.strip() for h in args.hours.split(",") if h.strip()]

    cat_root = exp_path / exp_run / "output" / "SMAP_EASEv2_M36_GLOBAL" / "cat" / "ens_avg"
    ana_root = exp_path / exp_run / "output" / "SMAP_EASEv2_M36_GLOBAL" / "ana" / "ens_avg"

    pattern = str(cat_root / "Y????" / "M??" / f"{exp_run}.tavg24_1d_lnd_Nt.*_1200z.nc4")
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No files found with pattern: {pattern}")

    for src_str in files:
        src = Path(src_str)
        # Extract date from filename: ...tavg24_1d_lnd_Nt.YYYYMMDD_1200z.nc4
        parts = src.name.split(".")
        if len(parts) < 3:
            raise ValueError(f"Unexpected filename format: {src.name}")
        date_token = parts[-2].split("_")[0]  # YYYYMMDD
        yyyy, mm, dd = date_token[:4], date_token[4:6], date_token[6:8]

        out_dir = out_root / exp_run_alias / yyyy / mm / dd
        ensure_dir(out_dir, args.dry_run)

        out_name = src.name.replace(exp_run, exp_run_alias)
        dst_nc = out_dir / out_name
        subset_nc(src, dst_nc, keep_vars, args.dry_run)

        # Copy 8 daily ObsFcstAna files
        for hh in hours:
            obs_name = f"{exp_run}.ens_avg.ldas_ObsFcstAna.{yyyy}{mm}{dd}_{hh}00z.bin"
            obs_src = ana_root / f"Y{yyyy}" / f"M{mm}" / obs_name
            obs_dst = out_dir / obs_name.replace(exp_run, exp_run_alias)
            if args.dry_run:
                print(f"[dry-run] copy {obs_src} -> {obs_dst}")
                continue
            if not obs_src.exists():
                raise FileNotFoundError(f"Missing ObsFcstAna file: {obs_src}")
            shutil.copy2(obs_src, obs_dst)


if __name__ == "__main__":
    main()
