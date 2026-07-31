"""
Independent (from-scratch) reader for the MATLAB write_seqbin_file.m
big-endian Fortran-sequential scaling-file format, used to cross-check
obs_scaling.seqbin.read_tb_scaling_file() rather than reusing it.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class MatlabScalingFile:
    asc_flag: int
    version: int
    start_time: tuple
    end_time: tuple
    n_grid: int
    n_angle: int
    angles: np.ndarray
    lon: np.ndarray
    lat: np.ndarray
    tile_id: np.ndarray
    data: np.ndarray  # (14, n_grid, n_angle)


def _rec(fp, dtype, count):
    tag = np.frombuffer(fp.read(4), dtype=">i4")[0]
    expected = np.dtype(dtype).itemsize * count
    assert tag == expected, f"record tag {tag} != expected {expected}"
    vals = np.frombuffer(fp.read(expected), dtype=dtype).copy()
    end_tag = np.frombuffer(fp.read(4), dtype=">i4")[0]
    assert end_tag == tag, f"mismatched tags {tag} {end_tag}"
    return vals


def read_matlab_scaling_file(path: str, n_fields: int = 14) -> MatlabScalingFile:
    with open(path, "rb") as fp:
        asc_flag, version = _rec(fp, ">i4", 2)
        st = _rec(fp, ">i4", 5)
        et = _rec(fp, ">i4", 5)
        n_grid, n_angle, n_subtile = _rec(fp, ">i4", 3)
        angles = _rec(fp, ">f4", n_angle)
        lon = _rec(fp, ">f4", n_grid)
        lat = _rec(fp, ">f4", n_grid)
        tile_id_cols = [_rec(fp, ">i4", n_grid) for _ in range(n_subtile)]
        tile_id = tile_id_cols[0]
        data = np.empty((n_fields, n_grid, n_angle), dtype=np.float32)
        for i in range(n_fields):
            for j in range(n_angle):
                data[i, :, j] = _rec(fp, ">f4", n_grid)
        trailing = fp.read()
        assert trailing == b"", f"unexpected trailing bytes: {len(trailing)}"
    return MatlabScalingFile(
        asc_flag=int(asc_flag), version=int(version),
        start_time=tuple(int(x) for x in st), end_time=tuple(int(x) for x in et),
        n_grid=int(n_grid), n_angle=int(n_angle),
        angles=angles, lon=lon, lat=lat, tile_id=tile_id, data=data,
    )
