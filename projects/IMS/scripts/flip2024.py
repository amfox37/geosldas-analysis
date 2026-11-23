"""Flip the 2024 IMS categorical field vertically (raw[::-1]) so it matches the 2000-2023 orientation.
lat/lon are left untouched (they describe the true grid; the 2024+ DATA array is what is upside down).
Usage: python3 -I flip2024.py <in.nc4> <out.nc4>"""
import sys, numpy as np, netCDF4 as nc
src, dst = sys.argv[1], sys.argv[2]
with nc.Dataset(src) as s, nc.Dataset(dst, "w") as d:
    for k, n in s.dimensions.items():
        d.createDimension(k, None if n.isunlimited() else len(n))
    for k, v in s.variables.items():
        o = d.createVariable(k, v.dtype, v.dimensions, zlib=True, fill_value=getattr(v, "_FillValue", None))
        o.setncatts({a: v.getncattr(a) for a in v.ncattrs() if a != "_FillValue"})
        if k == "ims_snowcover":
            for t in range(v.shape[0]):
                o[t] = np.asarray(v[t])[::-1]
        else:
            o[:] = v[:]
    d.setncatts({a: s.getncattr(a) for a in s.ncattrs()})
    d.history = "2024 ims_snowcover flipped vertically (raw[::-1]) to match 2000-2023 orientation; lat/lon unchanged"
