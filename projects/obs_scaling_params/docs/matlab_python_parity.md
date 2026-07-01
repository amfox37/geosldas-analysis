# MATLAB/Python Parity Notes

This project was migrated from `projects/matlab2python/obs_scaling_params`.
The MATLAB implementation remains in `matlab_reference/`.

## Local Comparison

Compared local files:

```text
/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/stats/M36_python_dedup_zscore_stats_2007_doy152_2024_doy151_W_75d_Nmin_20_sp_ALL_all_pentads.nc4
/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/stats/M36_dedup_zscore_stats_2007_doy152_2024_doy151_W_75d_Nmin_20_sp_ALL_all_pentads.nc4
```

Summary:

- Dimensions match exactly: `pentad=73`, `lon=1440`, `lat=720`.
- Coordinates/timing match exactly: `pentad`, `start_time`, `end_time`, `ll_lon`, `ll_lat`, `d_lon`, `d_lat`.
- Main 3D fields have about 9.4 million finite values.
- Python has 10 net extra finite points in the main 3D fields.
- Common finite points: 9,402,705.
- Finite-mask mismatch: 35 Python-only, 25 MATLAB-only.

Numeric differences on common finite points:

| Variable | Max abs diff | Mean abs diff |
|---|---:|---:|
| `o_mean` | 0.0299 | 8.1e-6 |
| `o_std` | 0.0515 | 5.0e-6 |
| `m_mean` | 0.0077 | 2.2e-6 |
| `m_std` | 0.0063 | 1.3e-6 |
| `n_data` | 3 observations | 0.025 observations |
| `m_min` | 0.0152 | 7.3e-8 |
| `m_max` | 0.0300 | 6.8e-7 |

Interpretation: the Python port is functionally faithful but not bit-for-bit
identical. The small residual differences are localized and consistent with
duplicate/edge handling and floating/order details, not an algorithm mismatch.

## MATLAB-Alignment Details Preserved in Python

The Python core keeps MATLAB behavior for:

- `t0_assim = t0_assim % dt_assim`
- treating `-9999` observation/model sentinels as missing before statistics
- rolling 75-day window accumulation
- combined-species handling
- de-duplication by rounded tile/species/lon/lat/obs/year key, keeping the
  later recent occurrence
- output fields and dimensions compatible with existing scaling files
