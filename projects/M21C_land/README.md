# M21C_land

Diagnostics for the M21C land DA configuration on the **cubed-sphere CF0360** tiling
(about 1.19 M land tiles), from the `M21C_testing` experiments on Discover
(`/gpfsm/dnb06/projects/p284/M21C_testing` = `/discover/nobackup/projects/land_da/M21C_testing`).

This is separate from `M21C_ls`, which covers the M36 Land Sweeper manuscript.

## Environment

Runs on Discover, where the experiment output lives. Use the GEOSpyD Python 3.13 environment:

```
/usr/local/other/GEOSpyD/25.3.1-0/2025-10-07/envs/py3.13/bin/python3
```

It has netCDF4, pandas, scipy, cartopy and jupyter, and needs no `g5_modules`.

## Notebooks

- `notebooks/cf0360_grid_and_tiles.ipynb`
  - Explains the three spatial layers: cubed-sphere grid cells (6 x 360 x 360), Pfafstetter catchments, and
    tiles (cell x catchment). Figures come from the BCS 30-arcsec rasters.
  - Shows how obs are assigned to tiles, and checks the rule against ObsFcstAna `tilenum`. For a cube-sphere
    tile space, GEOSldas uses a lat/lon "pert grid" (4N x 3N; `DE1440xPE1080` for C360). It picks the tile
    whose centre of mass is nearest (L1 distance in degrees) among the tiles whose centres fall in the obs's
    pert-grid cell. The obs is therefore not always inside its tile.
- `notebooks/cf0360_allsp_may2019_ofa_maps.ipynb`
  - ObsFcstAna maps and daily time series for the May 2019 all-species DA run
    (`M21C_test_CF0360_allsp_superob025ml_OL6sc_xc03125`) vs the 6-member open loop
    (`M21C_test_CF0360_OL6_superob025ml_clim`).
  - Each sensor is mapped on **its own obs grid**, so every coloured cell is one real obs footprint and grey
    gaps are real: EASE2 M36 for SMOS/SMAP/CYGNSS, 0.25 deg for H SAF super-obs, tile space for MODIS SCF.
    Robinson projection, north of 60S.
  - Obs are matched between runs, and the DA run's scaled obs are used as O for both (the OL ran with
    `scale = .false.`).
  - **Temporary:** uses `NMIN = 10` obs per cell instead of the repo default of 20, because one month is
    too short for 20. The override is specific to this notebook.
  - Caches go to `M21C_testing/analysis_M21C_land/`, outside the repo. Figures go to `figures/`
    (`*.png` is git-ignored).

- `notebooks/cf0360_allsp_may2019_states.ipynb`
  - Model states SFMC and RZMC for the same two runs. allsp vs OL6 uses the daily 0.5 deg `tavg24_2d_lnd_Nx` files,
    because OL6 wrote nothing else: daily global and latitude-band means, and maps of the monthly-mean and
    31 May differences.
  - allsp only: analysis increments (ANA - FCST from `inst3_1d_lndfcstana_Nt`) summed over the month in tile
    space, averaged to cube cells with `frac_cell` weights, and drawn on the true cube-cell shapes.

- `notebooks/cf0360_allsp_may2019_suspicious_increments.ipynb`
  - Flags excessive analysis increments per land tile and cycle, and maps the worst tile per cube cell.
    Groups: A = large compared with the forecast ensemble spread (z > 3 plus a minimum size; TSURF flags that come
    with a snow update are kept separate); B = physically large (|dW| > 100 mm, |dRZMC| > 0.05, |dSFMC| > 0.1,
    TC > 5 K); C = root zone pinned at saturation with collapsed spread for 5+ days after DA wetting;
    D = assimilated obs with |normalized innovation| > 5. Includes a top-events table with the obs assimilated nearby.

## Shared code used

- `projects/ascat_da/scripts/compare_thin_superob_ofa.py`: `species_groups()` maps the ObsFcstAna
  species index to a sensor group from the run's `LDASsa_SPECIAL_inputs_ensupd.nml`.
