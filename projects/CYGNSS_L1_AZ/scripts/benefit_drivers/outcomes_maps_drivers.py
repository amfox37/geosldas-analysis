#!/usr/bin/env python3
"""4-panel maps of the main drivers (SM-term skill, OL SM variability, elevation, L1 obs count) on AZ tile centres."""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs, cartopy.feature as cf
OUT = Path(__file__).resolve().parents[2] / 'output' / 'benefit_drivers'
d = pd.read_csv(OUT / 'analysis_table.csv')
m = d.in_l1_mask.astype(bool)
panels = [('op_rs_smterm', 'Operator SM-term skill (Spearman)', 'viridis'),
          ('sm_sfmc_std', 'OL surface SM variability (std, m3/m3)', 'viridis'),
          ('st_bcs_mean_elevation_m', 'Elevation (m)', 'terrain'),
          ('n_l1_obs', 'L1 obs count, 2020-2022', 'viridis')]
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from outcomes_common import tile_map
fig, axs = plt.subplots(2, 2, figsize=(11, 9), subplot_kw=dict(projection=ccrs.PlateCarree()))
for ax, (c, t, cm) in zip(axs.ravel(), panels):
    tile_map(ax, d, c, t, cmap=cm, sym=False)
fig.suptitle('Drivers of L1 DA benefit (grey = outside the L1 mask)')
fig.tight_layout(); fig.savefig(OUT / 'analysis_outcomes' / 'fig10_maps_drivers.png', dpi=130)
