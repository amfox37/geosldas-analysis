"""Plot ISMN soil_moisture station locations with data available in 2020, US Southwest."""
import warnings
warnings.filterwarnings('ignore')

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from ismn.interface import ISMN_Interface

# US Southwest bounding box (AZ, NM, NV, UT, CO, southern CA/west TX)
LON_MIN, LON_MAX = -125.0, -102.0
LAT_MIN, LAT_MAX = 31.0, 42.0

YEAR_START = pd.Timestamp('2020-01-01')
YEAR_END = pd.Timestamp('2020-12-31 23:59:59')

ismn_data = ISMN_Interface(Path('/Users/amfox/Desktop/ISMN_data'), parallel=False)
md = ismn_data.metadata

sm = md[md[('variable', 'val')] == 'soil_moisture'].copy()

overlaps_2020 = (sm[('timerange_from', 'val')] <= YEAR_END) & (sm[('timerange_to', 'val')] >= YEAR_START)
in_bbox = (
    (sm[('longitude', 'val')] >= LON_MIN) & (sm[('longitude', 'val')] <= LON_MAX) &
    (sm[('latitude', 'val')] >= LAT_MIN) & (sm[('latitude', 'val')] <= LAT_MAX)
)

sel = sm[overlaps_2020 & in_bbox]

stations = sel.groupby([('network', 'val'), ('station', 'val')]).agg(
    lon=(('longitude', 'val'), 'first'),
    lat=(('latitude', 'val'), 'first'),
).reset_index()
stations.columns = ['network', 'station', 'lon', 'lat']

print(f'{len(stations)} stations with soil_moisture data overlapping 2020 in US SW bbox')
print(stations['network'].value_counts())

fig, ax = plt.subplots(figsize=(9, 8), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='gray', linewidth=0.5)
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.add_feature(cfeature.BORDERS.with_scale('50m'), linewidth=0.7)

for network, grp in stations.groupby('network'):
    ax.scatter(grp['lon'], grp['lat'], s=20, transform=ccrs.PlateCarree(), label=network)

ax.set_title(f'ISMN soil_moisture sites with data in 2020 — US Southwest (n={len(stations)})')
ax.legend(fontsize=6, loc='lower left', ncol=2, framealpha=0.8)

out_path = Path(__file__).parent.parent / 'figures' / 'us_sw_sites_2020.png'
out_path.parent.mkdir(exist_ok=True)
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f'Saved: {out_path}')
