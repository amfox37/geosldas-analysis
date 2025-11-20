#!/usr/bin/env python
# coding: utf-8

import xarray as xr
import numpy as np

# Load the mask file
ds = xr.open_dataset('../test_data/clsm/subsurface_scattering_ASCAT_ERA5_Land.nc')

cold_mask = ds['cold_mask']
wet_mask = ds['wet_mask']
veg_mask = ds['veg_mask']
subsurface_mask = ds['subsurface_mask']
asc_lon = ds['lon']
asc_lat = ds['lat']

# Make a mask that combines the wet and subsurface masks
mask_combined = np.where(wet_mask == 1, 1, subsurface_mask)

dlon = 5.0  # Longitude grid spacing
dlat = 5.0  # Latitude grid spacing

ll_lon = -180  # Lower left corner longitude
ll_lat = -90  # Lower left corner latitude

# Create a regular grid
lon = np.arange(ll_lon, ll_lon + 360, dlon)
lat = np.arange(ll_lat, ll_lat + 180, dlat)
c_mask1 = np.empty([len(lat), len(lon)])
c_mask1.fill(np.nan)

loncnt = np.empty(len(lon)*len(lat))
latcnt = np.empty(len(lon)*len(lat))
maskcnt = np.empty(len(lon)*len(lat))
maskcnt.fill(np.nan)
cnt = 0

# Loop through the regular grid and assign the mask value. Use the nearest neighbour method with xarray
for i in range(len(lon)):
    # Print lon[i] to see progress if lon[i] is a multiple of 10
    print(lon[i])
    for j in range(len(lat)):
        
        distances = np.sqrt((asc_lon - lon[i])**2 + (asc_lat - lat[j])**2)
        closest_index = distances.argmin()
        # If distances at closest_index is more than 0.5 degrees, then the nearest neighbour is too far away so set to nan
        if distances[closest_index] > 0.14:
            c_mask1[j, i] = np.nan
        else:
            c_mask1[j, i] = mask_combined[closest_index]

        loncnt[cnt] = lon[i]
        latcnt[cnt] = lat[j]
        maskcnt[cnt] = c_mask1[j, i]

        cnt += 1

# Save the mask to a npz file
np.savez(f'ascat_combined_mask_{dlon}_degree.npz', lon=lon, lat=lat, c_mask1=c_mask1)

# Save the mask to a netCDF file
mask_out = np.where(np.isnan(c_mask1), -128, c_mask1)
mask_out = mask_out.astype(np.int8)
# mask_out = mask

ds = xr.Dataset({'mask': (['lat', 'lon'], mask_out)},
                coords={'lat': (['lat'], lat),
                        'lon': (['lon'], lon)})

# Add attributes to the 'mask' variable
ds['mask'].attrs['standard_name'] = 'subsurface_mask'
ds['mask'].attrs['long_name'] = 'Mask accounting for subsurface scattering'
ds['mask'].attrs['units'] = 'boolean'
ds['mask'].encoding['_FillValue'] = -128

ds.to_netcdf(f'ascat_combined_mask_{dlon}.nc')

# Create a new DataArray for ll_lon
ll_lon_da = xr.DataArray(ll_lon, name='ll_lon')
ll_lon_da.attrs['standard_name'] = 'longitude of lower left corner'
ll_lon_da.attrs['long_name'] = 'longitude of lower left corner'
ll_lon_da.attrs['units'] = 'degrees_east'
ll_lon_da.attrs['axis'] = 'X'

# Add ll_lon_da to the dataset
ds['ll_lon'] = ll_lon_da

# Write the dataset to a netCDF file
ds.to_netcdf(f'ascat_combined_mask_{dlon}.nc')

# Repeat for ll_lat
ll_lat_da = xr.DataArray(ll_lat, name='ll_lat')
ll_lat_da.attrs['standard_name'] = 'latitude of lower left corner'
ll_lat_da.attrs['long_name'] = 'latitude of lower left corner'
ll_lat_da.attrs['units'] = 'degrees_north'
ll_lat_da.attrs['axis'] = 'Y'

ds['ll_lat'] = ll_lat_da

ds.to_netcdf(f'ascat_combined_mask_{dlon}.nc')

# Create a new DataArray for dlon
dlon_da = xr.DataArray(dlon, name='d_lon')
dlon_da.attrs['standard_name'] = 'longitude grid spacing'
dlon_da.attrs['long_name'] = 'longitude grid spacing'
dlon_da.attrs['units'] = 'degrees'
dlon_da.attrs['axis'] = 'X'

# Add dlon_da to the dataset
ds['d_lon'] = dlon_da

# Write the dataset to a netCDF file
ds.to_netcdf(f'ascat_combined_mask_{dlon}.nc')

# Repeat for dlat
dlat_da = xr.DataArray(dlat, name='d_lat')
dlat_da.attrs['standard_name'] = 'latitude grid spacing'
dlat_da.attrs['long_name'] = 'latitude grid spacing'
dlat_da.attrs['units'] = 'degrees'
dlat_da.attrs['axis'] = 'Y'

ds['d_lat'] = dlat_da

ds.to_netcdf(f'ascat_combined_mask_{dlon}.nc')