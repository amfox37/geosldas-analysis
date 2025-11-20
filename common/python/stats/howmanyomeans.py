import xarray as xr

# Open both files
ds1 = xr.open_dataset('/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/CYGNSS_Experiments/OLv8_M36_cd/OLv8_M36_cd/output/SMAP_EASEv2_M36_GLOBAL/stats/z_score_clim_quarter_degree/M36_zscore_stats_2018_doy213_2024_doy181_W_75d_Nmin_20_sp_ALL_all_pentads.nc4')
ds2 = xr.open_dataset('/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/CYGNSS_Experiments/OLv8_M36_cd/OLv8_M36_cd/output/SMAP_EASEv2_M36_GLOBAL/stats/z_score_clim_quarter_degree/M36_zscore_stats_2018_doy213_2024_doy181_W_75d_Nmin_5_sp_ALL_all_pentads.nc4')

# Count non-NaN values in o_mean for each file
count1 = ds1['o_mean'].count().item()  # .item() to get as Python int
count2 = ds2['o_mean'].count().item()

print(f"Non-NaN o_mean values in Nmin_20: {count1}")
print(f"Non-NaN o_mean values in Nmin_5: {count2}")

# Compare the counts, give percent difference
percent_difference = ((count2 - count1) / count1) * 100 if count1 != 0 else float('inf')
print(f"Percent difference in counts: {percent_difference:.2f}%")