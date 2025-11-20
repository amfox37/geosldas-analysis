#!/usr/bin/env python3
"""
Simple script to plot observation maps from GEOSldas ObsFcstAna binary files

Usage:
    python plot_obs_maps.py

Modify the configuration section below to point to your experiment data.
"""

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
import os
import sys

# Add path to shared utilities
# sys.path.append('../../shared/python/')
from read_GEOSldas import read_ObsFcstAna, read_tilecoord, read_obs_param

# ============================================================================
# CONFIGURATION SECTION - MODIFY THESE PATHS FOR YOUR DATA
# ============================================================================

# Experiment configuration
EXPDIR = "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/"  # e.g., "/discover/nobackup/username/experiments/"
EXPID = "LS_OLv8_M36"                    # e.g., "SMAP_L4_SM_test"
DOMAIN = "SMAP_EASEv2_M36_GLOBAL"               # or whatever domain you're using

# Time period to plot (3-hourly files)
START_DATE = datetime(2020, 5, 1, 0, 0)         # Year, Month, Day, Hour, Minute
END_DATE = datetime(2020, 5, 3, 0, 0)           # Plot 2 days worth of data
DA_DT = 3 * 3600  # Data assimilation time step in seconds (3 hours)

# Observation parameter file (pick a representative date)
OBSPARAM_TIME = "20200101_0000"  # YYYYMMDD_HHMM format

# Output directory for plots
OUTPUT_DIR = "./obs_plots/"

# ============================================================================

def setup_paths():
    """Setup file paths based on configuration"""
    
    # Output directory for observation parameter and tile coordinate files
    rc_out_dir = f"{EXPDIR}/{EXPID}/output/{DOMAIN}/rc_out/"
    
    # Find the obsparam file (usually in a Y*/M*/ subdirectory)
    obsparam_file = None
    tilecoord_file = None
    
    # Look for files in rc_out subdirectories
    for root, dirs, files in os.walk(rc_out_dir):
        for file in files:
            if file.endswith(f".ldas_obsparam.{OBSPARAM_TIME}z.txt"):
                obsparam_file = os.path.join(root, file)
            if file.endswith(".ldas_tilecoord.bin"):
                tilecoord_file = os.path.join(root, file)
        if obsparam_file and tilecoord_file:
            break
    
    if not obsparam_file:
        print(f"Error: Could not find obsparam file for time {OBSPARAM_TIME}")
        print(f"Searched in: {rc_out_dir}")
        sys.exit(1)
        
    if not tilecoord_file:
        print(f"Error: Could not find tilecoord file")
        print(f"Searched in: {rc_out_dir}")
        sys.exit(1)
    
    print(f"Using obsparam file: {obsparam_file}")
    print(f"Using tilecoord file: {tilecoord_file}")
    
    return obsparam_file, tilecoord_file

def get_obs_file_path(date_time):
    """Generate the ObsFcstAna file path for a given datetime"""
    year_str = date_time.strftime('%Y')
    month_str = date_time.strftime('%m')
    datetime_str = date_time.strftime('%Y%m%d_%H%M')
    
    obs_file = f"{EXPDIR}/{EXPID}/output/{DOMAIN}/ana/ens_avg/Y{year_str}/M{month_str}/{EXPID}.ens_avg.ldas_ObsFcstAna.{datetime_str}z.bin"
    
    return obs_file

def create_map_plot(obs_data, tilecoord, species_info, date_time, variable='obs_obs'):
    """Create a map plot of observations for a specific species and time"""
    
    # Get unique species in the data
    unique_species = np.unique(obs_data['obs_species'])
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    for species_id in unique_species:
        # Skip if no data for this species
        species_mask = obs_data['obs_species'] == species_id
        if not np.any(species_mask):
            continue
            
        # Find species info
        species_name = f"Species_{int(species_id)}"
        for sp in species_info:
            if int(sp['species']) == species_id:
                species_name = sp['descr'].replace(' ', '_').replace('/', '_')
                break
        
        # Extract data for this species
        lons = obs_data['obs_lon'][species_mask]
        lats = obs_data['obs_lat'][species_mask]
        values = obs_data[variable][species_mask]
        
        # Remove NaN values
        valid_mask = ~np.isnan(values)
        lons = lons[valid_mask]
        lats = lats[valid_mask]
        values = values[valid_mask]
        
        if len(lons) == 0:
            print(f"No valid data for species {species_name} at {date_time}")
            continue
        
        # Create the plot
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        # Add map features
        ax.add_feature(cfeature.COASTLINE, alpha=0.8)
        ax.add_feature(cfeature.BORDERS, alpha=0.5)
        ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
        
        # Plot observations as scatter plot
        scatter = ax.scatter(lons, lats, c=values, s=1, alpha=0.7, 
                           transform=ccrs.PlateCarree(), cmap='viridis')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.05)
        
        # Set colorbar label based on variable
        var_labels = {
            'obs_obs': 'Observation Value',
            'obs_fcst': 'Forecast Value', 
            'obs_ana': 'Analysis Value',
            'obs_obsvar': 'Observation Variance',
            'obs_fcstvar': 'Forecast Variance',
            'obs_anavar': 'Analysis Variance'
        }
        cbar.set_label(var_labels.get(variable, variable))
        
        # Set extent to show global or regional view
        ax.set_global()
        
        # Add gridlines
        ax.gridlines(draw_labels=True, alpha=0.5)
        
        # Set title
        title = f"{species_name}\n{variable} - {date_time.strftime('%Y-%m-%d %H:%M')} UTC\nN_obs = {len(values)}"
        plt.title(title)
        
        # Save the plot
        time_str = date_time.strftime('%Y%m%d_%H%M')
        filename = f"{OUTPUT_DIR}/obs_map_{species_name}_{variable}_{time_str}.png"
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Saved: {filename}")

def plot_obs_minus_forecast(obs_data, tilecoord, species_info, date_time):
    """Create a map plot of observation minus forecast (O-F) innovations"""
    
    # Get unique species in the data  
    unique_species = np.unique(obs_data['obs_species'])
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    for species_id in unique_species:
        # Skip if no data for this species
        species_mask = obs_data['obs_species'] == species_id
        if not np.any(species_mask):
            continue
            
        # Find species info
        species_name = f"Species_{int(species_id)}"
        for sp in species_info:
            if int(sp['species']) == species_id:
                species_name = sp['descr'].replace(' ', '_').replace('/', '_')
                break
        
        # Extract data for this species
        lons = obs_data['obs_lon'][species_mask]
        lats = obs_data['obs_lat'][species_mask]
        obs_vals = obs_data['obs_obs'][species_mask]
        fcst_vals = obs_data['obs_fcst'][species_mask]
        
        # Compute O-F
        omf_vals = obs_vals - fcst_vals
        
        # Remove NaN values
        valid_mask = ~(np.isnan(omf_vals) | np.isnan(obs_vals) | np.isnan(fcst_vals))
        lons = lons[valid_mask]
        lats = lats[valid_mask]
        omf_vals = omf_vals[valid_mask]
        
        if len(lons) == 0:
            print(f"No valid O-F data for species {species_name} at {date_time}")
            continue
        
        # Create the plot
        fig = plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        # Add map features
        ax.add_feature(cfeature.COASTLINE, alpha=0.8)
        ax.add_feature(cfeature.BORDERS, alpha=0.5)
        ax.add_feature(cfeature.LAND, color='lightgray', alpha=0.5)
        ax.add_feature(cfeature.OCEAN, color='lightblue', alpha=0.3)
        
        # Use diverging colormap for O-F (centered on zero)
        vmax = np.percentile(np.abs(omf_vals), 95)  # Use 95th percentile to avoid outliers
        scatter = ax.scatter(lons, lats, c=omf_vals, s=20, alpha=0.7,
                           transform=ccrs.PlateCarree(), 
                           cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, pad=0.05)
        cbar.set_label('Observation - Forecast')
        
        # Set extent to show global or regional view
        ax.set_global()
        
        # Add gridlines
        ax.gridlines(draw_labels=True, alpha=0.5)
        
        # Set title with statistics
        mean_omf = np.mean(omf_vals)
        std_omf = np.std(omf_vals)
        title = f"{species_name} - Observation minus Forecast\n{date_time.strftime('%Y-%m-%d %H:%M')} UTC\n"
        title += f"N_obs = {len(omf_vals)}, Mean = {mean_omf:.3f}, Std = {std_omf:.3f}"
        plt.title(title)
        
        # Save the plot
        time_str = date_time.strftime('%Y%m%d_%H%M')
        filename = f"{OUTPUT_DIR}/obs_map_{species_name}_O-F_{time_str}.png"
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"Saved: {filename}")

def main():
    """Main function to process and plot observation data"""
    
    print("GEOSldas Observation Mapping Tool")
    print("=" * 50)
    
    # Setup file paths
    obsparam_file, tilecoord_file = setup_paths()
    
    # Read observation parameters and tile coordinates
    print("Reading observation parameters and tile coordinates...")
    obs_param = read_obs_param(obsparam_file)
    tilecoord = read_tilecoord(tilecoord_file)
    
    print(f"Found {len(obs_param)} observation species:")
    for i, sp in enumerate(obs_param):
        print(f"  {int(sp['species'])}: {sp['descr']}")
    
    # Process each time step
    current_time = START_DATE
    file_count = 0
    
    while current_time < END_DATE:
        print(f"\nProcessing time: {current_time}")
        
        # Get observation file path
        obs_file = get_obs_file_path(current_time)
        
        if os.path.exists(obs_file):
            print(f"Reading: {obs_file}")
            
            # Read observation data
            obs_data = read_ObsFcstAna(obs_file)
            
            # Check if we have observations
            if len(obs_data['obs_lon']) > 0:
                print(f"Found {len(obs_data['obs_lon'])} observations")
                
                # Create maps for different variables
                print("Creating observation maps...")
                create_map_plot(obs_data, tilecoord, obs_param, current_time, 'obs_obs')
                
                # print("Creating O-F innovation maps...")
                # plot_obs_minus_forecast(obs_data, tilecoord, obs_param, current_time)
                
                # Optionally create other plots
                # create_map_plot(obs_data, tilecoord, obs_param, current_time, 'obs_fcst')
                # create_map_plot(obs_data, tilecoord, obs_param, current_time, 'obs_ana')
                
                file_count += 1
            else:
                print("No observations found in file")
        else:
            print(f"File not found: {obs_file}")
        
        # Move to next time step
        current_time += timedelta(seconds=DA_DT)
    
    print(f"\nProcessed {file_count} observation files")
    print(f"Maps saved to: {OUTPUT_DIR}")
    print("Done!")

if __name__ == "__main__":
    main()