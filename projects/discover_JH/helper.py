import numpy as np
import struct
import os
from datetime import datetime, timedelta
from pyhdf.SD import SD, SDC


def read_ObsFcstAna(fname, isLDASsa=False):

    # Initialize outputs
    nodata = -9999

    date_time = {
        'year': nodata,
        'month': nodata,
        'day': nodata,
        'hour': nodata,
        'min': nodata,
        'sec': nodata,
        'dofyr': nodata,
        'pentad': nodata
    }

    obs_assim = []
    obs_species = []
    obs_tilenum = []
    obs_lon = []
    obs_lat = []
    obs_obs = []
    obs_obsvar = []
    obs_fcst = []
    obs_fcstvar = []
    obs_ana = []
    obs_anavar = []

    # Determine endianness
    machfmt = '>' if isLDASsa else '<'  # '>' for big-endian, '<' for little-endian

    if os.path.exists(fname):
        # print(f"reading from {fname}")

        with open(fname, 'rb') as ifp:
            # Read N_obs and time stamp entry
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            N_obs = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            year, month, day, hour, minute, second, dofyr, pentad = struct.unpack(f'{machfmt}8i', ifp.read(32))
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            date_time = {
                'year': year,
                'month': month,
                'day': day,
                'hour': hour,
                'min': minute,
                'sec': second,
                'dofyr': dofyr,
                'pentad': pentad
            }

            # Read observation assim flag
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            tmp_data = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}i').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            obs_assim = np.zeros(N_obs, dtype=int)
            obs_assim[tmp_data != 0] = 1

            # Read species information
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_species = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}i').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read tile number information
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_tilenum = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}i').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read longitude
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_lon = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read latitude
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_lat = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_obs = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_obsvar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space model forecast value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_fcst = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space model forecast variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_fcstvar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space analysis value
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_ana = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # Read observation-space analysis variance
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            obs_anavar = np.frombuffer(ifp.read(N_obs * 4), dtype=f'{machfmt}f').copy()
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

            # No-data check
            obs_obsvar[obs_obsvar == nodata] = np.nan
            obs_fcst[obs_fcst == nodata] = np.nan
            obs_fcstvar[obs_fcstvar == nodata] = np.nan
            obs_ana[obs_ana == nodata] = np.nan
            obs_anavar[obs_anavar == nodata] = np.nan

    else:
        print(f"file does not exist: {fname}")

    return {'date_time': date_time, 
            'obs_assim': obs_assim, 
            'obs_species': obs_species, 
            'obs_tilenum': obs_tilenum, 
            'obs_lon': obs_lon, 
            'obs_lat': obs_lat,
            'obs_obs': obs_obs, 
            'obs_obsvar': obs_obsvar, 
            'obs_fcst': obs_fcst, 
            'obs_fcstvar': obs_fcstvar, 
            'obs_ana': obs_ana, 
            'obs_anavar': obs_anavar}

def read_tilecoord(fname):
    int_precision = 'i'
    float_precision = 'f'

    # Determine endianness
    machfmt = '<'  # '>' for big-endian, '<' for little-endian

    print(f"reading from {fname}")

    tile_coord = {}

    with open(fname, 'rb') as ifp:
        fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
        tile_coord['N_tile'] = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
        fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

        Nt = tile_coord['N_tile']

        fields = ['tile_id', 'typ', 'pfaf', 'com_lon', 'com_lat', 'min_lon', 'max_lon',
                      'min_lat', 'max_lat', 'i_indg', 'j_indg', 'frac_cell', 'frac_pfaf',
                      'area', 'elev']

        for field in fields:
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]
            dtype = int_precision if field in ['tile_id', 'typ', 'pfaf', 'i_indg', 'j_indg'] else float_precision
            tile_coord[field] = np.frombuffer(ifp.read(Nt * 4), dtype=f'{machfmt}{dtype}')
            fortran_tag = struct.unpack(f'{machfmt}i', ifp.read(4))[0]

def get_tile_species_obs_values(ofa_data):
    """Get maximum observations per tile, the number of obs, and the number of obs > 0.9, split by species"""
    
    # Get unique values
    unique_tiles, tile_first_idx = np.unique(ofa_data['obs_tilenum'], return_index=True)
    unique_species = np.unique(ofa_data['obs_species'])
    
    # Get coordinates from first occurrences
    tile_lat = ofa_data['obs_lat'][tile_first_idx]
    tile_lon = ofa_data['obs_lon'][tile_first_idx]
    
    # Sort data for efficient processing
    sort_idx = np.argsort(ofa_data['obs_tilenum'])
    sorted_tiles = ofa_data['obs_tilenum'][sort_idx]
    sorted_species = ofa_data['obs_species'][sort_idx]
    sorted_obs = ofa_data['obs_obs'][sort_idx]
    
    # Initialize results dictionary
    results = {
        'tiles': unique_tiles,
        'lat': tile_lat,
        'lon': tile_lon,
        'max_values': {species: np.zeros(len(unique_tiles)) for species in unique_species},
        'num_obs': {species: np.zeros(len(unique_tiles)) for species in unique_species},
        'num_obs_gt_0.9': {species: np.zeros(len(unique_tiles)) for species in unique_species}
    }
    
    # Find split points for tiles
    tile_splits = np.searchsorted(sorted_tiles, unique_tiles)
    tile_splits = np.append(tile_splits, len(sorted_tiles))
    
    # Calculate max obs_obs for each tile and species
    for i in range(len(unique_tiles)):
        tile_data = sorted_obs[tile_splits[i]:tile_splits[i+1]]
        tile_species = sorted_species[tile_splits[i]:tile_splits[i+1]]
        
        for species in unique_species:
            species_mask = tile_species == species
            if np.any(species_mask):
                results['max_values'][species][i] = np.max(tile_data[species_mask])
                results['num_obs'][species][i] = np.sum(species_mask)
                results['num_obs_gt_0.9'][species][i] = np.sum(tile_data[species_mask] > 0.9)
    
    return results            


def read_modis_scf_hdf(fname, lon_min, lon_max, lat_min, lat_max, clear_index=20, snow_spatial=2):
    """Read MODIS Snow Cover Fraction from HDF4 file.
        
    # Usage over CONUS
    lon_min = -125.0
    lon_max = -66.0
    lat_min = 24.0
    lat_max = 50.0
    fname = "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/land_sweeper/LS_OLv8_M36/output/SMAP_EASEv2_M36_GLOBAL/ana/ens_avg/MYD10C1.A2005196.061.hdf"
    lon_out, lat_out, scf_out = read_modis_scf_hdf(fname, lon_min, lon_max, lat_min, lat_max)
"""
    
    # Constants
    CMG_N_lon = 7200
    CMG_N_lat = 3600
    CMG_ll_lon = -180.0
    CMG_ll_lat = -90.0
    CMG_ur_lon = 180.0
    CMG_ur_lat = 90.0
    CMG_dlon = 0.05
    CMG_dlat = 0.05
    
    # QC Parameters
    qc_snow_cover_max = 100
    qc_clear_index_min = clear_index
    qc_snow_spatial_max = snow_spatial
    
    # Calculate array indices for lat/lon bounds
    start_lon = int((lon_min - CMG_ll_lon)/CMG_dlon)
    start_lat = int((CMG_ur_lat - lat_max)/CMG_dlat)
    end_lon = int((lon_max - CMG_ll_lon)/CMG_dlon)
    end_lat = int((CMG_ur_lat - lat_min)/CMG_dlat)
    
    N_lon = end_lon - start_lon + 1
    N_lat = end_lat - start_lat + 1
    
    # Read HDF file
    hdf = SD(fname, SDC.READ)
    
    # Read datasets
    snow_cover = hdf.select('Day_CMG_Snow_Cover')[start_lat:end_lat+1, start_lon:end_lon+1]
    clear_index = hdf.select('Day_CMG_Clear_Index')[start_lat:end_lat+1, start_lon:end_lon+1]
    snow_spatial_qa = hdf.select('Snow_Spatial_QA')[start_lat:end_lat+1, start_lon:end_lon+1]
    
    # Generate lat/lon arrays
    lon_ind = np.arange(N_lon)
    lat_ind = np.arange(N_lat)
    
    lon_c = CMG_ll_lon + 0.5*CMG_dlon + (start_lon + lon_ind)*CMG_dlon
    lat_c = CMG_ur_lat - 0.5*CMG_dlat - (start_lat + lat_ind)*CMG_dlat
    
    # Apply QC and normalize SCF
    valid_mask = ((snow_cover <= qc_snow_cover_max) & 
                 (clear_index > qc_clear_index_min) & 
                 (snow_spatial_qa <= qc_snow_spatial_max))
    
    # Create output arrays
    lon_out = []
    lat_out = []
    scf_out = []
    
    for i in range(N_lon):
        for j in range(N_lat):
            if valid_mask[j,i]:
                scf = float(snow_cover[j,i])/float(clear_index[j,i])
                lon_out.append(lon_c[i])
                lat_out.append(lat_c[j])
                scf_out.append(scf)
    
    hdf.end()
    
    return np.array(lon_out), np.array(lat_out), np.array(scf_out)

def parse_modis_filename(filename):
    """Parse MODIS filename to get date
    
    # Example usage
    filename = "MOD10C1.A2005090.061.hdf"
    date = parse_modis_filename(filename)
    print(f"Date: {date.strftime('%Y-%m-%d')}")
"""
    # Extract date portion (assume fixed format)
    date_str = filename.split('.A')[1].split('.')[0]
    
    # Split into year and doy
    year = int(date_str[:4])
    doy_base = int(date_str[4:6])  # 19
    i = int(date_str[6:])          # 0-9
    
    # Combine base DOY and i
    doy = doy_base * 10 + i
    
    # Convert to datetime
    date = datetime(year, 1, 1) + timedelta(days=doy-1)
    
    return date
