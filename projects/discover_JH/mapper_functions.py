import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import xarray as xr
import re
import warnings

from shapely.errors import ShapelyDeprecationWarning

warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)

def format_number(num):
    if abs(num) < 0.01:
        # Use 4 decimal places if the number is less than 0.001
        return '{:.4f}'.format(num)
    elif abs(num) < 1.0:
        # Use 2 decimal places if the number is less than 1
        return '{:.2f}'.format(num)
    else:
        # Use 3 significant figures otherwise
        return '{:.3g}'.format(num)

def plot_global(array, saveflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = 'RdYlBu_r'
            cmap = plt.get_cmap('RdBu', 16)
        else:
            cmap = 'viridis' 
    
    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='lightgray')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    sc = ax.scatter(array[:, 1], array[:, 2],
                    c=array[:, 0], s=1, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
                    transform=ccrs.PlateCarree()) 
    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.12, fraction=0.04)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(units, fontsize=10)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=12)
    ax.text(0.45, -0.1,   'Longitude', fontsize=8, transform=ax.transAxes, ha='left')
    ax.text(-0.08, 0.4, 'Latitude', fontsize=8, transform=ax.transAxes, rotation='vertical', va='bottom')

    if saveflag:
        savename = plot_title+'.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 400)    

     # Show the plot
    plt.show()

def plot_global_tight(array, saveflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = 'RdYlBu_r'
            cmap = plt.get_cmap('RdBu_r', 20)
            #cmap = plt.get_cmap('viridis', 4)
        else:
            cmap = 'viridis'
            cmap = plt.get_cmap('viridis', 20)        
    
    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to North America
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=False,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='lightgrey')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    sc = ax.scatter(array[:, 1], array[:, 2],
                    c=array[:, 0], s=0.8, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
                    transform=ccrs.PlateCarree()) 
    # Set the colorbar properties
    
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04) #, format=mticker.FormatStrFormatter('%.3f'))
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(units, fontsize=12)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=16)
    # ax.text(0.45, -0.1,   'Longitude', fontsize=8, transform=ax.transAxes, ha='left')
    # ax.text(-0.08, 0.4, 'Latitude', fontsize=8, transform=ax.transAxes, rotation='vertical', va='bottom')

    if saveflag:
        savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 600, bbox_inches='tight')    

     # Show the plot
    plt.show()


##
def plot_na(array, saveflag=False, plot_title ='na_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Info for colorbar
    cmin, cmax, cmap = colorbar_info(array)
    
    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to North America
    ax.set_extent([-140, -50, 15, 60], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    ax.coastlines()
    ax.add_feature(cfeature.LAND, facecolor='lightgray')  # Set the land color to light gray
    ax.add_feature(cfeature.BORDERS)

    # scatter data
    sc = ax.scatter(array[:, 1], array[:, 2],
                    c=array[:, 0], s=3, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
                    transform=ccrs.PlateCarree()) 

    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.12, fraction=0.04)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(units, fontsize=10)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=12)
    ax.text(0.45, -0.1, 'Longitude', fontsize=8, transform=ax.transAxes, ha='left')
    ax.text(-0.08, 0.4, 'Latitude', fontsize=8, transform=ax.transAxes, rotation='vertical', va='bottom')

    if saveflag:
        savename = plot_title+'.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi=400)

    # Show the plot
    plt.show()


def colorbar_info(array):

    # Compute and print some stats for the data
    # -----------------------------------------
    stdev = np.nanstd(array[:,0])  # Standard deviation
    omean = np.nanmean(array[:, 0]) # Mean of the data
    datmi = np.nanmin(array[:, 0])  # Min of the data
    datma = np.nanmax(array[:, 0])  # Max of the data
    abmm = np.nanmax(np.abs(array[:, 0])) # Abs max of the data

    # Min max for colorbar
    # --------------------
    if np.nanmin(array[:, 0]) < 0:
        cmax = abmm
        cmin = abmm * -1
        cmap = 'RdBu'
    else:
        cmax = datma
        cmin = datmi
        cmap = 'viridis'

    return cmin, cmax, cmap

def plot_global_contour(lon2d, lat2d, field, saveflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None):

    # Check if cmin and cmax are None
    if cmin is None:
        cmin = np.nanmin(field)
    if cmax is None:
        cmax = np.nanmax(field)

    # Check if field has positive and negative values
    if cmin < 0:
    #    cmax = np.nanmax(np.abs(field))
    #    cmin = -cmax
        cmap = 'RdBu'
    else:
        cmap = 'viridis'

    levels = np.linspace(cmin,cmax,50)
    #levels = np.linspace(0, 1, 50)
    
    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='lightgray')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    sc = ax.contourf(lon2d, lat2d, field, transform=ccrs.PlateCarree(), cmap=cmap, levels=levels)
    
    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.12, fraction=0.04)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label(units, fontsize=10)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=12)
    ax.text(0.45, -0.1,   'Longitude', fontsize=8, transform=ax.transAxes, ha='left')
    ax.text(-0.08, 0.4, 'Latitude', fontsize=8, transform=ax.transAxes, rotation='vertical', va='bottom')

    if saveflag:
        savename = plot_title+'.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 400)    

     # Show the plot
    plt.show()

def plot_global_tight_pcm(array, saveflag=False, meanflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Read binary files and reshape to correct size
    # The number of rows and columns are in the file name
    lats = np.fromfile('../test_data/EASE2_M36km.lats.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    lons = np.fromfile('../test_data/EASE2_M36km.lons.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    
    ds = xr.open_dataset('DAv7_M36.inst3_1d_lndfcstana_Nt.20150901.nc4')
    lon = ds['lon']
    lat = ds['lat']
    n_tile = len(lat)
    
    # Convert to numpy array
    lon = lon.values
    lat = lat.values

    # Make an empty array with dimensions of the grid
    grid = np.full((406, 964), np.nan)
    lats_row = lats[:,1]
    lons_col = lons[1,:]

    # Fill the grid with the values from the dataset
    for i in range(n_tile):
    # Find the row and column of the grid
        row = np.abs(lats_row - lat[i]).argmin()
        col = np.abs(lons_col - lon[i]).argmin()

        # Check if row and col are within the valid range of indices for the grid array
        if row < grid.shape[0] and col < grid.shape[1]:
            grid[row, col] = -9998
        else:
            print('Row or column index out of bounds')

    # Put the array values onto the grid
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if row < grid.shape[0] and col < grid.shape[1]:
            # Check if the value in array[i, 0] is not a NaN
            if not np.isnan(array[i, 0]):
                grid[row, col] = array[i, 0]
        #else:
        #    print('Row or column index out of bounds')   
            # Use 4 decimal places if the number is less than 0.01
    # In the special case of having counts we want to have zeros not -9998. We can do this by setting all -9998 to 0 if plot_title contains 'Number' or 'Percent'
    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0    
    
    # Construct a text string containing the mean and +- standard deviation of non nan values in array[:,0]
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])


    textstr = 'Mean = {}±{} {}'.format(format_number(mean), format_number(std), units)

    if 'Relative $\Delta$ StdDev' in plot_title:
        textstr = 'Mean = {:.1f}±{:.1f} {}'.format(mean, std, units)
    
    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = plt.get_cmap('RdBu_r', 20).copy()
        else:
            cmap = plt.get_cmap('viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)

    cmap.set_under('lightgrey')        

    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to North America
    ax.set_extent([-180, 180, -60, 90], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=False,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='white')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    #sc = ax.scatter(array[:, 1], array[:, 2],
    #                c=array[:, 0], s=0.8, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
    #                transform=ccrs.PlateCarree())

    sc = ax.pcolormesh(lons, lats, grid, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cmin, vmax=cmax) 

    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04) #, format=mticker.FormatStrFormatter('%.3f'))
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('({})'.format(units), fontsize=12)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=18)

    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')

    if saveflag:
        savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 600, bbox_inches='tight')    

     # Show the plot
    plt.show()

def plot_aus_tight_pcm(array, saveflag=False, meanflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Read binary files and reshape to correct size
    # The number of rows and columns are in the file name
    lats = np.fromfile('../test_data/EASE2_M36km.lats.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    lons = np.fromfile('../test_data/EASE2_M36km.lons.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    
    ds = xr.open_dataset('DAv7_M36.inst3_1d_lndfcstana_Nt.20150901.nc4')
    lon = ds['lon']
    lat = ds['lat']
    n_tile = len(lat)
    
    # Convert to numpy array
    lon = lon.values
    lat = lat.values

    # Make an empty array with dimensions of the grid
    grid = np.full((406, 964), np.nan)
    lats_row = lats[:,1]
    lons_col = lons[1,:]

    # Fill the grid with the values from the dataset
    for i in range(n_tile):
    # Find the row and column of the grid
        row = np.abs(lats_row - lat[i]).argmin()
        col = np.abs(lons_col - lon[i]).argmin()

        # Check if row and col are within the valid range of indices for the grid array
        if row < grid.shape[0] and col < grid.shape[1]:
            grid[row, col] = -9998
        else:
            print('Row or column index out of bounds')

    # Put the array values onto the grid
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if row < grid.shape[0] and col < grid.shape[1]:
            # Check if the value in array[i, 0] is not a NaN
            if not np.isnan(array[i, 0]):
                grid[row, col] = array[i, 0]
        #else:
        #    print('Row or column index out of bounds')   

    # In the special case of having counts we want to have zeros not -9998. We can do this by setting all -9998 to 0 if plot_title contains 'Number' or 'Percent'
    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0    
    
    # Construct a text string containing the mean and +- standard deviation of non nan values in array[:,0]
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])

    def format_number(num):
        if abs(num) < 0.01:
            # Use 3 decimal places if the number is less than 0.001
            return '{:.4f}'.format(num)
        elif abs(num) < 1.0:
            # Use 2 decimal places if the number is less than 1
            return '{:.2f}'.format(num)
        else:
            # Use 3 significant figures otherwise
            return '{:.3g}'.format(num)

    textstr = 'Mean = {}±{} {}'.format(format_number(mean), format_number(std), units)

    if 'Relative $\Delta$ StdDev' in plot_title:
        textstr = 'Mean = {:.1f}±{:.1f} {}'.format(mean, std, units)
    
    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = plt.get_cmap('RdBu_r', 20).copy()
        else:
            cmap = plt.get_cmap('viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)

    cmap.set_under('lightgrey')        

    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to North America
    ax.set_extent([112, 154, -39, -10], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=False,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -90, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='white')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    #sc = ax.scatter(array[:, 1], array[:, 2],
    #                c=array[:, 0], s=0.8, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
    #                transform=ccrs.PlateCarree())

    sc = ax.pcolormesh(lons, lats, grid, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cmin, vmax=cmax) 

    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04) #, format=mticker.FormatStrFormatter('%.3f'))
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('({})'.format(units), fontsize=12)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=18)

    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')

    if saveflag:
        savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 600, bbox_inches='tight')    

     # Show the plot
    plt.show()


def plot_sa_tight_pcm(array, saveflag=False, meanflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Read binary files and reshape to correct size
    # The number of rows and columns are in the file name
    lats = np.fromfile('../test_data/EASE2_M36km.lats.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    lons = np.fromfile('../test_data/EASE2_M36km.lons.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    
    ds = xr.open_dataset('DAv7_M36.inst3_1d_lndfcstana_Nt.20150901.nc4')
    lon = ds['lon']
    lat = ds['lat']
    n_tile = len(lat)
    
    # Convert to numpy array
    lon = lon.values
    lat = lat.values

    # Make an empty array with dimensions of the grid
    grid = np.full((406, 964), np.nan)
    lats_row = lats[:,1]
    lons_col = lons[1,:]

    # Fill the grid with the values from the dataset
    for i in range(n_tile):
    # Find the row and column of the grid
        row = np.abs(lats_row - lat[i]).argmin()
        col = np.abs(lons_col - lon[i]).argmin()

        # Check if row and col are within the valid range of indices for the grid array
        if row < grid.shape[0] and col < grid.shape[1]:
            grid[row, col] = -9998
        else:
            print('Row or column index out of bounds')

    # Put the array values onto the grid
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if row < grid.shape[0] and col < grid.shape[1]:
            # Check if the value in array[i, 0] is not a NaN
            if not np.isnan(array[i, 0]):
                grid[row, col] = array[i, 0]
        #else:
        #    print('Row or column index out of bounds')   

    # In the special case of having counts we want to have zeros not -9998. We can do this by setting all -9998 to 0 if plot_title contains 'Number' or 'Percent'
    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0    
    
    # Construct a text string containing the mean and +- standard deviation of non nan values in array[:,0]
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])

    def format_number(num):
        if abs(num) < 0.01:
            # Use 3 decimal places if the number is less than 0.001
            return '{:.4f}'.format(num)
        elif abs(num) < 1.0:
            # Use 2 decimal places if the number is less than 1
            return '{:.2f}'.format(num)
        else:
            # Use 3 significant figures otherwise
            return '{:.3g}'.format(num)

    textstr = 'Mean = {}±{} {}'.format(format_number(mean), format_number(std), units)

    if 'Relative $\Delta$ StdDev' in plot_title:
        textstr = 'Mean = {:.1f}±{:.1f} {}'.format(mean, std, units)
    
    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = plt.get_cmap('RdBu_r', 20).copy()
        else:
            cmap = plt.get_cmap('viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)

    cmap.set_under('lightgrey')        

    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to North America
    ax.set_extent([0, 45, -36, -10], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=False,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    gl.xlocator = mticker.FixedLocator([-180, -135, -35, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='white')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    #sc = ax.scatter(array[:, 1], array[:, 2],
    #                c=array[:, 0], s=0.8, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
    #                transform=ccrs.PlateCarree())

    sc = ax.pcolormesh(lons, lats, grid, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cmin, vmax=cmax) 

    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04) #, format=mticker.FormatStrFormatter('%.3f'))
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('({})'.format(units), fontsize=12)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=18)

    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')

    if saveflag:
        savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 600, bbox_inches='tight')    

     # Show the plot
    plt.show()

# Plot the Tibetan Plateau
# -------------------------
# 
def plot_tb_tight_pcm(array, saveflag=False, meanflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Read binary files and reshape to correct size
    # The number of rows and columns are in the file name
    lats = np.fromfile('../test_data/EASE2_M36km.lats.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    lons = np.fromfile('../test_data/EASE2_M36km.lons.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    
    ds = xr.open_dataset('DAv7_M36.inst3_1d_lndfcstana_Nt.20150901.nc4')
    lon = ds['lon']
    lat = ds['lat']
    n_tile = len(lat)
    
    # Convert to numpy array
    lon = lon.values
    lat = lat.values

    # Make an empty array with dimensions of the grid
    grid = np.full((406, 964), np.nan)
    lats_row = lats[:,1]
    lons_col = lons[1,:]

    # Fill the grid with the values from the dataset
    for i in range(n_tile):
    # Find the row and column of the grid
        row = np.abs(lats_row - lat[i]).argmin()
        col = np.abs(lons_col - lon[i]).argmin()

        # Check if row and col are within the valid range of indices for the grid array
        if row < grid.shape[0] and col < grid.shape[1]:
            grid[row, col] = -9998
        else:
            print('Row or column index out of bounds')

    # Put the array values onto the grid
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if row < grid.shape[0] and col < grid.shape[1]:
            # Check if the value in array[i, 0] is not a NaN
            if not np.isnan(array[i, 0]):
                grid[row, col] = array[i, 0]
        #else:
        #    print('Row or column index out of bounds')   

    # In the special case of having counts we want to have zeros not -9998. We can do this by setting all -9998 to 0 if plot_title contains 'Number' or 'Percent'
    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0    
    
    # Construct a text string containing the mean and +- standard deviation of non nan values in array[:,0]
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])

    def format_number(num):
        if abs(num) < 0.01:
            # Use 3 decimal places if the number is less than 0.001
            return '{:.4f}'.format(num)
        elif abs(num) < 1.0:
            # Use 2 decimal places if the number is less than 1
            return '{:.2f}'.format(num)
        else:
            # Use 3 significant figures otherwise
            return '{:.3g}'.format(num)

    textstr = 'Mean = {}±{} {}'.format(format_number(mean), format_number(std), units)

    if 'Relative $\Delta$ StdDev' in plot_title:
        textstr = 'Mean = {:.1f}±{:.1f} {}'.format(mean, std, units)
    
    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = plt.get_cmap('RdBu_r', 20).copy()
        else:
            cmap = plt.get_cmap('viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)

    cmap.set_under('lightgrey')        

    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to Tibetan Plateau
    ax.set_extent([73, 105, 26, 40], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    #gl.xlocator = mticker.FixedLocator([-180, -135, -35, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='white')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    #ax.add_feature(cfeature.BORDERS)

    # scatter data
    #sc = ax.scatter(array[:, 1], array[:, 2],
    #                c=array[:, 0], s=0.8, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
    #                transform=ccrs.PlateCarree())

    sc = ax.pcolormesh(lons, lats, grid, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cmin, vmax=cmax) 

    # Plot a star at the center of the Tibetan Plateau
    ax.plot(84.0, 31.0, marker='*', color='red', markersize=10, transform=ccrs.PlateCarree())

    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04) #, format=mticker.FormatStrFormatter('%.3f'))
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('({})'.format(units), fontsize=12)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=18)

    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')

    if saveflag:
        savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 600, bbox_inches='tight')    

     # Show the plot
    plt.show()    


# Plot the Tibetan Plateau
# -------------------------
# 
def plot_NA_tight_pcm(array, saveflag=False, meanflag=False, plot_title ='global_plot', units='na', cmin=None, cmax=None, cmap=None):

    # Read binary files and reshape to correct size
    # The number of rows and columns are in the file name
    lats = np.fromfile('../test_data/EASE2_M36km.lats.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    lons = np.fromfile('../test_data/EASE2_M36km.lons.964x406x1.double', 
                          dtype=np.float64).reshape((406,964))
    
    ds = xr.open_dataset('DAv7_M36.inst3_1d_lndfcstana_Nt.20150901.nc4')
    lon = ds['lon']
    lat = ds['lat']
    n_tile = len(lat)
    
    # Convert to numpy array
    lon = lon.values
    lat = lat.values

    # Make an empty array with dimensions of the grid
    grid = np.full((406, 964), np.nan)
    lats_row = lats[:,1]
    lons_col = lons[1,:]

    # Fill the grid with the values from the dataset
    for i in range(n_tile):
    # Find the row and column of the grid
        row = np.abs(lats_row - lat[i]).argmin()
        col = np.abs(lons_col - lon[i]).argmin()

        # Check if row and col are within the valid range of indices for the grid array
        if row < grid.shape[0] and col < grid.shape[1]:
            grid[row, col] = -9998
        else:
            print('Row or column index out of bounds')

    # Put the array values onto the grid
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if row < grid.shape[0] and col < grid.shape[1]:
            # Check if the value in array[i, 0] is not a NaN
            if not np.isnan(array[i, 0]):
                grid[row, col] = array[i, 0]
        #else:
        #    print('Row or column index out of bounds')   

    # In the special case of having counts we want to have zeros not -9998. We can do this by setting all -9998 to 0 if plot_title contains 'Number' or 'Percent'
    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0    
    
    # Construct a text string containing the mean and +- standard deviation of non nan values in array[:,0]
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])

    def format_number(num):
        if abs(num) < 0.01:
            # Use 3 decimal places if the number is less than 0.001
            return '{:.4f}'.format(num)
        elif abs(num) < 1.0:
            # Use 2 decimal places if the number is less than 1
            return '{:.2f}'.format(num)
        else:
            # Use 3 significant figures otherwise
            return '{:.3g}'.format(num)

    textstr = 'Mean = {}±{} {}'.format(format_number(mean), format_number(std), units)

    if 'Relative $\Delta$ StdDev' in plot_title:
        textstr = 'Mean = {:.1f}±{:.1f} {}'.format(mean, std, units)
    
    # Check if cmin and cmax are None
    if cmin is None:
        # Info for colorbar
        cmin, cmax, cmap = colorbar_info(array)

    # Check if field has positive and negative values
    if cmap is None:
        if cmin < 0:
            cmap = plt.get_cmap('RdBu_r', 20).copy()
        else:
            cmap = plt.get_cmap('viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)

    cmap.set_under('lightgrey')        

    # Create the plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    # Set the extent to North America
    ax.set_extent([-125.0, -66.5, 24.5, 49.5], crs=ccrs.PlateCarree())

    # plot grid lines
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), draw_labels=True,
                          linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    #gl.xlocator = mticker.FixedLocator([-180, -135, -35, -45, 0, 45, 90, 135, 179.9])
    ax.tick_params(labelbottom=False, labeltop=False, labelleft=False, labelright=False)

    # ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='white')  # Set the land color to light gray
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)

    # scatter data
    #sc = ax.scatter(array[:, 1], array[:, 2],
    #                c=array[:, 0], s=0.8, linewidth=0, cmap=cmap, vmin=cmin, vmax=cmax,
    #                transform=ccrs.PlateCarree())

    sc = ax.pcolormesh(lons, lats, grid, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cmin, vmax=cmax) 

    # Plot a star at the center of the Tibetan Plateau
    ax.plot(84.0, 31.0, marker='*', color='red', markersize=10, transform=ccrs.PlateCarree())

    # Set the colorbar properties
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", pad=.05, fraction=0.04) #, format=mticker.FormatStrFormatter('%.3f'))
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('({})'.format(units), fontsize=12)

    # Set the axis and title labels
    plt.title(plot_title, fontsize=18)

    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')

    if saveflag:
        savename = re.sub('[^0-9a-zA-Z]+', '_', plot_title) + '.png'
        print(" Saving figure as", savename, "\n")
        plt.savefig(savename, dpi = 600, bbox_inches='tight')    

     # Show the plot
    plt.show()    

#################################################################################

def load_ease_grid(ease_path):
    """Load EASE grid data"""
    lats = np.fromfile(f'{ease_path}/EASE2_M36km.lats.964x406x1.double', 
                      dtype=np.float64).reshape((406,964))
    lons = np.fromfile(f'{ease_path}/EASE2_M36km.lons.964x406x1.double', 
                      dtype=np.float64).reshape((406,964))
    return lats, lons

def create_grid_mapping(array, lats, lons, lats_row, lons_col):
    """Map array data to grid"""
    grid = np.full((lats.shape[0], lons.shape[1]), -9998., dtype=np.float64)
    
    # Map array values to grid
    for i in range(len(array)):
        row = np.abs(lats_row - array[i, 2]).argmin()
        col = np.abs(lons_col - array[i, 1]).argmin()
        if row < grid.shape[0] and col < grid.shape[1]:
            if not np.isnan(array[i, 0]):
                grid[row, col] = array[i, 0]
    
    return grid

#################################################################################

def plot_region(array, lon_min, lon_max, lat_min, lat_max, 
                ease_path='../test_data',
                saveflag=False, 
                meanflag=False, 
                plot_title='regional_plot', 
                units='na', 
                cmin=None, 
                cmax=None, 
                cmap=None,
                output_dir='./plots',
                save_fmt='png',
                save_dpi=600,
                star_lon=None,
                star_lat=None):
    """
    Plot data for specified region using EASE grid.

    Call as:

    # Usage over CONUS
    lon_min = -125.0
    lon_max = -66.0
    lat_min = 24.0
    lat_max = 50.0
    
    # Plot the data
    plot_region(map_array, 
            lon_min, lon_max,
            lat_min, lat_max,
            meanflag=False,
            saveflag=False,
            units='SCF',
            plot_title=f'{expt_name} {start_date_str} - {end_date_str}:\n MODIS Snow Cover Fraction',
            cmin=0,
            cmax=1,
            cmap='viridis')
    
    Parameters
    ----------
    array : np.ndarray
        Array of shape (n,3) with values, lons, lats
    lon_min, lon_max : float
        Longitude bounds
    lat_min, lat_max : float
        Latitude bounds
    ease_path : str
        Path to EASE grid files
    """
    # Load EASE grid
    lats, lons = load_ease_grid(ease_path)
    lats_row = lats[:,1]
    lons_col = lons[1,:]
    
    # Create grid mapping
    grid = create_grid_mapping(array, lats, lons, lats_row, lons_col)
    
    # Handle counts/percentages
    if 'Number' in plot_title or 'Percent' in plot_title:
        grid[grid == -9998] = 0
    
    # Calculate statistics
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])
    textstr = format_stats(mean, std, units, plot_title)
    
    # Set up colormap
    if cmin is None or cmax is None:
        cmin, cmax, cmap = colorbar_info(array)
    if cmap is None:
        cmap = plt.get_cmap('RdBu_r' if cmin < 0 else 'viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)
    cmap.set_under('lightgrey')
    
    # Create plot
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    
    # Set region extent
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # Add map features
    setup_map_features(ax)
    
    # Plot data
    sc = ax.pcolormesh(lons, lats, grid, 
                       transform=ccrs.PlateCarree(), 
                       cmap=cmap, 
                       vmin=cmin, 
                       vmax=cmax)
    
    # Add star if star_lon/lat provided
    if star_lon is not None and star_lat is not None:
        ax.plot(star_lon, star_lat, 'r*', markersize=10, transform=ccrs.PlateCarree())
    
    # Add colorbar and labels
    setup_colorbar(sc, ax, cmin, cmax, units)
    plt.title(plot_title, fontsize=18)
    
    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')
    
    if saveflag:
        save_plot(plot_title)
    
    plt.show()
    
#################################################################################

def format_stats(mean, std, units, plot_title):
    """Format statistics string"""
    if 'Relative $\Delta$ StdDev' in plot_title:
        return f'Mean = {mean:.1f}±{std:.1f} {units}'
    
    def format_number(num):
        if abs(num) < 0.01:
            return f'{num:.4f}'
        elif abs(num) < 1.0:
            return f'{num:.2f}'
        return f'{num:.3g}'
    
    return f'Mean = {format_number(mean)}±{format_number(std)} {units}'

#################################################################################

def setup_map_features(ax):
    """Set up map features and gridlines"""
    gl = ax.gridlines(crs=ccrs.PlateCarree(central_longitude=0), 
                     draw_labels=True,
                     linewidth=1, 
                     color='gray', 
                     alpha=0.5, 
                     linestyle='-')
    gl.xlabel_style = {'size': 5, 'color': 'black'}
    gl.ylabel_style = {'size': 5, 'color': 'black'}
    ax.tick_params(labelbottom=False, labeltop=False, 
                  labelleft=False, labelright=False)
    ax.add_feature(cfeature.LAND, facecolor='white')
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)

def setup_colorbar(sc, ax, cmin, cmax, units):
    """Set up colorbar"""
    cbar = plt.colorbar(sc, ax=ax, orientation="horizontal", 
                       pad=.05, fraction=0.04)
    cbar.set_ticks(np.arange(cmin, cmax+0.000000001, (cmax-cmin)/4))
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(f'({units})', fontsize=12)
 

def save_plot(plot_title, output_dir='./plots', fmt='png', dpi=600):
    """
    Save plot to specified directory
    
    Parameters
    ----------
    plot_title : str
        Title of plot used for filename
    output_dir : str or Path
        Directory to save plots
    fmt : str
        File format (png, pdf, jpg)
    dpi : int
        Resolution for raster formats
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Clean filename
    clean_title = re.sub('[^0-9a-zA-Z]+', '_', plot_title)
    
    # Construct full save path
    savename = output_path / f"{clean_title}.{fmt}"
    print(f"Saving figure as {savename}")
    
    # Save with specified parameters
    plt.savefig(savename, dpi=dpi, bbox_inches='tight', format=fmt)

def colorbar_info(array):

    # Compute and print some stats for the data
    # -----------------------------------------
    stdev = np.nanstd(array[:,0])  # Standard deviation
    omean = np.nanmean(array[:, 0]) # Mean of the data
    datmi = np.nanmin(array[:, 0])  # Min of the data
    datma = np.nanmax(array[:, 0])  # Max of the data
    abmm = np.nanmax(np.abs(array[:, 0])) # Abs max of the data

    # Min max for colorbar
    # --------------------
    if np.nanmin(array[:, 0]) < 0:
        cmax = abmm
        cmin = abmm * -1
        cmap = 'RdBu'
    else:
        cmax = datma
        cmin = datmi
        cmap = 'viridis'

    return cmin, cmax, cmap    

#################################################################################

def plot_region_scatter(array, lon_min, lon_max, lat_min, lat_max, 
                grid_type='ease',
                ease_path='../test_data',
                saveflag=False, 
                meanflag=False, 
                plot_title='regional_plot', 
                units='na', 
                cmin=None, 
                cmax=None, 
                cmap=None,
                output_dir='./plots',
                save_fmt='png',
                save_dpi=600,
                star_lon=None,
                star_lat=None,
                point_size=6):
    """
    Plot data for specified region using either EASE or MODIS grid.
    
    Parameters
    ----------
    array : np.ndarray
        Array of shape (n,3) with values, lons, lats
    grid_type : str
        Either 'ease' or 'modis'
    """  
    
    # Calculate statistics
    mean = np.nanmean(array[:, 0])
    std = np.nanstd(array[:, 0])
    textstr = format_stats(mean, std, units, plot_title)

    # Extract values
    lons = array[:, 1]
    lats = array[:, 2]
    vals = array[:, 0]
    
    # Set up colormap
    if cmin is None or cmax is None:
        cmin, cmax, cmap = colorbar_info(array)
    if cmap is None:
        cmap = plt.get_cmap('RdBu_r' if cmin < 0 else 'viridis', 20).copy()
    else:
        cmap = plt.get_cmap(cmap)
    cmap.set_under('lightgrey')
    
    # Create plot
    fig = plt.figure(figsize=(15, 9))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    
    # Set region extent
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
    
    # Add map features
    setup_map_features(ax)
    
    # Plot data
    sc = ax.scatter(lons, lats, c=vals, s=point_size, cmap=cmap, 
                    vmin=cmin, vmax=cmax, edgecolor='none',
                    transform=ccrs.PlateCarree())
    
    # Add star if star_lon/lat provided
    if star_lon is not None and star_lat is not None:
        ax.plot(star_lon, star_lat, 'r*', markersize=10, transform=ccrs.PlateCarree())
    
    # Add colorbar and labels
    setup_colorbar(sc, ax, cmin, cmax, units)
    plt.title(plot_title, fontsize=18)
    
    if meanflag:
        ax.text(0.38, 0.05, textstr, fontsize=14, transform=ax.transAxes, ha='left')
    
    if saveflag:
        save_plot(plot_title)  

    # Overlay EASE grid lines (only if grid_type is 'modis')
    if grid_type == 'modis':
        ease_lat, ease_lon = load_grid(grid_type='ease', ease_path=ease_path)

        # Plot EASE latitude lines
        for i in range(0, ease_lat.shape[0], 1):  # every 1 rows
            ax.plot(ease_lon[i, :], ease_lat[i, :], color='grey', linewidth=0.5, alpha=0.6,
                    transform=ccrs.PlateCarree())

        # Plot EASE longitude lines
        for j in range(0, ease_lon.shape[1], 1):  # every 1 columns
            ax.plot(ease_lon[:, j], ease_lat[:, j], color='grey', linewidth=0.5, alpha=0.6,
                    transform=ccrs.PlateCarree())  
    
    plt.show()