# This script moves CYGNSS observations files from the "data" directory to directories with structure /CYGNSS_obs/YYYY/MM, using the timestamp in the filename.

import os
import shutil

# Define the path to the directory containing the CYGNSS observations files
data_dir = '/discover/nobackup/amfox/CYGNSS_obs/data'
output_dir = '/discover/nobackup/amfox/CYGNSS_obs'

# Get the list of files in the data directory that begin with "cyg.ddmi." and end with "soil-moisture-36km.a32.d33.nc"
files = [f for f in os.listdir(data_dir) if f.startswith('cyg.ddmi.s') and f.endswith('soil-moisture-36km.a32.d33.nc')]
print(f'Number of files: {len(files)}')

# Loop through the files

for file in files:
    print(f'File: {file}')
    # Get the timestamp from the filename for filenames like: cyg.ddmi.s20180801-030000-e20180801-210000.l3.grid-soil-moisture-36km.a32.d33.nc
    timestamp = file.split('ddmi.s')[1].split('-')[0]
    print(f'Timestamp: {timestamp}')    
    year = timestamp[:4]
    month = timestamp[4:6]
    day = timestamp[6:8]

    # Check if directory structure i.e. /CYGNSS_obs/Y2018/M08 exists, if not create it
    year_dir = os.path.join(output_dir, f'Y{year}')
    month_dir = os.path.join(year_dir, f'M{month}')
    if not os.path.exists(month_dir):
        os.makedirs(month_dir)

    # Move the file to the appropriate directory
    shutil.copy(os.path.join(data_dir, file), month_dir)

print('Files moved successfully!')

# Check that each year and month directory contains the correct number of daily files
# Sort the directories by year then month

for year in sorted(os.listdir(output_dir)):
    year_dir = os.path.join(output_dir, year)
    for month in sorted(os.listdir(year_dir)):
        month_dir = os.path.join(year_dir, month)
        files = os.listdir(month_dir)
        print(f'{year}/{month}: {len(files)} files')
        # Check that the number of files is correct
        # Number of days in the month
        if month in ['M01', 'M03', 'M05', 'M07', 'M08', 'M10', 'M12']:
            days = 31
        elif month in ['M04', 'M06', 'M09', 'M11']:
            days = 30
        elif month == 'M02':
            if int(year[1:]) % 4 == 0:
                days = 29
            else:
                days = 28
        if len(files) != days:
            print(f'Error: {len(files)} files found in {year}/{month}, expected {days} files')
        else:
            print('Number of files is correct')




