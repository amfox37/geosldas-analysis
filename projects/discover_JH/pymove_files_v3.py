import os
import re
import shutil

dir_path = "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/ASCAT_file_mover/Downloads_Nov_2024"  # Replace with the directory path containing the files
out_path = "/Users/amfox/Desktop/GEOSldas_diagnostics/test_data/ASCAT_file_mover"

# Create the satellite directories if they don't already exist
satellites = ['Metop-A', 'Metop-B', 'Metop-C']
for sat in satellites:
    sat_dir = os.path.join(out_path, sat)
    if not os.path.exists(sat_dir):
        os.makedirs(sat_dir)

# Find all the obs files starting with "W_XX-EUMETSAT-Darmstadt,SOUNDING+SATELLITE,METOP" in all subdirectories under the dir_path
# These files have names like "W_XX-EUMETSAT-Darmstadt,SOUNDING+SATELLITE,METOPB+ASCAT_C_EUMP_20230930011800_57245_eps_o_250_ssm_l2.bin" 

for root, dirs, files in os.walk(dir_path):

    for file in files:
        if file.startswith('W_XX-EUMETSAT-Darmstadt,SOUNDING+SATELLITE,METOP'):
            # Extract the satellite name from the file name
            sat_name = file.split('+')[1].split(',')[1]

            # Extract the time stamp from the file name
            timestamp = re.search(r'[0-9]{14}', file).group()
            year = timestamp[:4]
            month = timestamp[4:6]

            # Move the file to the appropriate directory based on the satellite name and time stamp
            if sat_name == 'METOPA':
                sat_dir = os.path.join(out_path, 'Metop-A')
            elif sat_name == 'METOPB':
                sat_dir = os.path.join(out_path, 'Metop-B')
            elif sat_name == 'METOPC':
                sat_dir = os.path.join(out_path, 'Metop-C')
            else:
                continue

            year_dir = os.path.join(sat_dir, 'Y'+year)
            if not os.path.exists(year_dir):
                os.makedirs(year_dir)
            month_dir = os.path.join(year_dir, 'M'+month)
            if not os.path.exists(month_dir):
                os.makedirs(month_dir)

            shutil.copy(os.path.join(root, file), os.path.join(month_dir, file))

            # Make symbolic link to the file in the same directory with a name such that W_XX-EUMETSAT-Darmstadt,SOUNDING+SATELLITE,METOPB+ASCAT_C_EUMP_20230930011800_57245_eps_o_250_ssm_l2.bin
            # becomes M01-ASCA-ASCSMO02-NA-5.0-20230930011800.bfr. M02- for METOPA M01- for METOPB, M03- for METOPC

            if sat_name == 'METOPA':
                sat_code = 'M02'
            elif sat_name == 'METOPB':
                sat_code = 'M01'
            elif sat_name == 'METOPC':
                sat_code = 'M03'

            symlink_name = sat_code + '-ASCA-ASCSMO02-NA-5.0-' + timestamp + '.bfr'

            # We need to remove the existing symlink if it exists (might have duplicate files)
            if os.path.exists(os.path.join(month_dir, symlink_name)):
                os.remove(os.path.join(month_dir, symlink_name))

            relative_path = os.path.relpath(os.path.join(month_dir, file), month_dir)
            os.symlink(relative_path, os.path.join(month_dir, symlink_name))

print("Files moved successfully!")
