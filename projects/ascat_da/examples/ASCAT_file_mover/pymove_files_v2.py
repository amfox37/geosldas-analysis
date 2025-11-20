import os
import re
import shutil
import tarfile
import gzip

dir_path = "/discover/nobackup/amfox/ASCAT_file_mover/archive.eumetsat.int/umarf-gwt/onlinedownload/amfox37/4859826"  # Replace with the directory path containing the files
out_path = "/discover/nobackup/amfox/ASCAT_file_mover"

# Create the satellite directories if they don't already exist
satellites = ['Metop-A', 'Metop-B', 'Metop-C']
for sat in satellites:
    sat_dir = os.path.join(out_path, sat)
    if not os.path.exists(sat_dir):
        os.makedirs(sat_dir)

# Extract any tarballs in the directory
for file in os.listdir(dir_path):
    if file.endswith('.tar'):
        tar = tarfile.open(os.path.join(dir_path, file), 'r')
        tar.extractall(path=dir_path)
        tar.close()

        for inner_file in os.listdir(dir_path):
            if inner_file.endswith('.gz'):
                print(f"Processing {inner_file}...")
                with gzip.open(os.path.join(dir_path, inner_file), 'rb') as f_in:
                    # create the output filename by removing the .gz extension
                    output_filename = os.path.splitext(inner_file)[0]
                    # open the output file for writing
                    with open(os.path.join(dir_path, output_filename), 'wb') as f_out:
                        # copy the contents of the gz file to the output file
                        shutil.copyfileobj(f_in, f_out)
                # remove the original compressed file
                os.remove(os.path.join(dir_path, inner_file))
                print(f"Finished processing {inner_file}.")


# Move the BFR files to the appropriate directories based on the satellite name in the file names
for file in os.listdir(dir_path):
    if file.endswith('.bfr'):
        # Extract the satellite name from the file name
        sat_name = file.split('-')[0]
        
        # Extract the time stamp from the file name
        timestamp = re.search(r'(?<=-)[0-9]{14}\.', file).group()[:-1]
        year = timestamp[:4]
        month = timestamp[4:6]
        
        # Move the file to the appropriate directory based on the satellite name and time stamp
        if sat_name == 'M01':
            sat_dir = os.path.join(out_path, 'Metop-B')
        elif sat_name == 'M02':
            sat_dir = os.path.join(out_path, 'Metop-A')
        elif sat_name == 'M03':
            sat_dir = os.path.join(out_path, 'Metop-C')
        else:
            continue
        
        year_dir = os.path.join(sat_dir, 'Y'+year)
        if not os.path.exists(year_dir):
            os.makedirs(year_dir)
        month_dir = os.path.join(year_dir, 'M'+month)
        if not os.path.exists(month_dir):
            os.makedirs(month_dir)
        
        shutil.move(os.path.join(dir_path, file), os.path.join(month_dir, file))
