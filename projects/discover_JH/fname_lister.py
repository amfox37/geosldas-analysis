#!/usr/bin/env python
# coding: utf-8

# In[5]:


import os
# import shutil

top_directory = "/home/amfox/smap/SMAP_Nature/ASCAT_EUMETSAT/Metop_C"  # Specify the top directory
out_directory = "/discover/nobackup/amfox/ASCAT_fname_lists" # Specify the top output directory
start_year = 2007  # Specify the start year
end_year = 2023  # Specify the end year

def parse_timestamp(timestamp):
    if len(timestamp) < 14:
        return None
    try:
        year = int(timestamp[:4])
        month = int(timestamp[4:6])
        day = int(timestamp[6:8])
        hour = int(timestamp[8:10])
        minute = int(timestamp[10:12])
        second = int(timestamp[12:14])
        return year, month, day, hour, minute, second
    except ValueError:
        return None

def create_directory(year, month, day):
    directory = os.path.join(out_directory, "Y{:04d}".format(year), "M{:02d}".format(month), "D{:02d}".format(day))
    os.makedirs(directory, exist_ok=True)
    return directory

def search_and_copy_files():
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            directory = os.path.join(top_directory, "Y{:04d}".format(year), "M{:02d}".format(month))
            if not os.path.isdir(directory):
                continue
            for filename in os.listdir(directory):
                if filename.startswith("M03-ASCA-ASCSMO02-NA-"):
                    timestamp = filename.split("-")[5]
                    file_year, file_month, file_day, _, _, _ = parse_timestamp(timestamp)
                    if file_year != year or file_month != month:
                        continue
                    destination_directory = create_directory(year, month, file_day)
                    destination_file = os.path.join(destination_directory, "M03-ASCA-ASCSMO02.txt")
                    with open(destination_file, "a") as f:
                        f.write(os.path.join(filename) + "\n")
                    #    f.write(os.path.join(directory, filename) + "\n")
                    #shutil.copy2(os.path.join(directory, filename), destination_directory)

search_and_copy_files()

