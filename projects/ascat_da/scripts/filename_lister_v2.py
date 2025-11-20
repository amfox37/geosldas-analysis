import os
import calendar

top_directory = "/discover/nobackup/projects/gmao/smap/SMAP_Nature/ASCAT_EUMETSAT/Metop_A"  # Specify the top directory
out_directory = "/discover/nobackup/amfox/ASCAT_fname_lists_v2"  # Specify the top output directory
start_year = 2007  # Specify the start year
end_year = 2024  # Specify the end year

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
            _, num_days = calendar.monthrange(year, month)  # Get the actual number of days in the month
            for day in range(1, num_days + 1):
                destination_directory = create_directory(year, month, day)
                destination_file = os.path.join(destination_directory, "M02-ASCA-ASCSMO02.txt")
                with open(destination_file, "w") as f:  # Create an empty file
                    pass
                directory = os.path.join(top_directory, "Y{:04d}".format(year), "M{:02d}".format(month))
                if os.path.isdir(directory):
                    for filename in sorted(os.listdir(directory)):
                        if filename.startswith("M02-ASCA-ASCSMO02-NA-"):
                            timestamp = filename.split("-")[5]
                            file_year, file_month, file_day, _, _, _ = parse_timestamp(timestamp)
                            if file_year != year or file_month != month or file_day != day:
                                continue
                            with open(destination_file, "a") as f:
                                f.write(os.path.join(filename) + "\n")

search_and_copy_files()
