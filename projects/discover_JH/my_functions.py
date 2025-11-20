import numpy as np
import os

from pybufrkit.decoder import generate_bufr_message
from pybufrkit.decoder import Decoder
from pybufrkit.dataquery import NodePathParser, DataQuerent

#####################################
def read_obsfcstana(path, file_name, printflag=False):

    # Precision of fortran tag
    int_precision = 'int32'
    # Precision of data in input file
    float_precision = 'float32'
    # Precision of data in input file
    logical_precision = 'int32'

    # Initialize outputs in case file does not exist or is empty
    nodata = -9999
    # date_time = {
    #     'year': nodata,
    #     'month': nodata,
    #     'day': nodata,
    #     'hour': nodata,
    #     'min': nodata,
    #     'sec': nodata,
    #     'dofyr': nodata,
    #     'pentad': nodata
    # }
    date_time = []
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

    # Determine machine format
    machfmt = 'b'


    # Get a list of files with a similar name in a directory
    file_ext = '.bin'
    # Use os.walk to recursively traverse all directories and subdirectories under the top-level directory
    files = [os.path.join(root, file) for root, dirs, files in os.walk(path) for file in files if file.startswith(file_name) and file.endswith(file_ext)]

    if printflag:
        print(files)

    # Open each file in turn
    mode = 'rb' if machfmt == 'b' else 'rl'

    for file in files:
        with open(os.path.join(path, file), mode) as ifp:
            if printflag:
                print ('Reading file ', file, '...')
            
            # Read N_obs and time stamp entry
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            N_obs = np.fromfile(ifp, int_precision, 1)
            N_obs = int(N_obs)
            year = np.fromfile(ifp, int_precision, 1)
            month = np.fromfile(ifp, int_precision, 1)
            day = np.fromfile(ifp, int_precision, 1)
            hour = np.fromfile(ifp, int_precision, 1)
            minute = np.fromfile(ifp, int_precision, 1)
            second = np.fromfile(ifp, int_precision, 1)
            dofyr = np.fromfile(ifp, int_precision, 1)
            pentad = np.fromfile(ifp, int_precision, 1)
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            # Populate date_time structure
            date_time_tmp = {
                'year': year,
                'month': month,
                'day': day,
                'hour': hour,
                'min': minute,
                'sec': second,
                'dofyr': dofyr,
                'pentad': pentad
            }
            date_time.append(date_time_tmp) 
            
            # Read observation assim flag
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            tmp_data = np.fromfile(ifp, logical_precision, N_obs)
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            tmp_data2 = np.zeros((N_obs, 1))
            indices = np.where(tmp_data != 0)[0]
            tmp_data2[indices] = 1
            obs_assim = np.append(obs_assim, tmp_data2)
            
            # Read species information
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_species = np.append(obs_species, np.fromfile(ifp, int_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            # Read tile number information
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_tilenum = np.append(obs_tilenum, np.fromfile(ifp, int_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read longitude
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_lon = np.append(obs_lon, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read latitude
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_lat = np.append(obs_lat, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            # Read observation value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_obs = np.append(obs_obs, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_obsvar = np.append(obs_obsvar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space model forecast value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_fcst = np.append(obs_fcst, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space model forecast variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_fcstvar = np.append(obs_fcstvar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space analysis value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_ana = np.append(obs_ana, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space analysis variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_anavar = np.append(obs_anavar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            
    # Close file
        ifp.close()
     
    # print('Total number of obs = ',len(obs_assim))

    return date_time, obs_species, obs_tilenum, obs_lon, obs_lat, obs_obs, obs_obsvar, obs_fcst, obs_fcstvar, obs_ana, obs_anavar 

#####################################
def read_obsfcstana_extend_datetime(path, file_name, printflag=False):
    
    # Precision of fortran tag
    int_precision = 'int32'
    # Precision of data in input file
    float_precision = 'float32'
    # Precision of data in input file
    logical_precision = 'int32'

    # Initialize outputs in case file does not exist or is empty
    nodata = -9999
    # date_time = {
    #     'year': nodata,
    #     'month': nodata,
    #     'day': nodata,
    #     'hour': nodata,
    #     'min': nodata,
    #     'sec': nodata,
    #     'dofyr': nodata,
    #     'pentad': nodata
    # }
    date_time = []
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

    # Determine machine format
    machfmt = 'b'


    # Get a list of files with a similar name in a directory
    file_ext = '.bin'
    # Use os.walk to recursively traverse all directories and subdirectories under the top-level directory
    files = [os.path.join(root, file) for root, dirs, files in os.walk(path) for file in files if file.startswith(file_name) and file.endswith(file_ext)]

    if printflag:
        print(files)

    # Open each file in turn
    mode = 'rb' if machfmt == 'b' else 'rl'

    for file in files:
        with open(os.path.join(path, file), mode) as ifp:
            if printflag:
                print ('Reading file ', file, '...')
            
            # Read N_obs and time stamp entry
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            N_obs = np.fromfile(ifp, int_precision, 1)
            N_obs = int(N_obs)
            year = np.fromfile(ifp, int_precision, 1)
            month = np.fromfile(ifp, int_precision, 1)
            day = np.fromfile(ifp, int_precision, 1)
            hour = np.fromfile(ifp, int_precision, 1)
            minute = np.fromfile(ifp, int_precision, 1)
            second = np.fromfile(ifp, int_precision, 1)
            dofyr = np.fromfile(ifp, int_precision, 1)
            pentad = np.fromfile(ifp, int_precision, 1)
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            # Populate date_time structure
            date_time_tmp = {
                'year': year,
                'month': month,
                'day': day,
                'hour': hour,
                'min': minute,
                'sec': second,
                'dofyr': dofyr,
                'pentad': pentad
            }
            # Extend date_time so it has the same length as everything else
            date_time.extend([date_time_tmp] * N_obs) 
            
            # Read observation assim flag
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            tmp_data = np.fromfile(ifp, logical_precision, N_obs)
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            tmp_data2 = np.zeros((N_obs, 1))
            indices = np.where(tmp_data != 0)[0]
            tmp_data2[indices] = 1
            obs_assim = np.append(obs_assim, tmp_data2)
            
            # Read species information
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_species = np.append(obs_species, np.fromfile(ifp, int_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            # Read tile number information
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_tilenum = np.append(obs_tilenum, np.fromfile(ifp, int_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read longitude
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_lon = np.append(obs_lon, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read latitude
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_lat = np.append(obs_lat, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            # Read observation value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_obs = np.append(obs_obs, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_obsvar = np.append(obs_obsvar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space model forecast value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_fcst = np.append(obs_fcst, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space model forecast variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_fcstvar = np.append(obs_fcstvar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space analysis value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_ana = np.append(obs_ana, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space analysis variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_anavar = np.append(obs_anavar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            
    # Close file
        ifp.close()
     
    # print('Total number of obs = ',len(obs_assim))

    return date_time, obs_species, obs_tilenum, obs_lon, obs_lat, obs_obs, obs_obsvar, obs_fcst, obs_fcstvar, obs_ana, obs_anavar 

#####################################
def read_obsfcstana_pentads(path, file_name, printflag=False):

    # Precision of fortran tag
    int_precision = 'int32'
    # Precision of data in input file
    float_precision = 'float32'
    # Precision of data in input file
    logical_precision = 'int32'

    # Initialize outputs in case file does not exist or is empty
    nodata = -9999
    # date_time = {
    #     'year': nodata,
    #     'month': nodata,
    #     'day': nodata,
    #     'hour': nodata,
    #     'min': nodata,
    #     'sec': nodata,
    #     'dofyr': nodata,
    #     'pentad': nodata
    # }
    date_time = []
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
    pentad2 = []

    # Determine machine format
    machfmt = 'b'


    # Get a list of files with a similar name in a directory
    file_ext = '.bin'
    files = [os.path.join(root, file) for root, dirs, files in os.walk(path) for file in files if file.startswith(file_name) and file.endswith(file_ext)]

    if printflag:
        print(files)

    # Open each file in turn
    mode = 'rb' if machfmt == 'b' else 'rl'

    for file in files:
        with open(os.path.join(path, file), mode) as ifp:
            if printflag:
                print ('Reading file ', file, '...')
            
            # Read N_obs and time stamp entry
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            N_obs = np.fromfile(ifp, int_precision, 1)
            N_obs = int(N_obs)
            year = np.fromfile(ifp, int_precision, 1)
            month = np.fromfile(ifp, int_precision, 1)
            day = np.fromfile(ifp, int_precision, 1)
            hour = np.fromfile(ifp, int_precision, 1)
            minute = np.fromfile(ifp, int_precision, 1)
            second = np.fromfile(ifp, int_precision, 1)
            dofyr = np.fromfile(ifp, int_precision, 1)
            pentad = np.fromfile(ifp, int_precision, 1)
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            # Populate date_time structure
            date_time_tmp = {
                'year': year,
                'month': month,
                'day': day,
                'hour': hour,
                'min': minute,
                'sec': second,
                'dofyr': dofyr,
                'pentad': pentad
            }
            date_time.append(date_time_tmp) 

            tmp_data2 = np.zeros((N_obs, 1))
            tmp_data2 = tmp_data2+pentad
            pentad2 = np.append(pentad2, tmp_data2)
            
            # Read observation assim flag
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            tmp_data = np.fromfile(ifp, logical_precision, N_obs)
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            tmp_data2 = np.zeros((N_obs, 1))
            indices = np.where(tmp_data != 0)[0]
            tmp_data2[indices] = 1
            obs_assim = np.append(obs_assim, tmp_data2)
            
            # Read species information
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_species = np.append(obs_species, np.fromfile(ifp, int_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            # Read tile number information
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_tilenum = np.append(obs_tilenum, np.fromfile(ifp, int_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read longitude
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_lon = np.append(obs_lon, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read latitude
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_lat = np.append(obs_lat, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            # Read observation value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_obs = np.append(obs_obs, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_obsvar = np.append(obs_obsvar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space model forecast value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_fcst = np.append(obs_fcst, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space model forecast variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_fcstvar = np.append(obs_fcstvar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space analysis value
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_ana = np.append(obs_ana, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)

            # Read observation-space analysis variance
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            obs_anavar = np.append(obs_anavar, np.fromfile(ifp, float_precision, N_obs))
            fortran_tag = np.fromfile(ifp, int_precision, 1)
            
            
    # Close file
        ifp.close()
     
    print('Total number of obs = ',len(obs_assim))

    return date_time, obs_species, obs_tilenum, obs_lon, obs_lat, obs_obs, obs_obsvar, obs_fcst, obs_fcstvar, obs_ana, obs_anavar, pentad2 

#####################################
def read_ascat_bufr(path, file_name, printflag=False):


    # Get a list of files with a similar name in a directory
    file_ext = '.bfr'
    files = [file for file in os.listdir(path) if file.startswith(file_name) and file.endswith(file_ext)]

    i = 1
    lat = []
    lon = []
    ssom = [] # Surface soil moisture
    tpcx = [] # Topographic complexity
    domo = [] # Direction of motion of moving observing platform
    smpf = [] # Soil moisture processing flag
    smcf = [] # Soil moisture correction flag
    alfr = [] # ASCAT land fraction
    iwfr = [] # Inundation And Wetland Fraction % 0-40-009
    snoc = [] # snow cover % 0-20-065
    flsf = [] # frozen land fraction %  0-40-008

    decoder = Decoder()

    # Open each file in turn
    for file in files:
        with open(os.path.join(path, file), 'rb') as fr:
            if printflag:
                print ('Reading file ', file, '...')
            
            # Loop for all messages
            for bufr_message in generate_bufr_message(decoder, fr.read()):
                #print ('Reading message number', i, '...')
                lat = np.append(lat, DataQuerent(NodePathParser()).query(bufr_message, '005001').all_values())
                lon = np.append(lon, DataQuerent(NodePathParser()).query(bufr_message, '006001').all_values())
                ssom = np.append(ssom, DataQuerent(NodePathParser()).query(bufr_message, '040001').all_values())
                # For QC-ing
                tpcx = np.append(tpcx, DataQuerent(NodePathParser()).query(bufr_message, '040010').all_values())
                domo = np.append(domo, DataQuerent(NodePathParser()).query(bufr_message, '001012').all_values())
                smpf = np.append(smpf, DataQuerent(NodePathParser()).query(bufr_message, '040006').all_values())
                smcf = np.append(smcf, DataQuerent(NodePathParser()).query(bufr_message, '040005').all_values())
                alfr = np.append(alfr, DataQuerent(NodePathParser()).query(bufr_message, '021166').all_values())
                iwfr = np.append(iwfr, DataQuerent(NodePathParser()).query(bufr_message, '040009').all_values())
                snoc = np.append(snoc, DataQuerent(NodePathParser()).query(bufr_message, '020065').all_values())
                flsf = np.append(flsf, DataQuerent(NodePathParser()).query(bufr_message, '040008').all_values())
                i = i + 1
        fr.close    
        
    print('Total number of obs = ',len(ssom))

    return lat, lon, ssom, tpcx, domo, smpf, smcf, alfr, iwfr, snoc, flsf


#####################################

def read_ascat_bufr_v2(path, file_name, printflag=False):
    file_ext = '.bfr'
    files = sorted([file for file in os.listdir(path) if file.startswith(file_name) and file.endswith(file_ext)])

    output = [[] for _ in range(11)]
    decoder = Decoder()

    for file in files:
        lat, lon, ssom, tpcx, domo, smpf, smcf, alfr, iwfr, snoc, flsf = [], [], [], [], [], [], [], [], [], [], []
        with open(os.path.join(path, file), 'rb') as fr:
            if printflag:
                print ('Reading file ', file, '...')
            
            for bufr_message in generate_bufr_message(decoder, fr.read()):
                lat.extend(DataQuerent(NodePathParser()).query(bufr_message, '005001').all_values())
                lon.extend(DataQuerent(NodePathParser()).query(bufr_message, '006001').all_values())
                ssom.extend(DataQuerent(NodePathParser()).query(bufr_message, '040001').all_values())
                tpcx.extend(DataQuerent(NodePathParser()).query(bufr_message, '040010').all_values())
                domo.extend(DataQuerent(NodePathParser()).query(bufr_message, '001012').all_values())
                smpf.extend(DataQuerent(NodePathParser()).query(bufr_message, '040006').all_values())
                smcf.extend(DataQuerent(NodePathParser()).query(bufr_message, '040005').all_values())
                alfr.extend(DataQuerent(NodePathParser()).query(bufr_message, '021166').all_values())
                iwfr.extend(DataQuerent(NodePathParser()).query(bufr_message, '040009').all_values())
                snoc.extend(DataQuerent(NodePathParser()).query(bufr_message, '020065').all_values())
                flsf.extend(DataQuerent(NodePathParser()).query(bufr_message, '040008').all_values())

        for i, lst in enumerate([lat, lon, ssom, tpcx, domo, smpf, smcf, alfr, iwfr, snoc, flsf]):
            output[i].extend(lst)

    output = [np.array(lst) for lst in output]

    print('Total number of files processed = ', len(files))

    return output

#####################################
def read_ascat_bufr_lat_lon(path, file_name, printflag=False):


    # Get a list of files with a similar name in a directory
    file_ext = '.bfr'
    files = [file for file in os.listdir(path) if file.startswith(file_name) and file.endswith(file_ext)]

    i = 1
    lat = []
    lon = []
    ssom = [] # Surface soil moisture

    decoder = Decoder()

    # Open each file in turn
    for file in files:
        with open(os.path.join(path, file), 'rb') as fr:
            if printflag:
                print ('Reading file ', file, '...')
            
            # Loop for all messages
            for bufr_message in generate_bufr_message(decoder, fr.read()):
                #print ('Reading message number', i, '...')
                lat = np.append(lat, DataQuerent(NodePathParser()).query(bufr_message, '005001').all_values())
                lon = np.append(lon, DataQuerent(NodePathParser()).query(bufr_message, '006001').all_values())
                ssom = np.append(ssom, DataQuerent(NodePathParser()).query(bufr_message, '040001').all_values())
                i = i + 1
        fr.close    
        
    print('Total number of obs = ',len(ssom))

    return lat, lon, ssom