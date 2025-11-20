import os
import eccodes as ecc
import numpy as np

def read_ascat_bufr_ec(path, file_name, printflag=False):
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

    # Open each file in turn
    for file in files:
        with open(os.path.join(path, file), 'rb') as f:
            if printflag:
                print ('Reading file ', file, '...')
            # Create a BUFR message iterator
            bufr = ecc.codes_bufr_new_from_file(f)
            while bufr is not None:
                #print ('Reading message number', i, '...')
                lat_array = ecc.codes_get_array(bufr, '005001')
                lat = np.append(lat, lat_array)
                lon_array = ecc.codes_get_array(bufr, '006001')
                lon = np.append(lon, lon_array)
                ssom_array = ecc.codes_get_array(bufr, '040001')
                ssom = np.append(ssom, ssom_array)
                tpcx_array = ecc.codes_get_array(bufr, '040010')
                tpcx = np.append(tpcx, tpcx_array)
                domo_array = ecc.codes_get_array(bufr, '001012')
                domo = np.append(domo, domo_array)
                smpf_array = ecc.codes_get_array(bufr, '040006')
                smpf = np.append(smpf, smpf_array)
                smcf_array = ecc.codes_get_array(bufr, '040005')
                smcf = np.append(smcf, smcf_array)
                alfr_array = ecc.codes_get_array(bufr, '021166')
                alfr = np.append(alfr, alfr_array)
                iwfr_array = ecc.codes_get_array(bufr, '040009')
                iwfr = np.append(iwfr, iwfr_array)
                snoc_array = ecc.codes_get_array(bufr, '020065')
                snoc = np.append(snoc, snoc_array)
                flsf_array = ecc.codes_get_array(bufr, '040008')
                flsf = np.append(flsf, flsf_array)
                i += 1
                ecc.codes_release(bufr)
                bufr = ecc.codes_bufr_new_from_file(f)
            f.close()    
        
    print('Total number of obs = ',len(ssom))

    return lat, lon, ssom, tpcx, domo, smpf, smcf, alfr, iwfr, snoc, flsf