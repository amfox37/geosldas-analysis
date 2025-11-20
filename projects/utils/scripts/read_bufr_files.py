#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np

from my_functions import read_ascat_bufr_lat_lon

path = '/discover/nobackup/amfox/ASCAT_file_mover/Metop-B/Y2023/M01'
file_name_start = 'M01-ASCA-ASCSMR02-NA-5.0-2023012'

printflag = True

lat, lon, ssom = read_ascat_bufr_lat_lon(path, file_name_start, printflag)


numobs = len(ssom)
obarray = np.empty([numobs, 3])

obarray[:, 0] = ssom
obarray[:, 1] = lon
obarray[:, 2] = lat

obarray[np.isnan(obarray)] = 120

# Save the array to a CSV file
np.savetxt('Metop_B_Y2023_M012_obarray.csv', obarray, delimiter=',')


# In[ ]:




