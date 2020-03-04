import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
import copy
IDs = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']



import sys
sys.path.insert(0,'../share_tools/')
from read_inference import read_fnu
#%%Translate flux into f_nu 

filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
               'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}

zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                  'F160W':25.9463,'ACS_F814W': 25.957}

fnu_list =  []
for i in range(len(IDs)):
    ID = IDs[i][2:]
    print(IDs[i][2:],":\t",filter_dict[ID], "\t",read_fnu(ID,count_n=[4, 4]))
#    if ID != "RXJ1131":
#        result = read_inf(ID, prop=0 , count_n=8)
#        weighted_value, rms_value = result[0][0], result[0][1]
#        uncer = rms_value/weighted_value
#        zp = zeropoint_dict[filter_dict[ID]]
#        mag = -2.5*np.log10(weighted_value) + zp
#        fnu = 10 ** ((mag-25)/(-2.5))
#        print(IDs[i][2:],":\t",filter_dict[ID], "\t",round(fnu,3), '+-',round(fnu*uncer,3))
#    else:
#        names = ['RXJ1131_bulge', 'RXJ1131_disk']
#        result = read_inf(ID, prop=0 , count_n=8)
#        for j in range(2):
#            weighted_value, rms_value = result[j][0], result[j][1]
#            uncer = rms_value/weighted_value
#            zp = zeropoint_dict[filter_dict[ID]]
#            mag = -2.5*np.log10(weighted_value) + zp
#            fnu = 10 ** ((mag-25)/(-2.5))
#            print(names[j],":\t",filter_dict[ID], "\t",round(fnu,3), '+-',round(fnu*uncer,3))    
