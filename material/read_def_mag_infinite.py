#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:50:31 2020

@author: Dartoon
"""
#!!! For the Lenstronomy version 1.5.0

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
from matplotlib.colors import LogNorm
import copy
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/my_code/py_tools/')
from flux_profile import cr_mask
import glob


# # =============================================================================
# # To reproduce
# # =============================================================================
# In lenstronomy version 1.3.0
# 1.Replae model_band_plot
# to change label size:
# 2.in plot._util.py
# change text_description()   


#%% Find the one with lowest Chisq
IDs = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']


#for i in range(len(IDs)):
#    print(i, ':', IDs[i])
#inp = int(input("Which source:\n"))
inp = 3
ID = IDs[inp][2:]
if ID == 'SDSS1206':
    folder_0 = '_singel_host_run_summary.pkl'
    folder_1 = 'singel_fit'
elif ID == 'HS2209' or ID == 'RXJ1131':
    folder_0 = '_3rd_run_summary.pkl'
    folder_1 = '3rd_fit'    
else:
    folder_0 = '_2nd_run_summary.pkl'
    folder_1 = '2nd_fit'
filename = '../share_tools/'+ ID+ folder_0
result = pickle.load(open(filename,'rb')) 
fit_values, chisq, labels, _ = result
chisq = [float(chisq[i]) for i in range(len(chisq))]

mag_result = []
for idx in range(len(chisq)):
    label = labels[idx]
    label = label.split(', ')
    filename = glob.glob('../'+IDs[inp]+'/model/{4}_PSFi_{0}/result_{1}_*_gammafix{2}_subg{3}.pkl'.format(label[-1], label[0], label[1][-3:], label[2][-1],folder_1) )
    readfile = filename[0] 
    result = pickle.load(open(readfile,'rb'))# ,encoding="latin1")
    fit_result, trans_result, kwargs_material, model_lists  = result
    #kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask = kwargs_material
    lens_model_list, source_model_list, lens_light_model_list, point_source_list = model_lists
    kwargs_result, chain_list = fit_result
    
    from lenstronomy.LightModel.light_model import LightModel
    light = LightModel(lens_light_model_list)
    fluxs = light.total_flux(kwargs_result['kwargs_lens_light'])
    
    filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
                   'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}
    
    zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                      'F160W':25.9463,'ACS_F814W': 25.957}
    zp = zeropoint_dict[filter_dict[ID]]
    f = np.sum(fluxs)
    if f <0:
        f = 100
    mag = -2.5*np.log10(np.sum(f))+zp
    mag_result.append(mag)

#%%   
count_n = 8 #len(labels) * 10 /100
chisq = [float(chisq[i]) for i in range(len(chisq))]
sort_chisq = np.argsort(np.asarray(chisq))    
Chisq_best = chisq[sort_chisq[0]]
Chisq_last= chisq[sort_chisq[count_n-1]]
inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
weight = np.zeros(len(chisq))
for i in sort_chisq[:count_n]:
    weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
weighted_value = np.sum(np.array(mag_result)*weight) / np.sum(weight)   
rms_value = np.sqrt(np.sum((np.array(mag_result)-weighted_value)**2*weight) / np.sum(weight)) 
print(weighted_value, rms_value)