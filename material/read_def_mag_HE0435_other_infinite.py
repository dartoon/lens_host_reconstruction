#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:50:31 2020

@author: Dartoon
"""


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

#%% Find the one with lowest Chisq
IDs = ['0_HE0435']


#for i in range(len(IDs)):
#    print(i, ':', IDs[i])
#0 = int(0ut("Which source:\n"))
ID = IDs[0][2:]
#file_160w = '/model/2nd_fit_PSFi_PSFrecons/result_PSF0_*_gammafix2.1_subg2.pkl'
filter = '814'
file = '/f{0}w_model/fit_PSFi_PSFrecons/result_PSF*_gammafix*.pkl'.format(filter)
#
mag_result = []
filename = glob.glob('../'+IDs[0]+file)
for idx in range(len(filename)-1):
    idx = idx+1
    readfile = filename[idx]
    result = pickle.load(open(readfile,'rb'))# ,encoding="latin1")
    fit_result, trans_result, kwargs_material, model_lists  = result
    kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask = kwargs_material
    lens_model_list, source_model_list, lens_light_model_list, point_source_list = model_lists
    kwargs_result, chain_list = fit_result
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]    
    lens_result = kwargs_result['kwargs_lens']
    lens_light_result = kwargs_result['kwargs_lens_light']
    source_result = kwargs_result['kwargs_source']
    ps_result = kwargs_result['kwargs_ps']
    from lenstronomy.LightModel.light_model import LightModel
    light = LightModel(lens_light_model_list)
    fluxs = light.total_flux(kwargs_result['kwargs_lens_light'])
    
    filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
                   'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}
    
    zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                      'F160W':25.9463,'ACS_F814W': 25.957}
    zp = zeropoint_dict['F{0}W'.format(filter)]
    f = np.sum(fluxs[0])
    if f <0:
        f = 0.0001
    mag = -2.5*np.log10(np.sum(f))+zp
    mag_result.append(mag)
print(np.mean(mag_result), np.std(mag_result))
