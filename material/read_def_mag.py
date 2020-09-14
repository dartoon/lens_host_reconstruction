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

from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import ImageData as Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel

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
filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
               'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}
zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                  'F160W':25.9463,'ACS_F814W': 25.957}

inp = 0
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
result = pickle.load(open(filename,'rb'))# ,encoding="latin1") 
fit_values, chisq, labels, _ = result
chisq = [float(chisq[i]) for i in range(len(chisq))]

def make_circle_msk(img,x=30,y=30, radius=5):
    """
    Create a mask for the image, at position x, y with radius as radius.
    """
    yi,xi = np.indices((len(img),len(img))).astype(np.float64)
    mask = [np.sqrt((xi-x)**2+(yi-y)**2)<radius]
    return mask[0]

fluxs_result = []
for idx in range(len(chisq)):
    label = labels[idx]
    label = label.split(', ')
    filename = glob.glob('../'+IDs[inp]+'/model/{4}_PSFi_{0}/result_{1}_*_gammafix{2}_subg{3}.pkl'.format(label[-1], label[0], label[1][-3:], label[2][-1],folder_1) )
    readfile = filename[0] 
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
    mcmc_new_list, labels_new = trans_result
    
    lens_model_class = LensModel(lens_model_list=lens_model_list)
    lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
    source_model_class = LightModel(light_model_list=source_model_list)
    point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])
    
    data_class = Data(**kwargs_data)
    kwargs_psf = kwargs_psf_updated
    psf_class = PSF(**kwargs_psf)
    kwargs_model = kwargs_model
    kwargs_result = kwargs_result
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]
    
    #%%
    #define a apeature:
    pixel_radius = 0.6/data_class.pixel_width

    fluxs = []
    if len(kwargs_result['kwargs_lens_light'])>1:
        for i in range(len(kwargs_result['kwargs_lens_light'])):
            kwargs_result_copy = copy.deepcopy(kwargs_result)
        #    #Cal Bulge and then disk:
            kwargs_result_copy['kwargs_lens_light'][1-i]['amp'] = 0
            imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lens_light_model_class,
                                            point_source_class=point_source_class, kwargs_numerics=kwargs_numerics)
            lens_profile = imageModel.image(kwargs_lens_light=kwargs_result_copy['kwargs_lens_light'], unconvolved=False)
            mask = make_circle_msk(lens_profile, x = len(lens_profile)/2, y=len(lens_profile)/2, radius = pixel_radius )
            fluxs.append(np.sum(lens_profile*mask) )
    imageModel = ImageModel(data_class, psf_class, lens_light_model_class=lens_light_model_class,
                                    point_source_class=point_source_class, kwargs_numerics=kwargs_numerics)
    lens_profile = imageModel.image(kwargs_ps=kwargs_ps, kwargs_lens_light=kwargs_result['kwargs_lens_light'], unconvolved=False )
#    lens_profile = data_class.data
    fluxs.append(np.sum(lens_profile*mask) )
    fluxs_result.append(fluxs)         

fluxs_result = np.array(fluxs_result)
zp = zeropoint_dict[filter_dict[ID]]
fluxs_result[fluxs_result<0] = 0.0001

mag_result = -2.5*np.log10(fluxs_result)+zp

#%%   
mag_result_copy = copy.deepcopy(mag_result)
read_i = 2  #0 bulge, 1 disk, 2 total
count_n = 8 #len(labels) * 10 /100
chisq = [float(chisq[i]) for i in range(len(chisq))]
sort_chisq = np.argsort(np.asarray(chisq))    
Chisq_best = chisq[sort_chisq[0]]
Chisq_last= chisq[sort_chisq[count_n-1]]
inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
weight = np.zeros(len(chisq))
for i in sort_chisq[:count_n]:
    weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
weighted_value = np.sum(np.array(mag_result_copy[:,read_i])*weight) / np.sum(weight)   
rms_value = np.sqrt(np.sum((np.array(mag_result_copy[:,read_i])-weighted_value)**2*weight) / np.sum(weight)) 
print "{0:.2f}".format(weighted_value), "+-" , " {0:.2f}".format(rms_value)
