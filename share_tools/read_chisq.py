#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 16:18:24 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import sys
import lenstronomy.Util.class_creator as class_creator
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools')
from flux_profile import cr_mask

def return_chisq(filename, lens_mask = None, fair_mask=True, fair_PSF =False):
    result = pickle.load(open(filename,'rb'))
    fit_result, trans_result, kwargs_material, model_lists  = result
    
    kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask_saved = kwargs_material
    lens_model_list, source_model_list, lens_light_model_list, point_source_list = model_lists

    kwargs_result, chain_list = fit_result
    
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]    
    kwargs_psf = kwargs_psf_updated

    if "kernel_point_source_init" in kwargs_psf.keys() and fair_PSF==True:    
        kwargs_psf['kernel_point_source'] = kwargs_psf['kernel_point_source_init']
    #Delete the PSF noise map
    if "psf_error_map" in kwargs_psf.keys() and fair_mask==True:    
        del kwargs_psf["psf_error_map"]    
    
    if lens_mask is None:
        lens_mask = lens_mask_saved
    
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]    
    imageModel = class_creator.create_im_sim(multi_band_list = multi_band_list, multi_band_type='multi-linear', kwargs_model = kwargs_model,
                                               bands_compute = [True] * len(multi_band_list),
                                               likelihood_mask_list=[lens_mask],
                                               band_index=0)
    logL = imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
    n_data = imageModel.num_data_evaluate    
    reduced_x2 = - logL * 2 / n_data    
    return reduced_x2

#print return_chisq('../0_HE0435/model/fit_PSFi_PSFrecons/result_PSF0_PSFrecons_gammafix1.9_subg2.pkl')
#print "bulge, disk flux, source plane:", host_flux0_total, host_flux1_total

#import corner
#fig = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
##fig.savefig('fig_PSF{0}_{2}{1}_corner.pdf'.format(psfno,subg,fname))
#plt.show()