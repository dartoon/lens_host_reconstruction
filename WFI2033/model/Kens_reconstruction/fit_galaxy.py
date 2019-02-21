#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 11:40:32 2019

@author: Dartoon

Fit as galaxy
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools')
from fit_qso import fit_galaxy

#galaxy_img = pyfits.getdata('flux-esmod_2band_composall_best_es002_sr.fits')
galaxy_img = pyfits.getdata('flux-esmod_2band_fid3_cov6_best_es002_sr.fits')
galaxy_msk = [galaxy_img!=0][0]

fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []
fixed_source.append({})  
kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0, 'center_y': 0})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

source_result, image_host, error_map, reduced_Chisq=fit_galaxy(galaxy_img, psf_ave=None, psf_std = None,
                                              background_rms=np.std(galaxy_img)/2,
                                              source_params=source_params, galaxy_msk = None,
                                              pix_sz = 0.048, no_MCMC =False,
                                              tag=None, deep_seed= True, pltshow=1, return_Chisq=True)
source_result = source_result[0]
import magVSamp
zp_assum = 25
mag = magVSamp.getMag(source_result['amp'],source_result,zp=25,deltaPix=0.08)
mag_flux = 10.**(-0.4*(mag-zp_assum))

print "Ken_image_flux:", galaxy_img.sum()
print "Sersic_total_flux:", mag_flux
print "Sersic_frame_flux", image_host[0].sum()

print "fit result:", source_result