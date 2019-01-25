#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 10:44:07 2019

@author: Dartoon

Test load text from PDF
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
from matplotlib.colors import LogNorm
import copy
import sys
sys.path.insert(0,'../../py_tools')
from flux_profile import cr_mask

from lenstronomy.ImSim.image_model import ImageModel

psf_i = raw_input("Which psf to load? 0, 1, 2, 3...:\n")

lens_result, source_result, lens_light_result, ps_result, cosmo_result,\
chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = pickle.load(open('fit_PSFi_QSO_mask/fit_result_PSF{0}_QSOmask_gammafix_subg2'.format(psf_i),'rb'))

print "lens_result:", lens_result,'\n'
print "source_light_result:", source_result,'\n'
print "lens_light_result:", lens_light_result,'\n'
print "point_source_result:", ps_result,'\n'



imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)
