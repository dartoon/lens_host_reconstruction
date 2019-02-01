#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 20:15:54 2019

@author: dxh

Read result from pickle
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
#PSFno = 
#subg= 2


#The values to pick up:
filename = "fit_PSFi_QSO_mask/result_PSF3_Dhost_Dlens_subg2"

result = pickle.load(open(filename,'rb'))

if len(result) == 2:
    fit_result, trans_result = result
elif len(result) ==3:
    fit_result, trans_result, kwargs_psf_updated = result