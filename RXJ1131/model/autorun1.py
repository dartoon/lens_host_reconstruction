#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 21:59:50 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os
path = os.getcwd()

psfno = 'avestd'

subg_list = [2,3]
for subg in subg_list:
#    print "psfno, subg:", psfno, subg
#    runfile(path+'/fit_PSFave_QSO_mask/2_modelling_QSOmask.py',
#            wdir=path+'/fit_PSFave_QSO_mask')
    
    print "psfno, subg:", psfno, subg
    runfile(path+'/fit_PSFave_PSFrecons/2_modelling_PSFrecons.py',
        wdir=path+'/fit_PSFave_PSFrecons')