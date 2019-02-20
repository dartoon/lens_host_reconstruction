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

psfnos = ['ave','avestd']
subg_list = [2,3]
for psfno in psfnos:
    for subg in subg_list:
        print "psfno, subg:", psfno, subg
        print "run PSF ave"
        runfile(path+'/fit_PSFave_PSFrecons/2_modelling_PSFrecons.py',
                wdir=path + '/fit_PSFave_PSFrecons')
        print "run PSF mask"
        runfile(path+'/fit_PSFave_QSOmask/2_modelling_QSOmask.py', 
                wdir=path+'/fit_PSFave_QSOmask')
