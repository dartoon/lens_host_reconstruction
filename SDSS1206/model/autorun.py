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

#psfnos = [0, 1]
psfnos = [3]
subg_list = [3]
for psfno in psfnos:
    for subg in subg_list:
        print "psfno, subg:", psfno, subg
        runfile(path+'/fit_PSFi_PSFrecons/2_modelling_PSFrecons.py',
                wdir=path + '/fit_PSFi_PSFrecons')
        runfile(path+'/fit_PSFi_QSOmask/2_modelling_QSOmask.py', 
                wdir=path+'/fit_PSFi_QSOmask')
