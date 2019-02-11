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

#psfnos = [0, 2]
#psfnos = [5, 6]
#psfnos = [7]
#subg_list = [3]
for psfno in psfnos:
    subg = 2
    print "psfno, subg:", psfno, subg
    runfile(path+'/2_modelling_PSFrecons_Dlenslight.py',
            wdir=path)
