#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:17:50 2019

@author: Dartoon

The pixel scale information
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

pixel_scale = {'HE0435' : 0.08, 'RXJ1131' : 0.05, 'WFI2033' : 0.08, 'SDSS1206': 0.08, 'HE1104': 0.08,
               'SDSS0246': 0.03, 'HS2209': 0.03, 'HE0047' : 0.03}

lens_redshift = {'HE0435' : 1.69, 'RXJ1131' : 0.65, 'WFI2033' : 1.66, 'SDSS1206': 1.79, 'HE1104': 2.32,
               'SDSS0246':  1.68, 'HS2209': 1.07, 'HE0047' : 1.66}

filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
               'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}

zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                  'F160W':25.9463,'ACS_F814W': 25.957}
