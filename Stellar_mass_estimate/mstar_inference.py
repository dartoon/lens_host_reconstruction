#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:25:07 2020

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

folders = ['0_HE0435', '1_RXJ1131_bulge', '1_RXJ1131_disk', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']

#%%
mstar_temp_dict = {}
for folder in folders:
    filename  = '../Stellar_mass_estimate/'+ folder + '/summary_1_PA00.fits'
    hdul_sum = pyfits.open(filename)
    name_sum = hdul_sum[1].columns
    table_sum = hdul_sum[1].data
    z_idx = [i for i in range(len(name_sum)) if 'zmc' in str(name_sum[i])][0]
#    
    filename  = '../Stellar_mass_estimate/'+ folder + '/SFH_1_PA00_param.fits'
    hdul = pyfits.open(filename)
    name = hdul[1].columns
    table = hdul[1].data
    mstar_idx = [i for i in range(len(name)) if 'Mstel' in str(name[i])][0]
    age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
#    print(folder,':\n', 'M_star:', round(table[1][mstar_idx],3), 'AGE:', str(round(10**table[1][age_idx],2))+' Gyrs', "redshift", table_sum[1][z_idx])
    mstar_temp_dict[folder] = round(table[1][mstar_idx],3)

systems = ['HE0435', 'RXJ1131', 'WFI2033', 'SDSS1206', 'HE1104', 'SDSS0246', 'HS2209', 'HE0047']
mstar_dict = {}
for ID in systems:
    if ID  == 'RXJ1131':
        mstar_dict[ID] = [mstar_temp_dict['1_RXJ1131_bulge'], mstar_temp_dict['1_RXJ1131_disk']]
    else:
        folder_ID = [folder for folder in folders if ID in folder][0]
        mstar_dict[ID] = mstar_temp_dict[folder_ID]

#%%     
mstar_dict = {'HE0435': 10.713,
 'RXJ1131': [10.264, 11.082],
 'WFI2033': 10.668,
 'SDSS1206': 10.779,
 'HE1104': 11.036,
 'SDSS0246': 10.686,
 'HS2209': 11.182,
 'HE0047': 10.913}