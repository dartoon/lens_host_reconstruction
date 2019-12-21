#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:41:16 2019

@author: Dartoon

Clean the PSF then saved, make it to be loaded easier for the future.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import copy
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
from flux_profile import profiles_compare, cr_mask_img
import glob

psfname_list = glob.glob('PSF?.fits') 
psfname_list.sort()
psfs = []
for i in range(len(psfname_list)):
    print "PSF:", psfname_list[i]
    psf_i = pyfits.getdata(psfname_list[i])
    print psf_i.sum()
    mask_files = glob.glob('PSF{0}_mask*.reg'.format(i))
    if mask_files != []:
        psf_mask = cr_mask_img(psf_i, mask_list=mask_files)
#        psf_mask = (1 - (psf_mask != 0)*1.)
        psf_clear = copy.deepcopy(psf_i)
        for i in range(len(psf_i)):
            for j in range(len(psf_i)):
                if psf_mask[i, j] == 0:
                    psf_clear[i, j] = (psf_i[-i-1, -j-1]*psf_mask[-i-1, -j-1])
        print "plot PSF after clear:"
        psf_i = psf_clear
    plt.imshow(psf_i, origin = 'low', norm=LogNorm())
    plt.show()
    psfs.append(psf_i)

fig = profiles_compare(psfs, scal_list=np.ones(len(psfs)),
                       prf_name_list=['PSF'+str(i) for i in range(len(psfs))],
                       gridspace = 'log', if_annuli=True)
fig.savefig("PSFs_profile.pdf")
fig.show()
    
#for i in range(len(psfs)):
#    psfs[i] = psfs[i]/np.sum(psfs[i])
#    pyfits.PrimaryHDU(psfs[i]).writeto('PSF{0}_clean.fits'.format(i),overwrite=True) 
#
#%%
import os   
if glob.glob('fit_PSFi_PSFrecons')  == []:
    os.mkdir('fit_PSFi_PSFrecons')
if glob.glob('fit_PSFi_QSOmask')  == []:
    os.mkdir('fit_PSFi_QSOmask')
