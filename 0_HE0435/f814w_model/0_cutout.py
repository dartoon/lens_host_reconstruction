#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 10:11:32 2019

@author: Dartoon

"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'../../../my_code/py_tools/')
from cut_image import sub_bkg, cut_center_bright, cut_image, save_loc_png
from mask_objects import mask_obj
from flux_profile import flux_profile, profiles_compare
ID = 'HE0435'

#==============================================================================
# Cut out the lens image
#==============================================================================
filename = '../data/HE0435_F814W_sci.fits'
fitsFile = pyfits.open(filename)

from my_astro_func import read_pixel_scale
pixscl = read_pixel_scale(filename)  #Should be 0.03 here, drizzed from 0.05
print "pixscl", pixscl
# img = fitsFile[0].data # check the back grounp
# img[img==0.] = np.nan
# img_sub, bkg_light = sub_bkg(img)
# pyfits.PrimaryHDU(img_sub).writeto('../data/HE0435_F814W_sci_sublight.fits',overwrite=True)
img_sub = pyfits.getdata('../data/HE0435_F814W_sci_sublight.fits')

#%%
center_QSO = np.array([3641-1,2288-1]) #c_psf_list[QSO_loc]
lens_image = cut_image(img_sub, center_QSO, 100)

plt.imshow(lens_image[98:-98,98:-98], norm=LogNorm(),origin='low')
plt.colorbar()
plt.show()
#lens_bkg_light = cut_image(bkg_light, center_QSO, 100)
#plt.imshow(lens_bkg_light, norm=LogNorm(),origin='low')
#plt.colorbar()
#plt.show()

#%%
#==============================================================================
# #Estimate the exposure time and error map:
#==============================================================================
from my_astro_func import read_fits_exp
exp = read_fits_exp(filename) #9712.0
mean_wht = exp * (pixscl/0.05)**2  #!!! Still need to confirm
wht =pyfits.getdata('../data/HE0435_F814W_wht.fits')
exp_map = exp * wht/mean_wht
#stdd = np.std(img_sub[mask ==False])  # TOO SMALL
#stdd = 0.01245   #!!!!Need to be change
stdd = 0.0123
print exp,np.median(exp_map)

lens_exp_map =  cut_image(exp_map, center_QSO, 100)
lens_rms =(lens_image/lens_exp_map+stdd**2)**0.5
plt.imshow(lens_rms, norm=LogNorm(),origin='low')
plt.colorbar()
plt.show()
#%%
#==============================================================================
# Cutout and compare the PSF
#==============================================================================
# from my_astro_func import auto_cut_PSFs
# PSFx, PSFy  = auto_cut_PSFs(img_sub)
#%%
psfs_pos = np.array([[1488, 3211], [1000, 3992]])
dist_psfs = (psfs_pos-center_QSO)[:,0]**2+(psfs_pos-center_QSO)[:,1]**2
psfs_pos = psfs_pos[dist_psfs.argsort()]
count = 0
PSF_gauss_centers, PSF_bright_centers=[],[]
PSFs_ = []
for i in range(len(psfs_pos)):
   print 'PSF',count
   PSF, PSF_center = cut_center_bright(image=img_sub, center=psfs_pos[i], radius=60,  return_center=True, plot=False, center_float=True)
   PSF_gauss_centers.append(PSF_center)
   _, PSF_br_center = cut_center_bright(image=img_sub, center=psfs_pos[i], radius=60, kernel = 'center_bright', return_center=True, plot=False)
   PSF_bright_centers.append(PSF_br_center)
   #Create and save PSFs mask
   _, PSF_mask = mask_obj(img=PSF, exp_sz=1.4)
   if len(PSF_mask) > 1:
       PSF_mask = np.sum(np.asarray(PSF_mask),axis=0)
   elif len(PSF_mask) == 1:
       PSF_mask = PSF_mask[0]
   PSF_mask = (1 - (PSF_mask != 0)*1.)
   PSFs_.append([PSF,PSF_center,PSF_mask])
   count += 1
center_match = (np.sum(abs(np.asarray(PSF_gauss_centers)-np.asarray(PSF_bright_centers)<0.7),axis = 1) == 2)
psf_pos_list = [PSF_gauss_centers[i] for i in range(len(PSFs_))] # if center_match[i]==True]
PSFs = [PSFs_[i] for i in range(len(PSFs_))] # if center_match[i]==True]
save_loc_png(img_sub,center_QSO,psf_pos_list, ID=ID, label='PSF', v_max = 3.2, v_min = 1.e-3, reg_ty = 'acs')
PSF_list = [PSFs[i][0] for i in range(len(PSFs))]
PSF_masks = [PSFs[i][2] for i in range(len(PSFs))]
fig = profiles_compare(PSF_list, scal_list=np.ones(len(PSFs)),
                      prf_name_list=['PSF'+str(i) for i in range(len(PSFs))],
                      gridspace = 'log', if_annuli=True)
##%%
## =============================================================================
## Test the bkg level for lensed QSO and PSF
## =============================================================================
#r_flux, r_grids, _=flux_profile(lens_image, [len(lens_image)/2]*2,
#                                      radius=len(lens_image)/2, grids=50, ifplot=True, fits_plot= True)
#
#for i in range(len(PSF_list)):
#    _, _, _ = flux_profile(PSF_list[i], [len(PSF_list[i])/2]*2,
#                            radius=len(PSF_list[i])/2, grids=50, ifplot=True, fits_plot= True)
#
#==============================================================================
# Save the lens fitting ingredients including:
#    lensing image and noise map
#    PSFs
#==============================================================================
pyfits.PrimaryHDU(lens_image).writeto('{0}_cutout.fits'.format(ID),overwrite=True) 
pyfits.PrimaryHDU(lens_rms).writeto('{0}_stdd.fits'.format(ID),overwrite=True) 
for i in range(len(PSFs)):
   pyfits.PrimaryHDU(PSF_list[i]).writeto('PSF{0}.fits'.format(i),overwrite=True) 
