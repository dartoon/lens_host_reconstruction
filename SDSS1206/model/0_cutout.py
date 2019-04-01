#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 10:11:32 2019

@author: Dartoon

1. Remove SDSS1206 background light.
2. Cutout SDSS1206 image.
3. Select potential PSFs and compare the SDSS1206 ones.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'../../py_tools')
from cut_image import sub_bkg, cut_center_bright, cut_image, save_loc_png
from mask_objects import mask_obj
from flux_profile import flux_profile, profiles_compare
ID = 'SDSS1206'

#==============================================================================
# Cut out the lens image
#==============================================================================
fitsFile = pyfits.open('../data/SDSS1206+4332_F160W_sci.fits')
img = fitsFile[0].data # check the back grounp
#img[img==0.] = np.nan
#img_sub, bkg_light = sub_bkg(img)
#pyfits.PrimaryHDU(img_sub).writeto('../data/SDSS1206+4332_sublight.fits',overwrite=True)
img_sub = pyfits.getdata('../data/SDSS1206+4332_sublight.fits')
center_QSO = np.array([1209,1348]) #c_psf_list[QSO_loc]

lens_image = cut_image(img_sub, center_QSO, 100)
plt.imshow(lens_image, norm=LogNorm(),origin='low')
plt.colorbar()
plt.show()
#lens_bkg_light = cut_image(bkg_light, center_QSO, 100)
#plt.imshow(lens_bkg_light, norm=LogNorm(),origin='low')
#plt.colorbar()
#plt.show()

#==============================================================================
# #Estimate the exposure time and error map:
#==============================================================================
wht = pyfits.getdata('../data/SDSS1206+4332_F160W_wht.fits')
exp = 8456.8892
mean_wht = exp * (0.08/0.12825)**2
exp_map = exp * wht/mean_wht
from photutils import make_source_mask
mask_0 = make_source_mask(img_sub, snr=2, npixels=5, dilate_size=11)
mask_1 = (np.isnan(img))
mask = mask_0 + mask_1
plt.imshow(img_sub * (1-mask*1), norm=LogNorm(), origin='low')
plt.show()
#stdd = np.std(img_sub[mask ==False])  # TOO SMALL
stdd = 0.005
print "stdd:",stdd
#filename = '{0}_cutout.fits'.format(ID)
#image = pyfits.open(filename)
lens_exp_map =  cut_image(exp_map, center_QSO, 100)
lens_rms =(lens_image/lens_exp_map+stdd**2)**0.5
plt.imshow(lens_rms, norm=LogNorm(),origin='low')
plt.colorbar()
plt.show()

#==============================================================================
# Cutout and compare the PSF
#==============================================================================
psfs_pos = np.array([[1112,2151]])
dist_psfs = (psfs_pos-center_QSO)[:,0]**2+(psfs_pos-center_QSO)[:,1]**2
psfs_pos = psfs_pos[dist_psfs.argsort()]
count = 0
PSF_gauss_centers, PSF_bright_centers=[],[]
PSFs_ = []
for i in range(len(psfs_pos)):
    print 'PSF',count
    PSF, PSF_center = cut_center_bright(image=img_sub, center=psfs_pos[i], radius=40,  return_center=True, plot=False, center_float=True)
    PSF_gauss_centers.append(PSF_center)
    _, PSF_br_center = cut_center_bright(image=img_sub, center=psfs_pos[i], radius=40, kernel = 'center_bright', return_center=True, plot=False)
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
save_loc_png(img_sub,center_QSO,psf_pos_list, ID=ID, label='PSF-' ,reg_ty = None)
PSF_list = [PSFs[i][0] for i in range(len(PSFs))]
PSF_masks = [PSFs[i][2] for i in range(len(PSFs))]
#fig = profiles_compare(PSF_list, scal_list=np.ones(len(PSFs)),
#                       prf_name_list=['PSF'+str(i) for i in range(len(PSFs))],
#                       gridspace = 'log', if_annuli=True)
## =============================================================================
## Test the bkg level for lensed QSO and PSF
## =============================================================================
#r_flux, r_grids, _=flux_profile(lens_image, [len(lens_image)/2]*2,
#                                      radius=len(lens_image)/2, grids=50, ifplot=True, fits_plot= True)
#
#for i in range(len(PSF_list)):
#    _, _, _ = flux_profile(PSF_list[i], [len(PSF_list[i])/2]*2,
#                            radius=len(PSF_list[i])/2, grids=50, ifplot=True, fits_plot= True)

##==============================================================================
## Save the lens fitting ingredients including:
##    lensing image and noise map
##    PSFs
##==============================================================================
#pyfits.PrimaryHDU(lens_image).writeto('SDSS1206_cutout.fits',overwrite=True) 
#pyfits.PrimaryHDU(lens_rms).writeto('SDSS1206_stdd.fits',overwrite=True) 
#for i in range(len(PSFs)):
#    pyfits.PrimaryHDU(PSF_list[i]).writeto('PSF{0}.fits'.format(i),overwrite=True) 
#
