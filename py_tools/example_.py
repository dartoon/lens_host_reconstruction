#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 14:04:02 2018

@author: Dartoon

Cut PSF and QSO for xxx
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'../../../py_tools')
from cut_image import cut_image, cut_center_bright, save_loc_png, grab_pos
import copy
import astropy.io.fits as pyfits

import os
path = os.getcwd()

ID = path.split('/')[-2]
fitsFile = pyfits.open('../astrodrz/final_drz.fits')

#img = fitsFile[1].data # check the back grounp
#from astropy.visualization import SqrtStretch
#from astropy.stats import SigmaClip
#from photutils import Background2D, SExtractorBackground  
#from astropy.visualization.mpl_normalize import ImageNormalize
#import matplotlib.pyplot as plt
#norm = ImageNormalize(stretch=SqrtStretch())         
#sigma_clip = SigmaClip(sigma=3., iters=10)
#bkg_estimator = SExtractorBackground()
#from photutils import make_source_mask
#mask_0 = make_source_mask(img, snr=2, npixels=5, dilate_size=11)
#mask_1 = (np.isnan(img))
#mask = mask_0 + mask_1
#bkg = Background2D(img, (50, 50), filter_size=(3, 3),
#                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
#                   mask=mask)
#from matplotlib.colors import LogNorm
#fig=plt.figure(figsize=(15,15))
#ax=fig.add_subplot(1,1,1)
#ax.imshow(img, norm=LogNorm(), origin='lower') 
##bkg.plot_meshes(outlines=True, color='#1f77b4')
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
#plt.show()  
#fig=plt.figure(figsize=(15,15))
#ax=fig.add_subplot(1,1,1)
#ax.imshow(mask, origin='lower') 
##bkg.plot_meshes(outlines=True, color='#1f77b4')
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
#plt.show()  
#
#back = bkg.background* ~mask_1
#fig=plt.figure(figsize=(15,15))
#ax=fig.add_subplot(1,1,1)
#ax.imshow(back, origin='lower', cmap='Greys_r')
#ax.xaxis.set_visible(False)
#ax.yaxis.set_visible(False)
#plt.show()
#
#img -= back              
#pyfits.PrimaryHDU(img).writeto('sub_coadd.fits',overwrite=True)
img = pyfits.getdata('sub_coadd.fits')

filename= '{0}.reg'.format(ID)
#c_psf_list, QSO_loc = grab_pos(filename,reg_ty = 'astrodrz_06', QSO_reg_return=True)
center_QSO = np.array([871,1163]) #c_psf_list[QSO_loc]

QSO, cut_center = cut_center_bright(image=img, center=center_QSO, radius=60, return_center=True, plot=False)
QSO_outer = cut_image(image=img, center=cut_center, radius=200)
pyfits.PrimaryHDU(QSO).writeto('{0}_cutout.fits'.format(ID),overwrite=True)
pyfits.PrimaryHDU(QSO_outer).writeto('{0}_cutout_outer.fits'.format(ID),overwrite=True)

wht = fitsFile[2].data # - (-0.002)  # check the back grounp
cut_wht = cut_image(image=wht, center=cut_center, radius=60)
pyfits.PrimaryHDU(cut_wht).writeto('wht_map.fits',overwrite=True)

PSFs = []
PSF_gauss_centers = []
PSF_bright_centers = []

count=0
psf_list = None
#psf_list = np.delete(c_psf_list, (QSO_loc), axis=0)
#dist = (psf_list-center_QSO)[:,0]**2+(psf_list-center_QSO)[:,1]**2
#psf_list = psf_list[dist.argsort()]
#for i in range(len(psf_list)):  
#    print 'PSF',i
#    PSF, PSF_center = cut_center_bright(image=img, center=psf_list[i], radius=60, return_center=True, plot=False)
#    PSFs.append([PSF, 1, PSF_center])
#    PSF_gauss_centers.append(PSF_center)
#    _, PSF_br_center = cut_center_bright(image=img, center=psf_list[i], radius=60, kernel = 'center_bright', return_center=True, plot=False)
#    PSF_bright_centers.append(PSF_br_center)
#    count += 1

#extra_psfs = None
extra_psfs = np.array([[2207,931],[2097,1329],[1209,1437]])
dist_extra = (extra_psfs-center_QSO)[:,0]**2+(extra_psfs-center_QSO)[:,1]**2
extra_psfs = extra_psfs[dist_extra.argsort()]
for i in range(len(extra_psfs)):
    print 'PSF',count
    PSF, PSF_center = cut_center_bright(image=img, center=extra_psfs[i], radius=60,  return_center=True, plot=False, center_float=True)
    PSFs.append([PSF,0, PSF_center])
    PSF_gauss_centers.append(PSF_center)
    _, PSF_br_center = cut_center_bright(image=img, center=extra_psfs[i], radius=60, kernel = 'center_bright', return_center=True, plot=False)
    PSF_bright_centers.append(PSF_br_center)
    count += 1

from mask_objects import mask_obj
print "QSO:"
a, QSO_mask = mask_obj(img=QSO, exp_sz=1.4)
if len(QSO_mask) > 1:
    QSO_mask = np.sum(np.asarray(QSO_mask),axis=0)
elif len(QSO_mask) == 1:
    QSO_mask = QSO_mask[0]
#print "QSO image:"
#plt.imshow((QSO_mask), origin='lower')
#plt.show()
QSO_mask = (1 - (QSO_mask != 0)*1.)

for i in range(len(PSFs)):
    print "PSF{0}:".format(i)
    _, PSF_mask = mask_obj(img=PSFs[i][0], exp_sz=1.4)
    if len(PSF_mask) > 1:
        PSF_mask = np.sum(np.asarray(PSF_mask),axis=0)
    elif len(PSF_mask) == 1:
        PSF_mask = PSF_mask[0]
#    print "PSF{0} image:".format(i)
#    plt.imshow(PSF_mask, origin='lower')
#    plt.show()
    PSF_mask = (1 - (PSF_mask != 0)*1.)
    PSFs[i].append(PSF_mask)

#center_match = (np.sum(abs(np.asarray(PSF_gauss_centers)-np.asarray(PSF_bright_centers)),axis = 1) == 0)
center_match = (np.sum(abs(np.asarray(PSF_gauss_centers)-np.asarray(PSF_bright_centers)<0.7),axis = 1) == 2)
PSF_c = len(PSFs[0][0])/2
bright_match = [(PSFs[i][0].max()-PSFs[i][0][PSF_c, PSF_c])/PSFs[i][0].max() < 0.1 for i in range(len(PSFs))]
PSFs_all = copy.deepcopy(PSFs)
PSFs=[]
for i in range(len(PSFs_all)):
    if center_match[i] == True and bright_match[i] == True:
        print i
        if np.sum(PSFs_all[i][0] * PSFs_all[i][3]) >100:
            PSFs.append(PSFs_all[i])

#del_list = [0]
#PSFs = [PSFs[i] for i in range(len(PSFs)) if i not in del_list]

#==============================================================================
# Compare the profile and derive the Average image
#==============================================================================
flux_list = []
for i in range(len(PSFs)):
    flux = np.sum(PSFs[i][0]*PSFs[i][3])
    print "tot_flux for PSF{0}".format(i), flux
    flux_list.append(flux)

#plot the first selection
if extra_psfs is None:
    save_loc_png(img,center_QSO,psf_list, ID=ID, label='ini' ,reg_ty = 'astrodrz_06')
else:
    save_loc_png(img,center_QSO,psf_list,extra_psfs, ID=ID, label='ini', reg_ty = 'astrodrz_06')

PSFs_familiy = [PSFs[i][1] for i in range(len(PSFs))]
if extra_psfs is None:
    loc_PSFs = psf_list
elif psf_list is None:
    loc_PSFs = extra_psfs
else:
    loc_PSFs = np.append(psf_list, extra_psfs, axis=0)
loc_ind_star =  [PSFs[i][2] for i in range(len(PSFs)) if PSFs[i][1]==1] #and flux_list[i]>100]
loc_like_star =  [PSFs[i][2] for i in range(len(PSFs)) if PSFs[i][1]==0] # and flux_list[i]>100]

if PSFs_familiy[-1] ==1:
    save_loc_png(img,center_QSO,loc_ind_star, ID=ID,reg_ty = 'astrodrz_06')
else:
    save_loc_png(img,center_QSO,loc_ind_star,loc_like_star, ID=ID,reg_ty = 'astrodrz_06')

PSF_list = [PSFs[i][0] for i in range(len(PSFs))]
PSF_masks = [PSFs[i][3] for i in range(len(PSFs))]
from flux_profile import QSO_psfs_compare
gridsp_l = ['log', None]
if_annuli_l = [False, True] 
for i in range(2):
    for j in range(2):
        plt_which_PSF = None
        plt_QSO = False
#        if i+j == 0:
#            plt_which_PSF = range(len(PSFs))
#            plt_QSO = True
        fig_psf_com = QSO_psfs_compare(QSO=QSO, QSO_msk=QSO_mask, psfs= PSF_list,
                                       plt_which_PSF=plt_which_PSF,
                                       PSF_mask_img=PSF_masks, grids=30,
                                       include_QSO=True, 
                                       plt_QSO = plt_QSO, norm_pix = 5.0,
                                       gridspace= gridsp_l[i], if_annuli=if_annuli_l[j])
#        fig_psf_com.savefig('PSFvsQSO{0}_{1}_{2}.pdf'.format(i,['xlog','xlin'][i],['circ','annu'][j]))
        if j==1:
            plt.show()
        else:
            plt.close()

import pickle
filename='{0}_PSFs_QSO'.format(ID)
datafile = open(filename, 'wb')
QSOs = [QSO,cut_center,QSO_outer]
pickle.dump([PSFs, QSOs], open(filename, 'wb'))
datafile.close()

#import pickle
#datafile = open('{0}_PSFs_QSO'.format(ID),'rb')
#PSFs, QSO=pickle.load(open('XID2202_PSFs_QSO','rb'))
#datafile.close()

