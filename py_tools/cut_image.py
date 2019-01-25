#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 11:45:45 2018

@author: Dartoon

The scripts for cut out the PSF candidate.
"""

import numpy as np
from regions import PixCoord, CirclePixelRegion 
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from flux_profile import pix_region
from matplotlib.colors import LogNorm
import copy 

def cut_image(image, center, radius):
    region = pix_region(center, radius=radius)
    cut = region.to_mask(mode='exact')
    cut_image = cut.cutout(image)
    return cut_image

def cut_center_bright(image, center, radius,kernel = 'gaussian', return_center=False,plot=True,center_float = False):
    """
    Auto cut the image, with with brightest center in the center.
        kernel: define the center.
    """
    temp_center = np.asarray(center)
#    print temp_center.astype(int)
    radius = radius
    img_test = cut_image(image=image, center=temp_center.astype(int), radius=radius)
    if kernel == 'gaussian':
        for i in range(3):
            frm_q = int(len(img_test)/2.5)
            from photutils import centroid_2dg
            test_center = frm_q + centroid_2dg(img_test[frm_q:-frm_q,frm_q:-frm_q])
            if i ==2 and plot==True :
                print test_center
                fig, ax = plt.subplots(1, 1)
                ax.imshow(img_test[frm_q:-frm_q,frm_q:-frm_q], origin='lower')
                marker = '+'
                ms, mew = 30, 2.
                plt.plot(test_center[0]-frm_q, test_center[1]-frm_q, color='b', marker=marker, ms=ms, mew=mew)
                plt.show()
                print 'test_center - radius', test_center, radius
            center_shift = np.array((test_center - radius))
            if i ==2 and plot==True :
                print 'center_shift',center_shift
            center = (center.astype(int) + np.round(center_shift))
            center_f = center.astype(int) + center_shift
            img_test = cut_image(image=image, center=center, radius=radius)
            if i ==2 and plot==True :
                plt_center = img_test[frm_q:-frm_q,frm_q:-frm_q].shape
                plt.plot(plt_center[0]/2, plt_center[1]/2, color='r', marker=marker, ms=ms, mew=mew)
                plt.imshow(img_test[frm_q:-frm_q,frm_q:-frm_q], origin='lower')
                plt.show()
        cut_c_b = img_test
    elif kernel == 'center_bright':
        frm_q = int(len(img_test)/2.5)
        test_center =  np.asarray(np.where(img_test == img_test[frm_q:-frm_q,frm_q:-frm_q].max()))[:,0]
        center_shift = np.array((test_center- radius))[::-1]
        center = (temp_center.astype(int) + np.round(center_shift))
        cut_c_b = cut_image(image=image, center=center, radius=radius)
        plt_center = img_test[frm_q:-frm_q,frm_q:-frm_q].shape
        if plot==True:
            marker = '+'
            ms, mew = 30, 2.
            plt.plot(plt_center[0]/2, plt_center[1]/2, color='c', marker='+', ms=30, mew=2)
            plt.imshow(cut_c_b[frm_q:-frm_q,frm_q:-frm_q], origin='lower')
            plt.show()
    if return_center==False:
        return cut_c_b
    elif return_center==True and center_float==False:
        return cut_c_b, center
    elif return_center==True and center_float==True:
        return cut_c_b, center_f
 
def save_loc_png(img, center_QSO, c_psf_list=None,extra_psfs=None,ID=None, label= None,
                 reg_ty=None, ifsave=True, label_shift_NO=(), shift_where=None):
    '''
    label shift_where: 1,2,3,4 --- up, right, down, left
    '''
    fig = plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    import copy, matplotlib
    my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
    my_cmap.set_bad('black')
    if reg_ty == None:
        vmax = 2.2
        vmin = 1.e-2
        QSO_box_size = 100
        PSF_box_size = 40
    cax=ax.imshow(img,origin='lower', cmap=my_cmap, norm=LogNorm(), vmin=vmin, vmax=vmax)
    QSO_reg = pix_region(center_QSO, radius= QSO_box_size)
    QSO_mask = QSO_reg.to_mask(mode='center')
    ax.text(center_QSO[0]-2*QSO_box_size, center_QSO[1]+1.5*QSO_box_size, 'QSO',color='white', fontsize=20)
    ax.add_patch(QSO_mask.bbox.as_artist(facecolor='none', edgecolor='white', linewidth=2))
    count=0
    count_shift = 0
    if label == None:
        name = 'PSF'
    elif label is not None:
        name = label
    if c_psf_list is not None:
        for i in range(len(c_psf_list)):
            PSF_reg = pix_region(c_psf_list[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_artist(facecolor='none', edgecolor='blue', linewidth=2))
            if count not in label_shift_NO:
                ax.text(c_psf_list[i][0]-2*PSF_box_size, c_psf_list[i][1]+2*PSF_box_size, '{1}{0}'.format(count,name),color='white', fontsize=15)
            else:
                if count in label_shift_NO:
                    shift = shift_label_index(shift_where[count_shift])
                    ax.text(c_psf_list[i][0]+shift[0]*PSF_box_size, c_psf_list[i][1]+shift[1]*PSF_box_size,
                            '{1}{0}'.format(count,name),color='white', fontsize=15)
                    count_shift += 1
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
#    plt.colorbar(cax)
    if extra_psfs is not None:
        for i in range(len(extra_psfs)):
            PSF_reg = pix_region(extra_psfs[i], radius= PSF_box_size)
            PSF_mask = PSF_reg.to_mask(mode='center')
            ax.add_patch(PSF_mask.bbox.as_artist(facecolor='none', edgecolor='yellow', linewidth=2))
            if count not in label_shift_NO:
                ax.text(extra_psfs[i][0]-2*PSF_box_size, extra_psfs[i][1]+2*PSF_box_size, '{1}{0}?'.format(count,name),color='white', fontsize=15)
            else:
                if count in label_shift_NO:
                    shift = shift_label_index(shift_where[count_shift])
                    print extra_psfs[i][0], shift[0]*PSF_box_size
                    ax.text(extra_psfs[i][0]+shift[0]*PSF_box_size, extra_psfs[i][1]+shift[1]*PSF_box_size,
                            '{1}{0}?'.format(count,name),color='white', fontsize=15)
                    count_shift += 1
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            count += 1
#    plt.colorbar(cax)
    if not ID == None:
        ax.text(len(img)*0.05, len(img)*0.8, ID,color='white', fontsize=30)
    if ifsave == True:
        fig.savefig('QSO_{0}_loc.pdf'.format(name))

def shift_label_index(s):
    if s==1:
        return np.array([-2, 2])
    elif s==2:
        return np.array([1, -1])
    elif s==3:
        return np.array([-2, -3])
    elif s==4:
        return np.array([-6, -1])    

def sub_bkg(img, plot=True):
    from astropy.stats import SigmaClip
    from photutils import Background2D, SExtractorBackground  
    sigma_clip = SigmaClip(sigma=3., iters=10)
    bkg_estimator = SExtractorBackground()
    from photutils import make_source_mask
    mask_0 = make_source_mask(img, snr=2, npixels=5, dilate_size=11)
    mask_1 = (np.isnan(img))
    mask = mask_0 + mask_1
    bkg = Background2D(img, (50, 50), filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       mask=mask)
    from matplotlib.colors import LogNorm
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    ax.imshow(img, norm=LogNorm(), origin='lower') 
    #bkg.plot_meshes(outlines=True, color='#1f77b4')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if plot:
        plt.show()  
    else:
        plt.close()
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    ax.imshow(mask, origin='lower') 
    #bkg.plot_meshes(outlines=True, color='#1f77b4')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if plot:
        plt.show()  
    else:
        plt.close()    
    back = bkg.background* ~mask_1
    fig=plt.figure(figsize=(15,15))
    ax=fig.add_subplot(1,1,1)
    ax.imshow(back, origin='lower', cmap='Greys_r')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    if plot:
        plt.show()  
    else:
        plt.close()
    return img-back, back    
