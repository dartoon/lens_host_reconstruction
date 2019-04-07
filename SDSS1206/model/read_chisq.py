#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 20 16:18:24 2019

@author: Dartoon
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
#from matplotlib.colors import LogNorm
#import copy
import sys
#import glob
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools')
from flux_profile import cr_mask
#from lenstronomy.ImSim.image_model import ImageModel
#from lenstronomy.Data.imaging_data import Data
#from lenstronomy.Data.psf import PSF
#from lenstronomy.PointSource.point_source import PointSource
#from lenstronomy.LensModel.lens_model import LensModel
#from lenstronomy.LightModel.light_model import LightModel

def reg_lens_model_list(lens_result):
    '''
    Recongnize the lens model as, e.g., 
    lens_model_list including type: ['SPEMD', 'SIS', 'SHEAR'].
    Given the lens_result
    '''
    lens_model_list = []
    for i in range(len(lens_result)):
        if 'theta_E' not in lens_result[i]:
            model = 'SHEAR'
        elif 'gamma' not in lens_result[i]:
            model = 'SIS'
        else:
            model = 'SPEMD'    
        lens_model_list.append(model)
    return lens_model_list
        
def reg_light_model_list(light_result):
    '''
    Recongnize the lens model as, e.g., 
    light_model_list = ['SERSIC_ELLIPSE','SERSIC'],
    Given the light_result
    '''
    light_model_list = []
    for i in range(len(light_result)):
        if 'e1' not in light_result[i]:
            model = 'SERSIC'
        else:
            model = 'SERSIC_ELLIPSE'    
        light_model_list.append(model)
    return light_model_list      

def return_chisq(filename, plt_show=0):
    import glob         
    lens_image_fname = glob.glob('./*_cutout.fits')
    lens_rms_fname = glob.glob('./*_stdd.fits')
    #subg = 2 #int(filename.split('subg')[1])
#    filename = 'fit_PSFi_PSFrecons/result_PSF0_2nd_UPversion_DhostDlens_fixG1c_SISG2_subg2'
    #filename = 'fit_PSFi_QSOmask/result_PSF0_UPversion_DhostDlens_fixG1c_SISG2_subg{0}'.format(subg)
    result = pickle.load(open(filename,'rb'))
    subg = int((filename.split('/')[1]).split('subg')[1][0])
    psfno = int((filename.split('/')[1]).split('PSF')[1][0])
    if len(result) == 2:
        fit_result, trans_result = result
    elif len(result) ==3:
        fit_result, trans_result, kwargs_psf_updated = result
    
    [lens_result, source_result, lens_light_result, ps_result,\
    cosmo_result,chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc] = fit_result
    [mcmc_new_list, labels_new] =  trans_result
    
    lens_image = pyfits.getdata(lens_image_fname[0])
    lens_rms = pyfits.getdata(lens_rms_fname[0])
    
    lens_mask = 1 - cr_mask(lens_image, 'lens_mask.reg')  #!!!
    obj_maks = cr_mask(lens_image, 'obj0_mask.reg')
    lens_mask = lens_mask * obj_maks
    fit_frame_size = 101  #!!!
    ct = (len(lens_image[2])-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]
    lens_image = lens_image[ct:-ct,ct:-ct]
    lens_rms = lens_rms[ct:-ct,ct:-ct]
    lens_mask = lens_mask[ct:-ct,ct:-ct]
    numPix = fit_frame_size
    deltaPix = 0.08  # pixel size in arcsec (area per pixel = deltaPix**2)  #!!!
    
    if len(result) == 3:
        psf = kwargs_psf_updated['kernel_point_source']
    else:
    #    psf = kwargs_psf_updated['kernel_point_source_init']
        folder_name = filename.split('/')[0].replace('QSOmask','PSFrecons')
        import glob
        filename_PSFave = glob.glob(folder_name+"/result*PSF{0}*".format(psfno))[0]
        psf = pickle.load(open(filename_PSFave,'rb'))[2]['kernel_point_source_init']
    
    if 'mask' in filename:
        from mask_objects import find_loc_max
        QSOx, QSOy =find_loc_max(lens_image)
        del QSOx[1]
        del QSOy[1]
        xy_index = np.indices((numPix,numPix))
        for i in range(len(QSOx)):
            if i == 0:
                areas = (np.sqrt((QSOx[i]-xy_index[1])**2+(QSOy[i]-xy_index[0])**2) <4 )  #!!! pixel numbers
            else:
                areas += (np.sqrt((QSOx[i]-xy_index[1])**2+(QSOy[i]-xy_index[0])**2) <4 )  #!!! pixel numbers
        lens_rms = lens_rms * (areas == 0) + 10**6 * (areas != 0)
    
    import lenstronomy.Util.simulation_util as sim_util
    kwargs_data = sim_util.data_configure_simple(numPix, deltaPix, inverse=True)
    kwargs_psf = sim_util.psf_configure_simple(psf_type='PIXEL', kernelsize=len(psf), deltaPix=deltaPix,kernel=psf)
    kwargs_data['image_data'] = lens_image
    kwargs_data['noise_map'] = lens_rms
    # =============================================================================
    # Setting up the fitting list, and guess parameter to set up the fitting parameter.
    # =============================================================================
    
    lens_model_list = reg_lens_model_list(lens_result) #['SPEMD', 'SIS', 'SHEAR','SIS']
    lens_light_model_list = reg_light_model_list(lens_light_result)  #['SERSIC_ELLIPSE','SERSIC_ELLIPSE','SERSIC']
    source_model_list = reg_light_model_list(source_result) #['SERSIC_ELLIPSE', 'SERSIC_ELLIPSE']
    point_source_list = ['LENSED_POSITION']
    from lenstronomy.Plots.output_plots import LensModelPlot
    kwargs_numerics = {'subgrid_res': subg, 'psf_subgrid': False}
    kwargs_numerics['mask'] = lens_mask    #Input the lens mask
    
#    data_class = Data(kwargs_data)
#    psf_class = PSF(kwargs_psf)
#    lens_model_class = LensModel(lens_model_list=lens_model_list)
#    lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
#    source_model_class = LightModel(light_model_list=source_model_list)
#    point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])
#    imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
#                                    lens_light_model_class,
#                                    point_source_class, kwargs_numerics=kwargs_numerics)
    kwargs_model = {'lens_model_list': lens_model_list,
                                   'source_light_model_list': source_model_list,
                                   'lens_light_model_list': lens_light_model_list,
                                   'point_source_model_list': point_source_list,
                                   'fixed_magnification_list': [False],
                                 }
    lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                                 lens_light_result, ps_result)
    if plt_show == 1:    
        f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)
        lensPlot.data_plot(ax=axes[0,0])
        lensPlot.model_plot(ax=axes[0,1])
        lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        lensPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=100)
        lensPlot.convergence_plot(ax=axes[1, 1], v_max=1)
        lensPlot.magnification_plot(ax=axes[1, 2])
        f.tight_layout()
        f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        plt.show()
        f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)
        lensPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True)
        lensPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True)
        lensPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True)
        lensPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
        lensPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True, unconvolved=True)
        lensPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        f.tight_layout()
        f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        plt.show()
    #    image_host_source0_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False,k=0)
    #    host_flux0_total = image_host_source0_plane.sum()
    #    image_host_source1_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False,k=1)
    #    host_flux1_total = image_host_source1_plane.sum()
    return lensPlot._reduced_x2
print return_chisq('fit_PSFi_PSFrecons/result_PSF0_2nd_UPversion_DhostDlens_fixG1c_SISG2_subg2')


#print "bulge, disk flux, source plane:", host_flux0_total, host_flux1_total

#import corner
#fig = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
##fig.savefig('fig_PSF{0}_{2}{1}_corner.pdf'.format(psfno,subg,fname))
#plt.show()