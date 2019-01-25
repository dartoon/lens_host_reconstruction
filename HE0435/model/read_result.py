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
from matplotlib.colors import LogNorm
import copy
import sys
sys.path.insert(0,'../../py_tools')
from flux_profile import cr_mask

psf_i = raw_input("Which psf to load? 0, 1, 2, 3...:\n")

lens_result, source_result, lens_light_result, ps_result, cosmo_result,\
chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = pickle.load(open('fit_result_PSF{0}'.format(psf_i),'rb'))

print "lens_result:", lens_result,'\n'
print "source_light_result:", source_result,'\n'
print "lens_light_result:", lens_light_result,'\n'
print "point_source_result:", ps_result,'\n'

import sys
sys.path.insert(0,'../../py_tools')
from mask_objects import detect_obj, mask_obj
from lenstronomy.SimulationAPI.simulations import Simulation
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel

ct = 30
lens_image = pyfits.getdata('HE0435_cutout.fits')
lens_rms = pyfits.getdata('HE0435_noisemap.fits')
#lens_mask_0 = cr_mask(lens_image, 'circle_mask.reg')
lens_image = lens_image[ct:-ct,ct:-ct]
lens_rms = lens_rms[ct:-ct,ct:-ct]

#lens_mask_0 = lens_mask_0[ct:-ct,ct:-ct]

_,lens_obj_mask  = mask_obj(img=lens_image, exp_sz=1)
lens_mask = (1-lens_obj_mask[0])*(1-lens_obj_mask[1]) #* (1 - lens_mask_0)

#print "plot lens image mask:"

psfno = psf_i

lens_rms = pyfits.getdata('HE0435_noisemap.fits')[ct:-ct,ct:-ct]
psf = pyfits.getdata('PSF{0}.fits'.format(psfno))
objs, psf_index = detect_obj(psf,snr=2.0,pltshow=0, exp_sz=1.7)
_, psf_mask = mask_obj(img=psf,snr=2.0,exp_sz=1.7, pltshow = 0)
if len(psf_mask) > 0:
    if len(psf_mask) > 1:
        psf_mask = np.sum(np.asarray(psf_mask),axis=0)
    elif len(psf_mask) == 1:
        psf_mask = psf_mask[0]
    psf_mask = (1 - (psf_mask != 0)*1.)
    plt.imshow(psf_mask, origin = 'low')
    plt.close()
    psf_clear = copy.deepcopy(psf)
    for i in range(len(psf)):
        for j in range(len(psf)):
            if psf_mask[i, j] == 0:
                psf_clear[i, j] = (psf[-i-1, -j-1]*psf_mask[-i-1, -j-1])
    print "plot PSF after clear:"
    plt.imshow(psf_clear, origin = 'low', norm=LogNorm())
    plt.close()
    psf = psf_clear
psf /= psf.sum()
#plt.imshow(lens_image, norm = LogNorm(), origin='lower')
#plt.close()

sigma_bkg = 0.01245   # inference from 0_cutout
exp_time = 10*4  # exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
numPix = len(lens_image)  # cutout pixel size
deltaPix = 0.08  # pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 0.16  # full width half max of PSF

SimAPI = Simulation()
kwargs_data = SimAPI.data_configure(numPix, deltaPix, exp_time, sigma_bkg)
data_class = Data(kwargs_data)
kwargs_psf = SimAPI.psf_configure(psf_type='PIXEL', fwhm=fwhm, kernelsize=len(psf), deltaPix=deltaPix,kernel=psf)
psf_class = PSF(kwargs_psf)

kwargs_data['image_data'] = lens_image
kwargs_data['noise_map'] = lens_rms
# =============================================================================
# Setting up the fitting list, and guess parameter to set up the fitting parameter.
# =============================================================================
lens_model_list = ['SPEMD', 'SHEAR']
lens_model_class = LensModel(lens_model_list=lens_model_list)
lens_light_model_list = ['SERSIC_ELLIPSE']
lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
source_model_list = ['SERSIC_ELLIPSE']
source_model_class = LightModel(light_model_list=source_model_list)
point_source_list = ['LENSED_POSITION']
point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])

from lenstronomy.Plots.output_plots import LensModelPlot
kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
kwargs_numerics['mask'] = lens_mask    #Input the lens mask

imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)

kwargs_model = {'lens_model_list': lens_model_list,
                               'source_light_model_list': source_model_list,
                               'lens_light_model_list': lens_light_model_list,
                               'point_source_model_list': point_source_list,
                               'fixed_magnification_list': [False],
                             }
lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                             lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")

f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

lensPlot.data_plot(ax=axes[0,0])
lensPlot.model_plot(ax=axes[0,1])
lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
lensPlot.source_plot(ax=axes[1, 0],convolution=False, deltaPix_source=0.01, numPix=100)
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