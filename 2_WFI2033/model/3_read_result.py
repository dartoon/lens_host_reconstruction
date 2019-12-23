#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 10:44:07 2019

@author: Dartoon

Test load text from PDF
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
from matplotlib.colors import LogNorm
import copy
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
from flux_profile import cr_mask


#readfile = 'fit_PSFi_PSFrecons/result_PSF0_PSFrecons_gammafix1.9_subg2.pkl'
readfile = 'fit_PSFi_QSOmask/result_PSF3_QSOmask_gammafix1.9_subg2.pkl'
result = pickle.load(open(readfile,'rb'))
fit_result, trans_result, kwargs_material, model_lists  = result

kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask = kwargs_material
lens_model_list, source_model_list, lens_light_model_list, point_source_list = model_lists

#%%
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import ImageData as Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.ImSim.image_linear_solve import ImageLinearFit
from lenstronomy.Plots import chain_plot
import lenstronomy.Plots.chain_plot as out_plot
#%%
kwargs_result, chain_list = fit_result

sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]    
lens_result = kwargs_result['kwargs_lens']
lens_light_result = kwargs_result['kwargs_lens_light']
source_result = kwargs_result['kwargs_source']
ps_result = kwargs_result['kwargs_ps']

mcmc_new_list, labels_new = trans_result

lens_model_class = LensModel(lens_model_list=lens_model_list)
lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
source_model_class = LightModel(light_model_list=source_model_list)
point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])

data_class = Data(**kwargs_data)
kwargs_psf = kwargs_psf_updated
#If want to use the original PSF to plot:
#kwargs_psf['kernel_point_source'] = kwargs_psf['kernel_point_source_init']
#del kwargs_psf["psf_error_map"]
psf_class = PSF(**kwargs_psf)
kwargs_model = kwargs_model
kwargs_result = kwargs_result
image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
multi_band_list = [image_band]

#%%
modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, 
                      arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
#param_class = fitting_seq.param_class
#print(param_class.num_param())
#print(chain_list)
for i in range(len(chain_list)):
    chain_plot.plot_chain_list(chain_list, i)

f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

modelPlot.data_plot(ax=axes[0,0])
modelPlot.model_plot(ax=axes[0,1])
modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=100)
modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
modelPlot.magnification_plot(ax=axes[1, 2])
f.tight_layout()
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
plt.show()

f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True)
modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True)
modelPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True)
modelPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
modelPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True, unconvolved=True)
modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
f.tight_layout()
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
plt.show()

imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)
ImageLinearFit
imageLinearFit = ImageLinearFit(data_class=data_class, psf_class=psf_class, lens_model_class = lens_model_class,
                                source_model_class = source_model_class,  point_source_class= point_source_class,
                                kwargs_numerics=kwargs_numerics) 
image_host_source0_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False)
host_flux0_total = image_host_source0_plane.sum()
#print "host flux, source plane:", host_flux0_total

if 'kernel_point_source_init' in kwargs_psf.keys():
    f, axes = out_plot.psf_iteration_compare(kwargs_psf_updated); f.show()
    plt.show()

import corner
fig = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
plt.show()
#
#
#n, num_param = np.shape(samples_mcmc)
#plot = corner.corner(samples_mcmc[:,:8], labels=param_mcmc[:8], show_titles=True)
#plot.show()
#plot = corner.corner(samples_mcmc[:,8:], labels=param_mcmc[8:], show_titles=True)
#plot.show()

#%%
import lenstronomy.Util.class_creator as class_creator
imageModel = class_creator.create_im_sim(multi_band_list = multi_band_list, multi_band_type='multi-linear', kwargs_model = kwargs_model,
                                           bands_compute = [True] * len(multi_band_list),
                                           likelihood_mask_list=[lens_mask],
                                           band_index=0)
logL = imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
n_data = imageModel.num_data_evaluate    
reduced_x2 = - logL * 2 / n_data