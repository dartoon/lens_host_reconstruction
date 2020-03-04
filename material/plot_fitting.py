#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 11:50:31 2020

@author: Dartoon
"""


import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
from matplotlib.colors import LogNorm
import copy
import sys
sys.path.insert(0,'/Users/Dartoon/Astro/Projects/my_code/py_tools/')
from flux_profile import cr_mask
import glob

#%% Find the one with lowest Chisq
IDs = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']


#for i in range(len(IDs)):
#    print(i, ':', IDs[i])
#inp = int(input("Which source:\n"))
inp = 3
ID = IDs[inp][2:]
if ID != 'SDSS1206':
    folder_0 = '_2nd_run_summary.pkl'
    folder_1 = '2nd_fit'
elif ID == 'SDSS1206':
    folder_0 = '_singel_host_run_summary.pkl'
    folder_1 = 'singel_fit'
filename = '../share_tools/'+ ID+ folder_0
result = pickle.load(open(filename,'rb') ,encoding="latin1") 
fit_values, chisq, labels, _ = result
chisq = [float(chisq[i]) for i in range(len(chisq))]
idx = np.where(chisq == np.min(chisq))[0][0]
label = labels[idx]
label = label.split(', ')

filename = glob.glob('../'+IDs[inp]+'/model/{4}_PSFi_{0}/result_{1}_*_gammafix{2}_subg{3}.pkl'.format(label[-1], label[0], label[1][-3:], label[2][-1],folder_1) )
#filename = glob.glob('../5_SDSS0246/model/2nd_fit_PSFi_PSFrecons/result_PSF0_*_gammafix2.1_subg2.pkl')
#filename = glob.glob('../5_SDSS0246/model/2nd_fit_PSFi_QSOmask/result_PSF0_QSOmask_gammafix2.1_subg2.pkl')
readfile = filename[0]

#%%
##
result = pickle.load(open(readfile,'rb') ,encoding="latin1")
fit_result, trans_result, kwargs_material, model_lists  = result
kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask = kwargs_material
lens_model_list, source_model_list, lens_light_model_list, point_source_list = model_lists
if label[-1] == 'PSFrecons':
    print("Change PSFs and ")
    filename = glob.glob('../'+IDs[inp]+'/model/{4}_PSFi_{0}/result_{1}_*_gammafix{2}_subg{3}.pkl'.format("QSOmask", label[0], label[1][-3:], label[2][-1],folder_1) )
    readfile = filename[0]
    result = pickle.load(open(readfile,'rb') ,encoding="latin1")
    _, trans_result, kwargs_material, model_lists  = result
    kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask = kwargs_material

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
#!!! For the Lenstronomy version 1.3.0
kwargs_result, chain_list = fit_result
if 'e1' in kwargs_result['kwargs_lens'][1].keys():
    g1 = kwargs_result['kwargs_lens'][1]['e1']
    g2 = kwargs_result['kwargs_lens'][1]['e2']
    del kwargs_result['kwargs_lens'][1]
    kwargs_result['kwargs_lens'].append({'gamma1': g1, 'gamma2': g2})

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
psf_class = PSF(**kwargs_psf)
kwargs_model = kwargs_model
kwargs_result = kwargs_result
image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
multi_band_list = [image_band]

#%%
modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, 
                      arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])

f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)
modelPlot.data_plot(ax=axes[0,0])
modelPlot.model_plot(ax=axes[0,1])
modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
modelPlot.subtract_from_data_plot(ax=axes[1, 0], #v_max=1, 
                                  point_source_add=True,lens_light_add=True, source_add=False, 
                                  text='Residual of lensed arcs')
modelPlot.decomposition_plot(ax=axes[1,1], text='Modelled source light', source_add=True)
#modelPlot.decomposition_plot(ax=axes[1,2], text='Source light', source_add=True, unconvolved=True)
#modelPlot.magnification_plot(ax=axes[1, 2])
delta_list = [0.02, 0.1, 0.02, 0.05, 0.06, 0.02, 0.02, 0.02, 0.02]
modelPlot.source_plot(ax=axes[1, 2], deltaPix_source=delta_list[inp], numPix=100, scale_size=0.5,
                      text='Reconstrucerd host, source plane',with_caustics=True)
#f.tight_layout()
f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
#plt.savefig('{0}_inference.png'.format(ID))
plt.show()

#f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)
#modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True)
#modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True)
#modelPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True)
#modelPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
#modelPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True, unconvolved=True)
#modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
#f.tight_layout()
#f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
#plt.show()

#imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
#                                lens_light_model_class,
#                                point_source_class, kwargs_numerics=kwargs_numerics)
#ImageLinearFit
#imageLinearFit = ImageLinearFit(data_class=data_class, psf_class=psf_class, lens_model_class = lens_model_class,
#                                source_model_class = source_model_class,  point_source_class= point_source_class,
#                                kwargs_numerics=kwargs_numerics) 
#image_host_source0_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False)
#host_flux0_total = image_host_source0_plane.sum()
##print "host flux, source plane:", host_flux0_total
#
#if 'kernel_point_source_init' in kwargs_psf.keys():
#    f, axes = out_plot.psf_iteration_compare(kwargs_psf_updated); f.show()
#    plt.show()
#
#
##%%
#import lenstronomy.Util.class_creator as class_creator
#imageModel = class_creator.create_im_sim(multi_band_list = multi_band_list, multi_band_type='multi-linear', kwargs_model = kwargs_model,
#                                           bands_compute = [True] * len(multi_band_list),
#                                           likelihood_mask_list=[lens_mask],
#                                           band_index=0)
#logL = imageModel.likelihood_data_given_model(source_marg=False, linear_prior=None, **kwargs_result)
#n_data = imageModel.num_data_evaluate    
#reduced_x2 = - logL * 2 / n_data