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

psf_i = 'PSF2'
subg = 'subg3'
fixgamma = '1.9'
#readfile = 'fit_PSFi_PSFrecons/result_PSF0_PSFrecons_gammafix1.9_subg2.pkl'
readfile = '2nd_fit_PSFi_PSFrecons/result_{0}_PSFrecons_gammafix{2}_{1}.pkl'.format(psf_i,subg, fixgamma)
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

psf_class = PSF(**kwargs_psf)
kwargs_model = kwargs_model
kwargs_result = kwargs_result
image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
multi_band_list = [image_band]

imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)
imageLinearFit = ImageLinearFit(data_class=data_class, psf_class=psf_class, lens_model_class = lens_model_class,
                                source_model_class = source_model_class,  point_source_class= point_source_class,
                                kwargs_numerics=kwargs_numerics) 
image_host_source0_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False)
host_flux0_total = image_host_source0_plane.sum()

modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, 
                      arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
#%%
import matplotlib
my_cmap = copy.copy(matplotlib.cm.get_cmap('gist_heat')) # copy the default cmap
my_cmap.set_bad('seashell') #linen
img = imageModel.image(kwargs_lens=lens_result, kwargs_source=source_result, kwargs_lens_light=lens_light_result, kwargs_ps=ps_result,
                 point_source_add=True,lens_light_add=True, source_add=False)
plt.imshow(img, origin='lower',norm=LogNorm(), cmap=my_cmap)
plt.show()
data = kwargs_data['image_data']
plt.imshow((data-img)*lens_mask, origin='lower',norm=LogNorm(), cmap=my_cmap)
plt.show()
print np.sum((data-img))

filt = 'f160w'
center_QSO = np.array([1318,1282])  #!!!!Need to be change
fitsFile = pyfits.open('../data/../data/WFI2033_sci.fits')
file_header = copy.deepcopy(fitsFile[0].header)
file_header['CRPIX1'] = file_header['CRPIX1']-center_QSO[0]+len(img)/2
file_header['CRPIX2'] = file_header['CRPIX2']-center_QSO[1]+len(img)/2
pyfits.PrimaryHDU((data-img)*lens_mask,header=file_header).writeto('lensed_arc_{3}_{0}_{1}_fixgamma{2}.fits'.format(psf_i,subg,fixgamma,filt),overwrite=True)
