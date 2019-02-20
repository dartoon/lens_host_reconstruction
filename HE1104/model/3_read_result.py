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

psf_i = 0 #raw_input("Which psf to load? 0, 1, 2, 3...:\n")

filename = 'fit_PSFi_QSOmask/result_PSF{0}_Shost_Dlenslight_nomask1_subg2'.format(psf_i)
#filename = 'fit_PSFi_PSFrecons/result_PSF{0}_Shost_Dlenslight_nomask1_subg2'.format(psf_i)
result = pickle.load(open(filename,'rb'))
subg = int(filename.split('subg')[1])


if len(result) == 2:
    fit_result, trans_result = result
elif len(result) ==3:
    fit_result, trans_result, kwargs_psf_updated = result

[lens_result, source_result, lens_light_result, ps_result,\
cosmo_result,chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc] = fit_result
 
[mcmc_new_list, labels_new] =  trans_result

print "lens_result:", lens_result,'\n'
print "source_light_result:", source_result,'\n'
print "lens_light_result:", lens_light_result,'\n'
print "point_source_result:", ps_result,'\n'
#=============================================================================
import sys
sys.path.insert(0,'../../py_tools')
from lenstronomy.SimulationAPI.simulations import Simulation
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel

lens_image = pyfits.getdata('WFI2033_cutout.fits')
lens_rms = pyfits.getdata('WFI2033_stdd.fits')
lens_mask = 1-cr_mask(lens_image, filename='lens_mask.reg')
mask2 = cr_mask(lens_image, filename='mask_2_bigger.reg')
lens_mask = lens_mask * mask2 #* mask2 * mask3

fit_frame_size = 101
ct = (len(lens_image[2])-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]

lens_image = lens_image[ct:-ct,ct:-ct]
lens_rms = lens_rms[ct:-ct,ct:-ct]
lens_mask = lens_mask[ct:-ct,ct:-ct]

#print "plot lens image mask:"
psfno = psf_i

if len(result) == 3:
    psf = kwargs_psf_updated['kernel_point_source']
else:
    folder_name = filename.split('/')[0].replace('QSOmask','PSFrecons')
    import glob
    filename_PSFave = glob.glob(folder_name+"/result*PSF{0}*".format(psf_i))[0]
    psf = pickle.load(open(filename_PSFave,'rb'))[2]['kernel_point_source_init']

numPix = len(lens_image)  # cutout pixel size
deltaPix = 0.08  # pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 0.16  # full width half max of PSF

if 'QSOmask' in filename:
    x_QSO = np.array([-0.59889414, -0.32321431, 0.94081951, 1.0572684])/deltaPix + numPix/2
    y_QSO = -np.array([-0.63639854, 1.47860106, -0.70573724,0.00098984])/deltaPix + numPix/2
    xy_index = np.indices((numPix,numPix))
    for i in range(len(x_QSO)):
        if i == 0:
            areas = (np.sqrt((x_QSO[i]-xy_index[0])**2+(y_QSO[i]-xy_index[1])**2) <3 )  # 3 piexls
        else:
            areas += (np.sqrt((x_QSO[i]-xy_index[0])**2+(y_QSO[i]-xy_index[1])**2) <3 )  # 3 piexls
    lens_rms = lens_rms * (areas == 0) + 10**6 * (areas != 0)

SimAPI = Simulation()
kwargs_data = SimAPI.data_configure(numPix, deltaPix)
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
lens_light_model_list = ['SERSIC_ELLIPSE'] * len(lens_light_result)
lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
source_model_list = ['SERSIC_ELLIPSE'] * len(source_result)
source_model_class = LightModel(light_model_list=source_model_list)
point_source_list = ['LENSED_POSITION']
point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])

from lenstronomy.Plots.output_plots import LensModelPlot
kwargs_numerics = {'subgrid_res': subg, 'psf_subgrid': False}
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
