#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 17:30:19 2019

@author: Dartoon

Modelling the RXJ1131

In this script, fix gamma as Ken's results. Use, psf mask.
Use PSF reconstruction.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle
import sys
sys.path.insert(0,'../../../py_tools')
from mask_objects import mask_obj, detect_obj
import copy
from psfs_average import psf_ave
from flux_profile import cr_mask

from lenstronomy.SimulationAPI.simulations import Simulation
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.Sampling.parameters import Param
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF


lens_image = pyfits.getdata('../WFI2033_cutout.fits')   #!!!need to change
lens_rms = pyfits.getdata('../WFI2033_stdd.fits')    #!!!need to change

lens_mask = 1-cr_mask(lens_image, filename='../lens_mask.reg')
mask2 = cr_mask(lens_image, filename='../mask_2_bigger.reg')
lens_mask = lens_mask * mask2 #* mask2 * mask3
plt.imshow(lens_mask, origin='low')
plt.show()

fit_frame_size = 101
ct = (len(lens_image[2])-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]

lens_image = lens_image[ct:-ct,ct:-ct]
lens_rms = lens_rms[ct:-ct,ct:-ct]
lens_mask = lens_mask[ct:-ct,ct:-ct]

sigma_bkg = 0.007   # inference from 0_cutout  #!!!need to change
exp_time = 26256.98 # exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
numPix = len(lens_image)  # cutout pixel size
deltaPix = 0.08  # pixel size in arcsec (area per pixel = deltaPix**2)  #!!!need to change
fwhm = 0.16  # full width half max of PSF

x_QSO = np.array([-0.59889414, -0.32321431, 0.94081951, 1.0572684])/deltaPix + numPix/2
y_QSO = -np.array([-0.63639854, 1.47860106, -0.70573724,0.00098984])/deltaPix + numPix/2

# =============================================================================
# Few things to set:
# =============================================================================
fix_gamma = 2.
fname = 'Shost_Dlenslight_nomask1_subg'
#psfno = 'ave'

psfname_list = ['../PSF{0}.fits'.format(i) for i in [0,1,2,6]]
psfs = []
for i in range(len(psfname_list)):
    psf_i = pyfits.getdata(psfname_list[i])
    _, psf_mask = mask_obj(img=psf_i,snr=2.6,exp_sz=1.7,pltshow=0)
    if 1 in psf_mask[0]:
        if len(psf_mask) > 1:
            psf_mask = np.sum(np.asarray(psf_mask),axis=0)
        elif len(psf_mask) == 1:
            psf_mask = psf_mask[0]
        psf_mask = (1 - (psf_mask != 0)*1.)
        psf_clear = copy.deepcopy(psf_i)
        for i in range(len(psf_i)):
            for j in range(len(psf_i)):
                if psf_mask[i, j] == 0:
                    psf_clear[i, j] = (psf_i[-i-1, -j-1]*psf_mask[-i-1, -j-1])
        print "plot PSF after clear:"
        plt.imshow(psf_clear, origin = 'low', norm=LogNorm())
        plt.show()
        psf_i = psf_clear
    psfs.append(psf_i)

psf_av, psf_std = psf_ave(psfs)
psf = psf_av/psf_av.sum()

SimAPI = Simulation()
kwargs_data = SimAPI.data_configure(numPix, deltaPix, exp_time, sigma_bkg)
kwargs_psf = SimAPI.psf_configure(psf_type='PIXEL', fwhm=fwhm, kernelsize=len(psf), deltaPix=deltaPix,kernel=psf)

if 'std' in psfno:
    print "including PSF error map"
    kwargs_psf['psf_error_map'] = psf_std

kwargs_data['image_data'] = lens_image
kwargs_data['noise_map'] = lens_rms
data_class = Data(kwargs_data)
# =============================================================================
# Setting up the fitting list, and guess parameter to set up the fitting parameter.
# =============================================================================
lens_model_list = ['SPEMD', 'SHEAR']
lens_model_class = LensModel(lens_model_list=lens_model_list)
lens_light_model_list = ['SERSIC_ELLIPSE'] * 3
lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
source_model_list = ['SERSIC_ELLIPSE']
source_model_class = LightModel(light_model_list=source_model_list)
point_source_list = ['LENSED_POSITION']
point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])

kwargs_numerics = {'subgrid_res': subg, 'psf_subgrid': False}
psf_class = PSF(kwargs_psf)
    
kwargs_numerics['mask'] = lens_mask    #Input the lens mask

imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)

kwargs_model = {'lens_model_list': lens_model_list,
                               'source_light_model_list': source_model_list,
                               'lens_light_model_list': lens_light_model_list,
                               'point_source_model_list': point_source_list,
                               'fixed_magnification_list': [False]
                             }

num_source_model = len(source_model_list)

kwargs_constraints = {'joint_lens_light_with_lens_light': [[0, 1, ['center_x', 'center_y']]],
                      'joint_source_with_point_source': [[0, 0]],
                              'num_point_source_list': [4],
                              }

kwargs_likelihood = {'check_bounds': True,
                             'force_no_add_image': False,
                             'source_marg': False,
                             'point_source_likelihood': False,
                             'position_uncertainty': 0.004,
                             'check_solver': True,
                             'solver_tolerance': 0.001
                             }

image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
multi_band_list = [image_band]

# initial guess of non-linear parameters, we chose different starting parameters than the truth #
kwargs_lens_init = [{'center_x': 0.0019689989121414583,
              'center_y': -0.011497666731504486,
              'e1': -0.03535314937396017,
              'e2': 0.18916890339488213,
              'gamma': 2.0,
              'theta_E': 1.1777960316138045},
             {'dec_0': 0,
              'e1': 0.0873129331489733,
              'e2': -0.014851825769766023,
              'ra_0': 0}]
kwargs_source_init = [{'R_sersic': 0.3513219841254192,
                  'center_x': 0.19058532038074247,
                  'center_y': 0.036543406152336075,
                  'e1': 0.03170624995752505,
                  'e2': -0.09535665666479456,
                  'n_sersic': 1.4280198296531827}]

kwargs_lens_light_init = [{'R_sersic': 0.2728900826606261,
                      'center_x': 0.04318525787908333,
                      'center_y': -0.005254803307251118,
                      'e1': -0.07135034593486825,
                      'e2': 0.15575122760293017,
                      'n_sersic': 4},
                     {'R_sersic': 1.793421489806939,
                      'center_x': 0.04318525787908333,
                      'center_y': -0.005254803307251118,
                      'e1': -0.11058462684580918,
                      'e2': 0.09759159851825082,
                      'n_sersic': 1},
                     {'R_sersic': 0.029065321208153066,
                      'center_x': 0.2804502700203146,
                      'center_y': 2.041699981873053,
                      'e1': -0.39005651113252693,
                      'e2': -0.1074667068129834,
                      'n_sersic': 0.6560768147213502}]

kwargs_ps_init = [{'dec_image': np.array([-0.5903102 , -0.32134738,  0.9434648 ,  1.0574707 ]),
  'ra_image': np.array([-0.63023983,  1.48375558, -0.70935364,  0.0019438 ])}]

# initial spread in parameter estimation #
kwargs_lens_sigma = [{'theta_E': 0.1, 'e1': 0.2, 'e2': 0.2, 'gamma': .1, 'center_x': 0.1, 'center_y': 0.1},
                     {'e1': 0.1, 'e2': 0.1}]
kwargs_source_sigma = [{'R_sersic': 0.2, 'n_sersic': .5, 'center_x': .1, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2}]
#kwargs_source_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
kwargs_lens_light_sigma = [{'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2}]
kwargs_lens_light_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
kwargs_lens_light_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
kwargs_ps_sigma = [{'ra_image': [0.02] * 4, 'dec_image': [0.02] * 4}]

# hard bound lower limit in parameter space #
kwargs_lower_lens = [{'theta_E': 0, 'e1': -0.5, 'e2': -0.5, 'gamma': 1.5, 'center_x': -10., 'center_y': -10},
#                     {'theta_E': 0, 'center_x': -10., 'center_y': -10},
                     {'e1': -0.2, 'e2': -0.2}]
kwargs_lower_source = [{'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10}]
#kwargs_lower_source.append({'R_sersic': 0.001,'n_sersic': 0.5,  'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
kwargs_lower_lens_light = [{'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10}]
kwargs_lower_lens_light.append({'R_sersic': 0.001,'n_sersic': 0.5,  'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
kwargs_lower_lens_light.append({'R_sersic': 0.001,'n_sersic': 0.5,  'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
kwargs_lower_ps = [{'ra_image': -10 * np.ones(4), 'dec_image': -10 * np.ones(4)}]

# hard bound upper limit in parameter space #
kwargs_upper_lens = [{'theta_E': 10, 'e1': 0.5, 'e2': 0.5, 'gamma': 2.5, 'center_x': 10., 'center_y': 10},
#                     {'theta_E': 10, 'center_x': 10., 'center_y': 10},
                     {'e1': 0.2, 'e2': 0.2}]
kwargs_upper_source = [{'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10}]
#kwargs_upper_source.append({'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10})
kwargs_upper_lens_light = [{'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10}]
kwargs_upper_lens_light.append({'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10})
kwargs_upper_lens_light.append({'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10})
kwargs_upper_ps = [{'ra_image': 10 * np.ones(4), 'dec_image': 10 * np.ones(4)}]

#Fix something:
fixed_source = [{}]
fixed_lens = [{'gamma' : fix_gamma}, {'ra_0': 0, 'dec_0': 0}]
fixed_lens_light = [{'n_sersic': 4},{'n_sersic': 1},{}]

lens_params = [kwargs_lens_init, kwargs_lens_sigma, fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma,
                     fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]   # Set as 1 as Sherry do.
ps_params = [kwargs_ps_init, kwargs_ps_sigma, [{}], kwargs_lower_ps, kwargs_upper_ps]

kwargs_params = {'lens_model': lens_params,
                'source_model': source_params,
                'lens_light_model': lens_light_params,
                'point_source_model': ps_params}

from lenstronomy.Workflow.fitting_sequence import FittingSequence
fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)

fitting_kwargs_list_0 = [
        {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 300,
         'n_iterations': 300}]  #first do the PSO
fitting_seq.fit_sequence(fitting_kwargs_list_0)

kwargs_psf_iter = {'stacking_method': 'median', 
                   'keep_error_map': True, 
                   'psf_symmetry': 1, 
                   'block_center_neighbour': 0.05}
fitting_kwargs_list_1 = [
        {'fitting_routine': 'psf_iteration', 'psf_iter_num': 100, 'psf_iter_factor': 0.2, 'kwargs_psf_iter': kwargs_psf_iter},
        {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100},
        {'fitting_routine': 'psf_iteration', 'psf_iter_num': 100, 'psf_iter_factor': 0.2, 'kwargs_psf_iter': kwargs_psf_iter},
        {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100},
        {'fitting_routine': 'MCMC', 'n_burn': 20, 'n_run': 20, 'walkerRatio': 10, 'mpi': False,
         'sigma_scale': .1},
    ]
lens_result, source_result, lens_light_result, ps_result, cosmo_result,\
chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list_1)

kwargs_data, kwargs_psf_updated, kwargs_numerics = fitting_seq.multi_band_list[0]
import lenstronomy.Plots.output_plots as out_plot
f, axes = out_plot.psf_iteration_compare(kwargs_psf_updated); f.show()
plt.savefig('fig_PSF{0}_PSFrecons_gammafix_comp.pdf'.format(psfno))
plt.show()

from flux_profile import profiles_compare
f_ = profiles_compare([kwargs_psf_updated['kernel_point_source_init'], kwargs_psf_updated['kernel_point_source']], scal_list=np.ones(2),
                       prf_name_list=['PSF', 'new PSF'],
                       gridspace = 'log', if_annuli=True)

##If to save the fitting reuslt as the pickle:
#filename='fit_result_PSF{0}_PSFrecons_subg2'.format(psfno)
#fit_result = [lens_result, source_result, lens_light_result, ps_result, kwargs_psf_updated]
#pickle.dump(fit_result, open(filename, 'wb'))
#If to load the fitting reuslt by the pickle:
#lens_result, source_result, lens_light_result, ps_result, cosmo_result,\
#chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc, kwargs_psf_updated = pickle.load(open(filename,'rb'))
#filename = 'fit_result_PSF{0}_{2}{1}'.format(psfno,subg, fname)
#lens_result, source_result, lens_light_result, ps_result, kwargs_psf_updated, host_flux_total = pickle.load(open(filename,'rb'))

##Update the PSF:
psf_class = PSF(kwargs_psf_updated)
#Update the imageModel with the new PSF.
imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)

from lenstronomy.Plots.output_plots import LensModelPlot

lensPlot = LensModelPlot(kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_result, source_result,
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
plt.savefig('fig_PSF{0}_{2}{1}_img0.pdf'.format(psfno,subg, fname))
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
plt.savefig('fig_PSF{0}_{2}{1}_img1.pdf'.format(psfno,subg,fname))
plt.show()
print(lens_result, source_result, lens_light_result, ps_result)
#print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
#print("parameters in order: ", param_mcmc)
#print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])

image_host_source0_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False,k=0)
host_flux_total = image_host_source0_plane.sum()
print "PSF:", psfno
print "host flux source plane:", host_flux_total

#If to save the fitting reuslt as the pickle:
#filename='fit_result_PSF{0}_{2}{1}'.format(psfno,subg, fname)
#fit_result = [lens_result, source_result, lens_light_result, ps_result, kwargs_psf_updated, host_flux_total]
#pickle.dump(fit_result, open(filename, 'wb'))

import corner
n, num_param = np.shape(samples_mcmc)
#plot = corner.corner(samples_mcmc[:,:8], labels=param_mcmc[:8], show_titles=True)
#plot = corner.corner(samples_mcmc[:,8:], labels=param_mcmc[8:], show_titles=True)

mcmc_new_list = []
param = Param(kwargs_model, kwargs_constraints, kwargs_fixed_lens = fixed_lens, kwargs_fixed_source=fixed_source, kwargs_fixed_lens_light=fixed_lens_light)
if param.num_param()[0] != len(samples_mcmc[0]):
    raise ValueError("The param shape is not as the samples_mcmc_ones.")
    
#labels_new = ["AGN flux in image plane", "Host flux image plance", "AGN flux in souce plane", "Host flux souce plane",  "Host Sersic", "Sersic Reff"]
labels_new = ["AGN flux in image plane", "Host flux image plance", "AGN flux in souce plane", "Host flux souce plane",  "Host Sersic", "Sersic Reff"] 
lensModel_simple = LensModel(lens_model_list)
for i in range(len(samples_mcmc)):
    kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo_out = param.args2kwargs(samples_mcmc[i])
    image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_lens=kwargs_lens_out, kwargs_source=kwargs_light_source_out, kwargs_lens_light=kwargs_light_lens_out, kwargs_ps=kwargs_ps_out)
    image_ps = imageModel.point_source(kwargs_ps_out, kwargs_lens=kwargs_lens_out, unconvolved=False)
    flux_quasar = np.sum(image_ps)
    image_host_source_plane = imageModel.source_surface_brightness(kwargs_light_source_out, kwargs_lens_out, de_lensed=True,unconvolved=False)
    host_flux_s = np.sum(image_host_source_plane)
    image_host_image_plane = imageModel.source_surface_brightness(kwargs_light_source_out, kwargs_lens_out, de_lensed=False,unconvolved=False)
    host_flux_i = np.sum(image_host_image_plane)
    x_image, y_image, AGN_fluxs = kwargs_ps_out[0]['ra_image'], kwargs_ps_out[0]['dec_image'], kwargs_ps_out[0]['point_amp']
    mag_macro = lensModel_simple.magnification(x_image, y_image, kwargs_lens_out)
    AGN_fluxi_s = AGN_fluxs/mag_macro
    AGN_fluxs_s = np.average(np.abs(AGN_fluxi_s))
#    print i, host_flux_s, host_flux_i
    mcmc_new_list.append([flux_quasar,host_flux_i, AGN_fluxs_s, host_flux_s, kwargs_light_source_out[0]['n_sersic'],kwargs_light_source_out[0]['R_sersic']])
    if i/1000 > (i-1)/1000 :
        print "finished translate:", i

fig = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
fig.savefig('fig_PSF{0}_{2}{1}_corner.pdf'.format(psfno,subg,fname))
plt.show()

picklename='result_PSF{0}_{2}{1}'.format(psfno,subg,fname)
fit_result = [lens_result, source_result, lens_light_result, ps_result, cosmo_result,chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc]
trans_result = [mcmc_new_list, labels_new]
pickle.dump([fit_result, trans_result, kwargs_psf_updated], open(picklename, 'wb'))

import os
os.system('say "your program of PSF{0} has finished"'.format(psfno))