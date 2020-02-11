#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 17:30:19 2019

@author: Dartoon

Modelling the HE0435

Trying to recover the H0LiCOW_IV's host flux ~47.

In this script, fix gamma as Ken's results. Use, psf mask.
Use PSF reconstruction.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle
import sys
sys.path.insert(0,'../../../../my_code/py_tools/')
#sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
from mask_objects import mask_obj
import copy, corner 
from psfs_average import psf_ave
from flux_profile import cr_mask
import glob

sys.path.insert(0,'../../../share_tools/')

import lenstronomy.Util.simulation_util as sim_util
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Data.psf import PSF
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Plots import chain_plot
from lenstronomy.Sampling.parameters import Param
from lenstronomy.ImSim.image_linear_solve import ImageLinearFit

import time

import os
path = os.getcwd()
ID = path.split('/')[-3][2:]

filename = '../../otherband_data/f555w_drc_sci.fits'
fitsFile = pyfits.open(filename)
from my_astro_func import read_pixel_scale
pixscl = read_pixel_scale(filename)

t1 = time.time()

lens_image = pyfits.getdata('../WFI2033_cutout.fits')   #!!!need to change
lens_rms = pyfits.getdata('../WFI2033_stdd.fits')    #!!!need to change

lens_mask = 1-cr_mask(lens_image, filename='../lens_mask.reg')
plt.imshow(lens_mask, origin='low')
plt.close()

framesize = 141
ct = (len(lens_image) - framesize)/2
lens_image = lens_image[ct:-ct,ct:-ct]
lens_rms = lens_rms[ct:-ct,ct:-ct]
lens_mask = lens_mask[ct:-ct,ct:-ct]

print "plot fitting image:"
plt.imshow(lens_image*lens_mask, origin='low', norm=LogNorm())
plt.close()


#%%
# =============================================================================
#  Fitting the QSO position
# =============================================================================
from mask_objects import find_loc_max
x, y =find_loc_max(lens_image)
center = len(lens_image)/2

plt.imshow(lens_image, norm = LogNorm(), origin='lower')
for i in range(len(x)):
    plt.plot(x[i], y[i], 'ro')
plt.close()

#from my_astro_func import read_pixel_scale
#pixscl = read_pixel_scale('../data/HE0047_F814W_sci.fits')
pix_sz = pixscl
x_s, y_s = [], []
for i in range(len(x)):
    x0, y0 =  -(float(x[i]) - center) * pix_sz , (float(y[i]) - center) * pix_sz
    x_s.append(x0)
    y_s.append(y0)

x_obj, y_obj =  -(float(52) - center) * pix_sz , (float(134) - center) * pix_sz    
# print x_obj, y_obj
#%%
#Things should be changed:
for fix_gamma in [2.1]:
    for psf_id in [0,1,2,3,4,5,6,7,8]:
        for subg in [2]:
            picklename='result_PSF{0}_PSFrecons_gammafix{1}_subg{2}.pkl'.format(psf_id, fix_gamma, subg)        
            if glob.glob(picklename) != []:
                continue
            print "psf_id:", psf_id, "fix_gamma:", fix_gamma, "subg:", subg
            psf = pyfits.getdata('../PSF{0}_clean.fits'.format(psf_id))
            numPix = len(lens_image)  # cutout pixel size
            deltaPix = pixscl
            kwargs_data = sim_util.data_configure_simple(numPix, deltaPix,inverse=True)
            kwargs_data['image_data'] = lens_image
            kwargs_data['noise_map'] = lens_rms
            data_class = ImageData(**kwargs_data)
            kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': psf, 'pixel_size': deltaPix}
            psf_class = PSF(**kwargs_psf)
            
            #%%This part would be different for lens to lens:
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
            
            kwargs_numerics = {'supersampling_factor': subg, 'supersampling_convolution': False}
            imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                            lens_light_model_class,
                                            point_source_class, kwargs_numerics=kwargs_numerics)
            
            kwargs_model = {'lens_model_list': lens_model_list,
                                           'source_light_model_list': source_model_list,
                                           'lens_light_model_list': lens_light_model_list,
                                           'point_source_model_list': point_source_list,
                                           'additional_images_list': [False],
                                           'fixed_magnification_list': [False],  # list of bools (same length as point_source_type_list). If True, magnification ratio of point sources is fixed to the one given by the lens model 
                                         }
            
            
            kwargs_constraints = {'joint_lens_light_with_lens_light': [[0, 1, ['center_x', 'center_y']]],
                                  'joint_source_with_point_source': [[0, 0]],
                                          'num_point_source_list': [4],
                                          }
            
    #        kwargs_likelihood = {'check_bounds': True,
    #                             'force_no_add_image': False,
    #                             'source_marg': False,
    #                             'position_uncertainty': 0.005,
    #                             'check_solver': True,
    #                             'solver_tolerance': 0.001,
    #                              'source_position_likelihood': True,
    #                              'image_likelihood_mask_list': [lens_mask]
    #                             }
            kwargs_likelihood = {'check_bounds': True,
                 'force_no_add_image': False,
                 'source_marg': False,
                 'image_position_uncertainty': 0.004,
                 'check_matched_source_position': True,
                 'source_position_tolerance': 0.001,
                 'source_position_sigma': 0.001,
                 'image_likelihood_mask_list': [lens_mask],
                 'check_positive_flux': True
                 }
                    
            
            image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
            multi_band_list = [image_band]
            kwargs_data_joint = {'multi_band_list': multi_band_list, 'multi_band_type': 'multi-linear'}
            
            kwargs_lens_init = [{'center_x': 0.0,
                               'center_y': -0.0,
                               'e1': 0.0,
                               'e2': -0.0,
                               'gamma': 1.95,
                               'theta_E': 1.4300127521568509},
                              {'dec_0': 0,
                               'gamma1': -0.0,
                               'gamma2': -0.0,
                               'ra_0': 0}]
            kwargs_source_init = [{'R_sersic': 0.4923466696533954,
                                   'center_x': 0.1551730435256952,
                                   'center_y': -0.0,
                                   'e1': -0.0,
                                   'e2': 0.0,
                                   'n_sersic': 1.4880395263822381}]
            
            kwargs_lens_light_init = [{'R_sersic': 0.6279035843511804,
                                       'center_x': 0.0,
                                       'center_y': -0.0,
                                       'e1': 0.0,
                                       'e2': 0.0,
                                       'n_sersic': 4},
                                      {'R_sersic': 2.796422092726421,
                                       'center_x': 0.0,
                                       'center_y': -0.0,
                                       'e1': 0.0,
                                       'e2': 0.0,
                                       'n_sersic': 1},
                                      {'R_sersic': 0.29485054460759136,
                                       'center_x': x_obj,
                                       'center_y': y_obj,
                                       'e1': 0.0,
                                       'e2': 0.0,
                                       'n_sersic': 1.7251055896870575}]
            
            kwargs_ps_init =    [{'dec_image': np.array(y_s),
                                   'ra_image': np.array(x_s)}]
                    
            # initial spread in parameter estimation #
            kwargs_lens_sigma = [{'theta_E': 0.1, 'e1': 0.2, 'e2': 0.2, 'gamma': .1, 'center_x': 0.1, 'center_y': 0.1},
                                 {'gamma1': 0.1, 'gamma2': 0.1}]
            kwargs_source_sigma = [{'R_sersic': 0.2, 'n_sersic': .5, 'center_x': .1, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2}]
            #kwargs_source_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
            kwargs_lens_light_sigma = [{'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2}]
            kwargs_lens_light_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
            kwargs_lens_light_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
            kwargs_ps_sigma = [{'ra_image': [0.02] * 4, 'dec_image': [0.02] * 4}]
            
            # hard bound lower limit in parameter space #
            kwargs_lower_lens = [{'theta_E': 0, 'e1': -0.5, 'e2': -0.5, 'gamma': 1.5, 'center_x': -10., 'center_y': -10},
            #                     {'theta_E': 0, 'center_x': -10., 'center_y': -10},
                                 {'gamma1': -0.2, 'gamma2': -0.2}]
            kwargs_lower_source = [{'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10}]
            #kwargs_lower_source.append({'R_sersic': 0.001,'n_sersic': 0.5,  'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
            kwargs_lower_lens_light = [{'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10}]
            kwargs_lower_lens_light.append({'R_sersic': 0.001,'n_sersic': 0.5,  'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
            kwargs_lower_lens_light.append({'R_sersic': 0.001,'n_sersic': 0.5,  'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10})
            kwargs_lower_ps = [{'ra_image': -10 * np.ones(4), 'dec_image': -10 * np.ones(4)}]
            
            # hard bound upper limit in parameter space #
            kwargs_upper_lens = [{'theta_E': 10, 'e1': 0.5, 'e2': 0.5, 'gamma': 2.5, 'center_x': 10., 'center_y': 10},
            #                     {'theta_E': 10, 'center_x': 10., 'center_y': 10},
                                 {'gamma1': 0.2, 'gamma2': 0.2}]
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
                                 fixed_lens_light, kwargs_lower_lens_light, kwargs_upper_lens_light]
            ps_params = [kwargs_ps_init, kwargs_ps_sigma, [{}], kwargs_lower_ps, kwargs_upper_ps]
            
            #%%
            kwargs_params = {'lens_model': lens_params,
                            'source_model': source_params,
                            'lens_light_model': lens_light_params,
                            'point_source_model': ps_params}
            
            fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
            
            fitting_kwargs_list_0 = [
                                    ['PSO', {'sigma_scale': 1., 'n_particles': 300, 'n_iterations': 300}]
                                    ]
            
            fitting_seq.fit_sequence(fitting_kwargs_list_0)
            
            kwargs_psf_iter = {'num_iter': 100, 'psf_iter_factor': 0.2,
                                'stacking_method': 'median', 
                               'keep_psf_error_map': False, 
                               'psf_symmetry': 1, 
                               'block_center_neighbour': 0.05}
            
            fitting_kwargs_list_1 = [
                                    ['psf_iteration', kwargs_psf_iter],
                                    ['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}],
                                    ['psf_iteration', kwargs_psf_iter],
                                    ['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}],
                                    ['MCMC', {'n_burn': 20, 'n_run': 20, 'walkerRatio': 4, 'sigma_scale': .1}]
                                    ]
            
            chain_list = fitting_seq.fit_sequence(fitting_kwargs_list_1)
            kwargs_result = fitting_seq.best_fit()
            
            sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]    
            lens_result = kwargs_result['kwargs_lens']
            lens_light_result = kwargs_result['kwargs_lens_light']
            source_result = kwargs_result['kwargs_source']
            ps_result = kwargs_result['kwargs_ps']
            
            kwargs_data, kwargs_psf_updated, kwargs_numerics = fitting_seq.multi_band_list[0]
            
            ##Update the PSF: Update the imageModel with the new PSF.
            psf_class = PSF(**kwargs_psf_updated)
            imageLinearFit = ImageLinearFit(data_class=data_class, psf_class=psf_class, lens_model_class = lens_model_class,
                                            source_model_class = source_model_class,  point_source_class= point_source_class,
                                            kwargs_numerics=kwargs_numerics)  
            imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                            lens_light_model_class,
                                            point_source_class, kwargs_numerics=kwargs_numerics)
            
            image_band = [kwargs_data, kwargs_psf_updated, kwargs_numerics]
            multi_band_list = [image_band]        
            
            modelPlot = ModelPlot(multi_band_list, kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat", likelihood_mask_list=[lens_mask])
            param_class = fitting_seq.param_class
            
            f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)
            modelPlot.data_plot(ax=axes[0,0])
            modelPlot.model_plot(ax=axes[0,1])
            modelPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
            modelPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=100)
            modelPlot.convergence_plot(ax=axes[1, 1], v_max=1)
            modelPlot.magnification_plot(ax=axes[1, 2])
            f.tight_layout()
            f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
            plt.savefig('fig_PSF{0}_PSFrecons_gammafix{1}_subg{2}_img0.pdf'.format(psf_id, fix_gamma, subg))
            plt.close()
            
            f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)
            modelPlot.decomposition_plot(ax=axes[0,0], text='Lens light', lens_light_add=True, unconvolved=True)
            modelPlot.decomposition_plot(ax=axes[1,0], text='Lens light convolved', lens_light_add=True)
            modelPlot.decomposition_plot(ax=axes[0,1], text='Source light', source_add=True, unconvolved=True)
            modelPlot.decomposition_plot(ax=axes[1,1], text='Source light convolved', source_add=True)
            modelPlot.decomposition_plot(ax=axes[0,2], text='All components', source_add=True, lens_light_add=True, unconvolved=True)
            modelPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
            f.tight_layout()
            f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
            plt.savefig('fig_PSF{0}_PSFrecons_gammafix{1}_subg{2}_img1.pdf'.format(psf_id, fix_gamma, subg))
            plt.close()
            mcmc_new_list = []
            param = Param(kwargs_model, kwargs_fixed_lens = fixed_lens, kwargs_fixed_lens_light=fixed_lens_light, **kwargs_constraints)
            if param.num_param()[0] != len(samples_mcmc[0]):
                raise ValueError("The param shape is not as the samples_mcmc_ones.")
                
            labels_new = ["AGN flux in image plane", "Host flux image plance", "AGN flux in souce plane", "Host flux souce plane",  "Host Sersic", "Sersic Reff"]
            lensModel_simple = LensModel(lens_model_list)
            for i in range(len(samples_mcmc)):
                param_out = param.args2kwargs(samples_mcmc[i])
                image_reconstructed, _, _, _ = imageLinearFit.image_linear_solve(kwargs_lens=param_out['kwargs_lens'], kwargs_source=param_out['kwargs_source'], kwargs_lens_light=param_out['kwargs_lens_light'], kwargs_ps=param_out['kwargs_ps'])
                image_ps = imageModel.point_source(param_out['kwargs_ps'], kwargs_lens=param_out['kwargs_lens'], unconvolved=False)
                flux_quasar = np.sum(image_ps)
                image_host_source_plane = imageModel.source_surface_brightness(param_out['kwargs_source'], param_out['kwargs_lens'], de_lensed=True,unconvolved=False)
                host_flux_s = np.sum(image_host_source_plane)
                image_host_image_plane = imageModel.source_surface_brightness(param_out['kwargs_source'], param_out['kwargs_lens'], de_lensed=False,unconvolved=False)
                host_flux_i = np.sum(image_host_image_plane)
                x_image, y_image, AGN_fluxs = param_out['kwargs_ps'][0]['ra_image'], param_out['kwargs_ps'][0]['dec_image'], param_out['kwargs_ps'][0]['point_amp']
                mag_macro = lensModel_simple.magnification(x_image, y_image, param_out['kwargs_lens'])
                AGN_fluxi_s = AGN_fluxs/mag_macro
                AGN_fluxs_s = np.average(np.abs(AGN_fluxi_s))
            #    print i, host_flux_s, host_flux_i
                mcmc_new_list.append([flux_quasar,host_flux_i, AGN_fluxs_s, host_flux_s, param_out['kwargs_source'][0]['n_sersic'],param_out['kwargs_source'][0]['R_sersic']])
                if i/1000 > (i-1)/1000 :
                    print "finished translate:", i
            
            fig = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
            fig.savefig('fig_PSF{0}_PSFrecons_gammafix{1}_subg{2}_corner.pdf'.format(psf_id, fix_gamma, subg))
            plt.close()
            
            fit_result = [kwargs_result, chain_list]
            trans_result = [mcmc_new_list, labels_new]
            pickle.dump([fit_result, trans_result, 
                         [kwargs_data, kwargs_psf_updated, kwargs_numerics, kwargs_model, lens_mask],
                         [lens_model_list, source_model_list, lens_light_model_list, point_source_list]], open(picklename, 'wb'))
            t2 = time.time()
            print "time cost (s):", round((t2-t1),3)

#import os
#os.system('say "your program of PSF{0} has finished"'.format(psf_id))
