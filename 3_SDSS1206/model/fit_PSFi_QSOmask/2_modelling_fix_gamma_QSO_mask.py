#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 17:30:19 2019

@author: Dartoon

Modelling the HE0435

Trying to recover the H0LiCOW_IV's host flux ~47.

The QSO center noise level is boost to inf.
psf_error_map not taken.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle
import sys
#sys.path.insert(0,'../../../../my_code/py_tools/')
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools/')
from mask_objects import mask_obj, find_loc_max
import copy, corner 
from psfs_average import psf_ave
from flux_profile import cr_mask
import glob

sys.path.insert(0,'../../../share_tools/')
from lens_information import pixel_scale

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

t1 = time.time()

lens_image = pyfits.getdata('../{0}_cutout.fits'.format(ID))   #!!!need to change
lens_rms = pyfits.getdata('../{0}_stdd.fits'.format(ID))    #!!!need to change

lens_mask = 1 - cr_mask(lens_image, '../lens_mask.reg')
obj_maks = cr_mask(lens_image, '../obj0_mask.reg')
lens_mask = lens_mask * obj_maks
plt.imshow(lens_mask, norm = LogNorm(), origin='lower')
plt.close()

framesize = 101
ct = (len(lens_image) - framesize)/2
lens_image = lens_image[ct:-ct,ct:-ct]
lens_rms = lens_rms[ct:-ct,ct:-ct]
lens_mask = lens_mask[ct:-ct,ct:-ct]


if glob.glob('{0}_stdd_boost.fits'.format(ID)) == []: 
    print "re-run the rms boost process"
    x_s, y_s = find_loc_max(lens_image)
    plt.imshow(lens_image, origin='low', norm=LogNorm())
    for i in range(len(x_s)):
        p_x, p_y = x_s[i], y_s[i]
        plt.text(p_x, p_y, 'obj{0}'.format(i), fontsize = 12)
    plt.show()
    qso_ids = [0,2]
    QSO_pos= []
    for i in qso_ids:
        QSO_pos.append([float(x_s[i]-0.5), float(y_s[i]-0.5)])
    xy_index = np.indices((len(lens_image),len(lens_image)))
    for i in range(len(QSO_pos)):
        if i == 0:
            areas = (np.sqrt((QSO_pos[i][1]-xy_index[0])**2+(QSO_pos[i][0]-xy_index[1])**2) <3. )  # five piexls
        else:
            areas += (np.sqrt((QSO_pos[i][1]-xy_index[0])**2+(QSO_pos[i][0]-xy_index[1])**2) <3. )  # five piexls
    lens_rms = lens_rms * (areas == 0) + 10**6 * (areas != 0)
    plt.imshow(lens_rms, origin='low', norm=LogNorm(), vmax = 0.2)
    plt.show()
    pyfits.PrimaryHDU(lens_rms).writeto('{0}_stdd_boost.fits'.format(ID),overwrite=False) 

lens_rms = pyfits.getdata('{0}_stdd_boost.fits'.format(ID))

#print "plot fitting image:"
#plt.imshow(lens_image*lens_mask, origin='low', norm=LogNorm())
#plt.close()

id_ind = int(sys.argv[1])
fix_gamma = [1.9, 2.0, 2.1][id_ind]
for psf_id in [0,1,2]:
    for subg in [3,2]:
        picklename='result_PSF{0}_QSOmask_gammafix{1}_subg{2}.pkl'.format(psf_id, fix_gamma, subg)        
        if glob.glob(picklename) != []:
            continue
        print "psf_id:", psf_id, "fix_gamma:", fix_gamma, "subg:", subg
        psf = pyfits.getdata('../PSF{0}_clean.fits'.format(psf_id))
        numPix = len(lens_image)  # cutout pixel size
        deltaPix = pixel_scale[ID]
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
        lens_model_list = ['SPEMD', 'SIS', 'SHEAR','SIS']
        lens_model_class = LensModel(lens_model_list=lens_model_list)
        lens_light_model_list = ['SERSIC_ELLIPSE','SERSIC_ELLIPSE','SERSIC']
        lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
        source_model_list = ['SERSIC_ELLIPSE', 'SERSIC_ELLIPSE']
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
        
        
        kwargs_constraints = {'joint_lens_light_with_lens_light': [[0,1,['center_x', 'center_y']]],
                              'joint_source_with_point_source': [[0, 0]],
                                      'num_point_source_list': [2],
                                      'joint_lens_with_light': [[2,1,['center_x', 'center_y']]],  #Light 2, Mass 1 (the small satellite)
                                      }

        kwargs_likelihood = {'check_bounds': True,
                             'force_no_add_image': False,
                             'source_marg': False,
                             'position_uncertainty': 0.005,
                             'check_solver': True,
                             'solver_tolerance': 0.001,
                              'source_position_likelihood': True,
                              'image_likelihood_mask_list': [lens_mask]
                             }
        
        image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
        multi_band_list = [image_band]
        kwargs_data_joint = {'multi_band_list': multi_band_list, 'multi_band_type': 'multi-linear'}
        
        kwargs_lens_init = [{'center_x': 0.025294402098616758,
                      'center_y': -0.08729829325047747,
                      'e1': 0.13444626263283052,
                      'e2': 0.0814560912560452,
                      'gamma': 2,
                      'theta_E': 1.3464756286523458},
                     {'center_x': -0.1417700702867299,
                      'center_y': 0.9217736735669627,
                      'theta_E': 0.22848433131987803},
                     {'dec_0': 0,
                      'e1': 0.009890975967993336,
                      'e2': 0.11343461385823016,
                      'ra_0': 0},
                      {'center_x': -1.6,
                      'center_y': 2.5,
                      'theta_E': 0.1248433131987803}
                      ]
        kwargs_source_init = [{'R_sersic': 0.022705023884574294,
                              'center_x': 0.025951627602455007,
                              'center_y': -0.18860759511690517,
                              'e1': -0.3952447566966397,
                              'e2': -0.15070924390456186,
                              'n_sersic': 4},
                               {'R_sersic': 0.25649613868935184,
                              'center_x': 0.001640436298278629,
                              'center_y': -0.19081466313798995,
                              'e1': -0.012136714690905138,
                              'e2': -0.12915566177270657,
                              'n_sersic': 1}]
        kwargs_lens_light_init = [{'R_sersic': 0.6359566788792801,
                              'center_x': 0.012154579012256295,
                              'center_y': 0.0036468079277633602,
                              'e1': 0.09966388423356236,
                              'e2': -0.08879198307974209,
                              'n_sersic': 4},
                             {'R_sersic': 0.9846407252383856,
                              'center_x': -0.0033623504780080655,
                              'center_y': 0.04104484433700866,
                              'e1': 0.24374681250138644,
                              'e2': -0.28042640127143054,
                              'n_sersic': 1},
                             {'R_sersic': 0.08463894365889412,
                              'center_x': -0.1417700702867299,
                              'center_y': 0.9217736735669627,
                              'n_sersic': 1}]
        kwargs_ps_init =[{'dec_image': np.array([ 1.20957114, -1.81481053]),
                          'ra_image': np.array([-0.51918733, -0.56006618])}]
                
        # initial spread in parameter estimation #
        kwargs_lens_sigma = [{'theta_E': 0.1, 'e1': 0.2, 'e2': 0.2, 'gamma': .1, 'center_x': 0.1, 'center_y': 0.1},
                             {'theta_E': 0.1, 'center_x': 0.1, 'center_y': 0.1},
                             {'e1': 0.1, 'e2': 0.1},
                             {'theta_E': 0.1, 'center_x': 0.1, 'center_y': 0.1}]
        kwargs_source_sigma = [{'R_sersic': 0.2, 'n_sersic': .5, 'center_x': .1, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2},
                               {'R_sersic': 0.2, 'n_sersic': .5, 'center_x': .1, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2}]
        kwargs_lens_light_sigma = [{'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2},
                                   {'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2},
                                   {'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1}]
        #kwargs_lens_light_sigma.append({'R_sersic': 0.2, 'center_x': .1, 'n_sersic': .5, 'center_y': 0.1, 'e1': 0.2, 'e2': 0.2})
        kwargs_ps_sigma = [{'ra_image': [0.02] * 2, 'dec_image': [0.02] * 2}]
        
        # hard bound lower limit in parameter space #
        kwargs_lower_lens = [{'theta_E': 0, 'e1': -0.5, 'e2': -0.5, 'gamma': 1.5, 'center_x': -10., 'center_y': -10},
                             {'theta_E': 0, 'center_x': -10., 'center_y': -10},
                             {'e1': -0.2, 'e2': -0.2},
                             {'theta_E': 0, 'center_x': -10., 'center_y': -10}]
        kwargs_lower_source = [{'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10},
                               {'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10}]
        kwargs_lower_lens_light = [{'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10},
                                   {'R_sersic': 0.001, 'n_sersic': 0.5, 'e1': -0.5, 'e2': -0.5, 'center_x': -10, 'center_y': -10},
                                   {'R_sersic': 0.001, 'n_sersic': 0.5, 'center_x': -10, 'center_y': -10}]
        kwargs_lower_ps = [{'ra_image': -10 * np.ones(2), 'dec_image': -10 * np.ones(2)}]
        
        # hard bound upper limit in parameter space #
        kwargs_upper_lens = [{'theta_E': 10, 'e1': 0.5, 'e2': 0.5, 'gamma': 2.5, 'center_x': 10., 'center_y': 10},
                             {'theta_E': 10, 'center_x': 10., 'center_y': 10},
                             {'e1': 0.2, 'e2': 0.2},
                             {'theta_E': 10, 'center_x': 10., 'center_y': 10}]
        kwargs_upper_source = [{'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10}
                                ,{'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10}]
        kwargs_upper_lens_light = [{'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10},
                                   {'R_sersic': 10, 'n_sersic': 5., 'e1': 0.5, 'e2': 0.5, 'center_x': 10, 'center_y': 10},
                                   {'R_sersic': 10, 'n_sersic': 5., 'center_x': 10, 'center_y': 10}]
        kwargs_upper_ps = [{'ra_image': 10 * np.ones(2), 'dec_image': 10 * np.ones(2)}]
        
        #Fix something:
        fixed_source = [{'n_sersic': 4},{'n_sersic': 1}]
        fixed_lens = [{'gamma': fix_gamma}, {}, {'ra_0': 0, 'dec_0': 0},{'center_x': -1.6, 'center_y': 2.5}]
        fixed_lens_light = [{'n_sersic': 4},{'n_sersic': 1},{'n_sersic': 1}]
                
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
        
        fitting_kwargs_list = [
                                ['PSO', {'sigma_scale': 1., 'n_particles': 300, 'n_iterations': 300}],
                                ['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}],
                                ['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 100}],                                
                                ['MCMC', {'n_burn': 20, 'n_run': 20, 'walkerRatio': 4, 'sigma_scale': .1}]
                                ]
        chain_list = fitting_seq.fit_sequence(fitting_kwargs_list)
        kwargs_result = fitting_seq.best_fit()
        
        sampler_type, samples_mcmc, param_mcmc, dist_mcmc  = chain_list[-1]    
        lens_result = kwargs_result['kwargs_lens']
        lens_light_result = kwargs_result['kwargs_lens_light']
        source_result = kwargs_result['kwargs_source']
        ps_result = kwargs_result['kwargs_ps']
        
        imageLinearFit = ImageLinearFit(data_class=data_class, psf_class=psf_class, lens_model_class = lens_model_class,
                                        source_model_class = source_model_class,  point_source_class= point_source_class,
                                        kwargs_numerics=kwargs_numerics)  
        
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
        plt.savefig('fig_PSF{0}_QSOmask_gammafix{1}_subg{2}_img0.pdf'.format(psf_id, fix_gamma, subg))
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
        plt.savefig('fig_PSF{0}_QSOmask_gammafix{1}_subg{2}_img1.pdf'.format(psf_id, fix_gamma, subg))
        plt.close()
        mcmc_new_list = []
        param = Param(kwargs_model, kwargs_fixed_lens = fixed_lens, kwargs_fixed_lens_light=fixed_lens_light, **kwargs_constraints)
        if param.num_param()[0] != len(samples_mcmc[0]):
            raise ValueError("The param shape is not as the samples_mcmc_ones.")
            
        labels_new = ["AGN flux in image plane", "Bulge flux image plance",  "Disk flux image plance",
               "AGN flux in souce plane", "Bulge flux souce plane", "Disk flux souce plane", "Bulge Reff", "Disk Reff"]  
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
        fig.savefig('fig_PSF{0}_QSOmask_gammafix{1}_subg{2}_corner.pdf'.format(psf_id, fix_gamma, subg))
        plt.close()
        
        fit_result = [kwargs_result, chain_list]
        trans_result = [mcmc_new_list, labels_new]
        pickle.dump([fit_result, trans_result, 
                     [kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_mask],
                     [lens_model_list, source_model_list, lens_light_model_list, point_source_list]], open(picklename, 'wb'))
        t2 = time.time()
        print "time cost (s):", round((t2-t1),3)

#import os
#os.system('say "your program of PSF{0} has finished"'.format(psf_id))
