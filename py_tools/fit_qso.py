#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 20:23:36 2018

@author: Dartoon

A func for fitting the QSO with a input PSF. The PSF std is optional.
"""
from matplotlib.pylab import plt
import numpy as np
import corner

def fit_qso(QSO_im, psf_ave, psf_std=None, source_params=None,ps_param=None, background_rms=0.04, pix_sz = 0.168,
            exp_time = 300., fix_n=None, image_plot = True, corner_plot=True,
            flux_ratio_plot=True, deep_seed = False, fixcenter = False, QSO_msk=None, QSO_std=None,
            tag = None, no_MCMC= False,pltshow = 1):
    '''
    A quick fit for the QSO image with (so far) single sersice + one PSF. The input psf noise is optional.
    
    Parameter
    --------
        QSO_im: An array of the QSO image.
        psf_ave: The psf image.
        psf_std: The psf noise, optional.
        source_params: The prior for the source. Default is given.
        background_rms: default as 0.04
        exp_time: default at 2400.
        deep_seed: if Ture, more mcmc steps will be performed.
        tag: The name tag for save the plot
            
    Return
    --------
        Will output the fitted image (Set image_plot = True), the corner_plot and the flux_ratio_plot.
        source_result, ps_result, image_ps, image_host
    
    To do
    --------
        
    '''
    # data specifics need to set up based on the data situation
    background_rms = background_rms  #  background noise per pixel (Gaussian)
    exp_time = exp_time  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = len(QSO_im)  #  cutout pixel size
    deltaPix = pix_sz
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size = len(psf_ave)
    kernel = psf_ave
    
    if psf_std is not None:
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False, 'psf_error_map': True}     #Turn on the PSF error map
    else: 
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
    
    if source_params is None:
        # here are the options for the host galaxy fitting
        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        
        # Disk component, as modelled by an elliptical Sersic profile
        if fix_n == None:
            fixed_source.append({})  # we fix the Sersic index to n=1 (exponential)
            kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -10, 'center_y': -10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 10, 'center_y': 10})
        elif fix_n is not None:
            fixed_source.append({'n_sersic': fix_n})
            kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': fix_n, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.001, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': fix_n, 'center_x': -10, 'center_y': -10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3, 'n_sersic': fix_n, 'center_x': 10, 'center_y': 10})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    else:
        source_params = source_params
    
    if ps_param is None:
        center_x = 0.0
        center_y = 0.0
        point_amp = QSO_im.sum()/2.
        fixed_ps = [{}]
        kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
        kwargs_ps_init = kwargs_ps
        kwargs_ps_sigma = [{'ra_image': [0.05], 'dec_image': [0.05]}]
        kwargs_lower_ps = [{'ra_image': [-0.6], 'dec_image': [-0.6]}]
        kwargs_upper_ps = [{'ra_image': [0.6], 'dec_image': [0.6]}]
        ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
    else:
        ps_param = ps_param
    
    kwargs_params = {'source_model': source_params,
                     'point_source_model': ps_param}
    
    #==============================================================================
    #Doing the QSO fitting 
    #==============================================================================
    from lenstronomy.SimulationAPI.simulations import Simulation
    from lenstronomy.Data.imaging_data import Data
    SimAPI = Simulation()
    kwargs_data = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms)
    data_class = Data(kwargs_data)
    kwargs_psf =  SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=kernel_size, deltaPix=deltaPix, kernel=kernel)
    from lenstronomy.Data.psf import PSF
    psf_class = PSF(kwargs_psf)
    data_class.update_data(QSO_im)
    
    from lenstronomy.LightModel.light_model import LightModel
    light_model_list = ['SERSIC_ELLIPSE'] * len(source_params[0])
    lightModel = LightModel(light_model_list=light_model_list)
    from lenstronomy.PointSource.point_source import PointSource
    point_source_list = ['UNLENSED'] * len(ps_param[0])
    pointSource = PointSource(point_source_type_list=point_source_list)
    
    from lenstronomy.ImSim.image_model import ImageModel
    if QSO_msk is not None:
        kwargs_numerics['mask'] = QSO_msk
    
    kwargs_model = { 'source_light_model_list': light_model_list,
                    'point_source_model_list': point_source_list
                    }
    # numerical options and fitting sequences
    
    if fixcenter == False:
        kwargs_constraints = {'num_point_source_list': [1] * len(ps_param[0])
                              }
    elif fixcenter == True:
        kwargs_constraints = {'joint_source_with_point_source': [[0, 0]],
                              'num_point_source_list': [1]
                              }
    
    kwargs_likelihood = {'check_bounds': True,  #Set the bonds, if exceed, reutrn "penalty"
                         'source_marg': False,  #In likelihood_module.LikelihoodModule -- whether to fully invert the covariance matrix for marginalization
                          'check_positive_flux': True,        
                         }
    kwargs_data['image_data'] = QSO_im
    if QSO_std is not None:
        kwargs_data['noise_map'] = QSO_std
    
    if psf_std is not None:
        kwargs_psf['psf_error_map'] = psf_std
    
    
    imageModel = ImageModel(data_class, psf_class, source_model_class=lightModel,
                                    point_source_class=pointSource, kwargs_numerics=kwargs_numerics)
                  
    image_band = [kwargs_data, kwargs_psf, kwargs_numerics]
    multi_band_list = [image_band]
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
#    mpi = False  # MPI possible, but not supported through that notebook.
    # The Params for the fitting. kwargs_init: initial input. kwargs_sigma: The parameter uncertainty. kwargs_fixed: fixed parameters;
    #kwargs_lower,kwargs_upper: Lower and upper limits.

    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    if deep_seed == False:
        fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 0.8, 'n_particles': 80,
             'n_iterations': 60},
            {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}
            ]
    elif deep_seed == True:
         fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': 150,
             'n_iterations': 150},
            {'fitting_routine': 'MCMC', 'n_burn': 50, 'n_run': 100, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}
            ]
    if no_MCMC == True:
        fitting_kwargs_list = [fitting_kwargs_list[0],
                               ]        
    
    import time
    start_time = time.time()
    lens_result, source_result, lens_light_result, ps_result, cosmo_temp, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
    end_time = time.time()
    print(end_time - start_time, 'total time needed for computation')
    print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    # this is the linear inversion. The kwargs will be updated afterwards
    image_reconstructed, error_map, _, _ = imageModel.image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)
    image_ps = imageModel.point_source(ps_result)
    image_host = []
    for i in range(len(source_result)):
        image_host.append(imageModel.source_surface_brightness(source_result, k=i))
    # let's plot the output of the PSO minimizer
    from lenstronomy.Plots.output_plots import LensModelPlot
    lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                             lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")
    
    if image_plot:
        f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
        lensPlot.data_plot(ax=axes[0,0], text="Data")
        lensPlot.model_plot(ax=axes[0,1])
        lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
        
        lensPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
        lensPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
        lensPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
        
        lensPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
        lensPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
        lensPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
        
        f.tight_layout()
        #f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
        if tag is not None:
            f.savefig('{0}_fitted_image.pdf'.format(tag))
        if pltshow == 0:
            plt.close()
        else:
            plt.show()
        
    if corner_plot==True and no_MCMC==False:
        # here the (non-converged) MCMC chain of the non-linear parameters
        if not samples_mcmc == []:
           n, num_param = np.shape(samples_mcmc)
           plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
           if tag is not None:
               plot.savefig('{0}_para_corner.pdf'.format(tag))
           if pltshow == 0:
               plt.close()
           else:
               plt.show()
           
        
    if flux_ratio_plot==True and no_MCMC==False:
        from lenstronomy.Workflow.parameters import Param
        param = Param(kwargs_model, kwargs_constraints, kwargs_fixed_source=source_params[2], kwargs_fixed_ps=fixed_ps)
        mcmc_new_list = []
        labels_new = [r"Quasar flux", r"host_flux", r"source_x", r"source_y"]
        
        # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
        for i in range(len(samples_mcmc)/10):
            kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo = param.getParams(samples_mcmc[i+ len(samples_mcmc)/10*9])
            image_reconstructed, _, _, _ = imageModel.image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)
            
            image_ps = imageModel.point_source(kwargs_ps_out)
            flux_quasar = np.sum(image_ps)
            image_disk = imageModel.source_surface_brightness(kwargs_light_source_out, k=0)
            flux_disk = np.sum(image_disk)
            source_x = kwargs_ps_out[0]['ra_image']
            source_y = kwargs_ps_out[0]['dec_image']
            #    image_buldge = imageModel.source_surface_brightness(kwargs_light_source_out, k=1)
            #    flux_buldge = np.sum(image_buldge)
#            kwargs_ps_out
            if flux_disk>0:
                mcmc_new_list.append([flux_quasar, flux_disk, source_x, source_y])
        plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
        if tag is not None:
            plot.savefig('{0}_HOSTvsQSO_corner.pdf'.format(tag))
        if pltshow == 0:
            plt.close()
        else:
            plt.show()
    if QSO_std is None:
        return source_result, ps_result, image_ps, image_host, np.sqrt(data_class.C_D+np.abs(error_map))
    else:
        return source_result, ps_result, image_ps, image_host, np.sqrt(QSO_std**2+np.abs(error_map)) #error_map=0



def fit_qso_multiband(QSO_im_list, psf_ave_list, psf_std_list=None, source_params=None,ps_param=None,
                      background_rms_list=[0.04]*5, pix_sz = 0.168,
                      exp_time = 300., fix_n=None, image_plot = True, corner_plot=True,
                      flux_ratio_plot=True, deep_seed = False, fixcenter = False, QSO_msk_list=None,
                      QSO_std_list=None, tag = None, no_MCMC= False, pltshow = 1, new_band_seq=None):
    '''
    A quick fit for the QSO image with (so far) single sersice + one PSF. The input psf noise is optional.
    
    Parameter
    --------
        QSO_im: An array of the QSO image.
        psf_ave: The psf image.
        psf_std: The psf noise, optional.
        source_params: The prior for the source. Default is given.
        background_rms: default as 0.04
        exp_time: default at 2400.
        deep_seed: if Ture, more mcmc steps will be performed.
        tag: The name tag for save the plot
            
    Return
    --------
        Will output the fitted image (Set image_plot = True), the corner_plot and the flux_ratio_plot.
        source_result, ps_result, image_ps, image_host
    
    To do
    --------
        
    '''
    # data specifics need to set up based on the data situation
    background_rms_list = background_rms_list  #  background noise per pixel (Gaussian)
    exp_time = exp_time  #  exposure time (arbitrary units, flux per pixel is in units #photons/exp_time unit)
    numPix = len(QSO_im_list[0])  #  cutout pixel size
    deltaPix = pix_sz
    fwhm = 0.1  # full width half max of PSF (only valid when psf_type='gaussian')
    psf_type = 'PIXEL'  # 'gaussian', 'pixel', 'NONE'
    kernel_size_list = [len(psf_ave_list[i]) for i in range(len(psf_ave_list))]
    kernel_list = psf_ave_list
    if new_band_seq == None:
        new_band_seq= range(len(QSO_im_list))
    
    if psf_std_list is not None:
        kwargs_numerics_list = [{'subgrid_res': 1, 'psf_subgrid': False, 'psf_error_map': True}] * len(QSO_im_list)     #Turn on the PSF error map
    else: 
        kwargs_numerics_list = [{'subgrid_res': 1, 'psf_subgrid': False}] * len(QSO_im_list)
    
    if source_params is None:
        # here are the options for the host galaxy fitting
        fixed_source = []
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        
        # Disk component, as modelled by an elliptical Sersic profile
        if fix_n == None:
            fixed_source.append({})  # we fix the Sersic index to n=1 (exponential)
            kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': 0.3, 'center_x': -10, 'center_y': -10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3., 'n_sersic': 7., 'center_x': 10, 'center_y': 10})
        elif fix_n is not None:
            fixed_source.append({'n_sersic': fix_n})
            kwargs_source_init.append({'R_sersic': 0.3, 'n_sersic': fix_n, 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
            kwargs_source_sigma.append({'n_sersic': 0.001, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
            kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.1, 'n_sersic': fix_n, 'center_x': -10, 'center_y': -10})
            kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 3, 'n_sersic': fix_n, 'center_x': 10, 'center_y': 10})
        source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]
    else:
        source_params = source_params
    
    if ps_param is None:
        center_x = 0.0
        center_y = 0.0
        point_amp = QSO_im_list[0].sum()/2.
        fixed_ps = [{}]
        kwargs_ps = [{'ra_image': [center_x], 'dec_image': [center_y], 'point_amp': [point_amp]}]
        kwargs_ps_init = kwargs_ps
        kwargs_ps_sigma = [{'ra_image': [0.01], 'dec_image': [0.01]}]
        kwargs_lower_ps = [{'ra_image': [-10], 'dec_image': [-10]}]
        kwargs_upper_ps = [{'ra_image': [10], 'dec_image': [10]}]
        ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]
    else:
        ps_param = ps_param
    
    kwargs_params = {'source_model': source_params,
                     'point_source_model': ps_param}
    
    #==============================================================================
    #Doing the QSO fitting 
    #==============================================================================
    from lenstronomy.SimulationAPI.simulations import Simulation
    from lenstronomy.Data.imaging_data import Data
    from lenstronomy.Data.psf import PSF
    SimAPI = Simulation()
    kwargs_data_list, data_class_list = [], []
    for i in range(len(QSO_im_list)):
        kwargs_data_i = SimAPI.data_configure(numPix, deltaPix, exp_time, background_rms_list[i])
        kwargs_data_list.append(kwargs_data_i)
        data_class_list.append(Data(kwargs_data_i))
    kwargs_psf_list = []
    psf_class_list = []
    for i in range(len(QSO_im_list)):
        kwargs_psf_i =  SimAPI.psf_configure(psf_type=psf_type, fwhm=fwhm, kernelsize=kernel_size_list[i],
                                           deltaPix=deltaPix, kernel=kernel_list[i])
        kwargs_psf_list.append(kwargs_psf_i)
        psf_class_list.append(PSF(kwargs_psf_i))
        data_class_list[i].update_data(QSO_im_list[i])
    
    from lenstronomy.LightModel.light_model import LightModel
    light_model_list = ['SERSIC_ELLIPSE'] * len(source_params[0])
    lightModel = LightModel(light_model_list=light_model_list)
    from lenstronomy.PointSource.point_source import PointSource
    point_source_list = ['UNLENSED']
    pointSource = PointSource(point_source_type_list=point_source_list)
    
    from lenstronomy.ImSim.image_model import ImageModel
    imageModel_list = []
    for i in range(len(QSO_im_list)):
        kwargs_data_list[i]['image_data'] = QSO_im_list[i]
        if QSO_msk_list is not None:
            kwargs_numerics_list[i]['mask'] = QSO_msk_list[i]
        if QSO_std_list is not None:
            kwargs_data_list[i]['noise_map'] = QSO_std_list[i]
        if psf_std_list is not None:
            kwargs_psf_list[i]['psf_error_map'] = psf_std_list[i]
    
    image_band_list = []
    for i in range(len(QSO_im_list)):
        imageModel_list.append(ImageModel(data_class_list[i], psf_class_list[i], source_model_class=lightModel,
                                        point_source_class=pointSource, kwargs_numerics=kwargs_numerics_list[i]))
                  
        
        image_band_list.append([kwargs_data_list[i], kwargs_psf_list[i], kwargs_numerics_list[i]])
    multi_band_list = [image_band_list[i] for i in range(len(QSO_im_list))]
    
    # numerical options and fitting sequences
    
    kwargs_model = { 'source_light_model_list': light_model_list,
                    'point_source_model_list': point_source_list
                    }
    
    if fixcenter == False:
        kwargs_constraints = {'num_point_source_list': [1]
                              }
    elif fixcenter == True:
        kwargs_constraints = {'joint_source_with_point_source': [[0, 0]],
                              'num_point_source_list': [1]
                              }
    
    kwargs_likelihood = {'check_bounds': True,  #Set the bonds, if exceed, reutrn "penalty"
                         'source_marg': False,  #In likelihood_module.LikelihoodModule -- whether to fully invert the covariance matrix for marginalization
                          'check_positive_flux': True,        
                         }
    
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
#    mpi = False  # MPI possible, but not supported through that notebook.
    # The Params for the fitting. kwargs_init: initial input. kwargs_sigma: The parameter uncertainty. kwargs_fixed: fixed parameters;
    #kwargs_lower,kwargs_upper: Lower and upper limits.

    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)
    
    if deep_seed == False:
        fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 0.8, 'n_particles': 80,
            'n_iterations': 60, 'compute_bands': [True]+[False]*(len(QSO_im_list)-1)},
            {'fitting_routine': 'align_images', 'n_particles': 10, 'n_iterations': 10,
            'compute_bands': [False]+[True]*(len(QSO_im_list)-1)},
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 0.8, 'n_particles': 100,
            'n_iterations': 200, 'compute_bands': [True]*len(QSO_im_list)},
            {'fitting_routine': 'MCMC', 'n_burn': 10, 'n_run': 20, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}                
            ]
    elif deep_seed == True:
         fitting_kwargs_list = [
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 0.8, 'n_particles': 150,
            'n_iterations': 60, 'compute_bands': [True]+[False]*(len(QSO_im_list)-1)},
            {'fitting_routine': 'align_images', 'n_particles': 20, 'n_iterations': 20,
            'compute_bands': [False]+[True]*(len(QSO_im_list)-1)},
            {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 0.8, 'n_particles': 150,
            'n_iterations': 200, 'compute_bands': [True]*len(QSO_im_list)},
            {'fitting_routine': 'MCMC', 'n_burn': 50, 'n_run': 100, 'walkerRatio': 50, 'mpi': False,   ##Inputs  to CosmoHammer:
               #n_particles - particleCount; n_burn - burninIterations; n_run: sampleIterations (n_burn and n_run usually the same.); walkerRatio: walkersRatio.
            'sigma_scale': .1}             
            ]
    if no_MCMC == True:
        del fitting_kwargs_list[-1]
    
    import time
    start_time = time.time()
    lens_result, source_result, lens_light_result, ps_result, cosmo_temp, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(fitting_kwargs_list)
    end_time = time.time()
    print(end_time - start_time, 'total time needed for computation')
    print('============ CONGRATULATION, YOUR JOB WAS SUCCESSFUL ================ ')
    source_result_list, ps_result_list = [], []
    image_reconstructed_list, error_map_list, image_ps_list, image_host_list, shift_RADEC_list=[], [], [], [],[]
    for k in range(len(QSO_im_list)):
    # this is the linear inversion. The kwargs will be updated afterwards
        image_reconstructed_k, error_map_k, _, _ = imageModel_list[k].image_linear_solve(kwargs_source=source_result, kwargs_ps=ps_result)
        print "source_result", 'for', "k", source_result
        [kwargs_data_k, kwargs_psf_k, kwargs_numerics_k] = fitting_seq.multi_band_list[k]
        data_class_k = Data(kwargs_data_k)
        psf_class_k = PSF(kwargs_psf_k)
        imageModel_k = ImageModel(data_class_k, psf_class_k, source_model_class=lightModel,
                                point_source_class=pointSource, kwargs_numerics=kwargs_numerics_list[k])
        
        image_ps_k = imageModel_k.point_source(ps_result)
        image_host_k = []
        for i in range(len(source_result)):
            image_host_k.append(imageModel_list[k].source_surface_brightness(source_result, k=i))
        # let's plot the output of the PSO minimizer
        from lenstronomy.Plots.output_plots import LensModelPlot
        lensPlot = LensModelPlot(kwargs_data_list[k], kwargs_psf_list[k], kwargs_numerics_list[k], kwargs_model, lens_result, source_result,
                                 lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat")
        
        image_reconstructed_list.append(image_reconstructed_k)
        source_result_list.append(source_result)
        ps_result_list.append(ps_result)
        error_map_list.append(error_map_k)
        image_ps_list.append(image_ps_k)
        image_host_list.append(image_host_k)
        if 'ra_shift' in fitting_seq.multi_band_list[k][0].keys():
            shift_RADEC_list.append([fitting_seq.multi_band_list[k][0]['ra_shift'], fitting_seq.multi_band_list[k][0]['dec_shift']])
        else:
            shift_RADEC_list.append([0,0])
        if image_plot:
            f, axes = plt.subplots(3, 3, figsize=(16, 16), sharex=False, sharey=False)
            lensPlot.data_plot(ax=axes[0,0], text="Data")
            lensPlot.model_plot(ax=axes[0,1])
            lensPlot.normalized_residual_plot(ax=axes[0,2], v_min=-6, v_max=6)
            
            lensPlot.decomposition_plot(ax=axes[1,0], text='Host galaxy', source_add=True, unconvolved=True)
            lensPlot.decomposition_plot(ax=axes[1,1], text='Host galaxy convolved', source_add=True)
            lensPlot.decomposition_plot(ax=axes[1,2], text='All components convolved', source_add=True, lens_light_add=True, point_source_add=True)
            
            lensPlot.subtract_from_data_plot(ax=axes[2,0], text='Data - Point Source', point_source_add=True)
            lensPlot.subtract_from_data_plot(ax=axes[2,1], text='Data - host galaxy', source_add=True)
            lensPlot.subtract_from_data_plot(ax=axes[2,2], text='Data - host galaxy - Point Source', source_add=True, point_source_add=True)
            f.tight_layout()
            if tag is not None:
                f.savefig('{0}_fitted_image_band{1}.pdf'.format(tag,new_band_seq[k]))
            if pltshow == 0:
                plt.close()
            else:
                plt.show()
            
            if corner_plot==True and no_MCMC==False and k ==0:
                # here the (non-converged) MCMC chain of the non-linear parameters
                if not samples_mcmc == []:
                   n, num_param = np.shape(samples_mcmc)
                   plot = corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
                   if tag is not None:
                       plot.savefig('{0}_para_corner.pdf'.format(tag))
                   if pltshow == 0:
                       plt.close()
                   else:
                       plt.show()
            if flux_ratio_plot==True and no_MCMC==False:
                from lenstronomy.Workflow.parameters import Param
                param = Param(kwargs_model, kwargs_constraints, kwargs_fixed_source=source_params[2], kwargs_fixed_ps=fixed_ps)
                mcmc_new_list = []
                labels_new = [r"Quasar flux", r"host_flux", r"source_x", r"source_y"]
                # transform the parameter position of the MCMC chain in a lenstronomy convention with keyword arguments #
                for i in range(len(samples_mcmc)/10):
                    kwargs_lens_out, kwargs_light_source_out, kwargs_light_lens_out, kwargs_ps_out, kwargs_cosmo = param.getParams(samples_mcmc[i+ len(samples_mcmc)/10*9])
                    image_reconstructed, _, _, _ = imageModel_list[k].image_linear_solve(kwargs_source=kwargs_light_source_out, kwargs_ps=kwargs_ps_out)
                    
                    image_ps = imageModel_list[k].point_source(kwargs_ps_out)
                    flux_quasar = np.sum(image_ps)
                    image_disk = imageModel_list[k].source_surface_brightness(kwargs_light_source_out, k=0)
                    flux_disk = np.sum(image_disk)
                    source_x = kwargs_ps_out[0]['ra_image']
                    source_y = kwargs_ps_out[0]['dec_image']
                    #    image_buldge = imageModel.source_surface_brightness(kwargs_light_source_out, k=1)
                    #    flux_buldge = np.sum(image_buldge)
                    #kwargs_ps_out
                    if flux_disk>0:
                        mcmc_new_list.append([flux_quasar, flux_disk, source_x, source_y])
                plot = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
                if tag is not None:
                    plot.savefig('{0}_HOSTvsQSO_corner_band{1}.pdf'.format(tag,new_band_seq[k]))
                if pltshow == 0:
                    plt.close()
                else:
                    plt.show()
    errp_list = []
    for k in range(len(QSO_im_list)):
        if QSO_std_list is None:
            errp_list.append(np.sqrt(data_class_list[k].C_D+np.abs(error_map_list[k])))
        else:
            errp_list.append(np.sqrt(QSO_std_list[k]**2+np.abs(error_map_list[k])))
            
    return source_result_list, ps_result_list, image_ps_list, image_host_list, errp_list, shift_RADEC_list, fitting_seq     #fitting_seq.multi_band_list

    