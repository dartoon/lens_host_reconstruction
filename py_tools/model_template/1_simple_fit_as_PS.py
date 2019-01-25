#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 22:34:48 2019

@author: Dartoon

Fit the double as two Point source + a Sersic.
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import sys
sys.path.insert(0,'../../py_tools')
from fit_qso import fit_qso
lens_image = pyfits.getdata('HE1104_cutout.fits')
lens_rms = pyfits.getdata('HE1104_noisemap.fits')
psf = pyfits.getdata('PSF1.fits')
psf /= psf.sum()

fit_frame_size = 101
ct = (len(lens_image[2])-fit_frame_size)/2     # If want to cut to 61, QSO_im[ct:-ct,ct:-ct]

lens_image = lens_image[ct:-ct,ct:-ct]
lens_rms = lens_rms[ct:-ct,ct:-ct]
plt.imshow(lens_image, norm = LogNorm(), origin='lower')
plt.show()

# =============================================================================
#  The direct measure is that:
#      deflect position: 101, 101, (0,0)
#      lensed_QSO1 pos: 113, 107, (12,6), intensity 2382
#      lensed_QSO2 pos: 76, 90, (-25,-11) intensity 834
#      image distance = np.sqrt((12+25)**2 + (11+6)**2) = 40.71 pixels = 3.25
#      Einstein radius ~ 1.629''
# =============================================================================
fixed_source = []
kwargs_source_init = []
kwargs_source_sigma = []
kwargs_lower_source = []
kwargs_upper_source = []      
fixed_source.append({})  
kwargs_source_init.append({'R_sersic': 1, 'n_sersic': 2., 'e1': 0., 'e2': 0., 'center_x': 0., 'center_y': 0.})
kwargs_source_sigma.append({'n_sersic': 0.5, 'R_sersic': 0.5, 'e1': 0.1, 'e2': 0.1, 'center_x': 0.1, 'center_y': 0.1})
kwargs_lower_source.append({'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.3, 'center_x': -0.5, 'center_y': -0.5})
kwargs_upper_source.append({'e1': 0.5, 'e2': 0.5, 'R_sersic': 10., 'n_sersic': 7., 'center_x': 0.5, 'center_y': 0.5})

fixed_ps = []
kwargs_ps_init = []
kwargs_ps_sigma = []
kwargs_lower_ps = []
kwargs_upper_ps = []
fixed_ps.append({})
kwargs_ps_init.append({'ra_image': [-0.96], 'dec_image': [0.48], 'point_amp': [2382]})
kwargs_ps_sigma.append({'ra_image': [0.02], 'dec_image': [0.02]})
kwargs_lower_ps.append({'ra_image': [-1.26], 'dec_image': [0.18]})
kwargs_upper_ps.append({'ra_image': [-0.66], 'dec_image': [0.78]})

fixed_ps.append({})
kwargs_ps_init.append({'ra_image': [2.0], 'dec_image': [-0.88], 'point_amp': [834]})
kwargs_ps_sigma.append({'ra_image': [0.02], 'dec_image': [0.02]})
kwargs_lower_ps.append({'ra_image': [1.7], 'dec_image': [-1.18]})
kwargs_upper_ps.append({'ra_image': [2.3], 'dec_image': [-0.58]})
source_params = [kwargs_source_init, kwargs_source_sigma, fixed_source, kwargs_lower_source, kwargs_upper_source]

ps_param = [kwargs_ps_init, kwargs_ps_sigma, fixed_ps, kwargs_lower_ps, kwargs_upper_ps]            
source_result, ps_result, image_ps, image_host, error_map=fit_qso(lens_image, psf_ave=psf, psf_std = None, ps_param = ps_param,
                                                                  background_rms=0.00793274,
                                                                  source_params=source_params, QSO_msk = None, fixcenter=False,
                                                                  pix_sz = 0.08, no_MCMC =True,
                                                                  QSO_std =lens_rms, tag=None, deep_seed= True, pltshow = 1)   

# =============================================================================
# Derive the simple result
# =============================================================================
#The fitting result:
#[{'dec_image': array([ 0.52940031]),
#  'point_amp': array([ 2928.75990164]),
#  'ra_image': array([-0.97788635])},
# {'dec_image': array([-0.8035575]),
#  'point_amp': array([ 940.39042383]),
#  'ra_image': array([ 1.9284998])}]  
# now we fit the point sources with a simpler model (only SPEMD and SHEAR without the NFWs)
import copy
import numpy as np
for i in range(len(ps_result)):
    if i ==0:
        x_image, y_image = ps_result[i]['ra_image'], ps_result[i]['dec_image']
    else:
        x_image = np.append(x_image, ps_result[i]['ra_image'])
        y_image = np.append(y_image, ps_result[i]['dec_image'])
# import the lens model class 
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.lens_model_extensions import LensModelExtensions
# import the lens equation solver class (finding image plane positions of a source position)
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
# import lens model solver with 2 image positions constrains
from lenstronomy.LensModel.Solver.solver2point import Solver2Point
lens_model_list_simple = ['SPEMD', 'SHEAR']
lensModelSimple = LensModel(lens_model_list_simple)
solver2Point_simple = Solver2Point(lensModel=lensModelSimple, solver_type='ELLIPSE')

# options for Solver2Point are:
#    'CENTER': solves for 'center_x', 'center_y' parameters of the first lens model
#    'ELLIPSE': solves for 'e1', 'e2' of the first lens  (can also be shear)
#    'SHAPELETS': solves for shapelet coefficients c01, c10
#    'THETA_E_PHI: solves for Einstein radius of first lens model and shear angle of second model


kwargs_lens_init = [{'theta_E': 1.524, 'gamma': 2., 'e1': 0, 'e2': 0, 'center_x': 0., 'center_y': 0}, {'e1': 0.00, 'e2': -0.0}]
kwargs_fit_simple, precision = solver2Point_simple.constraint_lensmodel(x_pos=x_image, y_pos=y_image,
                                                                        kwargs_list=kwargs_lens_init, xtol=1.49012e-10)
print("the fitted macro-model parameters are: ", kwargs_fit_simple)

# check whether this simpler solution obeys the lens equation
lensModel_simple = LensModel(lens_model_list_simple)
beta_x_new, beta_y_new = lensModel_simple.ray_shooting(x_image, y_image, kwargs_fit_simple)

# we can now set a new estimate of the source position
beta_x_new = np.mean(beta_x_new)
beta_y_new = np.mean(beta_y_new)
print("The relative position in the source plane (should match) is: ", beta_x_new, beta_y_new)
# and solve for the new image positions (which should come very close to the true ones)
lensEquationSolver_new = LensEquationSolver(lensModel=lensModel_simple)
x_image_new, y_image_new = lensEquationSolver_new.image_position_from_source(kwargs_lens=kwargs_fit_simple,
                            sourcePos_x=beta_x_new, sourcePos_y=beta_y_new, min_distance=0.01, search_window=5,
                            precision_limit=10**(-10), num_iter_max=100)

#print "The fitted souce pos:", beta_x_new, beta_y_new
print "The input image:", x_image, y_image
print "The reconstruct image:", x_image_new, y_image_new

mag_macro = lensModel_simple.magnification(x_image, y_image, kwargs_fit_simple)
#print "mag use input pos:", mag_macro
#mag_macro = lensModel_simple.magnification(x_image_new, y_image_new, kwargs_fit_simple)
#print "mag use recons pos:", mag_macro

print "The Mag ratio:", mag_macro[0]/mag_macro[1] 
print "The Desired ratio:", 2382/834.
# =============================================================================
# Plot the fit image 
# =============================================================================
from mpl_toolkits.axes_grid1 import make_axes_locatable
from lenstronomy.SimulationAPI.simulations import Simulation
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.LightModel.light_model import LightModel
from lenstronomy.Sampling.parameters import Param
# =============================================================================
# Setting up the fitting list, and guess parameter to set up the fitting parameter.
# =============================================================================
ct = 30
lens_image = pyfits.getdata('HE1104_cutout.fits')[ct:-ct,ct:-ct]
lens_rms = pyfits.getdata('HE1104_noisemap.fits')[ct:-ct,ct:-ct]
psf = pyfits.getdata('PSF1.fits')

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
# =============================================================================
# Setting up the fitting list, and guess parameter to set up the fitting parameter.
# =============================================================================
lens_model_list = ['SPEMD', 'SHEAR']
kwargs_shear = kwargs_fit_simple[1]  # gamma_ext: shear strength, psi_ext: shear angel (in radian)
kwargs_spemd = kwargs_fit_simple[0]
kwargs_lens = [kwargs_spemd, kwargs_shear]
lens_model_class = LensModel(lens_model_list=lens_model_list)

lens_light_model_list = ['SERSIC_ELLIPSE']
kwargs_sersic = source_result[0]
kwargs_lens_light = [kwargs_sersic]
lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
# 'SERSIC_ELLIPSE': elliptical Sersic profile

source_model_list = ['SERSIC_ELLIPSE']
ra_source, dec_source = beta_x_new, beta_y_new
kwargs_sersic_ellipse = {'amp': 50., 'R_sersic': .2, 'n_sersic': 3, 'center_x': ra_source,
                         'center_y': dec_source, 'e1': -0.1, 'e2': 0.01}
kwargs_source = [kwargs_sersic_ellipse]
source_model_class = LightModel(light_model_list=source_model_list)

lensEquationSolver = LensEquationSolver(lens_model_class)
x_image, y_image = lensEquationSolver.findBrightImage(ra_source, dec_source, kwargs_lens, numImages=2,
                                                      min_distance=deltaPix, search_window=numPix * deltaPix)
mag = lens_model_class.magnification(x_image, y_image, kwargs=kwargs_lens)
kwargs_ps = [{'ra_image': x_image, 'dec_image': y_image,
                           'point_amp': np.abs(mag)*1}]  # quasar point source position in the source plane and intrinsic brightness
point_source_list = ['LENSED_POSITION']
point_source_class = PointSource(point_source_type_list=point_source_list, fixed_magnification_list=[False])
kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False}
imageModel = ImageModel(data_class, psf_class, lens_model_class, source_model_class,
                                lens_light_model_class,
                                point_source_class, kwargs_numerics=kwargs_numerics)
image_sim = SimAPI.simulate(imageModel, kwargs_lens, kwargs_source,
                                         kwargs_lens_light, kwargs_ps)
kwargs_data['image_data'] = image_sim
data_class.update_data(image_sim)

plt.matshow(image_sim, norm = LogNorm(), origin='lower')
plt.show()

# =============================================================================
# Print the well fitted parameter for the fitting in the next step.
# =============================================================================
print "lensed x, y position:", x_image, y_image
print "fitted lens light:", kwargs_sersic
print "Source position, ra_source, dec_source:", ra_source, dec_source
print "Lensing mass: SPEMD {0}, Shear {1}".format(kwargs_spemd, kwargs_shear)

