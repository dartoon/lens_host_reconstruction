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
sys.path.insert(0,'/Users/Dartoon/Astro/my_code/py_tools')
from flux_profile import cr_mask

psf_i = 0 #raw_input("Which psf to load? 0, 1, 2, 3...:\n")

subg = 3 #int(filename.split('subg')[1])
filename = 'fit_PSFi_PSFrecons/result_PSF{0}_Shost_Slens_subg{1}'.format(psf_i, subg)
#filename = 'fit_PSFi_QSOmask/result_PSF{0}_Shost_Slens_subg{1}'.format(psf_i, subg)
result = pickle.load(open(filename,'rb'))

print filename

if len(result) == 2:
    fit_result, trans_result = result
elif len(result) ==3:
    fit_result, trans_result, kwargs_psf_updated = result

[lens_result, source_result, lens_light_result, ps_result,\
cosmo_result,chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc] = fit_result
 
[mcmc_new_list, labels_new] =  trans_result
#[imageModel, lensPlot, param] = model_sets

print "lens_result:", lens_result,'\n'
print "source_light_result:", source_result,'\n'
print "lens_light_result:", lens_light_result,'\n'
print "point_source_result:", ps_result,'\n'


#%%
#=============================================================================
from lenstronomy.SimulationAPI.simulations import Simulation
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LightModel.light_model import LightModel

lens_image = pyfits.getdata('SDSS0246_cutout.fits')
lens_rms = pyfits.getdata('SDSS0246_stdd.fits')

lens_mask = 1-cr_mask(lens_image, filename='lens_mask.reg')
mask1 = cr_mask(lens_image, filename='obj0_mask.reg')
lens_mask = lens_mask * mask1 #* mask2

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
#    psf = kwargs_psf_updated['kernel_point_source_init']
    folder_name = filename.split('/')[0].replace('QSOmask','PSFrecons')
    import glob
    filename_PSFave = glob.glob(folder_name+"/result*PSF{0}*".format(psf_i))[0]
    psf = pickle.load(open(filename_PSFave,'rb'))[2]['kernel_point_source_init']

numPix = len(lens_image)  # cutout pixel size
deltaPix = 0.03 # pixel size in arcsec (area per pixel = deltaPix**2)
fwhm = 0.16  # full width half max of PSF

if 'mask' in filename:
    from mask_objects import find_loc_max
    QSOx, QSOy =find_loc_max(lens_image)
    del QSOx[1]
    del QSOy[1]
    xy_index = np.indices((numPix,numPix))
    for i in range(len(QSOx)):
        if i == 0:
            areas = (np.sqrt((QSOx[i]-xy_index[1])**2+(QSOy[i]-xy_index[0])**2) <4 )  # 3 piexls
        else:
            areas += (np.sqrt((QSOx[i]-xy_index[1])**2+(QSOy[i]-xy_index[0])**2) <4 )  # 3 piexls
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
lens_light_model_list = ['SERSIC_ELLIPSE']
lens_light_model_class = LightModel(light_model_list=lens_light_model_list)
source_model_list = ['SERSIC_ELLIPSE']
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
lensPlot.source_plot(ax=axes[1, 0], deltaPix_source=0.01, numPix=100)
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

image_host_source0_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False,k=0)
host_flux0_total = image_host_source0_plane.sum()
image_host_source1_plane = imageModel.source_surface_brightness(source_result, lens_result, de_lensed=True,unconvolved=False,k=1)
host_flux1_total = image_host_source1_plane.sum()
#print "bulge, disk flux, source plane:", host_flux0_total, host_flux1_total

import corner
fig = corner.corner(mcmc_new_list, labels=labels_new, show_titles=True)
#fig.savefig('fig_PSF{0}_{2}{1}_corner.pdf'.format(psfno,subg,fname))
plt.show()