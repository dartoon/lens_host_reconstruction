#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 10:13:22 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from specutils import Spectrum1D
import glob

import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

#%%To estimate the AB mag for a given stellar template
#filename  = 'SED_files/summary_*_PA00.fits'
#filename = glob.glob(filename)[0]
#hdul_sum = pyfits.open(filename)
#name_sum = hdul_sum[1].columns
#table_sum = hdul_sum[1].data
#z_idx = [i for i in range(len(name_sum)) if 'zmc' in str(name_sum[i])][0]
#mel_idx = [i for i in range(len(name_sum)) if 'logZsun' in str(name_sum[i])][0]
#z = table_sum[1][z_idx] #Redshift
#mel = table_sum[1][mel_idx] #Metallicity 
filename_p  = '../2_WFI2033_color/SFH_*_PA00_param.fits'
filename_p = glob.glob(filename_p)[0]
hdul = pyfits.open(filename_p)
table = hdul[1].data
name = hdul[1].columns
age_idx = [i for i in range(len(name)) if 'T_MW' in str(name[i])][0]
age = str(round(10**table[1][age_idx],3)) # 'AGE:' Gyr
z_idx = [i for i in range(len(name)) if 'zmc' in str(name[i])][0]
z = str(round(table[1][z_idx],2))
mel_idx = [i for i in range(len(name)) if 'Z_MW' in str(name[i])][0]
mel = str(round(10**table[1][mel_idx],3)) # 'Mel:' Z*/Z_sun

#%%
spec1d = pyfits.open("/Users/Dartoon/Astro/Packages/gsf/gsf/example/templates/gsf_spec_2033.fits")  
name_spec = spec1d[1].columns
table_spec = spec1d[1].data
array_spec = np.zeros((len(table_spec), 7))
for i in range(len(table_spec)):
    array_spec[i, :] = table_spec[i]

#%%
plt.figure(figsize=(10, 6))
array_spec[:,2] =  array_spec[:,2]/ array_spec[:,2].max() * 2.65
plt.plot(array_spec[:,0]/10000., array_spec[:,2])

f160w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/205.fil')
f160w_fil[:,2] = f160w_fil[:,2]/f160w_fil[:,2].max()* array_spec[:,2].max()/6
plt.plot(f160w_fil[60:,1]/10000., f160w_fil[60:,2], label='F160W response')

f555w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/235.fil')
f555w_fil[:,2] = f555w_fil[:,2]/f555w_fil[:,2].max()* array_spec[:,2].max()/6
plt.plot(f555w_fil[:,1]/10000., f555w_fil[:,2], label='F555W response')

f814w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/239.fil')
f814w_fil[:,2] = f814w_fil[:,2]/f814w_fil[:,2].max()* array_spec[:,2].max()/6
plt.plot(f814w_fil[:,1]/10000., f814w_fil[:,2], label='F814W response')

f140w_fil = np.loadtxt('/Users/Dartoon/Astro/Packages/gsf/gsf/example/filter/204.fil')
f140w_fil[:,2] = f140w_fil[:,2]/f140w_fil[:,2].max()* array_spec[:,2].max()/6
plt.plot(f140w_fil[:,1]/10000., f140w_fil[:,2], label='F140W response')


fnu = np.array([316.14340, 329.82858, 17.06660, 44.70776])
lam = np.array([13971.05, 15418.99, 5373.10, 8084.25])
yerr_l = np.array([155.291/562.887, 155.291/562.887, 155.291/562.887, 155.291/562.887])
mag = -2.5*np.log10(fnu) + 25
flambda = 10**(-0.4*(mag+2.402+5.0*np.log10(lam))) * 10**17
lam = lam/10000.
plt.scatter(lam, flambda, c='r', zorder =100)
plt.errorbar(lam, flambda, yerr=yerr_l*flambda,fmt='.',color='gray',markersize=1, zorder =90)


#plt.plot(sov_jwst_f144w_fil[:,0]/10000., sov_jwst_f144w_fil[:,1], label='NIRCam F444W response curve')
plt.xlim(0.25, 3)
xmin, xmax, ymin, ymax = plt.axis()
plt.text( (xmax-xmin)*0.7, (ymax-ymin)*0.6, 'stellar population with:', fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.53, 'z={0}'.format(z), fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.45, 'age={0:.3f} Gyr'.format(float(0.625)), fontsize=17)
plt.text( (xmax-xmin)*0.8, (ymax-ymin)*0.37, 'Z$_*$/Z$_\odot$={0}'.format(mel), fontsize=17)
plt.legend(prop={'size':15}, ncol=2, loc = 1)
plt.tick_params(labelsize=15)
plt.xlabel("um",fontsize=27)
plt.ylabel(r"f$_\lambda$ 10$^{\rm -17}$ erg/s/cm$^2$/$\AA$",fontsize=27)
#plt.yticks([])
plt.savefig('WFI2033_host_color.pdf')