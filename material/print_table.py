#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:41:40 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt

import glob
import sys
sys.path.insert(0,'../share_tools/')

from read_inference import read_inf
from lens_information import zeropoint_dict, filter_dict, lens_redshift

IDs = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']

full_name = ['HE0435$-$1223', 'RXJ1131$-$1231', 'WFI2033$-$4723', 'SDSS1206$+$4332', 
			 'HE1104$-$1805', 'SDSS0246$-$0825', 'HS2209$+$1914', 'HE0047$-$1756']


Camera = ['WFC3-IR', 'ACS', 'WFC3-IR', 'WFC3-IR', 'WFC3-IR', 'WFC3-UVIS', 'WFC3-UVIS', 'WFC3-UVIS']

Filter = ['F160W', 'F814W', 'F160W', 'F160W', 'F160W', 'F814W','F814W','F814W']

path = "/Volumes/Seagate_Expansion_Drive/Projects_backup/lens_host_reconstruction/"
exp_list = [9340.4, 1980.0, 26257.0, 8456.9, 14698.3, 8481.0, 4542.0+9696.0, 9712.0]

resolution = ['$0\\farcs{}08$', '$0\\farcs{}05$', '$0\\farcs{}08$', '$0\\farcs{}08$',
			  '$0\\farcs{}08$', '$0\\farcs{}03$', '$0\\farcs{}03$', '$0\\farcs{}03$']

#%%Table 1 :
print("Object ID & $z$ & camera & Filter & exposure & pixel scale \\\\")
print(" & & & & time (s) & (drizzled) \\\\ \\hline")
for i in range(len(IDs)):
#    files = glob.glob(path+IDs[i]+'/data/*sci.fits')
#    if len(files)>1:
#        print('Some thing is wrong!!! in:', IDs[i])
#    datafile = files[0]
#    fitsFile = pyfits.open(datafile)
#    exp_list.append(round(fitsFile[0].header['EXPTIME'],1))
    print(full_name[i], '&',lens_redshift[IDs[i][2:]]  , '&', Camera[i], '&', Filter[i] , '&', round(exp_list[i]), '&', resolution[i], "\\\\")

#%%Table 2
#Pring BH information
sys.path.insert(0,'../MBH_estimator/')
from est_MBH import BL_info, cal_MBH
print("Object ID & Line(s) used & FWHM & log($L_\\lambda$) & $\\log$\\mbh \\\\")
print("& & (\\kms) & ($[\\rm erg~s^[-1]]$) & (M$_{\\odot}$) \\\\ \\hline")
for ID in ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '7_HE0047']:
	ID = ID[2:]
	line_infor = BL_info[ID]
#	for key in BL_info[ID].keys():
	if len(BL_info[ID].keys()) == 1:
		key = list(BL_info[ID].keys())[0]
		print(ID, '&', '\\{0}'.format(key), '&', '{0:.0f}'.format(BL_info[ID][key]['FWHM']), '&',
		'{0:.2f}'.format(np.mean(BL_info[ID][key]['L'])), '&', '{0:.2f}'.format(cal_MBH(key, np.mean(BL_info[ID][key]['L']), BL_info[ID][key]['FWHM'])) ,"\\\\" )
	elif len(BL_info[ID].keys()) == 2:
		key1, key2 = list(BL_info[ID].keys())
		print(ID, '&', '\\{0}/\\{1}'.format(key1, key2), '&', '{0:.0f}/{1:.0f}'.format(BL_info[ID][key1]['FWHM'],BL_info[ID][key2]['FWHM']), '&',
		'{0:.2f}/{1:.2f}'.format(np.mean(BL_info[ID][key1]['L']),np.mean(BL_info[ID][key2]['L'])), 
		'&', '{0:.2f}/{1:.2f}'.format(cal_MBH(key1, np.mean(BL_info[ID][key1]['L']), BL_info[ID][key1]['FWHM']), cal_MBH(key2, np.mean(BL_info[ID][key2]['L']), BL_info[ID][key2]['FWHM'])), "\\\\" )
		

#%%Table 3:    
#Get the infernece from pickle: 
print("Object ID & Magnitude & Host-Total Flux Ratio & Reff & \\sersic\\ $n$ & adopted AGE & $\\log (M_{*}$)  \\\\")
print(" & (AB system) & ($\\%$) & (arcsec) & & (Gyr) & (M$_{\\odot}$) \\\\ \\hline")

for i in range(len(IDs)):
	ID = IDs[i][2:]
	host_flux_list = read_inf(IDs[i][2:],prop=0)
	AGN_flux_list = read_inf(IDs[i][2:],prop=3)
	Reff_list  = read_inf(IDs[i][2:],prop=1)
	n_list  = read_inf(IDs[i][2:],prop=2)
	zp =  zeropoint_dict[filter_dict[ID]]
	
	if ID != 'RXJ1131':
		idx = 0 
		sed_name = IDs[i]
		sed_filename  = '../Stellar_mass_estimate/'+ sed_name + '/SFH_*_PA00_param.fits'
		sed_filename = glob.glob(sed_filename)[0]
		hdul = pyfits.open(sed_filename)
		name = hdul[1].columns
		table = hdul[1].data	
		tab_name = hdul[1].columns
		mstar_idx = [ii for ii in range(len(tab_name)) if 'Mstel' in str(tab_name[ii])][0]
		m_star =  10**table[1][mstar_idx]
		age_idx = [ii for ii in range(len(tab_name)) if 'T_MW' in str(tab_name[ii])][0]
		age = 10**table[1][age_idx]
		
		host_flux_hml = [host_flux_list[idx][0]+ host_flux_list[idx][1], host_flux_list[idx][0], host_flux_list[idx][0] - host_flux_list[idx][1]]
		m_star_lmh = [ np.log10(m_star*host_flux_hml[2]/host_flux_hml[1]), np.log10(m_star), np.log10(m_star*host_flux_hml[0]/host_flux_hml[1]) ]
		m_star_mul = [ m_star_lmh[1], (m_star_lmh[2] - m_star_lmh[1]), (m_star_lmh[1] - m_star_lmh[0])]
		
		mag = [-2.5* np.log10(flux) + zp for flux in host_flux_hml]
		mag_mul = [(mag[1]), (mag[2]-mag[1]), (mag[0]-mag[1])]  #median, upper, lower
		flux_r_mD = [host_flux_list[idx][0]/(host_flux_list[idx][0]+AGN_flux_list[0][0])*100, host_flux_list[idx][1]/(host_flux_list[idx][0]+AGN_flux_list[0][0])*100]
		Reff_mD = [Reff_list[idx][0], Reff_list[idx][1]]
		n_mD = [n_list[idx][0], n_list[idx][1]]
		print(ID, '&','${0:.3f}\substack[+{1:.3f}\\\\{2:.3f}]$'.format(mag_mul[0], mag_mul[1], mag_mul[2]), '&' ,
				'${0:.3f}\pm{1:.3f}$'.format(flux_r_mD[0], flux_r_mD[1]), '&' ,
				'${0:.3f}\pm{1:.3f}$'.format(Reff_mD[0], Reff_mD[1]), '&' ,
				'${0:.3f}\pm{1:.3f}$'.format(n_mD[0], n_mD[1]), '&', 
				'${0:.3f}$'.format(age),
				'&','${0:.2f}\substack[+{1:.2f}\\\\-{2:.2f}]$'.format(m_star_mul[0], m_star_mul[1], m_star_mul[2]),
				'\\\\')
	if ID == 'RXJ1131':
		sed_name_l = ['1_RXJ1131_bulge', '1_RXJ1131_disk']
		flux_r_mD = [(host_flux_list[0][0]+host_flux_list[1][0])/(host_flux_list[0][0]+host_flux_list[1][0]+AGN_flux_list[0][0])*100, 
					   np.sqrt(host_flux_list[0][1]**2+host_flux_list[1][1]**2)/(host_flux_list[0][0]+host_flux_list[1][0]+AGN_flux_list[0][0])*100]

		n = [4, 1]
		ID_name = ['RXJ1131$_[\\rm bulge]$', 'RXJ1131$_[\\rm disk]$']
		for j in range(2):
			idx = j
			sed_name = sed_name_l[j]
			sed_filename  = '../Stellar_mass_estimate/'+ sed_name + '/SFH_*_PA00_param.fits'
			sed_filename = glob.glob(sed_filename)[0]
			hdul = pyfits.open(sed_filename)
			name = hdul[1].columns
			table = hdul[1].data	
			tab_name = hdul[1].columns
			mstar_idx = [ii for ii in range(len(tab_name)) if 'Mstel' in str(tab_name[ii])][0]
			m_star =  round(table[1][mstar_idx],3)
			age_idx = [ii for ii in range(len(tab_name)) if 'T_MW' in str(tab_name[ii])][0]
			age = round(10**table[1][age_idx],3)

			host_flux_hml = [host_flux_list[idx][0]+ host_flux_list[idx][1], host_flux_list[idx][0], host_flux_list[idx][0] - host_flux_list[idx][1]]
			m_star_lmh = [ np.log10(m_star*host_flux_hml[2]/host_flux_hml[1]), np.log10(m_star), np.log10(m_star*host_flux_hml[0]/host_flux_hml[1]) ]
			m_star_mul = [ m_star_lmh[1], (m_star_lmh[2] - m_star_lmh[1]), (m_star_lmh[1] - m_star_lmh[0])]

			
			host_flux_hml = [host_flux_list[idx][0]+ host_flux_list[idx][1], host_flux_list[idx][0], host_flux_list[idx][0] - host_flux_list[idx][1]]
			mag = [-2.5* np.log10(flux) + zp for flux in host_flux_hml]
			mag_mul = [(mag[1]), (mag[2]-mag[1]), (mag[0]-mag[1])]
			Reff_mD = [Reff_list[idx][0], Reff_list[idx][1]]
			flux_r_mD = [(host_flux_list[idx][0])/(host_flux_list[0][0]+host_flux_list[1][0]+AGN_flux_list[0][0])*100, 
						   (host_flux_list[idx][1])/(host_flux_list[0][0]+host_flux_list[1][0]+AGN_flux_list[0][0])*100]
			
			print(ID_name[j], '&','${0:.3f}\substack[+{1:.3f}\\\\{2:.3f}]$'.format(mag_mul[0], mag_mul[1], mag_mul[2]), '&' ,
					'${0:.3f}\pm{1:.3f}$'.format(flux_r_mD[0], flux_r_mD[1]), '&' ,
					'${0:.3f}\pm{1:.3f}$'.format(Reff_mD[0], Reff_mD[1]), '&' ,
					'fix to {0}'.format(n[j]), '&', 
					'${0:.3f}$'.format(age),
					'&','${0:.2f}\substack[+{1:.2f}\\\\-{2:.2f}]$'.format(m_star_mul[0], m_star_mul[1], m_star_mul[2]),
					'\\\\')		
