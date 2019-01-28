#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 15:58:15 2019

@author: Dartoon

Summarize the fitting result
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle

PSFno_ = [0,1,2,3,4] 
subg_ = [2,3]

#The values to pick up:
pick = 3
pick_names =  ['AGN flux in image plane',
   'Host flux image plance',
   'AGN flux in souce plane',
   'Host flux souce plane',
   'Host Sersic',
   'Sersic Reff']


fit_type = ['QSO_mask', 'PSF_recon']
folder_i = ['fit_PSFi_QSO_mask/', 'fit_PSFi_PSFrecons/']
fit_value_l, fit_value_m, fit_value_h, = [], [], []
Reff = []
labels = []
for PSFno in PSFno_:
    for subg in subg_:
        filename = ['result_PSF{0}_QSOmask_gammafix_subg{1}'.format(PSFno,subg),
                    'result_PSF{0}_PSFrecons_gammafix_subg{1}'.format(PSFno,subg)]
        for i in range(len(filename)):
            readfile = folder_i[i]+filename[i]
            result = pickle.load(open(readfile,'rb'))
            if len(result) == 2:
                fit_result, trans_result = result
            elif len(result) ==3:
                fit_result, trans_result, kwargs_psf_updated = result
            fit_value = np.asarray(trans_result[0])[:,pick]
            fit_value_l.append(np.percentile(fit_value,16,axis=0))
            fit_value_m.append(np.percentile(fit_value,50,axis=0))
            fit_value_h.append(np.percentile(fit_value,84,axis=0))
            labels.append("PSF{0}, subg{1}, {2}".format(PSFno, subg, fit_type[i]))

folder_ave = ['fit_PSFave_QSO_mask/', 'fit_PSFave_PSFrecons/']
for PSFtyp in ['','std']:
    for subg in subg_:
        filename = ['result_PSFave{0}_QSOmask_gammafix_subg{1}'.format(PSFtyp,subg),
                    'result_PSFave{0}_PSFrecons_gammafix_subg{1}'.format(PSFtyp,subg)]
        for i in range(len(filename)):
            readfile = folder_ave[i]+filename[i]
            result = pickle.load(open(readfile,'rb'))
            if len(result) == 2:
                fit_result, trans_result = result
            elif len(result) ==3:
                fit_result, trans_result, kwargs_psf_updated = result
            fit_value = np.asarray(trans_result[0])[:,pick]
            fit_value_l.append(np.percentile(fit_value,16,axis=0))
            fit_value_m.append(np.percentile(fit_value,50,axis=0))
            fit_value_h.append(np.percentile(fit_value,84,axis=0))
            labels.append("PSFave{0}, subg{1}, {2}".format(PSFtyp, subg, fit_type[i]))

plt.figure(figsize=(10, 8))
bars = [labels[i*4].split(',')[0] for i in range(len(labels)/4)]
length = len(fit_value_l)/4
x_pos = np.arange(length) + 0.5
fmt =['--','--','-','-']
for ftype in range(4):
    index = [i*4+ftype for i in range(length)]
    flux_source = [fit_value_m[i] for i in index]
    asm_error = [[fit_value_m[i]-fit_value_l[i] for i in index],
                        [fit_value_h[i]-fit_value_m[i] for i in index]]
    plt.errorbar(x_pos, flux_source, yerr=asm_error, fmt =fmt[ftype],
                 label="{0}, {1}".format(labels[ftype].split(',')[1],labels[ftype].split(',')[2]))
plt.legend(numpoints=1,ncol=2,loc=2,prop={'size':18})
plt.title('The fitting result with different PSF assumption.')
plt.ylabel(pick_names[pick])
plt.ylim(0,np.max(fit_value_m)*1.3)
plt.xticks(x_pos, bars)
plt.show()



