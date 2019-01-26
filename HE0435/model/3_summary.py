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

fit_type = ['QSO_mask', 'PSF_recon']
folder = ['fit_PSFi_QSO_mask/', 'fit_PSFi_PSFrecons/']
host_flux_source_l, host_flux_source_m, host_flux_source_h, = [], [], []
labels = []
for PSFno in PSFno_:
    for subg in subg_:
        filename = ['result_PSF{0}_QSOmask_gammafix_subg{1}'.format(PSFno,subg),
                    'result_PSF{0}_PSFrecons_gammafix_subg{1}'.format(PSFno,subg)]
        for i in range(len(filename)):
            readfile = folder[i]+filename[i]
            result = pickle.load(open(readfile,'rb'))
            if len(result) == 2:
                fit_result, trans_result = result
            elif len(result) ==3:
                fit_result, trans_result, kwargs_psf_updated = result
            host_flux_source = np.asarray(trans_result[0])[:,3]
            host_flux_source_l.append(np.percentile(host_flux_source,16,axis=0))
            host_flux_source_m.append(np.percentile(host_flux_source,50,axis=0))
            host_flux_source_h.append(np.percentile(host_flux_source,84,axis=0))
            labels.append("PSF{0}, subg{1}, {2}".format(PSFno, subg, fit_type[i]))

plt.figure(figsize=(10, 8))
bars = ['PSF{0}'.format(i) for i in PSFno_]
x_pos = np.arange(5) + 0.5
for ftype in range(4):
    index = [i*4+ftype for i in range(5)]
    flux_source = [host_flux_source_m[i] for i in index]
    asm_error = [[host_flux_source_m[i]-host_flux_source_l[i] for i in index],
                        [host_flux_source_h[i]-host_flux_source_m[i] for i in index]]
    plt.errorbar(x_pos, flux_source, yerr=asm_error,
                 label="{0}, {1}".format(labels[ftype].split(',')[1],labels[ftype].split(',')[2]))
plt.legend(numpoints=1,ncol=2,loc=2,prop={'size':18})
plt.title('My title')
plt.ylabel('Host flux values')
plt.ylim(0,160)
plt.xticks(x_pos, bars)
plt.show()
            

