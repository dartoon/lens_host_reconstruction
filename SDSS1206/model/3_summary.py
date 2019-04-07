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
import glob
from read_chisq import return_chisq 

PSFno_ = [0]
subg_ = [2,3]

#The values to pick up:
pick = 4
#pick_names =  ["AGN flux in image plane", "Host flux image plance",
#               "AGN flux in souce plane", "Host flux souce plane", 
#               "Host Sersic", "Sersic Reff"] 
labels = []
fit_value_l, fit_value_m, fit_value_h, = [], [], []
chisq = []

folder_i = ['fit_PSFi_QSOmask', 'fit_PSFi_PSFrecons']
for folder in folder_i:
    filenames = glob.glob('{0}/result_PSF?_fit_again_UPversion*'.format(folder))
    filenames.sort()
    for filename in filenames:
        PSFtyp = filename.split('/')[1].split('_')[1]
        if int(PSFtyp.split('PSF')[1]) in PSFno_:
            print filename
            subg = filename[(filename).find('subg')+len('subg')]
            fit_type = filename.split('/')[0].split('_')[-1]
            result = pickle.load(open(filename,'rb'))
            chisq.append(repr(round(return_chisq(filename),3)))
            if len(result) == 2:
                fit_result, trans_result = result
            elif len(result) ==3:
                fit_result, trans_result, kwargs_psf_updated = result
            fit_value = np.asarray(trans_result[0])[:,pick]
            fit_value_l.append(np.percentile(fit_value,16,axis=0))
            fit_value_m.append(np.percentile(fit_value,50,axis=0))
            fit_value_h.append(np.percentile(fit_value,84,axis=0))
            labels.append("{0}, subg{1}, {2}".format(PSFtyp, subg, fit_type))

#folder_ave = ['fit_PSFave_QSOmask', 'fit_PSFave_PSFrecons']
#for folder in folder_ave:
#    filenames = glob.glob('{0}/result*'.format(folder))
#    filenames.sort()
#    for filename in filenames:
##        print filename
#        PSFtyp = filename.split('/')[1].split('_')[1]
#        subg = filename[(filename).find('subg')+len('subg')]
#        fit_type = filename.split('/')[0].split('_')[-1]
#        result = pickle.load(open(filename,'rb'))
#        if len(result) == 2:
#            fit_result, trans_result = result
#        elif len(result) ==3:
#            fit_result, trans_result, kwargs_psf_updated = result
#        fit_value = np.asarray(trans_result[0])[:,pick]
#        fit_value_l.append(np.percentile(fit_value,16,axis=0))
#        fit_value_m.append(np.percentile(fit_value,50,axis=0))
#        fit_value_h.append(np.percentile(fit_value,84,axis=0))
#        labels.append("{0}, subg{1}, {2}".format(PSFtyp, subg, fit_type))
#%%
pick_names = trans_result[1]
fit_value_l = [x for _,x in sorted(zip(labels,fit_value_l))]
fit_value_m = [x for _,x in sorted(zip(labels,fit_value_m))]
fit_value_h = [x for _,x in sorted(zip(labels,fit_value_h))]
chisq = [x for _,x in sorted(zip(labels,chisq))]

labels = [x for _,x in sorted(zip(labels,labels))]  # Sort the label at the last

#bars = [labels[i].split(',')[0] for i in range(len(labels))]
#bars = sorted(set(bars), key=bars.index)

if len(labels)/4. !=  int(len(labels)/4.):
    raise ValueError("The labels is not 4 times")
else:
    length = len(labels)/4
    index0 = [i*4 for i in range(length)]
    label0 = [labels[i].split(',')[0] for i in index0]
    index1 = [i*4+3 for i in range(length)]
    label1 = [labels[i].split(',')[0] for i in index1]
    if label0 != label1:
        print label0, label1
        raise ValueError("The labels is wrong")
    else:
        bars = label0
    
plt.figure(figsize=(10, 8))
length = len(fit_value_l)/4
x_pos = np.arange(length) + 0.5
fmt =['--','--','-','-']
for ftype in range(4):
    index = [i*4+ftype for i in range(length)]
    value_result = [fit_value_m[i] for i in index]
    asm_error = [[fit_value_m[i]-fit_value_l[i] for i in index],
                        [fit_value_h[i]-fit_value_m[i] for i in index]]
    plt.errorbar(x_pos+0.01*ftype, value_result, yerr=asm_error, fmt =fmt[ftype],
                 label="{0}, {1}".format(labels[ftype].split(',')[1],labels[ftype].split(',')[2]))
    plt.text((x_pos+0.01*ftype)[0], value_result[0]*9.2/9, chisq[ftype],fontsize=15)
##If want to put horizontal line:
#fill_ref = True
#if fill_ref:
#    if pick == 3 or pick == 3-len(pick_names):
#        ref_value_l,ref_value, ref_value_h = 33.1, (43.8+33.1)/2, 43.8
#    elif pick == 4 or pick == 4-len(pick_names): # Sersic index
#        ref_value_l,ref_value, ref_value_h = 1.3, 1.45, 1.6
#    elif pick == 5 or pick == 5-len(pick_names):  # Sersic Reff
#        ref_value_l,ref_value, ref_value_h = 0.53, (0.53+0.6)/2 , 0.60
#    xs = np.linspace(x_pos[0], x_pos[-1])
#    ys_l = xs*0 + ref_value_l
#    ys_m = xs*0 + ref_value
#    ys_h = xs*0 + ref_value_h
#    plt.plot(xs,ys_l, xs, ys_m, xs,ys_h,color = 'gray')
#    plt.fill_between(xs, ys_l, ys_h, facecolor='gray', alpha = 0.2)
plt.legend(numpoints=1,ncol=2,loc=2,prop={'size':18})
plt.title('The fitting result with different PSF assumption.')
plt.ylabel(pick_names[pick],fontsize=15)
plt.ylim(np.min(fit_value_l)/1.3,np.max(fit_value_m)*1.3)
plt.tick_params(labelsize=15)
plt.xticks(x_pos, bars)
plt.show()


