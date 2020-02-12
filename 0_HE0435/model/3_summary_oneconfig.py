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
import sys
import copy
from adjustText import adjust_text   # avoid the overlapping while ploting
sys.path.insert(0,'../../share_tools/')
from read_chisq import return_chisq 

pkl_folderfiles = glob.glob("fit_PSFi_PSFrecons/result_PSF*.pkl")
pkl_files = [pkl_folderfiles[i].split('/')[1] for i in range(len(pkl_folderfiles))]
psf_num = np.max([int(pkl_files[i].split('_PSF')[1]) for i in range(len(pkl_files))]) + 1

if float(len(pkl_files)/6) != len(pkl_files)/6. :
    raise ValueError("The number of the pickle result is the Multiples of six")

import os
path = os.getcwd()
ID = path.split('/')[-2]
filt = 'f160w'
#%%
sort_label_type=['subg_2, PSFrecons']

PSFno_ = range(psf_num)
#The values to pick up:
filename = pkl_folderfiles[0]
result = pickle.load(open(filename,'rb'))
fit_result, trans_result, kwargs_material, model_lists  = result
pick_names = trans_result[-1]

_, _, _, _, lens_mask = kwargs_material

fair_PSF = True

labels = []
fit_value_l_list, fit_value_m_list, fit_value_h_list, = [], [], []
chisq = []
fixgamma_list = []

folder_i = ['fit_PSFi_PSFrecons']
for folder in folder_i:
#    if folder == 'fit_PSFi_QSOmask':
#        filenames = glob.glob('{0}/result_PSF?_S*'.format(folder))
#    elif 'fit_PSFi_PSFrecons':
    filenames = glob.glob('{0}/result_PSF?_*subg{1}.pkl'.format(folder, sort_label_type[0][5]))        
    filenames.sort()
    for filename in filenames:
        PSFtyp = filename.split('/')[1].split('_')[1]
        if int(PSFtyp.split('PSF')[1]) in PSFno_:
            subg = filename[(filename).find('subg')+len('subg')]
            fit_type = filename.split('/')[0].split('_')[-1]
            fixgamma = (filename.split('gammafix')[1]).split('_subg')[0]
            if fixgamma not in fixgamma_list:
                fixgamma_list.append(fixgamma)
            print filename
            result = pickle.load(open(filename,'rb'))
            fit_result, trans_result, kwargs_material, model_lists  = result
            fit_value_l_i, fit_value_m_i, fit_value_h_i, = [], [], []
            for pick in range(len(pick_names)):
                fit_value = np.asarray(trans_result[0])[:,pick]
                fit_value_l_i.append(np.percentile(fit_value,16,axis=0))
                fit_value_m_i.append(np.percentile(fit_value,50,axis=0))
                fit_value_h_i.append(np.percentile(fit_value,84,axis=0))
            fit_value_l_list.append(fit_value_l_i)
            fit_value_m_list.append(fit_value_m_i)
            fit_value_h_list.append(fit_value_h_i)
            labels.append("{0}, fixgamma_{1}, subg_{2}, {3}".format(PSFtyp, fixgamma, subg, fit_type))
            chisq.append(repr(round(return_chisq(filename, lens_mask= lens_mask, fair_mask=True, fair_PSF = fair_PSF),3)))   #Fair_mask = True means same rms

labels_orig = copy.deepcopy(labels)

fit_value_l_list = [x for _,x in sorted(zip(labels,fit_value_l_list))]  #!!!
fit_value_m_list = [x for _,x in sorted(zip(labels,fit_value_m_list))]
fit_value_h_list = [x for _,x in sorted(zip(labels,fit_value_h_list))]
chisq = [x for _,x in sorted(zip(labels,chisq))]
labels = [x for _,x in sorted(zip(labels,labels))]  # Sort the label at the last

#%%
        
for i in range(len(pick_names)):
    print i, ':', pick_names[i]
print pick_names
pick = input('select the names to plot?:\n')

fit_value_l = [fit_value_l_list[i][pick] for i in range(len(labels))]
fit_value_m = [fit_value_m_list[i][pick] for i in range(len(labels))]
fit_value_h = [fit_value_h_list[i][pick] for i in range(len(labels))]

bars = [labels[i].split(',')[0] for i in range(len(labels))]
bars = sorted(set(bars), key=bars.index)

lines_num = 3 #gamma grids
index0 = [i*3 for i in range(len(bars))]
label0 = [labels[i].split(',')[0] for i in index0]
index1 = [i*3+2 for i in range(len(bars))]
label1 = [labels[i].split(',')[0] for i in index1]

# #count_n =len(labels) * 10 /100
# chisq = [float(chisq[i]) for i in range(len(chisq))]
# sort_chisq = np.argsort(np.asarray(chisq))    
# Chisq_best = chisq[sort_chisq[0]]
# Chisq_last= chisq[sort_chisq[-1]]
# weight = np.zeros(len(chisq))
# for i in range(len(sort_chisq)):
#     weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best))
# weighted_value = np.sum(np.array(fit_value_m)*weight) / np.sum(weight)
# rms_value = np.sqrt(np.sum((np.array(fit_value_m)-weighted_value)**2*weight) / np.sum(weight))
# print "weighted_value, rms_value:", weighted_value, rms_value

#Plot it out    
x_pos = np.arange(len(bars)) + 0.5
fmt =['-']
color = ['green', 'blue', 'red']
gamma_grid = [1.9, 2.0, 2.1]
texts = []
fig = plt.figure(figsize=(12, 10))

for ftype in range(lines_num):
    index = [i for i in range(len(labels)) if 'fixgamma_{0}'.format(gamma_grid[ftype]) in labels[i]]
    ### index = [index_all[(i*lines_num+ftype)] for i in range(len(bars))]
    for i in range(len(index)):
        if labels[index[i]][6:] != labels[index[0]][6:]:
            raise ValueError("The labels is wrong for some reason")
    value_result = [fit_value_m[i] for i in index]
    asm_error = [[fit_value_m[i]-fit_value_l[i] for i in index],
                        [fit_value_h[i]-fit_value_m[i] for i in index]]
    plt.errorbar(x_pos+0.01*ftype, value_result, yerr=asm_error, color=color[ftype],
                  label="{0}, {1}".format(labels[ftype].split(',')[1],labels[ftype].split(',')[3]))
    for i in range(len(index)):
        texts.append(plt.text((x_pos+0.01*ftype)[i], value_result[i]*9.2/9, chisq[index[i]],fontsize=15, color=color[ftype]) )
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))            
#     xs = np.linspace(x_pos[0], x_pos[-1])
# #    count_n = 8 #len(labels) * 10 /100
# #    chisq = [float(chisq[i]) for i in range(len(chisq))]
# #    sort_chisq = np.argsort(np.asarray(chisq))    
# #    Chisq_best = chisq[sort_chisq[0]]
# #    Chisq_last= chisq[sort_chisq[count_n-1]]
# #    inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
# #    weight = np.zeros(len(chisq))
# #    for i in sort_chisq[:count_n]:
# #        weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
# #    weighted_value = np.sum(np.array(fit_value_m)*weight) / np.sum(weight)
# #    rms_value = np.sqrt(np.sum((np.array(fit_value_m)-weighted_value)**2*weight) / np.sum(weight))
# #    print weighted_value, rms_value
# #    plt.plot(xs, xs*0+(weighted_value - rms_value), xs, xs*0+weighted_value, xs, xs*0+(weighted_value + rms_value),color = 'red')
# #    plt.fill_between(xs, weighted_value - rms_value, weighted_value + rms_value, facecolor='red', alpha = 0.1)    
plt.legend(numpoints=1,ncol=2,loc=2,prop={'size':18})
plt.title('The fitting result for {0}, filter {1}'.format(ID[2:], filt),fontsize=25)
plt.ylabel(pick_names[pick],fontsize=25)
plt.ylim(np.min(fit_value_l)/1.3,np.max(fit_value_m)*1.3)
plt.tick_params(labelsize=15)
plt.xticks(x_pos, bars)
plt.savefig('fig_result_pick{0}_oneconfig.pdf'.format(pick,fixgamma))
plt.show()
    
