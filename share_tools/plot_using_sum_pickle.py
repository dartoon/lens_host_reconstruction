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
import copy
from adjustText import adjust_text   # avoid the overlapping while ploting
from read_chisq import return_chisq 

IDs = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']
folder_type = '2nd_fit_'
# IDs = ['3_SDSS1206']
# folder_type = 'singel_fit_'
#IDs = ['6_HS2209']
#folder_type = '3rd_fit_'
#%%
for i in range(len(IDs)):
    print(i, ':', IDs[i])
inp = int(input("Which source:\n"))
ID = IDs[inp][2:]
filename = ID+'_' +folder_type[:4] + 'run_summary.pkl'
# filename = ID+'_singel_host_run_summary.pkl'
result = pickle.load(open(filename,'rb') ,encoding="latin1") 
fit_values, chisq, labels, pick_names = result
fit_value_l_list, fit_value_m_list, fit_value_h_list = fit_values

sort_label_type=['subg_2, PSFrecons', ', subg_2, QSOmask', 'subg_3, PSFrecons',', subg_3, QSOmask']
#if labels[0][-17:]!='subg_2, PSFrecons' or labels[1][-15:]!='subg_2, QSOmask' or labels[0][-17:]!='subg_3, PSFrecons' or labels[1][-15:]!='subg_3, QSOmask'
for i in range(len(labels)):
    if labels[i][-17:] != sort_label_type[i%4]:
        raise ValueError("The labels is wrong for some reason")
        
for i in range(len(pick_names)):
    print(i, ':', pick_names[i])
print(pick_names)
pick = int(input('select the names to plot?:\n'))

fit_value_l = [fit_value_l_list[i][pick] for i in range(len(labels))]
fit_value_m = [fit_value_m_list[i][pick] for i in range(len(labels))]
fit_value_h = [fit_value_h_list[i][pick] for i in range(len(labels))]

zp = 24. #!!!
if_trans_mag = False
if if_trans_mag and pick_names[pick]=='Host flux souce plane':
    fit_value_l = -2.5*np.log10(fit_value_h)+zp
    fit_value_m = -2.5*np.log10(fit_value_m)+zp
    fit_value_h = -2.5*np.log10(fit_value_l)+zp

bars = [labels[i].split(',')[0] for i in range(len(labels))]
bars = sorted(set(bars), key=bars.index)

index0 = [i*12 for i in range(len(bars))]
label0 = [labels[i].split(',')[0] for i in index0]
index1 = [i*12+11 for i in range(len(bars))]
label1 = [labels[i].split(',')[0] for i in index1]
if label0 != label1 or bars!= label0:
    print(label0, label1)
    raise ValueError("The labels is wrong for some reason")
color = ['green', 'orange']

count_n = 8#int(len(labels) * 20 /100)
chisq = [float(chisq[i]) for i in range(len(chisq))]
sort_chisq = np.argsort(np.asarray(chisq))    
Chisq_best = chisq[sort_chisq[0]]
#Chisq_last= chisq[sort_chisq[-1]]
Chisq_last= chisq[sort_chisq[count_n-1]]
weight = np.zeros(len(chisq))
inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
for i in sort_chisq[:count_n]:
#for i in range(len(sort_chisq)):
#    weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best))
    weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
weighted_value = np.sum(np.array(fit_value_m)*weight) / np.sum(weight)
rms_value = np.sqrt(np.sum((np.array(fit_value_m)-weighted_value)**2*weight) / np.sum(weight))
print("weighted_value, rms_value:", weighted_value, rms_value, 'used sets:', count_n)

def plt_result(fixgamma, chisq=chisq):
    #Plot it out    
    fig = plt.figure(figsize=(12, 10))
    x_pos = np.arange(len(bars)) + 0.5
    fmt =['--','--','-','-']
    lines_num = 4 #subg 2-3, PSFrecon and QSO mask
    texts = []
    for ftype in range(lines_num):
        index_all = [i for i in range(len(labels)) if 'fixgamma_{0}'.format(fixgamma) in labels[i]]
        index = [index_all[(i*lines_num+ftype)] for i in range(len(bars))]
        for i in range(len(index)):
            if labels[index[i]][6:] != labels[index[0]][6:]:
                raise ValueError("The labels is wrong for some reason")
#        print index
#        label_type = [labels[i] for i in index]
        value_result = [fit_value_m[i] for i in index]
        asm_error = [[fit_value_m[i]-fit_value_l[i] for i in index],
                            [fit_value_h[i]-fit_value_m[i] for i in index]]
        plt.errorbar(x_pos+0.01*ftype, value_result, yerr=asm_error, fmt =fmt[ftype], color=color[ftype%2],
                     label="{0}, {1}".format(labels[ftype].split(',')[2],labels[ftype].split(',')[3]))
        for i in range(len(index)):
            texts.append(plt.text((x_pos+0.01*ftype)[i], value_result[i], chisq[index[i]],fontsize=15) )
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))            
    
    xs = np.linspace(x_pos[0], x_pos[-1])
    ys_l = xs*0 + weighted_value-rms_value
    ys_m = xs*0 + weighted_value
    ys_h = xs*0 + weighted_value+rms_value
    plt.plot(xs,ys_l, xs, ys_m, xs,ys_h,color = 'red')
    plt.fill_between(xs, ys_l, ys_h, facecolor='red', alpha = 0.06)

#    ##If want to put horizontal line:
#    fill_ref = False
#    if pick == 3 or pick == 3-len(pick_names):  #
#        ref_value_l,ref_value, ref_value_h = 42.317, 47.7 , 53.76
#        fill_ref = True
#    elif pick == 4 or pick == 4-len(pick_names):
#        ref_value_l,ref_value, ref_value_h = 3.94 -0.14, 3.94 ,3.94 + 0.14
#        fill_ref = True
#    elif pick == 5 or pick == 5-len(pick_names):
#        ref_value_l,ref_value, ref_value_h = 0.82 - 0.14, 0.82, 0.82 + 0.14
#        fill_ref = True
#    ys_l = xs*0 + ref_value_l
#    ys_m = xs*0 + ref_value
#    ys_h = xs*0 + ref_value_h
#    if if_trans_mag and pick_names[pick]=='Host flux souce plane':
#        ys_l = -2.5*np.log10(ys_h)+zp
#        ys_m = -2.5*np.log10(ys_m)+zp
#        ys_h = -2.5*np.log10(ys_l)+zp
#    if fill_ref ==True:
#        plt.plot(xs,ys_l, xs, ys_m, xs,ys_h,color = 'gray')
#        plt.fill_between(xs, ys_l, ys_h, facecolor='gray', alpha = 0.2)
    plt.legend(numpoints=1,ncol=2,loc=2,prop={'size':18})
    plt.title('The fitting result for {0}, fix gamma = {1}'.format(ID[2:], fixgamma),fontsize=25)
    plt.ylabel(pick_names[pick],fontsize=25)
    
    plt.ylim(np.min(fit_value_l)/1.3,np.max(fit_value_m)*1.3)
    if if_trans_mag and pick_names[pick]=='Host flux souce plane':
        plt.ylim(np.min(fit_value_l)/1.05,np.max(fit_value_m)*1.05)
        plt.ylabel('Magnitude',fontsize=25)
    plt.tick_params(labelsize=15)
    plt.xticks(x_pos, bars)
#    plt.savefig('fig_result_pick{0}_fixgamma{1}.pdf'.format(pick,fixgamma))
#    plt.close()
    
    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    return image

import imageio
kwargs_write = {'fps':1.0, 'quantizer':'nq'}
imageio.mimsave('fig_result_{0}_PSF.gif'.format(pick_names[pick]), [plt_result(fix_gamma) for fix_gamma in  ['1.9', '2.0', '2.1']], fps=1)

#

