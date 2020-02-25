#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 11:44:45 2020

@author: Dartoon
"""

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
import copy

def read_inf(ID, prop=0 , count_n=8, rod = 3):
	if ID == "SDSS1206" :
		filename = ID+'_singel_host_run_summary.pkl'
	else:
		filename = ID+'_2nd_run_summary.pkl'
	result = pickle.load(open('../share_tools/'+filename,'rb'),encoding="latin1")
	fit_values, chisq, labels, pick_names = result
	sort_label_type=['subg_2, PSFrecons', ', subg_2, QSOmask', 'subg_3, PSFrecons',', subg_3, QSOmask']
	fit_value_l_list, fit_value_m_list, fit_value_h_list = fit_values
	for i in range(len(labels)):
		if labels[i][-17:] != sort_label_type[i%4]:
			raise ValueError("The labels is wrong for some reason")
	if prop == 0 :  #Flux
		if 'Disk' in pick_names[2]:
			picks = [4,5]
		else:
			picks = [3]
	elif prop == 1: #Reff
		if 'Disk' in pick_names[2]:
			picks = [6,7]
		else:
			picks = [5]
	elif prop == 2: #Sersic n
		if 'Disk' in pick_names[2]:
			picks = []
		else:
			picks = [4]
	elif prop == 3: #AGN flux
		if 'Disk' in pick_names[2]:
			picks = [3]
		else:
			picks = [2]
			
#	print(ID, pick_names)
	result = []
	for pick in picks:
		fit_value_m = [fit_value_m_list[i][pick] for i in range(len(labels))]
		count_n = count_n
		
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
		result.append([round(weighted_value,rod), round(rms_value,rod)])
	return result

#print(read_inf('HE0435', prop = 2))
#print(read_inf('RXJ1131', prop = 2))