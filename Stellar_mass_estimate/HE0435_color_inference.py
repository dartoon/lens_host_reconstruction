import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
import copy
import sys
sys.path.insert(0,'../share_tools/')
from read_chisq import return_chisq 

weighted_value_list, rms_value_list = [], []

ID_name = '0_HE0435'
ID = 'HE0435'
folder_list = ['model/2nd_fit_PSFi_*','f555w_model/fit_PSFi_PSFrecons/', 'f814w_model/fit_PSFi_PSFrecons/']

#%%
labels = []
fit_value_l_list, fit_value_m_list, fit_value_h_list, = [], [], []
chisq = []
fair_PSF = True
for i in range(len(folder_list)):
	filenames = glob.glob('../'+ID_name+'/'+folder_list[i]+'/result_PSF?_*subg?.pkl')   
	filenames.sort()
	chisq = []
	fit_value_l, fit_value_m, fit_value_h = [], [], []
	for filename in filenames:
		fit_type = filename.split('/')[0].split('_')[-1]
# 		print(filename)
		result = pickle.load(open(filename,'rb') ,encoding="latin1") 
		fit_result, trans_result, kwargs_material, model_lists  = result
		pick_names = trans_result[-1]
		pick = 1
		fit_value = np.asarray(trans_result[0])[:,pick]
		fit_value_l.append(np.percentile(fit_value,16,axis=0))
		fit_value_m.append(np.percentile(fit_value,50,axis=0))
		fit_value_h.append(np.percentile(fit_value,84,axis=0))
		chisq.append(repr(round(return_chisq(filename, lens_mask= None, fair_mask=True, fair_PSF = fair_PSF),3)))   #Fair_mask = True means same rms
	count_n = len(chisq)
	if count_n>10:
		count_n = 8
	chisq = [float(chisq[i]) for i in range(len(chisq))]
	sort_chisq = np.argsort(np.asarray(chisq))
	Chisq_best = chisq[sort_chisq[0]]
 	#Chisq_last= chisq[sort_chisq[-1]]
	Chisq_last= chisq[sort_chisq[count_n-1]]
	weight = np.zeros(len(chisq))
	inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
	
	for i in sort_chisq[:count_n]:
		weight[i] = np.exp(-1/2. * (chisq[i]-Chisq_best)/(Chisq_best* inf_alp))
	weighted_value = np.sum(np.array(fit_value_m)*weight) / np.sum(weight)
	rms_value = np.sqrt(np.sum((np.array(fit_value_m)-weighted_value)**2*weight) / np.sum(weight))
	weighted_value_list.append(weighted_value)
	rms_value_list.append(rms_value)

filter_list = ['F160W', 'ACS_F555W', 'ACS_F814W']
zeropoint_dict = { 'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                  'F160W':25.9463,'ACS_F814W': 25.957, 'ACS_F555W': 25.735}

print(ID+" fnu inference in image plane:")
for i in range(len(filter_list)):
	uncer =  rms_value_list[i]/weighted_value_list[i]
	zp = zeropoint_dict[filter_list[i]]
	mag = -2.5*np.log10(weighted_value_list[i]) + zp
	fnu = 10 ** ((mag-25)/(-2.5))
#	print(filter_list[i], round(mag,3)) 
	print(filter_list[i], round(fnu,3), '+-', round(fnu*uncer,3))
	
#%%After runing the fit:
#The inferred stellar mass:
M_star_boost = 12.04
mag = 1100.7300575792065/70.10881988478579
M_star = np.log10(10**12.04/mag)
print(M_star)