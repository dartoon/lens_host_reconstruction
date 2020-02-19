import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
import copy
IDs = ['0_HE0435', '1_RXJ1131', '1_RXJ1131_disk', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']

weighted_value_list, rms_value_list = [], []

prop = 0

#%%
for i in range(len(IDs)):
	ID = IDs[i][2:]
	if ID == "SDSS1206" :
		filename = ID+'_singel_host_run_summary.pkl'
	else:
		filename = ID+'_2nd_run_summary.pkl'
	disk_note = 0
	if ID == "RXJ1131_disk" :
		ID = 'RXJ1131'
		filename = ID+'_2nd_run_summary.pkl'
		disk_note = 1
	# filename = ID+'_2nd_run_summary.pkl'	
	result = pickle.load(open("../share_tools/"+filename,'rb'),encoding="latin1") 
	fit_values, chisq, labels, pick_names = result
	sort_label_type=['subg_2, PSFrecons', ', subg_2, QSOmask', 'subg_3, PSFrecons',', subg_3, QSOmask']
	fit_value_l_list, fit_value_m_list, fit_value_h_list = fit_values
	for i in range(len(labels)):
		if labels[i][-17:] != sort_label_type[i%4]:
			raise ValueError("The labels is wrong for some reason")
	if 'Disk' in pick_names[2]:
		pick = 4
	else:
		pick = 3
	if disk_note == 1:
		pick = 5
	# print(ID, pick_names)
	fit_value_l = [fit_value_l_list[i][pick] for i in range(len(labels))]
	fit_value_m = [fit_value_m_list[i][pick] for i in range(len(labels))]
	fit_value_h = [fit_value_h_list[i][pick] for i in range(len(labels))]
	count_n = 8
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
	weighted_value_list.append(weighted_value)
	rms_value_list.append(rms_value)
	# print pick_names[pick]
	# print "weighted_value, rms_value:", round(weighted_value,2), round(rms_value,2)

fig = plt.figure(figsize=(12, 10))
x_pos = np.arange(len(IDs))
for i in range(len(IDs)):
 	plt.errorbar(x_pos[i], weighted_value_list[i], yerr=rms_value_list[i], fmt = '*')#, color=color[ftype%2],
    # plt.legend(numpoints=1,ncol=2,loc=2,prop={'size':18})
plt.title('The overall fitting results',fontsize=25)
plt.ylabel(pick_names[pick],fontsize=25)
plt.tick_params(labelsize=15)
plt.xticks(x_pos, [IDs[i][2:] for i in range(len(IDs))])
plt.close()

#%%Translate flux into f_nu 

filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
               'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}

zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                  'F160W':25.9463,'ACS_F814W': 25.957}

fnu_list =  []
for i in range(len(IDs)):
	ID = IDs[i][2:]
	if ID == "RXJ1131_disk" :
		ID = 'RXJ1131'	
		# print(filter_dict[ID])
	uncer = rms_value_list[i]/weighted_value_list[i]
	zp = zeropoint_dict[filter_dict[ID]]
	mag = -2.5*np.log10(weighted_value_list[i]) + zp
	fnu = 10 ** ((mag-25)/(-2.5))
	print(IDs[i][2:],":\t",filter_dict[ID], "\t",round(fnu,3), '+-',round(fnu*uncer,3))
    # fnu_list.append(fnu)
    
 