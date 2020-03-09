#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 15:58:15 2019

@author: Dartoon

Summarize the fitting result

Based on python 2
"""
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import pickle
import glob
import copy
#from adjustText import adjust_text   # avoid the overlapping while ploting
from read_chisq import return_chisq 

#IDs = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']
#folder_type = '2nd_fit_'
# IDs = ['3_SDSS1206']
# folder_type = 'singel_fit_'
IDs = ['1_RXJ1131']
folder_type = '3rd_fit_'
for i in range(len(IDs)):
    ID = IDs[i][2:]
    root = '../'+ IDs[i]+ '/model/'
    pkl_folderfiles = glob.glob(root+folder_type+"PSFi_PSFrecons/result_PSF*.pkl")
    pkl_files = [pkl_folderfiles[i].split('/')[-1] for i in range(len(pkl_folderfiles))]
    psf_num = np.max([int(pkl_files[i].split('_PSF')[1]) for i in range(len(pkl_files))]) + 1
   
    if float(len(pkl_files)/6) != len(pkl_files)/6. :
        raise ValueError("The number of the pickle result is the Multiples of six")
   
    PSFno_ = range(psf_num)
    #The values to pick up:
    filename = pkl_folderfiles[0]
    result = pickle.load(open(filename,'rb'))
    fit_result, trans_result, kwargs_material, model_lists  = result
    pick_names = trans_result[-1]
    _, _, _, _, lens_mask = kwargs_material
    stdd_file = glob.glob(root+'fit_PSFi_QSOmask/*_stdd_boost.fits')
    if len(stdd_file)>1:
        raise ValueError("More than one *_stdd_boost.fits exits") 
    rms_boost = pyfits.getdata(stdd_file[0])
    lens_mask[rms_boost == rms_boost.max()] = 0
    plt.imshow(lens_mask,origin='low')
    plt.show()
   
    fair_PSF = True
    labels = []
    fit_value_l_list, fit_value_m_list, fit_value_h_list, = [], [], []
    chisq = []
    fixgamma_list = []
   
    folder_i = [folder_type+'PSFi_QSOmask', folder_type+'PSFi_PSFrecons']
    for folder in folder_i:
        filenames = glob.glob(root+'{0}/result_PSF?_*'.format(folder))        
        filenames.sort()
        for filename in filenames:
            PSFtyp = filename.split('/')[-1].split('_')[1]
            if int(PSFtyp.split('PSF')[1]) in PSFno_:
                subg = filename[(filename).find('subg')+len('subg')]
                fit_type = filename.split('/')[-2].split('_')[-1]
                fixgamma = (filename.split('gammafix')[1]).split('_subg')[0]
                if fixgamma not in fixgamma_list:
                    fixgamma_list.append(fixgamma)
                print(filename)
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
                chisq.append(repr(round(return_chisq(filename, lens_mask= lens_mask, 
                                                     fair_mask=True, fair_PSF = fair_PSF, py_version = 2),3)))   
    labels_orig = copy.deepcopy(labels)
    fit_value_l_list = [x for _,x in sorted(zip(labels,fit_value_l_list))]  #!!!
    fit_value_m_list = [x for _,x in sorted(zip(labels,fit_value_m_list))]
    fit_value_h_list = [x for _,x in sorted(zip(labels,fit_value_h_list))]
    chisq = [x for _,x in sorted(zip(labels,chisq))]
    labels = [x for _,x in sorted(zip(labels,labels))]  # Sort the label at the last
    picklename = ID+'_' +folder_type[:4] + 'run_summary.pkl'
    to_be_saved = [[fit_value_l_list, fit_value_m_list, fit_value_h_list],chisq, labels, pick_names]
    pickle.dump(to_be_saved, open(picklename, 'wb'))

