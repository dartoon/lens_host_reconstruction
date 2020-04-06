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

def read_inf(ID, prop=0 , count_n=[4, 4], rod = 3):
    '''
    Return the photometry of the infernece from the summarized pickle file.
    
    Parameter
    --------
        count_n: if a singel number, combine QSO_mask and PSF_recon. If list, combining them separately.
    Return
    --------
        A sth sth
    '''
    if ID == "SDSS1206" :
        filename = ID+'_singel_host_run_summary.pkl'
    elif ID == 'HS2209' or ID == 'RXJ1131':
        filename = ID+'_3rd_run_summary.pkl'
    else:
        filename = ID+'_2nd_run_summary.pkl'
    result = pickle.load(open('/Users/Dartoon/Astro/Projects/lens_host_reconstruction/share_tools/'+filename,'rb'),encoding="latin1")
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
    elif prop == 4 :  #Flux image plane
        if 'Disk' in pick_names[2]:
            picks = [1,2]
        else:
            picks = [1]            
#    print(ID, pick_names)
    result = []
#    print(filename)
    for pick in picks:
        if isinstance(count_n, list) == False:
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
        elif isinstance(count_n, list):
            _key = ['QSOmask', 'PSFrecons']
            weight = np.zeros(len(labels))
            for ii in range(2):
                count = count_n[ii]
                fit_value_m = [fit_value_m_list[i][pick] for i in range(len(labels))]
                chisq_ig = [float(chisq[i]) for i in range(len(chisq))]
                for i in range(len(chisq)):
                    if _key[ii] in labels[i]:
                        chisq_ig[i] = 10**4
                sort_chisq = np.argsort(np.asarray(chisq_ig))    
                Chisq_best = chisq_ig[sort_chisq[0]]
                Chisq_last= chisq_ig[sort_chisq[count-1]]
                inf_alp = (Chisq_last-Chisq_best) / (2*2.* Chisq_best)
                for i in sort_chisq[:count]:
                    weight[i] = np.exp(-1/2. * (chisq_ig[i]-Chisq_best)/(Chisq_best* inf_alp))
#            print(weight)
            weighted_value = np.sum(np.array(fit_value_m)*weight) / np.sum(weight)
            rms_value = np.sqrt(np.sum((np.array(fit_value_m)-weighted_value)**2*weight) / np.sum(weight))
            result.append([round(weighted_value,rod), round(rms_value,rod)])            
    return result
#result = read_inf('HE0435', prop = 0)
#print(read_inf('RXJ1131', prop = 4, count_n=[4, 4]))
#print(read_inf('HS2209', prop = 1, count_n=8))
#print(read_inf('HS2209',  count_n=8))    

filter_dict = {'HE0435': 'F160W', 'RXJ1131': 'ACS_F814W', 'WFI2033': 'F160W', 'SDSS1206': 'F160W', 
               'HE1104': 'F160W', 'SDSS0246': 'F814W', 'HS2209': 'F814W', 'HE0047': 'F814W'}

zeropoint_dict = {'F814W' :25.110, 'F555W': 25.807,'F125W': 26.2303,'F140W': 26.4524, 
                  'F160W':25.9463,'ACS_F814W': 25.957}

def read_fnu(ID, count_n=[4, 4]):
    result = []
    if ID != "RXJ1131":
        _result = read_inf(ID, prop=0 , count_n=count_n, rod = 5)
        weighted_value, rms_value = _result[0][0], _result[0][1]
        uncer = rms_value/weighted_value
        zp = zeropoint_dict[filter_dict[ID]]
        mag = -2.5*np.log10(weighted_value) + zp
        fnu = 10 ** ((mag-25)/(-2.5))
        result.append((fnu, fnu*uncer))
    else: #For RXJ1131, return bulge then disk.
        _result = read_inf(ID, prop=0 , count_n=count_n, rod = 5)
        for j in range(2):
            weighted_value, rms_value = _result[j][0], _result[j][1]
            uncer = rms_value/weighted_value
            zp = zeropoint_dict[filter_dict[ID]]
            mag = -2.5*np.log10(weighted_value) + zp
            fnu = 10 ** ((mag-25)/(-2.5))
            result.append((fnu, fnu*uncer))
    return result
#result = read_fnu('HE0435',  count_n=[4, 4])
#print(read_fnu('HE0435',  count_n=[4, 4]))



def read_mstar(ID, count_n=[4, 4], rod = 5):
    result = []
    org_fnu = {'HE0435' :  [29.32598834079476],
        'RXJ1131' :  [14.808995717960412, 189.5940277253136],
        'WFI2033' :  [27.447888815355586],
        'SDSS1206' :  [30.566894649185073],
        'HE1104' :  [30.933900240069743],
        'SDSS0246' :  [3.7138818625670735],
        'HS2209' :  [70.16711652527341],
        'HE0047' :  [6.476528070125638]}
    #Correct based on new :
    new_fnw = read_fnu(ID, count_n=count_n)
    if ID != "RXJ1131": 
        #Read org mstar
        sed_name = '*' + ID
        sed_filename  = '../Stellar_mass_estimate/'+ sed_name + '/SFH_*_PA00_param.fits'
        sed_filename = glob.glob(sed_filename)[0]
        hdul = pyfits.open(sed_filename)
        table = hdul[1].data    
        tab_name = hdul[1].columns
        mstar_idx = [ii for ii in range(len(tab_name)) if 'Mstel' in str(tab_name[ii])][0]
        m_star_org =  10**table[1][mstar_idx]        
#        print(m_star_org)
        m_star = [m_star_org/org_fnu[ID][0] * new_fnw[0][0], m_star_org/ org_fnu[ID][0] * new_fnw[0][1]]
        result.append(m_star)
    else:
        sed_name_l = ['1_RXJ1131_bulge', '1_RXJ1131_disk']
        for j in range(2):
            sed_name = sed_name_l[j]
            sed_filename  = '../Stellar_mass_estimate/'+ sed_name + '/SFH_*_PA00_param.fits'
            sed_filename = glob.glob(sed_filename)[0]
            hdul = pyfits.open(sed_filename)
            table = hdul[1].data    
            tab_name = hdul[1].columns
            mstar_idx = [ii for ii in range(len(tab_name)) if 'Mstel' in str(tab_name[ii])][0]
            m_star_org =  10**table[1][mstar_idx]        
            m_star = [m_star_org/ org_fnu[ID][j] * new_fnw[j][0], m_star_org/ org_fnu[ID][j] * new_fnw[j][0] * (new_fnw[j][1]/new_fnw[j][0])]
            result.append(m_star)
    return result
#result = read_mstar('HE0435',  count_n=8)
#result = read_mstar('RXJ1131',  count_n=8)