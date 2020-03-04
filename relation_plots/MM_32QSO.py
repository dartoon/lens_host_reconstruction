#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:40:00 2020

@author: Dartoon
"""

import numpy as np
np.set_printoptions(precision=4)
from matplotlib import colors
import matplotlib as mat
mat.rcParams['font.family'] = 'STIXGeneral'

import matplotlib as mpl
mpl.rc('image', cmap='jet')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
############### with evolution-corrected or not? #################
import sys
#from dmag import pass_dmag
from adjustText import adjust_text   # avoid the overlapping while ploting

plt.figure(figsize=(11.5,12))
########input 25 local by Bennert++2011 ############
bloc = np.array([[ 0.054 , 10.12  ,  0.24  ,  7.436 ,  0.4   ],
       [ 0.043 , 10.95  ,  0.23  ,  8.006 ,  0.4   ],
       [ 0.076 , 10.33  ,  0.22  ,  7.436 ,  0.4   ],
       [ 0.041 , 10.38  ,  0.23  ,  7.166 ,  0.4   ],
       [ 0.051 , 10.5   ,  0.23  ,  8.516 ,  0.4   ],
       [ 0.0524, 10.32  ,  0.23  ,  6.976 ,  0.4   ],
       [ 0.0475,  9.83  ,  0.24  ,  7.666 ,  0.4   ],
       [ 0.055 , 10.41  ,  0.23  ,  7.836 ,  0.4   ],
       [ 0.0355, 10.33  ,  0.22  ,  7.766 ,  0.4   ],
       [ 0.021 , 10.2   ,  0.22  ,  7.366 ,  0.4   ],
       [ 0.038 , 10.26  ,  0.22  ,  7.356 ,  0.4   ],
       [ 0.0229,  9.94  ,  0.24  ,  7.496 ,  0.4   ],
       [ 0.047 , 10.14  ,  0.22  ,  7.646 ,  0.4   ],
       [ 0.0559,  9.65  ,  0.22  ,  6.976 ,  0.4   ],
       [ 0.0501, 10.11  ,  0.23  ,  7.956 ,  0.4   ],
       [ 0.0541, 10.04  ,  0.23  ,  7.076 ,  0.4   ],
       [ 0.0558, 10.73  ,  0.21  ,  7.506 ,  0.4   ],
       [ 0.0365, 10.3   ,  0.24  ,  7.476 ,  0.4   ],
       [ 0.0304, 10.24  ,  0.24  ,  7.826 ,  0.4   ],
       [ 0.0481,  9.92  ,  0.22  ,  7.276 ,  0.4   ],
       [ 0.0483, 10.    ,  0.23  ,  7.626 ,  0.4   ],
       [ 0.0465,  9.82  ,  0.23  ,  7.546 ,  0.4   ],
       [ 0.0532,  9.95  ,  0.23  ,  7.776 ,  0.4   ],
       [ 0.0585, 10.33  ,  0.22  ,  7.216 ,  0.4   ],
       [ 0.0409, 10.33  ,  0.22  ,  7.34  ,  0.4   ]])
plt.errorbar(bloc[:,1],bloc[:,3], xerr=bloc[:,2] ,yerr=bloc[:,4],fmt='.',color='gray',markersize=15)
#plt.plot(bloc[:,1],bloc[:,3], '.',color='gray',markersize=15)
Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)

########input 30 local by Haring 04 ############
hloc = np.array([[ 0.        , 11.77815125,  0.18      ,  9.47712125,  0.21298437],
       [ 0.        , 10.36172784,  0.18      ,  7.14612804,  0.14261786],
       [ 0.        , 10.83250891,  0.18      ,  8.        ,  0.10205999],
       [ 0.        , 11.5563025 ,  0.18      ,  8.63346846,  0.21305862],
       [ 0.        , 11.5563025 ,  0.18      ,  8.71600334,  0.32388976],
       [ 0.        , 11.74818803,  0.18      ,  8.72427587,  0.31358181],
       [ 0.        , 11.462398  ,  0.18      ,  8.51851394,  0.20838037],
       [ 0.        ,  9.79239169,  0.18      ,  7.14612804,  0.18220849],
       [ 0.        , 11.11394335,  0.18      ,  7.56820172,  0.16526173],
       [ 0.        , 11.462398  ,  0.18      ,  9.39794001,  0.39967027],
       [ 0.        , 10.56820172,  0.18      ,  7.65321251,  0.18655821],
       [ 0.        ,  8.90308999,  0.18      ,  6.39794001,  0.34062062],
       [ 0.        , 10.83884909,  0.18      ,  7.64345268,  0.469419  ],
       [ 0.        , 10.88081359,  0.18      ,  7.14612804,  0.22796598],
       [ 0.        , 11.07918125,  0.18      ,  9.        ,  0.23856063],
       [ 0.        , 10.83250891,  0.18      ,  8.32221929,  0.2835412 ],
       [ 0.        , 10.20411998,  0.18      ,  8.        ,  0.61649806],
       [ 0.        , 10.30103   ,  0.18      ,  7.20411998,  0.43571349],
       [ 0.        , 10.98677173,  0.18      ,  8.2787536 ,  0.25971825],
       [ 0.        , 11.11394335,  0.18      ,  8.49136169,  0.20586415],
       [ 0.        , 10.07918125,  0.18      ,  8.47712125,  0.2480033 ],
       [ 0.        , 10.96378783,  0.18      ,  8.04139269,  0.56734929],
       [ 0.        , 10.64345268,  0.18      ,  7.74818803,  0.43388101],
       [ 0.        , 11.43136376,  0.18      ,  9.        ,  0.33450339],
       [ 0.        , 11.69019608,  0.18      ,  9.30103   ,  0.04845501],
       [ 0.        , 11.04139269,  0.18      ,  8.23044892,  0.62621233],
       [ 0.        , 10.56820172,  0.18      ,  8.38021124,  0.10653741],
       [ 0.        , 10.17609126,  0.18      ,  7.11394335,  0.18446512],
       [ 0.        ,  9.84509804,  0.18      ,  6.54406804,  0.14739052],
       [ 0.        , 10.04139269,  0.18      ,  6.56820172,  0.15706652]])
# +  np.log10(Haring04[:,4] * 10 ** Haring04[:,5]))/2 #Sigma LgMBH
plt.errorbar(hloc[:,1],hloc[:,3], xerr=hloc[:,2] ,yerr=hloc[:,4],fmt='.',color='black',markersize=15)
#plt.plot(hloc[:,1],hloc[:,3], '.',color='black',markersize=15)
Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)

#############################################################
###################fitting all togethera with MCMC#########################
x=np.append(bloc[:,1], hloc[:,1])
y=np.append(bloc[:,3], hloc[:,3])
yerr=(np.append(bloc[:,2], hloc[:,2])**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # 0.2 is the uncertainty level for the L_R
def lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
#    yerr=(m*(np.append(bloc[:,2], hloc[:,2]))**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # right error level
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [0.93027905, -1.95536508, 0.35], args=(x, y, yerr))
m_ml, b_ml,sint_ml= result["x"]

xp = np.array([5, 13])
#plt.plot(xp, m_ml*xp+b_ml, 'r-')
def lnprior(theta):
    m, b, sint	 = theta
    if -5.0 < m < 5 and -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
ndim, nwalkers = 3, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 1000)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
xl = np.linspace(5, 13, 100)
m, b, sint =np.percentile(samples, 50,axis=0)
plt.plot(xl, m_ml*xl+b_ml, color="k", linewidth=4.0,zorder=-0.5)

def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]
m=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
for i in range(100):
    posi=np.random.uniform(16,84)
    m=np.percentile(samples,posi,axis=0)[0]
    b=samples[:,1][samples[:,0]==find_n(samples[:,0],m)][0]   #may find out many numbers
    plt.plot(xl, m*xl+b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1000)

plt.text(9.3, 6.24, r"log(M$_{\rm BH}/$10$^{7}$M$_{\odot}$)=%s+%slog(M$_*/$10$^{10}$M$_{\odot}$)"%(round(b_ml+m_ml*10-7,2),round(m_ml,2)),color='blue',fontsize=25)

#%%
inp_SS13 = 1
ss = np.array([[ 0.717, 10.66 ,  6.936,  1.   ],
       [ 1.065,  9.8  ,  6.836,  0.5  ],
       [ 0.96 , 10.49 ,  7.116,  0.23 ],
       [ 0.97 , 10.29 ,  7.986,  1.   ],
       [ 0.544, 11.01 ,  8.376,  1.   ],
       [ 1.044, 10.33 ,  7.656,  0.35 ],
       [ 0.675, 10.96 ,  7.596,  1.   ],
       [ 0.569, 10.55 ,  7.636,  0.4  ],
       [ 0.737, 10.66 ,  8.896,  0.32 ],
       [ 0.664, 10.19 ,  6.806,  0.41 ],
       [ 0.837, 10.39 ,  8.106,  0.33 ],
       [ 0.74 , 10.76 ,  7.776,  0.24 ],
       [ 0.733, 10.88 ,  7.656,  0.67 ],
       [ 0.622, 10.88 ,  7.376,  0.49 ],
       [ 1.034, 10.91 ,  7.866,  0.77 ],
       [ 0.841, 11.34 ,  8.406,  0.7  ]])
if inp_SS13 ==1:
    plt.scatter(ss[:,1],ss[:,2],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')
  
inp_b11= 1
b11 = np.array([[ 1.227, 10.58 ,  8.87 ,  9.83 ],
       [ 1.9  , 10.64 ,  9.17 , 10.64 ],
       [ 1.22 , 10.54 ,  8.24 , 10.54 ],
       [ 1.031, 10.78 ,  7.85 ,  9.53 ],
       [ 1.617, 10.61 ,  8.08 , 10.61 ],
       [ 1.615, 10.45 ,  8.3  , 10.45 ],
       [ 1.037,  9.62 ,  7.75 ,  9.62 ],
       [ 1.218, 10.71 ,  8.37 , 10.71 ],
       [ 1.371, 10.9  ,  8.27 ,  9.99 ],
       [ 1.021, 10.96 ,  8.35 ,  9.29 ],
       [ 1.45 , 10.74 ,  8.77 , 10.74 ]])
if inp_b11 ==1:
    plt.scatter(b11[:,1],b11[:,2],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')

inp_Cis= 1
cis11 = np.array([[ 0.73,  7.72, 10.3 ],
       [ 0.34,  8.29, 11.23],
       [ 0.34,  8.08, 10.65],
       [ 0.34,  8.39, 11.02],
       [ 0.35,  7.39, 10.54],
       [ 0.34,  8.66, 11.14],
       [ 0.38,  7.77, 10.68],
       [ 0.35,  7.24, 10.95],
       [ 0.85,  8.29, 11.07],
       [ 0.7 ,  8.15, 11.17],
       [ 0.44,  7.79, 10.57],
       [ 0.35,  7.59, 10.47],
       [ 0.37,  8.58, 10.57],
       [ 0.77,  8.49, 10.86],
       [ 0.73,  8.03, 11.01],
       [ 0.83,  8.07, 10.81],
       [ 0.52,  8.01, 10.54],
       [ 0.73,  7.41, 10.36],
       [ 0.36,  8.07, 11.28],
       [ 0.55,  7.75, 11.08],
       [ 0.69,  7.91, 10.66],
       [ 0.53,  8.22, 10.84],
       [ 0.62,  7.35, 10.53],
       [ 0.67,  7.73, 10.75],
       [ 0.79,  8.24, 10.53],
       [ 0.52,  8.38, 11.15],
       [ 0.37,  7.7 , 10.48],
       [ 0.55,  8.61, 11.2 ],
       [ 0.63,  7.5 , 10.73],
       [ 0.82,  7.82, 10.77],
       [ 0.66,  8.19, 11.03],
       [ 0.38,  8.25, 11.08]])
if inp_Cis ==1:
    plt.scatter(cis11[:,2],cis11[:,1],c='darkseagreen',marker="^",s=180,zorder=100, edgecolors='white')

# inp_Knud= 0
# Knud = np.array([[ 1.841, 46.58 ,  8.55 ,  6.56 , 10.69 ],
#        [ 2.039, 47.15 ,  8.68 , 24.63 , 10.81 ],
#        [ 2.086, 46.76 ,  8.67 , 10.22 , 10.8  ],
#        [ 1.96 , 46.92 ,  8.65 , 14.33 , 10.77 ],
#        [ 2.064, 46.88 ,  8.67 , 13.11 , 10.79 ],
#        [ 1.86 , 46.72 ,  8.6  ,  9.15 , 10.73 ],
#        [ 1.879, 46.61 ,  8.62 ,  7.21 , 10.75 ],
#        [ 1.953, 46.89 ,  8.67 , 13.74 , 10.79 ],
#        [ 1.894, 46.96 ,  8.6  , 15.88 , 10.73 ],
#        [ 1.807, 46.66 ,  8.64 ,  7.92 , 10.77 ],
#        [ 1.846, 46.59 ,  8.63 ,  6.78 , 10.76 ],
#        [ 1.873, 46.85 ,  8.58 , 12.38 , 10.72 ],
#        [ 2.023, 46.55 ,  8.58 ,  6.19 , 10.72 ],
#        [ 1.934, 46.74 ,  8.64 ,  9.77 , 10.77 ],
#        [ 1.96 , 46.83 ,  8.56 , 11.95 , 10.69 ],
#        [ 1.943, 46.71 ,  8.68 ,  8.96 , 10.81 ],
#        [ 1.966, 46.66 ,  8.69 ,  8.14 , 10.81 ],
#        [ 2.145, 46.83 ,  8.61 , 11.67 , 10.74 ],
#        [ 1.984, 46.95 ,  8.62 , 15.9  , 10.75 ],
#        [ 1.914, 47.02 ,  8.51 , 18.22 , 10.65 ],
#        [ 2.004, 46.7  ,  8.63 ,  8.67 , 10.76 ]])
# if inp_Knud ==1:
#     plt.scatter(Knud[:,4], Knud[:,2],c='blue',marker="o",s=180,zorder=101, edgecolors='white')

tx, ty = 11.85, 7.
plt.text(tx, ty, "intermediate\n  sample\nuncertainties",  fontsize=20)
plt.errorbar(tx+0.2,ty-0.08, xerr=0.2, yerr=0.4, color='darkseagreen',ecolor='black', fmt='^',zorder=-500,markersize=10)
#%%
#==============================================================================
# My new inference
#==============================================================================
from scipy.integrate import quad
# from load_result import load_host_p, load_MBH, load_err
def EZ(z,om):
      return 1/np.sqrt(om*(1+z)**3+(1-om))
def EE(z):
      return quad(EZ, 0, z , args=(om))[0]
vec_EE=np.vectorize(EE)

h0=70.             #km/s/Mpc
om=0.3
c=299790.        #speed of light [km/s]

ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'CDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','XID2138','XID2202',\
'XID2396', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
MB_ID = ['CDFS-1', 'CID543','CID70',  'SXDS-X735', 'CDFS-229', 'ECDFS-321', 'CID1174',\
'CID216', 'CID237','CID3242','CID3570','CID452', 'CID454',\
'CID50','CID607','LID1273', 'LID1538','LID360','SXDS-X1136',\
'SXDS-X50', 'SXDS-X717','SXDS-X763','SXDS-X969','LID1820','LID1622',\
'LID1878', 'CID206', 'ECDFS-358', 'CDFS-724', 'CID597','CID1281','CID255']
zs = np.array([1.63 , 1.301, 1.667, 1.447, 1.326, 1.57 , 1.552, 1.567, 1.618,
       1.532, 1.244, 1.407, 1.478, 1.239, 1.294, 1.617, 1.527, 1.579,
       1.325, 1.411, 1.276, 1.412, 1.585, 1.551, 1.516, 1.6  , 1.483,
       1.626, 1.337, 1.272, 1.445, 1.664])
host_n = np.array([4.8064, 0.4933, 3.616 , 2.0273, 0.4632, 2.2601, 1.8625, 6.1593,
       4.7361, 6.1424, 0.7278, 1.4017, 0.6239, 3.2167, 3.418 , 1.2204,
       2.8195, 0.7517, 1.9862, 1.6483, 5.58  , 2.3938, 2.0529, 1.2396,
       3.9565, 0.7716, 3.0717, 1.6708, 1.6177, 1.7523, 3.1462, 4.2099])
Mstar = np.array([10.2864, 10.4345, 10.5518, 10.7849, 10.6283, 11.1009, 10.6314,
       10.6296, 10.7542, 10.747 , 10.7134, 10.8578, 10.6966, 10.801 ,
       10.7502, 10.8897, 10.7078, 10.6594, 10.4861, 10.5363, 10.5055,
        9.6791, 10.6088, 10.4758, 10.7369, 10.6927, 10.4496, 10.7348,
        9.7876, 10.4582,  9.9995, 10.6496])
MBs = np.array([7.8378, 8.1868, 8.2689, 8.4403, 8.1997, 8.457 , 7.9915, 7.8539,
       8.2891, 8.4494, 7.8876, 8.2953, 8.3114, 8.4229, 8.5263, 8.5606,
       8.4694, 8.4497, 8.3326, 7.9417, 8.2039, 8.4646, 8.1956, 8.5519,
       8.4632, 8.455 , 8.5278, 8.1218, 8.0246, 7.8128, 7.7488, 8.2727])
Mstar_err = np.array([[-0.15,  0.19], [-0.15,  0.19], [-0.14,  0.16], [-0.14,  0.17], 
                      [-0.11,  0.12], [-0.2 ,  0.3 ], [-0.15,  0.18], [-0.1 ,  0.1 ], 
                      [-0.13,  0.14], [-0.15,  0.17], [-0.1 ,  0.1 ], [-0.1 ,  0.1 ], 
                      [-0.1 ,  0.1 ], [-0.21,  0.36], [-0.18,  0.25], [-0.12,  0.13], 
                      [-0.12,  0.13], [-0.11,  0.11], [-0.13,  0.14], [-0.13,  0.15], 
                      [-0.12,  0.12], [-0.24,  0.48], [-0.17,  0.23], [-0.12,  0.12], 
                      [-0.14,  0.16], [-0.19,  0.28], [-0.25,  0.53], [-0.14,  0.16], 
                      [-0.18,  0.25], [-0.18,  0.24], [-0.15,  0.18], [-0.15,  0.18]])
plt.scatter(Mstar,MBs,c='tomato',s=280,marker="*",zorder=100, edgecolors='k')
plt.errorbar(Mstar,MBs, xerr=[np.abs(Mstar_err)[:,0], np.abs(Mstar_err)[:,1]], yerr=0.4, color='blue',ecolor='orange', fmt='.',zorder=-500,markersize=1)

#%%
#folders = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']
h0licow = ['HE0435', 'RXJ1131', 'WFI2033', 'SDSS1206', 'HE1104', 'SDSS0246', 'HE0047']#, 'HS2209']
sys.path.insert(0,'../Stellar_mass_estimate')
#from mstar_inference import mstar_dict
sys.path.insert(0,'../share_tools')
from read_inference import read_mstar
sys.path.insert(0,'../MBH_estimator')
from est_MBH import MBH_dic
texts = []
for ID in h0licow:
    Mstar = np.log10(read_mstar(ID, count_n=[4, 4])[0][0])
    MBs = MBH_dic[ID]
    plt.scatter(Mstar, MBs, c='blue',s=580,marker=".",zorder=100, edgecolors='k')
    texts.append( plt.text(Mstar, MBs, ID, fontsize=15, zorder=200) )
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red')) 
    
#%%    

#==============================================================================
# The colorbar label setting up
#==============================================================================

plt.title(r"M$_{\rm BH}-$M$_*$ relation",fontsize=35)
plt.xlabel(r"log(M$_*$/M$_{\odot})$",fontsize=35)
plt.ylabel(r'log(M$_{\rm BH}$/M$_{\odot}$)',fontsize=35)
plt.xlim(9,12.5)
plt.ylim(6.0,10)
plt.grid(linestyle='--')
plt.tick_params(labelsize=25)

Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)
SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=13)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=16,markeredgecolor='k')

plt.legend([Bkc,Hkc,SS13,new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"intermediate redshift AGNs",
"Ding+2020 32 QSO"\
],scatterpoints=1,numpoints=1,loc=2,prop={'size':20},ncol=2)
plt.savefig("MBH-Mstar.pdf")
plt.show()