#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 12:54:19 2020

@author: Dartoon
"""
import numpy as np
np.set_printoptions(precision=4)
import matplotlib.pyplot as plt
import matplotlib as mat
import matplotlib.lines as mlines
from matplotlib import colors
mat.rcParams['font.family'] = 'STIXGeneral'
host=plt.figure(figsize=(14.5,12))
ax=host.add_subplot(111)   #to get the log(1+z) and z label

import matplotlib as mpl
mpl.rc('image', cmap='jet')
import sys
from adjustText import adjust_text   # avoid the overlapping while ploting

# sys.path.insert(0,'../py_tools')

########## input local data ####
#==============================================================================
# The seleting for dm and host_total and dmag are in this local
#==============================================================================

#from dmag import pass_dmag

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
########input 30 local by Haring 04 ############
hloc = np.array([[3.7484e-03, 1.1778e+01, 1.8000e-01, 9.4771e+00, 2.1298e-01],
       [3.4930e-03, 1.0362e+01, 1.8000e-01, 7.1461e+00, 1.4262e-01],
       [2.4703e-03, 1.0833e+01, 1.8000e-01, 8.0000e+00, 1.0206e-01],
       [4.2822e-03, 1.1556e+01, 1.8000e-01, 8.6335e+00, 2.1306e-01],
       [7.3369e-03, 1.1556e+01, 1.8000e-01, 8.7160e+00, 3.2389e-01],
       [2.4297e-02, 1.1748e+01, 1.8000e-01, 8.7243e+00, 3.1358e-01],
       [1.3564e-02, 1.1462e+01, 1.8000e-01, 8.5185e+00, 2.0838e-01],
       [3.6091e-03, 9.7924e+00, 1.8000e-01, 7.1461e+00, 1.8221e-01],
       [5.6030e-03, 1.1114e+01, 1.8000e-01, 7.5682e+00, 1.6526e-01],
       [6.7825e-03, 1.1462e+01, 1.8000e-01, 9.3979e+00, 3.9967e-01],
       [1.7743e-04, 1.0568e+01, 1.8000e-01, 7.6532e+00, 1.8656e-01],
       [1.8910e-04, 8.9031e+00, 1.8000e-01, 6.3979e+00, 3.4062e-01],
       [2.6564e-03, 1.0839e+01, 1.8000e-01, 7.6435e+00, 4.6942e-01],
       [5.3251e-03, 1.0881e+01, 1.8000e-01, 7.1461e+00, 2.2797e-01],
       [2.2610e-03, 1.1079e+01, 1.8000e-01, 9.0000e+00, 2.3856e-01],
       [4.8618e-03, 1.0833e+01, 1.8000e-01, 8.3222e+00, 2.8354e-01],
       [2.6099e-03, 1.0204e+01, 1.8000e-01, 8.0000e+00, 6.1650e-01],
       [2.7029e-03, 1.0301e+01, 1.8000e-01, 7.2041e+00, 4.3571e-01],
       [5.3483e-03, 1.0987e+01, 1.8000e-01, 8.2788e+00, 2.5972e-01],
       [6.0889e-03, 1.1114e+01, 1.8000e-01, 8.4914e+00, 2.0586e-01],
       [3.5627e-03, 1.0079e+01, 1.8000e-01, 8.4771e+00, 2.4800e-01],
       [3.6556e-03, 1.0964e+01, 1.8000e-01, 8.0414e+00, 5.6735e-01],
       [3.4930e-03, 1.0643e+01, 1.8000e-01, 7.7482e+00, 4.3388e-01],
       [2.2842e-03, 1.1431e+01, 1.8000e-01, 9.0000e+00, 3.3450e-01],
       [3.9109e-03, 1.1690e+01, 1.8000e-01, 9.3010e+00, 4.8455e-02],
       [2.7262e-03, 1.1041e+01, 1.8000e-01, 8.2304e+00, 6.2621e-01],
       [6.0195e-03, 1.0568e+01, 1.8000e-01, 8.3802e+00, 1.0654e-01],
       [5.3483e-03, 1.0176e+01, 1.8000e-01, 7.1139e+00, 1.8447e-01],
       [3.0748e-03, 9.8451e+00, 1.8000e-01, 6.5441e+00, 1.4739e-01],
       [2.3349e-05, 1.0041e+01, 1.8000e-01, 6.5682e+00, 1.5707e-01]])

zs =  np.array([3.7484e-03, 3.4930e-03, 2.4703e-03, 4.2822e-03, 7.3369e-03,
       2.4297e-02, 1.3564e-02, 3.6091e-03, 5.6030e-03, 6.7825e-03,
       1.7743e-04, 1.8910e-04, 2.6564e-03, 5.3251e-03, 2.2610e-03,
       4.8618e-03, 2.6099e-03, 2.7029e-03, 5.3483e-03, 6.0889e-03,
       3.5627e-03, 3.6556e-03, 3.4930e-03, 2.2842e-03, 3.9109e-03,
       2.7262e-03, 6.0195e-03, 5.3483e-03, 3.0748e-03, 2.3349e-05])
hloc[:,0] = zs

#############################################################
###################fitting with MCMC#########################
x=np.append(bloc[:,1], hloc[:,1])
y=np.append(bloc[:,3], hloc[:,3])
yerr=(np.append(bloc[:,2], hloc[:,2])**2+np.append(bloc[:,4], hloc[:,4])**2)**0.5  # 0.2 is the uncertainty level for the L_R
def lnlike(theta, x, y, yerr):
    m, b, sint= theta
    model = m * x + b
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [1.036, -1.947, 0.3], args=(x, y, yerr))
m_ml, b_ml,sint_ml= result["x"]
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
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

m_mid, b_mid, sint_mid =np.percentile(samples, 50,axis=0)

######################
style = 1#input("0 as SS13, 1 as delta(logMBH):\n")
if style == 0:
    plt.plot(np.log10(hloc[:,0]+1),
             10**(hloc[:,3]-hloc[:,1]), '.',color='black',markersize=10)
    
    plt.plot(np.log10(bloc[:,0]+1),
                 10**(bloc[:,3]-bloc[:,1]),'.',color='gray',markersize=10)
elif style ==1:
    xl = np.linspace(-0.9, 13, 100)
    plt.errorbar(np.log10(hloc[:,0]+1),
                 hloc[:,3]-(m_ml*hloc[:,1]+b_ml),yerr=(hloc[:,2]**2 + hloc[:,4]**2)**0.5 ,fmt='.',color='black',markersize=10)
    plt.errorbar(np.log10(bloc[:,0]+1),
                 bloc[:,3]-(m_ml*bloc[:,1]+b_ml),yerr=(bloc[:,2]**2 + bloc[:,4]**2)**0.5 ,fmt='.',color='gray',markersize=10)
    ty=xl*0
    ty1=xl*0+np.std(y-(m_ml*x+b_ml))
    ty2=xl*0-np.std(y-(m_ml*x+b_ml))
    plt.fill_between(xl,ty1,ty2,color='linen',zorder=-50)

Bkc=mlines.Line2D([], [], color='gray', ls='', marker='.', markersize=15)
Hkc=mlines.Line2D([], [], color='black', ls='', marker='.', markersize=15)
######################

# f0 ='data/SS13_MM.txt'
# ss = np.loadtxt(f0)[:,1:]  #0 redshift; 1 M*; 2 BH mass;

# f1 ='data/B11_MM.txt'
# b11 = np.loadtxt(f1)[:,1:]  #0 redshift; 1 M*; 2 BH mass;

# f2 = 'data/Cisternas_data.txt'
# cis11 = np.loadtxt(f2)  #0 redshift;

# f3 = 'data/high_edd_agn.txt'
# Knud = np.loadtxt(f3)[:,2:]  # 0 redshift; 1 L_bol; 2 M_BH; 3 M_acc; 4 M_*
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

#%%
#==============================================================================
# My new inference
#==============================================================================
# from load_result import load_host_p, load_MBH, load_err
# from load_result import load_zs, load_n
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
yerr_highz = [((m_ml*Mstar_err[:,0])**2+0.4**2)**0.5, ((m_ml*Mstar_err[:,1])**2+0.4**2)**0.5]

#plt.scatter(Mstar,MBs,c=zs,s=880,marker="*",zorder=100, vmin=0.3, vmax=2, edgecolors='k')

#plt.errorbar(np.log10(1+zs),MBs-(m_ml*lumi_s+b_ml),yerr=(0.4**2+0.2**2)**0.5,fmt='x',color='royalblue',markersize=28,zorder=250)#, mec='k')
# if style ==0:
#     plt.scatter(np.log10(1+ss[:,0]), 10** (ss[:,2]-ss[:,1]), c='darkseagreen',
#                 s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')
    
#     plt.scatter(np.log10(1+b11[:,0]), 10** (b11[:,2]-b11[:,1]), c='darkseagreen',
#                 s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')    
    
#     plt.scatter(np.log10(1+cis11[:,0]),10** (cis11[:,1]-cis11[:,2]), c='darkseagreen',
#                 s=180,marker="^", zorder=100,vmin=0.3, vmax=2, edgecolors='white')
    
#     plt.scatter(np.log10(1+zs), 10**( MBs- Mstar),c='tomato',
#                 s=580,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')

#     plt.scatter(np.log10(1+zs[9]), np.median(np.concatenate([10**(hloc[:,3]-hloc[:,1]),10**(bloc[:,3]-bloc[:,1])])),
#                 facecolors='none', s=280,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='blue', linewidth=6, alpha=0.5)       
#     plt.scatter(np.log10(1+zs[9]), np.median(np.concatenate([10**(hloc[:,3]-hloc[:,1]+0.21),10**(bloc[:,3]-bloc[:,1]+0.21)])),
#                 facecolors='none', s=280,marker="o",zorder=900, vmin=0.3, vmax=5, edgecolors='blue', linewidth=6, alpha=0.5)      
#     plt.arrow(np.log10(1+zs[9]),np.median(np.concatenate([10**(hloc[:,3]-hloc[:,1]),10**(bloc[:,3]-bloc[:,1])])), 0, +0.0009,
#               zorder=800, head_length= 0.001374/8,head_width= 0.005,fc='k',ec='k')
    
# if style ==1:
plt.errorbar(np.log10(1+ss[:,0]),ss[:,2]-(m_ml*ss[:,1]+b_ml),yerr=(0.4**2+(m_ml*0.2)**2)**0.5,fmt='^',color='darkseagreen',markersize=9)
plt.errorbar(np.log10(1+b11[:,0]),b11[:,2]-(m_ml*b11[:,1]+b_ml),yerr=(0.4**2+(m_ml*0.2)**2)**0.5,fmt='^',color='darkseagreen',markersize=9)  

plt.errorbar(np.log10(1+cis11[:,0]),cis11[:,1]-(m_ml*cis11[:,2]+b_ml),yerr=(0.4**2+(m_ml*0.35)**2)**0.5,fmt='^',color='darkseagreen',markersize=9) 

#    plt.errorbar(np.log10(1+Knud[:,0]),Knud[:,2]-(m_ml*Knud[:,4]+b_ml),yerr=0,fmt='o',color='blue',markersize=9) 


plt.scatter(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),c='tomato',
            s=280,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
plt.errorbar(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),
             yerr= yerr_highz,
             color='tomato',ecolor='orange', fmt='.',markersize=1)    

plt.scatter(np.log10(1+zs[9]),0,facecolors='none',
            s=280,marker="o",zorder=900, edgecolors='blue', linewidth=6, alpha=0.5)    
plt.scatter(np.log10(1+zs[9]),0+0.21,facecolors='none',
            s=280,marker="o",zorder=900, edgecolors='blue', linewidth=6, alpha=0.5)    
plt.arrow(np.log10(1+zs[9]),0, 0, +0.19, zorder=800, head_length=0.05,head_width=0.005,fc='k',ec='k')

#####fit the evolution##########
################################
z_ss, y_ss = ss[:,0], ss[:,2]-(m_ml*ss[:,1]+b_ml)

z_b11, y_b11 = b11[:,0], b11[:,2]-(m_ml*b11[:,1]+b_ml)

z_cis, y_cis = cis11[:,0], cis11[:,1]-(m_ml*cis11[:,2]+b_ml)

z_cosmos, y_cosmos = zs, MBs-(m_ml*Mstar+b_ml)
yerr_hz = (yerr_highz[0]+ yerr_highz[1])/2
                              
z=np.concatenate((z_ss, z_b11, z_cis, z_cosmos),axis=0)
y=np.concatenate((y_ss, y_b11, y_cis, y_cosmos),axis=0)
yerr_imd= np.zeros(len(z_ss)+len(z_b11))+(0.4**2+(m_ml*0.2)**2)**0.5   # the error for the fitting
yerr_cis = np.zeros(len(z_cis)) + (0.4**2+(m_ml*0.35)**2)**0.5 
yerr = np.concatenate((yerr_imd,yerr_cis, yerr_hz),axis=0)

#if consider 32 AGN only:
z=z_cosmos
y=y_cosmos
yerr = yerr_hz    

yerr = np.sqrt(yerr**2 + sint_ml**2)
#### fit with emcee ###############
x=np.log10(1+z)
y=y

def lnlike(theta, x, y, yerr):
    b, sint= theta
    model = b*x
    sigma2 = (yerr**2 + sint**2)
    if sint>=0 :
      return -0.5*(np.sum((y-model)**2/sigma2)+np.sum(np.log(2*np.pi*sigma2)))
    else:
      return -np.inf

import scipy.optimize as op
nll = lambda *args: -lnlike(*args)
result = op.minimize(nll, [1.8, 0.3], args=(x, y, yerr))
b_ml_offset,_= result["x"]

xp = np.array([5, 13])
def lnprior(theta):
    b, sint	 = theta
    if -10 < b < 10.0 and 0 < sint < 10:
        return 0.0
    return -np.inf
def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)
ndim, nwalkers = 2, 100
pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
import emcee
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
sampler.run_mcmc(pos, 500)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

b_ml_offset, _ =np.percentile(samples, 50,axis=0)
#print "lnlike=",lnlike(theta=[b_ml_offset, sint_mid],x=x, y=y, yerr=yerr)
xl = np.linspace(0, 5, 100)
plt.plot(xl, xl*0+xl*b_ml_offset, color="red", linewidth=4.0,zorder=0)

plt.plot(xl, xl*0, color="black", linewidth=2.0,zorder=0)
def find_n(array,value):           #get the corresponding b for a given m 
    idx= (np.abs(array-value)).argmin()
    return array[idx]
b=np.percentile(samples,50,axis=0)[0]
#print samples[:,1][samples[:,0]==find_n(samples[:,0],m)]
for i in range(100):
    posi=np.random.uniform(16,84)
    b=np.percentile(samples,posi,axis=0)[0]    
    #print b
    plt.plot(xl, xl*0+xl*b, color="lightgray", alpha=0.2,linewidth=7.0,zorder=-1)
value=round(b_ml_offset,2)
#####################
value,sig=round(b_ml_offset,2),round((np.percentile(samples,84,axis=0)[0]-np.percentile(samples,16,axis=0)[0])/2,2)

#%%
#folders = ['0_HE0435', '1_RXJ1131', '2_WFI2033', '3_SDSS1206', '4_HE1104', '5_SDSS0246', '6_HS2209', '7_HE0047']
h0licow = ['HE0435', 'RXJ1131', 'WFI2033', 'SDSS1206', 'HE1104', 'SDSS0246', 'HE0047']#, 'HS2209']
sys.path.insert(0,'../Stellar_mass_estimate')
from mstar_inference import mstar_dict
sys.path.insert(0,'../MBH_estimator')
from est_MBH import MBH_dic
sys.path.insert(0,'../share_tools')
from lens_information import lens_redshift
from read_inference import read_mstar

texts = []
for ID in h0licow:
    Mstar = mstar_dict[ID]
    MBs = MBH_dic[ID]
    if ID == 'RXJ1131':
        Mstar = np.log10(read_mstar(ID, count_n=[4, 4])[0][0])
#    plt.scatter(Mstar, MBs[0], c='blue',s=580,marker=".",zorder=100, edgecolors='k')
    plt.scatter(np.log10(1+lens_redshift[ID]),MBs-(m_ml*Mstar+b_ml),c='blue',
            s=280,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
    plt.text(np.log10(1+lens_redshift[ID]),MBs-(m_ml*Mstar+b_ml), ID, fontsize=15, zorder=200) 
#adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red')) 
  

#plt.scatter(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),c='tomato',
#            s=280,marker="*",zorder=300, vmin=0.3, vmax=5, edgecolors='k')
#plt.errorbar(np.log10(1+zs),MBs-(m_ml*Mstar+b_ml),
#             yerr= yerr_highz,
#             color='tomato',ecolor='orange', fmt='.',markersize=1)    

#%% Where loop ends
plt.xlabel(r"log(1+z)",fontsize=45)
new_sample = mlines.Line2D([], [], color='tomato', ls='', marker='*', markersize=20,markeredgecolor='k')

plt.xticks(np.arange(-0.1,1,0.1))
xl=-0.01
xh=np.log10(1+2.5)
if style ==0:
    ax.set_yscale('log')
    plt.axis([xl,xh,0,0.5])
    plt.ylabel(r"M$_{\rm BH}$/M$_*$",fontsize=45)
if style ==1:
    plt.yticks(np.arange(-5.5,6,0.5))
    plt.axis([xl,xh,-2.0,3.5])
    plt.ylim([-2.0,3.5])
    plt.ylabel(r"$\Delta$logM$_{\rm BH}$ (vs M$_*$)",fontsize=45)
plt.grid()
plt.tick_params(labelsize=35)

ax2=ax.twiny()
tticks=np.array([10**xl-1,0.5,1,1.5,2,10**xh-1])
ax2.set_xticks([np.log(t+1) for t in tticks])  # for the entire scale
ax2.set_xticklabels([0,0.5,1,1.5,2,2.5])  # 0 actuall is corresponds to 10**-0.01-1
ax2.set_xlabel('z',fontsize=45)
plt.tick_params(labelsize=35)

SS13 = mlines.Line2D([], [], color='darkseagreen', ls='', marker='^', markersize=8)

plt.legend([Bkc, Hkc, SS13, new_sample],[
'Local by Bennert+11',\
"Local by H&R",
"Intermediate redshift AGNs",
"Ding+2020 32 QSO"
],scatterpoints=1,numpoints=1,loc=2,prop={'size':28},ncol=2,handletextpad=0)
plt.savefig("MBH-Mstar_vz.pdf")
plt.show()