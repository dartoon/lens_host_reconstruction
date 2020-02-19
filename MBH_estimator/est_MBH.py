import numpy as np
MBH_by_Domini =  {'HE0435': [8.76, 0.44], 'RXJ1131': [7.9, 0.6], 'WFI2033': [8.63, 0.35], 'SDSS1206': [8.62, 0.24], 
               'HE1104': [9.37, 0.5], 'SDSS0246': [8.59, 0.36], 'HE0047': [8.83, 0.23]}#,'HS2209': [None, None],}

#HE1104 From Peng 06: Arxiv: 0603248 #For Civ: F_A*2.99792*10**5/(1+z)/1549 F_A = 103.0 , z = 2.32 (Table 3)

#SDSS1206: See Birrer J1206 paper and Shen 2011. Civ don't have board line
#For L, the L_Mgii are translated to L_3100 based on the equation (11) in Shen2011
#http://vizier.cfa.harvard.edu/viz-bin/VizieR-5?-ref=VIZ5af3dce413049&-out.add=.&-source=J/ApJS/194/45/catalog&recno=54910

#For Dominique's sample
#For bolometric correction. 
#Lbol for Mgii should be -np.log10(5.15)
#Lbol for Hb should be -np.log10(9.6), 
#Lbol for Civ should be -np.log10(3.81), 

BL_info =  {'HE0435': {'Mgii': {'L': [45.7 -np.log10(5.15), 46.0 -np.log10(5.15)], 'FWHM':  4930.0}},
                 'RXJ1131': {'Mgii': {'L': [45.0-np.log10(5.15)], 'FWHM':  5630.0},
                             'Hb': {'L': [45.0-np.log10(9.6)], 'FWHM':  4545.0}},  #Should -np.loq10(9.6)
                 'WFI2033':  {'Mgii': {'L': [46.0-np.log10(5.15),  45.8-np.log10(5.15)], 'FWHM':   3960.0}},
                 'SDSS1206': {'Mgii': {'L': [(44.23-2.22)/0.909 -np.log10(4.)], 'FWHM':   5632.4}},
                 'HE1104':  {'Civ': {'L': [46.18], 'FWHM':   6004.0}}, 
                 'SDSS0246': {'Mgii': {'L': [46.0-np.log10(5.15), 45.8-np.log10(5.15)], 'FWHM':   3700.0}},
                 'HE0047': {'Mgii': {'L': [46.3-np.log10(5.15), 46.3-np.log10(5.15)], 'FWHM':   4145.0}}}
#                 'HS2209': [None, None]} 

def cal_MBH(line, L, FWHM):
    if line == 'Hb':
#        a=np.log10(5.9)+6         ## =6.7853   
#        b=0.69
        a = 6.91
        b = 0.5
        logM=a + b*(L-44.) +2 * np.log10(FWHM/1000.)
    if line == 'Mgii':
        """
        #Dominique recipe. (note L - 5.15.) for his calbration:
        a, b = 7.13, 0.5 # 1.51 * FWHM_term
        logM=a + b*(L-44.) + 1.51 * np.log10(FWHM/1000.)
        #Shen 2011 recipe. L - 5.15 for correction.:
        a, b = 6.74, 0.62
        logM=a + b*(L-44.) + 2.0 * np.log10(FWHM/1000.)
        """
        a, b = 6.623, 0.47  #Final adopted receipe
        logM=a + b*(L-44.) +2 * np.log10(FWHM/1000.)
#        a, b = 7.13, 0.5 #Dominique;s paper
#        logM=a + b*(L-44.) + 1.51 * np.log10(FWHM/1000.)
#        a, b = 0.74, 0.62 #Shen 2011 recipe
#        logM=a + b*(L-44.) + 2.0 * np.log10(FWHM)
    if line == 'Civ':
#        a=np.log10(4.5)+6         ## =6.6532   
#        b=0.53
        a, b = 6.322, 0.53 #H0LiCOW VII
        logM=a + b*(L-44.) + 2.0 * np.log10(FWHM/1000.)
    return logM
        
def est_MBH(ID):
    logMBH = []
    for key in BL_info[ID].keys():
        logMBH.append( cal_MBH(key, np.mean(BL_info[ID][key]['L']), BL_info[ID][key]['FWHM']) )
    return logMBH

MBH_dic = {}
for key in MBH_by_Domini.keys():
    MBH_dic[key] = np.mean(est_MBH(key))
    print(key, est_MBH(key), MBH_by_Domini[key][0])
    
#MBH_dic =     {'HE0435': 8.543644440905101,
# 'RXJ1131': 8.347036339665747,
# 'WFI2033': 8.376840974201665,
# 'SDSS1206': 8.232486979634622,
# 'HE1104': 8.854081367142467,
# 'SDSS0246': 8.31785405048463,
# 'HE0047': 8.604499672123223}