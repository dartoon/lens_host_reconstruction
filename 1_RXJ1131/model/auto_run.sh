#!/bin/bash
fixg_i=2
#target="2_WFI2033"
cd ./2nd_fit_PSFi_QSOmask/
#echo "2nd_fit_PSFi_QSOmask/"
python 2_modelling_fix_gamma_QSO_mask.py ${fixg_i}
#cd ../2nd_fit_PSFi_PSFrecons/
#echo "2nd_fit_PSFi_PSFrecons/"
#python 2_modelling_fix_gamma_PSF_recons.py ${fixg_i}

