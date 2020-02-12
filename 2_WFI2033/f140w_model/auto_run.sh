#!/bin/bash
fixg_i=2
#target="2_WFI2033"
cd ./fit_PSFi_QSOmask/
echo "fit_PSFi_QSOmask/"
python 2_modelling_fix_gamma_QSO_mask.py ${fixg_i}
cd ../fit_PSFi_PSFrecons/
echo "fit_PSFi_PSFrecons/"
python 2_modelling_fix_gamma_PSF_recons.py ${fixg_i}

