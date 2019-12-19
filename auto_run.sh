#!/bin/bash
i=2
target="2_WFI2033"
cd $target/model/fit_PSFi_QSOmask/
echo "$target/model/fit_PSFi_QSOmask/"
python 2_modelling_fix_gamma_QSO_mask.py ${i}
cd ../fit_PSFi_PSFrecons/
echo "$target/model/fit_PSFi_PSFrecons/"
python 2_modelling_fix_gamma_PSF_recons.py ${i}

