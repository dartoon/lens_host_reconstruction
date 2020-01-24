#!/bin/bash
fixg_i=2
for target in 2_WFI2033 3_SDSS1206 4_HE1104 5_SDSS0246 6_HS2209 7_HE0047;do        #the number of parameter
#target="2_WFI2033"
	cd /Users/Dartoon/Astro/HOLICOW/lens_host_reconstruction
	cd $target/model/2nd_fit_PSFi_QSOmask/
	echo "$target/model/2nd_fit_PSFi_QSOmask/"
	python 2_modelling_fix_gamma_QSO_mask.py ${fixg_i}
	cd ../2nd_fit_PSFi_PSFrecons/
	echo "$target/model/2nd_fit_PSFi_PSFrecons/"
	python 2_modelling_fix_gamma_PSF_recons.py ${fixg_i}
done
