import numpy as np
import astropy.io.fits as pyfits
import os
import string
import glob
#files='./esmod_3band_*'
files=glob.glob("esmod_*")
for i in range(len(files)):
  filename=files[i]
  d = pyfits.open(filename)[0].data.copy()
  arc=pyfits.open(filename,mode = 'update')
  frame = len(d)
  dx= (arc[0].header[12]-arc[0].header[11])/frame
  dy= (arc[0].header[14]-arc[0].header[13])/frame
  sq= dx*dy/(0.08**2)
  print np.sqrt(sq)
  arc.flush()
  d=d*sq
#  pyfits.PrimaryHDU(d).writeto('flux-{0}'.format(filename),overwrite=True)
