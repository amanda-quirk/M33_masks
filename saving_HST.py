from astropy.io import fits
import numpy as np 
import h5py 

#HST data-- using the old i filter; this has already been dust corrected
#HST_hdu = fits.open('/Users/amandaquirk/Documents/M33/data/m33-merged-F475W-F814W.fits', memmap=True)
#HST_data = HST_hdu[1].data 
#HST_RA = HST_data['ra'] #degrees
# HST_Dec = HST_data['dec'] #degrees
# F814W = HST_data['f814w_vega']

# output_file = h5py.File('HST_data_for_isolation.hdf5', 'w')
# output_file.create_dataset("RAs", data=HST_RA)
# output_file.create_dataset("Decs", data=HST_Dec)
# output_file.create_dataset("mag", data=F814W)
# output_file.close() 

input_file = h5py.File('HST_julia_isolated.hdf5', 'r')
tag = input_file["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated 
print(len(tag))
