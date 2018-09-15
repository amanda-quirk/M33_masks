from astropy.io import fits
import h5py 
import numpy as np 

#HST data
# HST_hdu = fits.open('/Volumes/Enterprise/HST_M33_Data/merged-6filt-bricks23.fits', memmap=True)
# HST_data = HST_hdu[1].data 
# HST_RA = HST_data['ra'] #degrees
# HST_Dec = HST_data['dec'] #degrees
# F814W = HST_data['f814w_vega']

RA, Dec, i = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/CFHT_data.txt', unpack=True)

# f = h5py.File('/Users/amandaquirk/Desktop/m33_initsample.hdf5', 'r')
# dset = f['df']
# dset0 = dset['block0_values']
# HST_RA = dset0[:, 90]
# HST_Dec = dset0[:, 5]
# F814W = dset0[:, 86]

output_file = h5py.File('/Users/amandaquirk/Documents/M33/Data/CFHT_data_for_isolation.hdf5', 'w')
output_file.create_dataset("RAs", data=RA)
output_file.create_dataset("Decs", data=Dec)
output_file.create_dataset("mag", data=i)
output_file.close() 

# input_file = h5py.File('/Users/amandaquirk/Documents/M33/Data/HST_julia_isolated.hdf5', 'r')
# tag = input_file["isolation_tag"][...] #[...] loads the data; 1 means it is NOT isolated 
# print(len(tag))
