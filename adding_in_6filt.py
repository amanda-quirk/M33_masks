'''
adding all the photometry to the original 2019 target list BECAUSE PAST AMANDA DIDN'T DO IT
'''
import numpy as np
from astropy.io import fits 

#grab the target list
ref_data = np.genfromtxt('/Users/amandaquirk/Documents/M33/Data/target_list_RGB_2019.in', dtype=None, names='ID, ras, decs, frame, F814W, filter1, priority, list_assignment, selection_flag')
ID = ref_data['ID']
ID = [a.decode("utf-8") for a in ID]

#break the IDs down to get the index of the row in original data file
RGB_ind = []
CFHT_ind = []
for name in ID:
	if name.startswith('RGB') == True:
		RGB_ind.append(int(name[4:]))
	else:
		CFHT_ind.append(int(name[4:]))

#grab the photometry
HST_hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/fall18_refine.fits')
HST_data = HST_hdu[1].data 
HST_F475W = HST_data['F475W_VEGA']
HST_F275W = HST_data['F275W_VEGA']
HST_F336W = HST_data['F336W_VEGA']
HST_F110W = HST_data['F110W_VEGA']
HST_F160W = HST_data['F160W_VEGA']

hdu = fits.open('/Users/amandaquirk/Documents/M33/data/M33.GI.matchcat.extcorr.fits')
data = hdu[1].data
CFHT_g_mag = data['MAG1_AUTO_0']
CFHT_i_mag = data['MAG2_AUTO_0']

#get the photometry for the data that's on the target list
HST_F475W = HST_F475W[RGB_ind]
HST_F275W = HST_F275W[RGB_ind]
HST_F336W = HST_F336W[RGB_ind]
HST_F110W = HST_F110W[RGB_ind]
HST_F160W = HST_F160W[RGB_ind]

CFHT_g_mag = CFHT_g_mag[CFHT_ind]
CFHT_i_mag = CFHT_i_mag[CFHT_ind]
CFHT_mag_extra = np.ones(len(CFHT_ind)) * np.nan

all_mags4 = np.concatenate((HST_F160W, CFHT_mag_extra), axis = None)
all_mags5 = np.concatenate((HST_F275W, CFHT_mag_extra), axis = None)
all_mags6 = np.concatenate((HST_F336W, CFHT_mag_extra), axis = None)
all_mags3 = np.concatenate((HST_F110W, CFHT_mag_extra), axis = None)
all_mags2 = np.concatenate((HST_F475W, CFHT_i_mag), axis = None)

#save ID, ras, decs, frame, F814W, filter1, F475W, filter2, priority, list_assignment, selection_flag, F110W, F160W, F275W, F336W
np.savetxt('/Users/amandaquirk/Documents/M33/Data/target_list_RGB_2019_6filt.in', np.c_[ID, ref_data['ras'], ref_data['decs'], ref_data['frame'], ref_data['F814W'], ref_data['filter1'], all_mags2, ref_data['filter1'], ref_data['priority'], ref_data['list_assignment'], ref_data['selection_flag'], all_mags3, all_mags4, all_mags5, all_mags6], fmt="%-s", delimiter='\t', header='ID, ras, decs, frame, F814W, filter1, F475W, filter2, priority, list_assignment, selection_flag, F110W, F160W, F275W, F336W') 