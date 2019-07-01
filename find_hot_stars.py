'''
purpose: find stars labeled MS and BHeB for Yong's observing run 
will return 3 lists: all of the host stars from the HST photometry, all of the hot stars that made it onto our masks, all of the hot stars on our masks that were given zqual of 2
'''

import numpy as np 
from astropy.io import fits
import pandas as pd 

#full list =======================================================================================
# HST_hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/fall18_refine.fits')
# HST_data = HST_hdu[1].data 
# HST_RA = HST_data['RA']
# HST_Dec = HST_data['Dec']
# HST_F814W = HST_data['F814W_VEGA']
# HST_F475W = HST_data['F475W_VEGA']
# HST_F814W_crowd = HST_data['F814CROWDMAG']
# HST_F475W_crowd = HST_data['F475CROWDMAG']
# HST_ID = HST_data['TARGTYPE']

# print('All the damn data read in!')

# all_hot_ra = []
# all_hot_dec = []
# all_hot_814 = []
# all_hot_475 = []
# all_hot_814_crowd = []
# all_hot_475_crowd = []
# all_hot_ID = []

# for i in range(len(HST_RA)):
# 	if HST_ID[i] == 'MS' or HST_ID[i] == ('BHeB_bright'):
# 		all_hot_ra.append(HST_RA[i])
# 		all_hot_dec.append(HST_Dec[i])
# 		all_hot_814.append(HST_F814W[i])
# 		all_hot_475.append(HST_F475W[i])
# 		all_hot_814_crowd.append(HST_F814W_crowd)
# 		all_hot_475_crowd.append(HST_F475W_crowd)
# 		all_hot_ID.append('{}_{}'.format(HST_ID[i], i))
# 	if i % 10000 == 0:
		# print(i / len(HST_RA))

#np.savetxt('/Users/amandaquirk/Desktop/all_host_stars.text', np.c_[all_hot_ID, all_hot_ra, all_hot_dec, all_hot_475, all_hot_814, all_hot_475_crowd, all_hot_814_crowd])


ref_data = np.genfromtxt('/Users/amandaquirk/Documents/M33/Data/all_target_list.in', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, magnitude2, filter2, priority, list_assignment, selection_flag, HST_F110W, HST_F160W, HST_F275W, HST_F336W')

ref_ID = ref_data['ID']
ref_ID = [a.decode("utf-8") for a in ref_ID]
ra = ref_data['ras']
ra = np.array([a.decode("utf-8") for a in ra])
dec = ref_data['decs']
dec = np.array([a.decode("utf-8") for a in dec])
mag1 = ref_data['magnitude1']
mag2 = ref_data['magnitude2']

print('Read in all target list data')

#mask out everything but young stars selected from HST
ref_ID = pd.Series(ref_ID)
HST_selected = (ref_ID.str.startswith('MS') == True) | (ref_ID.str.startswith('BHeB') == True) 

# HST_selected = np.array((HST_selected))
ra = ra[HST_selected]
dec = dec[HST_selected]
mag1 = mag1[HST_selected]
mag2 = mag2[HST_selected]
ref_ID = ref_ID[HST_selected]

print('There are {} HST young stars'.format(sum(HST_selected)))

#full list ==================================================================================
#only using ones from my target list, since these are the ones that passed the isolation tests
np.savetxt('/Users/amandaquirk/Desktop/all_hot_stars.txt', np.c_[ref_ID, ra, dec, mag2, mag1], fmt="%-s", delimiter='\t', header='ID, RA, Dec, F475W, F814W')

#list of the ones reduced with zqual 2 ========================================================
def import_zspec_data(maskname):
	hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	data = hdu[1].data
	zqual = data['ZQUALITY']
	zID = data['OBJNAME']

	return zqual, zID

zqual_a1, ID_a1 = import_zspec_data('A1M33P')
zqual_b1, ID_b1 = import_zspec_data('B1M33P')
zqual_c1, ID_c1 = import_zspec_data('C1M33P')

all_IDs = list(ID_a1) + list(ID_b1) + list(ID_c1)
all_zquals = np.array(list(zqual_a1) + list(zqual_b1) + list(zqual_c1))

#look at zqual 2 hot stars
all_IDs = pd.Series(all_IDs)
bad_hot_star = ((all_IDs.str.startswith('MS') == True) | (all_IDs.str.startswith('BHeB') == True)) & (all_zquals == 2)
all_zquals = all_zquals[bad_hot_star]
all_IDs = np.array(all_IDs)
all_IDs = all_IDs[bad_hot_star]

zqual2_ra = []
zqual2_dec = []
zqual2_475 = []
zqual2_814 = []
zqual2_ID = []
ref_ID = list(ref_ID)
for i in range(len(all_IDs)):
	N = ref_ID.index(all_IDs[i]) #index of where the observed star is in our target list
	zqual2_ra.append(ra[N])
	zqual2_dec.append(dec[N])
	zqual2_475.append(mag2[N])
	zqual2_814.append(mag1[N])
	zqual2_ID.append(ref_ID[N])

np.savetxt('/Users/amandaquirk/Desktop/zqual2_hot_stars.txt', np.c_[zqual2_ID, zqual2_ra, zqual2_dec, zqual2_475, zqual2_814], fmt="%-s", delimiter='\t', header='ID, RA, Dec, F475W, F814W')

#list of all the ones on our mask =============================================================
def read_mask_data(maskname):
	ID = np.genfromtxt('/Users/amandaquirk/Documents/M33/Masks/Mask_Designs/{}.out'.format(maskname), usecols=(0,), dtype='str', max_rows=201) #only want to look at the stars that are actually on our list so want to cut it off early
	ID = pd.Series(ID)
	hot_stars = (ID.str.startswith('MS') == True) | (ID.str.startswith('BHeB') == True) #eliminate CHFT and KG targets
	return ID[hot_stars]

#read in all the data
mask_ID_a1 = read_mask_data('A1M33P')
mask_ID_a2 = read_mask_data('A2M33P')
mask_ID_b1 = read_mask_data('B1M33P')
mask_ID_b2 = read_mask_data('B2M33P')
mask_ID_c1 = read_mask_data('C1M33P')
mask_ID_c2 = read_mask_data('C2M33P')
mask_ID_d1 = read_mask_data('D1M33P')
mask_ID_d2 = read_mask_data('D2M33P')
mask_ID_e1 = read_mask_data('E1M33P')
mask_ID_e2 = read_mask_data('E2M33P')
mask_ID_k1 = read_mask_data('K1M33P')

print('Read in the mask data')

all_mask_ID = list(mask_ID_a1) + list(mask_ID_a2) + list(mask_ID_b1) + list(mask_ID_b2) + list(mask_ID_c1) + list(mask_ID_c2) + list(mask_ID_d1) + list(mask_ID_d2) + list(mask_ID_e1) + list(mask_ID_e2) + list(mask_ID_k1)

mask_ra = []
mask_dec = []
mask_475 = []
mask_814 = []
mask_ID = []
for i in range(len(all_mask_ID)):
	N = ref_ID.index(all_mask_ID[i]) #index of where the observed star is in our target list
	mask_ra.append(ra[N])
	mask_dec.append(dec[N])
	mask_475.append(mag2[N])
	mask_814.append(mag1[N])
	mask_ID.append(ref_ID[N])

np.savetxt('/Users/amandaquirk/Desktop/all_masks_hot_stars.txt', np.c_[mask_ID, mask_ra, mask_dec, mask_475, mask_814], fmt="%-s", delimiter='\t', header='ID, RA, Dec, F475W, F814W')
	

