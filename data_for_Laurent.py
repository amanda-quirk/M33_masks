from astropy.io import fits
import numpy as np 

'''
put together all of the M33 data to send to Laurent; mostly just copy full_dataset.py
for now, file format: ID, RA, Dec, mag, filter name -- for 2016 data, there will be no photometry
ignores serendip stars
'''

def import_data(maskname, zspeced = True): #import spectroscopy data from zspec files
	if zspeced == True:
		hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	else:
		hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/zspecs/mgzresult.{}.fits'.format(maskname))
	data = hdu[1].data
	RA = data['RA']
	Dec = data['Dec'] 
	ID = data['OBJNAME'] 
	return ID, RA, Dec 

#all of our masks ========================================================================================================================================
masks_2018 = ['A1M33P', 'B1M33P', 'C1M33P', 'D1M33P', 'E1M33P'] #photometry from M33/Data/all_target_list.in
masks_2018_not_done = ['A2M33P', 'B2M33P', 'C2M33P', 'D2M33P', 'E2M33P', 'K1M33P'] #photometry from M33/Data/all_target_list.in
masks_2016 = ['M33D2A', 'M33D2B', 'M33D3A', 'M33D3B', 'M33D3D', 'M33D4A', 'M33D4B'] #don't have photometry
masks_2019 = ['D1M33R', 'D2M33R', 'E1M33R', 'E2M33R'] #E1M33R photometry from target_list_RGB_2019.in and rest photometry from M33/Masks/2019b/Dec/final/target_list_RGB_2019_expanded.in

#read in all of the mask data ============================================================================================================================
ras_19 = []
ids_19 = []
decs_19 = []
for mask in masks_2019:
	ids, ras, decs = import_data(mask, zspeced=True)
	ids_19 = list(ids_19) + list(ids)
	ras_19 = list(ras_19) + list(ras)
	decs_19 = list(decs_19) + list(decs)

ras_18 = []
ids_18 = []
decs_18 = []
for mask in masks_2018:
	ids, ras, decs = import_data(mask, zspeced=True)
	ids_18 = list(ids_18) + list(ids)
	ras_18 = list(ras_18) + list(ras)
	decs_18 = list(decs_18) + list(decs)
for mask in masks_2018_not_done:
	ids, ras, decs = import_data(mask, zspeced=False)
	ids_18 = list(ids_18) + list(ids)
	ras_18 = list(ras_18) + list(ras)
	decs_18 = list(decs_18) + list(decs)

ras_16 = []
ids_16 = []
decs_16 = []
for mask in masks_2016:
	ids, ras, decs = import_data(mask, zspeced=True)
	ids_16 = list(ids_16) + list(ids)
	ras_16 = list(ras_16) + list(ras)
	decs_16 = list(decs_16) + list(decs)

#remove duplicate IDs =========================================================================================================================================
unique_19, indices_19 = np.unique(ids_19, return_index=True) #returns an array of SORTED unique ID values and then the indices
ids_19 = np.array(ids_19)[indices_19] #grabbing the data for unique stars
ras_19 = np.array(ras_19)[indices_19]
decs_19 = np.array(decs_19)[indices_19]

unique_18, indices_18 = np.unique(ids_18, return_index=True) #returns an array of SORTED unique ID values and then the indices
ids_18 = np.array(ids_18)[indices_18] #grabbing the data for unique stars
ras_18 = np.array(ras_18)[indices_18]
decs_18 = np.array(decs_18)[indices_18]

unique_16, indices_16 = np.unique(ids_16, return_index=True) #returns an array of SORTED unique ID values and then the indices
ids_16 = np.array(ids_16)[indices_16] #grabbing the data for unique stars
ras_16 = np.array(ras_16)[indices_16]
decs_16 = np.array(decs_16)[indices_16]

#combining arrays of data that will be match to photometry ===================================================================================================
ids_phot = list(ids_19) + list(ids_18)
ras_phot = list(ras_19) + list(ras_18)
decs_phot = list(decs_19) + list(decs_18)

#data from target list (RA, Dec, and photometry) ==============================================================================================================
def import_phot(path):
	ref_data = np.genfromtxt('/Users/amandaquirk/Documents/M33/{}'.format(path), dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, magnitude2, filter2, priority, list_assignment, selection_flag, HST_F110W, HST_F160W, HST_F275W, HST_F336W')
	#pull out the fields that Laurent will need / I want for organizational purposes later
	ref_ID = ref_data['ID']
	ref_ID = [a.decode("utf-8") for a in ref_ID]
	ra = ref_data['ras']
	ra = [a.decode("utf-8") for a in ra]
	dec = ref_data['decs']
	dec = [a.decode("utf-8") for a in dec]
	mag2 = ref_data['magnitude2']

	return ref_ID, ra, dec, mag2

#bringing in the photometry from the 2018 and 2019 data
paths = ['Data/all_target_list.in', 'Data/target_list_RGB_2019.in', 'Masks/2019b/Dec/final/target_list_RGB_2019_expanded.in']

#the arrays below will contain many duplicates but there's no need to remove them
ref_IDs = []
ref_ras = []
ref_dec = []
ref_mag2 = []
for path in paths:
	ids, ras, decs, mags = import_phot(path)
	ref_IDs = list(ref_IDs) + list(ids)
	ref_ras = list(ref_ras) + list(ras)
	ref_dec = list(ref_dec) + list(decs)
	ref_mag2 = list(ref_mag2) + list(mags)

#match photometry to the zspec files ========================================================================================================================
nonserendip_g_filter = []
nonserendip_ID = []
nonserendip_ra = []
nonserendip_dec = []
filter_flag = []
#loop through so that if ID starts with CFHT, give appropriate filter (g) and if not, give other appropriate filter (F475W) 
for i in range(len(ids_phot)):
	if ids_phot[i].startswith('serendip') == False: 
		N = ref_IDs.index(ids_phot[i])
		nonserendip_g_filter.append(ref_mag2[N])
		nonserendip_ra.append(ref_ras[N])
		nonserendip_dec.append(ref_dec[N])
		nonserendip_ID.append(ref_IDs[N])
		if ids_phot[i].startswith('CFHT') == True:
			filter_flag.append('g')
		else:
			filter_flag.append('F475W')

#make nan list for 2016 data 
g_filter_16 = [np.nan for i in ras_16]
filter_flag_16 = [np.nan for i in ras_16]

#combining the data =========================================================================================================================================
all_ids = list(nonserendip_ID) + list(ids_16)
all_ras = list(nonserendip_ra) + list(ras_16)
all_decs = list(nonserendip_dec) + list(decs_16)
all_g = list(nonserendip_g_filter) + list(g_filter_16)
all_filters = list(filter_flag) + list(filter_flag_16)

#save it all to a file =====================================================================================================================================
np.savetxt('/Users/amandaquirk/Desktop/deimos_m33_targets.txt', np.c_[all_ids, all_ras, all_decs, all_g, all_filters], header='ID, ra, dec, mag, filter', fmt='%s')








