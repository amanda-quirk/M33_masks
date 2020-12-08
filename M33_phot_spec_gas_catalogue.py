'''
Combines the photometry and spectroscopy for 2018b and 2019b M33 data. Also adds in all of the gas data from Laurent. Deals with duplicates. Corrects all velocities to heliocentric; does nothing with serendips

output: catalogue with IDs, RA, Dec, Z, V_corr w/o aband, w/ aband, V_err, Aband, MJD, Zqual, F275W, F336W, F475W, F814W, F110W, F160W, age tag, paired HI, paired CO, paired Halpha   
'''

from astropy.io import fits
import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy import units as u
from deprojecting_rotation_funcs import correct_vel

#import zspec data ==================================================================================================================================
def import_data(maskname): #import spectroscopy data from zspec files
	hdu = fits.open('/Volumes/Titan/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	data = hdu[1].data
	z = data['Z'] #redshift
	error = data['Z_ERR'] #error on the velocity
	zqual = data['ZQUALITY'] #eventually only want to use 1, 3, or 4
	aband = data['ABAND']
	ID = data['OBJNAME']
	time = data['MJD']
	ra = data['ra']
	dec = data['dec']
	mask = data['maskname']
	return z, error, zqual, aband, ID, time, ra, dec, mask 

#the 2016 masks
masks_2016 = ['M33D2A', 'M33D2B', 'M33D3A', 'M33D3B', 'M33D3D', 'M33D4A', 'M33D4B', 'M33D4C', 'M33MA1', 'M33MA2']

#the 2018 and 2019 masks
masks_2018 = ['A1M33P', 'B1M33P', 'C1M33P', 'D1M33P', 'E1M33P','A2M33P', 'B2M33P', 'E2M33P', 'K1M33P', 'C2M33P', 'D2M33P'] #photometry from M33/Data/all_target_list.in; not all of these masks have been zspeced yet so some are commented out
masks_2019 = ['D1M33R', 'D2M33R', 'E1M33R', 'E2M33R'] #E1M33R photometry from target_list_RGB_2019.in and rest photometry from M33/Masks/2019b/Dec/final/target_list_RGB_2019_expanded.in
mask_list = masks_2016 + masks_2018 + masks_2019

#read in the data into one list for each parameter
all_ids = []
all_ras = []
all_decs = []
all_zs = []
all_errs = []
all_zquals = []
all_abands = []
all_times = []
all_masks = []
for name in mask_list:
	print(name)
	zs, errors, zquals, abands, IDs, times, ras, decs, masks = import_data(name)
	all_ids = all_ids + list(IDs)
	all_ras = all_ras + list(ras)
	all_decs = all_decs + list(decs)
	all_zs = all_zs + list(zs)
	all_errs = all_errs + list(errors)
	all_zquals = all_zquals + list(zquals)
	all_abands = all_abands + list(abands)
	all_times = all_times + list(times)
	all_masks = all_masks + list(masks)

#shift the HST astrometry for the 2016b masks according to what Karrie did when making the masks -- this puts them on the CFHT reference frame
corrected_ras = np.zeros_like(all_ras)
corrected_decs = np.zeros_like(all_decs)
for i in range(len(all_ids)):
	if len(all_ids[i]) > 6: #not 2016 CFHT
		if all_ids[i].startswith('5'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - 0.043 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec + 0.37 * u.arcsecond
		elif all_ids[i].startswith('6'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra + 0.0255 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.18 * u.arcsecond
		elif all_ids[i].startswith('7'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - 0.0285 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.235 * u.arcsecond
		elif all_ids[i].startswith('8'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - -0.035 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec + 0.09 * u.arcsecond
		elif all_ids[i].startswith('9'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra + 0.1 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.28 * u.arcsecond
		elif all_ids[i].startswith('10'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra + 0.036 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.27 * u.arcsecond
		elif all_ids[i].startswith('11'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - 0.029 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec + 0.09 * u.arcsecond
		elif all_ids[i].startswith('12'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra + 0.022 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.37 * u.arcsecond
		elif all_ids[i].startswith('13'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - 0.043 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.28 * u.arcsecond
		elif all_ids[i].startswith('14'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - 0.0285 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec + 0.00 * u.arcsecond
		elif all_ids[i].startswith('15'):
			coord = SkyCoord(all_ras[i], all_decs[i], unit=(u.hourangle,u.deg))
			ra_shift = coord.ra - 0.05 * u.arcsecond
			corrected_ras[i] = ra_shift.hourangle
			corrected_decs[i] = coord.dec - 0.185 * u.arcsecond
		else: #2018 or 2019 data, needs no offset
			corrected_ras[i] = all_ras[i]
			corrected_decs[i] = all_decs[i]
	else: #2016 CFHT data, needs no offset
		corrected_ras[i] = all_ras[i]
		corrected_decs[i] = all_decs[i]

all_ras = corrected_ras #there will be different formats in here but all are the same unit so ok
all_decs = corrected_decs #there will be different formats in here but all are the same unit so ok

#deal with duplicate stars ==========================================================================================================================
#first, let's remove all serendips
not_serendip = (np.array(all_ids) != 'serendip1') & (np.array(all_ids) != 'serendip2')
all_ids = np.array(all_ids)[not_serendip]
all_ras = np.array(all_ras)[not_serendip]
all_decs = np.array(all_decs)[not_serendip]
all_zs = np.array(all_zs)[not_serendip]
all_errs = np.array(all_errs)[not_serendip] * const.c.to(u.km/u.s) #converting to km/s
all_zquals = np.array(all_zquals)[not_serendip]
all_abands = np.array(all_abands)[not_serendip]
all_times = np.array(all_times)[not_serendip]
all_masks = np.array(all_masks)[not_serendip]

#now, let's figure out which stars were observed once and those observed more than once 
dup_ind = [] #list of arrays
first_dup = [] #the index of the first entry to use for parameters not observation dependent
unique_ind = [] #none of these stars have repeated observations
names = [] #don't want to double count our duplicates lol 
for i in range(len(all_ids)):
	if all_ids[i] not in names: #don't want to go through the other entries of a duplicate after it's been counted the first time 
		N = np.where(all_ids == all_ids[i])
		if len(N[0]) > 1:
			dup_ind.append(N[0])
			first_dup.append(N[0][0])
		elif len(N[0]) == 1:
			unique_ind.append(N[0][0])
		else:
			print('yikes')
		names.append(all_ids[i])

#data for all of the stars observed once
unique_ras = all_ras[unique_ind]
unique_decs = all_decs[unique_ind]
unique_zs = all_zs[unique_ind]
unique_errs = all_errs[unique_ind]
unique_zquals = all_zquals[unique_ind]
unique_abands = all_abands[unique_ind]
unique_times = all_times[unique_ind]
unique_masks = all_masks[unique_ind]
unique_IDs = all_ids[unique_ind]

#data that is observationally independent for the repeated stars 
dup_ras = all_ras[first_dup]
dup_decs = all_decs[first_dup]
dup_IDs = all_ids[first_dup]

#duplicates -- if the z qualities are similar, we want to take the weighted average of the corrected velocities; if they aren't similar, we take the data from the high zquality 
#first, correct the velocities
#real quick do the easy case
unique_vels = correct_vel(unique_ras, unique_decs, unique_times, unique_zs, unique_abands)
unique_vels_aband = correct_vel(unique_ras, unique_decs, unique_times, unique_zs, unique_abands, apply_aband=True)
print('Done correcting unique velocities')

#now with the duplicates -- going to go through them individually 
def calc_weights(err):
        return 1 / (err**2)

def normed_weight(w):
        sum_weights=sum(w)
        return w / sum_weights

#function does the weighted meean
def weighted_mean(data,norm_w):
	return sum(data * norm_w)

def weighted_mean_error(err):
	addition = sum(calc_weights(err))
	return np.sqrt(1 / addition)

dup_z = []
dup_err = []
dup_zqual = []
dup_aband = []
dup_time = []
dup_vel = []
dup_vel_aband = []
dup_mask = [] #if I take the weighted average of velocities, I just combine the mask names together
count_zquals_close = 0
count_zquals_not_close_2 = 0
count_zquals_not_close_3_1 = 0
count_zquals_not_close_3_2 = 0
for i in range(len(dup_ind)):
	#evaluate zqualities first
	zquals = all_zquals[dup_ind[i]]
	std = np.std(zquals)
	maskname = ''
	if std < 0.7: #zquals are close enough that we will take the weighted avg of the velocities
		vels = []
		vels_aband = []
		for ind in dup_ind[i]: #go through each entry of a duplicate
			vels.append(correct_vel(all_ras[ind], all_decs[ind], all_times[ind], all_zs[ind], all_abands[ind]))
			vels_aband.append(correct_vel(all_ras[ind], all_decs[ind], all_times[ind], all_zs[ind], all_abands[ind], apply_aband=True))
			maskname += all_masks[ind] 
		weights = calc_weights(all_errs[dup_ind[i]])
		weights_norm = normed_weight(weights)
		dup_vel.append(weighted_mean(vels, weights_norm))
		dup_vel_aband.append(weighted_mean(vels_aband, weights_norm))
		dup_zqual.append(max(zquals))
		dup_time.append('weight_avg_dup')
		dup_aband.append(weighted_mean(all_abands[dup_ind[i]], weights_norm)) #changed this and the next line from weight_avg_dup
		dup_z.append(weighted_mean(all_zs[dup_ind[i]], weights_norm))
		dup_err.append(weighted_mean_error(all_errs[dup_ind[i]]))
		dup_mask.append(maskname)
		count_zquals_close += 1
	else: #taking the data from the max zquality entry since duplicate entries are of much worse quality
		if len(zquals) == 2:
			where_best = np.argmax(zquals) 
			best_ind = dup_ind[i][where_best]
			dup_z.append(all_zs[best_ind])
			dup_err.append(all_errs[best_ind])
			dup_zqual.append(all_zquals[best_ind])
			dup_aband.append(all_abands[best_ind])
			dup_time.append(all_times[best_ind])
			dup_mask.append(all_masks[best_ind])
			dup_vel.append(correct_vel(all_ras[best_ind], all_decs[best_ind], all_times[best_ind], all_zs[best_ind], all_abands[best_ind]))
			dup_vel_aband.append(correct_vel(all_ras[best_ind], all_decs[best_ind], all_times[best_ind], all_zs[best_ind], all_abands[best_ind], apply_aband=True))
			count_zquals_not_close_2 += 1
		else: #there are 3 or more observations of the same star, take data from the highest zqual
			best = np.max(zquals)
			where_best = np.argmax(zquals)
			N = np.where(zquals == best)
			if len(N[0]) == 1: #only one best zqual, just take that data
				best_ind = dup_ind[i][where_best]
				dup_z.append(all_zs[best_ind])
				dup_err.append(all_errs[best_ind])
				dup_zqual.append(all_zquals[best_ind])
				dup_aband.append(all_abands[best_ind])
				dup_time.append(all_times[best_ind])
				dup_mask.append(all_masks[best_ind])
				dup_vel.append(correct_vel(all_ras[best_ind], all_decs[best_ind], all_times[best_ind], all_zs[best_ind], all_abands[best_ind]))
				dup_vel_aband.append(correct_vel(all_ras[best_ind], all_decs[best_ind], all_times[best_ind], all_zs[best_ind], all_abands[best_ind], apply_aband=True))
				count_zquals_not_close_3_1 += 1
			else: #there are multiple highest zqual entries so thake their weighted avg
				vels = []
				vels_aband = []
				real_inds = []
				maskname = ''
				for ind in N[0]:
					real_ind = dup_ind[i][ind]
					real_inds.append(dup_ind[i][ind])
					vels.append(correct_vel(all_ras[real_ind], all_decs[real_ind], all_times[real_ind], all_zs[real_ind], all_abands[real_ind]))
					vels_aband.append(correct_vel(all_ras[real_ind], all_decs[real_ind], all_times[real_ind], all_zs[real_ind], all_abands[real_ind], apply_aband=True))
					maskname += all_masks[real_ind]
				weights = calc_weights(all_errs[real_inds])
				weights_norm = normed_weight(weights)
				dup_vel.append(weighted_mean(vels, weights_norm))
				dup_vel_aband.append(weighted_mean(vels_aband, weights_norm))
				dup_zqual.append(max(zquals))
				dup_time.append('weight_avg_dup')
				dup_aband.append(weighted_mean(all_abands[real_inds], weights_norm)) #changed this and the next line from 'weight_avg_dup'
				dup_z.append(weighted_mean(all_zs[real_inds], weights_norm))
				dup_err.append(weighted_mean_error(all_errs[real_inds]))
				dup_mask.append(maskname)
				count_zquals_not_close_3_2 += 1
print('Done with dup velocities', count_zquals_close, count_zquals_not_close_2, count_zquals_not_close_3_1, count_zquals_not_close_3_2)

#combing arrays for dups and nondups ======================================================================================================
sorted_ras = list(unique_ras) + list(dup_ras)
sorted_decs = list(unique_decs) + list(dup_decs)
sorted_zs = list(unique_zs) + list(dup_z)
sorted_errs = list(unique_errs) + list(dup_err)
sorted_zquals = list(unique_zquals) + list(dup_zqual)
sorted_abands = list(unique_abands) + list(dup_aband)
sorted_times = list(unique_times) + list(dup_time)
sorted_vels = list(unique_vels) + list(dup_vel)
sorted_vels_aband = list(unique_vels_aband) + list(dup_vel_aband)
sorted_ids = list(unique_IDs) + list(dup_IDs)
sorted_masks = list(unique_masks) + list(dup_mask)

#need to remove the 0 that is in front of some of the IDs
fixed_ids = []
for name in sorted_ids:
	if name.startswith('0') == True:
		fixed_ids.append(name[1:])
	else:
		fixed_ids.append(name)
sorted_ids = fixed_ids

#pulling in the photometry from the target lists ==========================================================================================
def import_phot(path):
	if path.endswith('input.fits'):
		h = fits.open('/Volumes/Titan/M33/Data/{}'.format(path), memmap = True)
		ref_data = h[1].data
		ref_ID = ref_data['OBJNAME']
		ref_ID = [int(a) for a in ref_ID] #converting from string to int to get rid of extra space
		ref_ID = [str(a) for a in ref_ID] #converting it back -.- so that it matches sorted_ids
		mag3 = ref_data['MAG1_ACS'] #F475W OR F606W
		mag4 = ref_data['MAG2_ACS'] #F814W
		return ref_ID, np.zeros(len(mag3)), np.zeros(len(mag3)), mag3, mag4, np.zeros(len(mag3)), np.zeros(len(mag3))
	elif path.endswith('matchcat.fits'): #extcorr i/o matchcat
		h = fits.open('/Volumes/Titan/M33/Data/{}'.format(path), memmap = True)
		ref_data = h[1].data
		ref_ID = ref_data['OBJNO1'] 
		ref_ID = [int(a) for a in ref_ID] #converting from string to int to get rid of extra space
		ref_ID = [str(a) for a in ref_ID] #converting it back -.- so that it matches sorted_ids
		mag3 = ref_data['MAG1_AUTO'] #g #_0 for the extcorr file
		mag4 = ref_data['MAG2_AUTO'] #i #_0 for the extcorr file 
		return ref_ID, np.zeros(len(mag3)), np.zeros(len(mag3)), mag3, mag4, np.zeros(len(mag3)), np.zeros(len(mag3))
	else:
		ref_data = np.genfromtxt('/Volumes/Titan/M33/Data/{}'.format(path), dtype=None, names='ID, ras, decs, frame, F814W, filter1, 	F475W, filter2, priority, list_assignment, selection_flag, F110W, F160W, F275W, F336W')
		#pull out the fields that Laurent will need / I want for organizational purposes later
		ref_ID = ref_data['ID']
		ref_ID = [a.decode("utf-8") for a in ref_ID]
		mag1 = ref_data['F275W']
		mag2 = ref_data['F336W']
		mag3 = ref_data['F475W']
		mag4 = ref_data['F814W']
		mag5 = ref_data['F110W']
		mag6 = ref_data['F160W']
		return ref_ID, mag1, mag2, mag3, mag4, mag5, mag6

#bringing in the photometry from the 2018 and 2019 data
print('Reading in photometry')
paths = ['all_target_list.in', 'target_list_RGB_2019_6filt.in', 'target_list_RGB_2019_expanded_6filt.in', 'HST_all_list3.input.fits', 'M33.GI.matchcat.fits', 'HST_15b_all_nolist3dups.input.fits', 'HST_16b_all_nolist3dups.input.fits'] #need to swap out the current M33 file for M33.GI.matchcat.fits once I get it from Karrie

#the arrays below will contain many duplicates but there's no need to remove them
ref_IDs = []
F275Ws = []
F336Ws = []
F475Ws = []
F814Ws = []
F110Ws = []
F160Ws = []
for path in paths:
	ref_ids, mag1s, mag2s, mag3s, mag4s, mag5s, mag6s = import_phot(path)
	ref_IDs = ref_IDs + list(ref_ids)
	F275Ws = F275Ws + list(mag1s)
	F336Ws = F336Ws + list(mag2s)
	F475Ws = F475Ws + list(mag3s)
	F814Ws = F814Ws + list(mag4s)
	F110Ws = F110Ws + list(mag5s)
	F160Ws = F160Ws + list(mag6s)
print('Photometry loaded')

#match photometry to the zspec files =======================================================================================================
age_tag = []
sorted_F275W = []
sorted_F336W = []
sorted_F475W = []
sorted_F814W = []
sorted_F110W = []
sorted_F160W = []
unmatched_ID = []
for i in range(len(sorted_ids)):
	try:
		N = ref_IDs.index(sorted_ids[i])
		sorted_F275W.append(F275Ws[N]) 
		sorted_F336W.append(F336Ws[N]) 
		sorted_F475W.append(F475Ws[N]) 
		sorted_F814W.append(F814Ws[N]) 
		sorted_F110W.append(F110Ws[N]) 
		sorted_F160W.append(F160Ws[N])
		#adding an age tag for organizational purposes
		if sorted_ids[i].startswith('CFHT'):
			age_tag.append('CFHT')
		elif sorted_ids[i].startswith('AGB'):
			age_tag.append('AGB')
		elif sorted_ids[i].startswith('BHeB'):
			age_tag.append('BHeB')
		elif sorted_ids[i].startswith('RHeB'):
			age_tag.append('RHeB') 
		elif sorted_ids[i].startswith('RGB'):
			age_tag.append('RGB')
		elif sorted_ids[i].startswith('MS'):
			age_tag.append('MS')
		elif sorted_ids[i].startswith('HeB'):
			age_tag.append('HeB')
		elif sorted_ids[i].startswith('KG'):
			age_tag.append('Xray')
		else: #2016 data is just numbers. Sources with 7 or 8 digit ID numbers are from the HST catalogs. Sources with < 7 digit ID 	numbers are from the CFHT catalog 
			age_tag.append('unknown')
	except ValueError:
		unmatched_ID.append(sorted_ids[i])
print('The number of unmatched stars is', len(unmatched_ID))
print('Sorted by ages')

#add in the gas ============================================================================================================================
#THIS WILL BE INCOMPLETE UNTIL I GET THE NEW FILE FROM LAURENT
gas_data = np.genfromtxt('/Volumes/Titan/M33/Data/velocity-m33-hi-co-ha.ascii', dtype=None, names='ID, ras, decs, mag, filter, HI, CO, Ha')
gas_IDs = gas_data['ID']
gas_IDs = [a.decode("utf-8") for a in gas_IDs]
gas_HI = gas_data['HI'] #LSR frame
gas_CO = gas_data['CO'] #LSR frame
gas_Ha = gas_data['Ha'] #heliocentric frame

#match to stars
sorted_HI = []
sorted_CO = []
sorted_Ha = []
missing_gas_count = 0
for i in range(len(sorted_ids)):
	if sorted_ids[i] in gas_IDs:
 		N = gas_IDs.index(sorted_ids[i])
 		if gas_HI[N] < 900:
 			sorted_HI.append(gas_HI[N])
 		elif gas_HI[N] > 900: #Laurent put in 999 for values he doesn't have a measurement 
 			sorted_HI.append(np.nan)
 		if gas_CO[N] < 900:
 			sorted_CO.append(gas_CO[N])
 		elif gas_CO[N] > 900:
 			sorted_CO.append(np.nan)
 		if gas_Ha[N] < 900:
 			sorted_Ha.append(gas_Ha[N])
 		elif gas_Ha[N] > 900:
 			sorted_Ha.append(np.nan)
	else:
		missing_gas_count += 1
		sorted_HI.append(np.nan)
		sorted_CO.append(np.nan)
		sorted_Ha.append(np.nan)
print('Number of stars missing gas data is', missing_gas_count, 'out of', len(sorted_ids)) #hopefully this becomes 0 once all the data is in

#correct the reference frame for HI and CO -- LSR to heliocentric; this is a rough estimate for now
def LSR_to_helio(v, ra, dec): #I am still sketpcial of this formula
	v_corr = np.zeros(len(v))
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
	for i in range(len(v)):
		if v[i] == np.nan:
			v_corr[i] = v[i]
		else:
			#v_corr[i] = v[i] + 0.59 
			l = sc[i].galactic.l.value
			b = sc[i].galactic.b.value
			v_corr[i] = v[i] - 9 * np.cos(l) * np.cos(b) - 12 * np.sin(l) * np.cos(b) - 7 * np.sin(b)
	return v_corr #km/s

sorted_corrected_HI = LSR_to_helio(sorted_HI, sorted_ras, sorted_decs)
sorted_corrected_CO = LSR_to_helio(sorted_CO, sorted_ras, sorted_decs)

#saving the catalogue!! =====================================================================================================================
sorted_errs = [a.value for a in sorted_errs] #strip the unit

np.savetxt('/Volumes/Titan/M33/Data/M33_2018b_phot_spec.txt', np.c_[sorted_ids, sorted_ras, sorted_decs, sorted_F275W, sorted_F336W, sorted_F475W, sorted_F814W, sorted_F110W, sorted_F160W, sorted_zs, sorted_vels, sorted_vels_aband, sorted_errs, sorted_zquals, sorted_abands, sorted_times, sorted_masks, age_tag, sorted_corrected_HI, sorted_corrected_CO, sorted_Ha], fmt='%s', delimiter='\t', header='ID, RA, Dec, F275W, F336W, F475W/g/F606W, F814W/i, F110W, F160W, redshift, heliocorrected vel (km/s), helio+aband corrected vel (km/s), velocity error (km/s), zquality, A band, MJD, mask name, age tag, HI (km/s), CO (km/s), Halpha (km/s)') 

#np.savetxt('/Users/amandaquirk/Desktop/M33_all_coords.txt', np.c_[sorted_ids, sorted_ras, sorted_decs], fmt='%s', delimiter='\t', header='ID, RA, Dec') 













