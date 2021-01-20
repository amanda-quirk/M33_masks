import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits 
import matplotlib.pyplot as plt 

#read in the data
def get_median_spec(mask):
	data = np.genfromtxt('/data/d/aquirk/2020b_reductions/qual_check/{}_out_selected.txt'.format(mask), names='id, ra, dec, 	frame, mag, fillter, priority, flag1, flag2', dtype=None)
	
	targ_id = data['id']
	ra = data['ra']
	ra = np.array([a.decode("utf-8") for a in ra])
	dec = data['dec']
	dec = np.array([a.decode("utf-8") for a in dec])
	mag = data['mag']
	slit_idx = np.linspace(0, len(targ_id), len(targ_id) + 1)
	print('Read in the data')
	
	#read in the data from the spec1d files
	print('Doing the median spec stuff')
	median_spec = []
	median_snr = []
	final_mag = []
	final_id = []
	final_ra = []
	final_dec = []
	chips = []
	for i in range(len(targ_id)):
		#print(targ_id[i])
		if slit_idx[i] < 10:
			filepath = '/data/d/aquirk/2020b_reductions/{}/spec1d.{}.00{}.{}.fits.gz'.format(mask, mask, int(slit_idx[i]), targ_id[i]) 
		elif (slit_idx[i] > 9) & (slit_idx[i] < 100):
			filepath = '/data/d/aquirk/2020b_reductions/{}/spec1d.{}.0{}.{}.fits.gz'.format(mask, mask, int(slit_idx[i]), targ_id[i]) 
		else:
			filepath = '/data/d/aquirk/2020b_reductions/{}/spec1d.{}.{}.{}.fits.gz'.format(mask, mask, int(slit_idx[i]), targ_id[i])  
		try:	
			with fits.open(filepath) as fits_data:
				chips.append(fits_data[1].header['CHIPNO'])
				blue_data = fits_data[1].data
				red_data = fits_data[2].data
				data_specs = blue_data.field('SPEC').tolist() + red_data.field('SPEC').tolist()
				data_lambdas = blue_data.field('LAMBDA').tolist() + red_data.field('LAMBDA').tolist()
				data_ivars = blue_data.field('IVAR').tolist() + red_data.field('IVAR').tolist()
				data_spec = [item for sublist in data_specs for item in sublist]
				data_lambda = [item for sublist in data_lambdas for item in sublist]
				data_ivar = [item for sublist in data_ivars for item in sublist]

				median_computation_spec = []
				median_computation_lambda = []
				SNR = []
				for k, j in enumerate(data_spec):
					if data_lambda[k] >= 6600 and data_lambda[k] <= 9000:
						median_computation_spec.append(data_spec[k])
						median_computation_lambda.append(data_lambda[k])
						SNR.append(data_spec[k] * np.sqrt(data_ivar[k]))
					else:
						continue
			median_spec.append(np.median(median_computation_spec))
			median_snr.append(np.median(SNR))
			final_mag.append(mag[i])
			final_id.append(int(targ_id[i]))
			final_ra.append(ra[i])
			final_dec.append(dec[i])
		except IOError:
			print(filepath)
			print('Slit skipped')
	print('Done!')
	
	print('Converting Coords')
	coords = SkyCoord(ra=final_ra, dec=final_dec, unit=(u.hourangle, u.deg))
	
	print('Saving Data')
	np.savetxt('median_spec_data_{}.txt'.format(mask), np.c_[final_id, coords.ra.deg, coords.dec.deg, median_spec, median_snr, final_mag, chips], header='id, ra (deg), dec (deg), median spec, median snr, mag, chip_num')
	return np.array(median_spec), np.array(median_snr), np.array(final_mag), np.array(chips)

def make_plots(mask, clouds):
	#read in data from function above
	data = np.genfromtxt('median_spec_data_{}.txt'.format(mask), names='final_id, ra, dec, median_spec, median_snr, final_mag, chip_num', dtype=None)
	mag = data['final_mag']
	spec = data['median_spec']
	snr = data['median_snr']
	chip = data['chip_num']
	#spec, snr, mag, chip = get_median_spec(mask)
	chip5 = np.array(chip == 1, dtype=bool) 

	#calculate the median offset from x = y
	y = 26 - 2.5 * np.log10(spec)
	offset = mag - y  
	print(mask, '=', np.nanmedian(offset))

	if clouds == True:
		plt.scatter(mag, y, c='b')
		plt.scatter(mag[chip5], y[chip5], c='r', label='CCD5')
	if clouds == False:
		plt.scatter(mag, y, edgecolors='b', facecolors='none')
		plt.scatter(mag[chip5], y[chip5], edgecolors='r', facecolors='none', label='CCD5')

	x = np.linspace(20, 24, 10)
	plt.plot(x,x, c='k')
	plt.ylabel(r'$26.5 -2.5 \times \log{\rm (Median\ Continuum)}$')
	plt.xlabel('i')
	plt.ylim(19, 26)
	plt.xlim(24, 20)
	plt.legend()
	plt.title('Mask = {}'.format(mask))
	plt.savefig('median_spec_mag_{}.png'.format(mask))
	plt.close()

	if clouds == True:
		plt.scatter(mag, np.log10(snr), c='b')
		plt.scatter(mag[chip5], np.log10(snr[chip5]), c='r', label='CCD5')
	if clouds == False:
		plt.scatter(mag, np.log10(snr), edgecolors='b', facecolors='none')
		plt.scatter(mag[chip5], np.log10(snr[chip5]), edgecolors='r', facecolors='none', label='CCD5')
	plt.ylabel(r'$\rm log(Median\ SNR)$')
	plt.xlabel('i')
	plt.xlim(24, 20)
	plt.ylim(-2.5, 1.2)
	plt.legend()
	plt.title('Mask = {}'.format(mask))
	plt.savefig('median_snr_mag_{}.png'.format(mask))
	plt.close()

#checking masks for trends with declination =============================================================================
def dec_offset(mask, pa, clouds):
	data = np.genfromtxt('median_spec_data_{}.txt'.format(mask), names='final_id, ra, dec, median_spec, median_snr, final_mag, chip_num', dtype=None)
	mag = data['final_mag']
	spec = data['median_spec']
	dec = data['dec']
	ra = data['ra']
	
	#calculate the offset from x = y
	y = 26 - 2.5 * np.log10(spec)
	offset = mag - y
	
	#sort by declinations 
	sorted_dec = np.array(sorted(dec))
	sorted_ra = np.array(sorted(ra))
	if pa == 90:
		sorted_offest = offset[np.argsort(ra)]
		#go by 15 and get median of dec and of offset; skipping the first and last element so this is easier
		med_ra = []
		med_offset = []
		x = 0
		while (x + 1) * 15 < len(mag):
			med_ra.append(np.nanmedian(sorted_ra[x * 15: (x+1) * 15]))
			med_offset.append(np.nanmedian(sorted_offest[x * 15: (x+1) * 15]))
			x= x + 1

		med_ra.append(np.nanmedian(sorted_ra[x * 15:]))
		med_offset.append(np.nanmedian(sorted_offest[x * 15:]))

		#plotting
		if clouds == True:
			plt.scatter(med_ra, med_offset, c='b')
		if clouds == False:
			plt.scatter(med_ra, med_offset, facecolors='none', edgecolors='b')
		plt.xlabel('Median RA')
		plt.gca().invert_xaxis()
		plt.ylim(-1.6, 1.1)
		plt.ylabel('Median Vertical Offset')
		plt.title('Mask = {}'.format(mask))
		plt.savefig('median_offset_ra_{}.png'.format(mask))
		plt.close()
	
	else:	
		sorted_offest = offset[np.argsort(dec)]
		#go by 15 and get median of dec and of offset; skipping the first and last element so this is easier
		med_dec = []
		med_offset = []
		x = 0
		while (x + 1) * 15 < len(mag):
			med_dec.append(np.nanmedian(sorted_dec[x * 15: (x+1) * 15]))
			med_offset.append(np.nanmedian(sorted_offest[x * 15: (x+1) * 15]))
			x= x + 1
	
		med_dec.append(np.nanmedian(sorted_dec[x * 15:]))
		med_offset.append(np.nanmedian(sorted_offest[x * 15:]))

		#plotting
		if clouds == True:
			plt.scatter(med_dec, med_offset, c='b')
		if clouds == False:
			plt.scatter(med_dec, med_offset, facecolors='none', edgecolors='b')
		plt.xlabel('Median Declination')
		plt.ylim(-1.6, 1.1)
		plt.ylabel('Median Vertical Offset')
		plt.title('Mask = {}'.format(mask))
		plt.savefig('median_offset_dec_{}.png'.format(mask))
		plt.close()

masks = ['pTN1b', 'pTN1a', 'pTN3', 'pTN4', 'pTS2', 'pTN2a', 'pTN5', 'pTE1', 'pTS1', 'pTN2b', 'pTS3']
pas = [22.5, 22.5, 90, 90, 90, 22.5, 90, 90, 22.5, 22.5, 90]
weather = [False, False, False, True, True, True, False, True, True, False, True]
for i in range(len(masks)):
	#print(masks[i])
	#make_plots(masks[i], weather[i])
	dec_offset(masks[i], pas[i], weather[i])



