from astropy.io import fits
import numpy as np 
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as const
from astropy import units as u

'''
this adds together the photomerty and the spectroscopy for the 2018B masks

**right now this only pulls HST selected stars for AY9***
'''

def import_data(maskname): #import spectroscopy data from zspec files
	hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	data = hdu[1].data
	z = data['Z'] #redshift
	error = data['Z_ERR'] #error on the velocity
	zqual = data['ZQUALITY'] #eventually only want to use 1, 3, or 4
	aband = data['ABAND']
	ID = data['OBJNAME']
	time = data['MJD']
	return z, error, zqual, aband, ID, time   

z_a1, error_a1, zqual_a1, aband_a1, ID_a1, time_a1 = import_data('A1M33P')
z_b1, error_b1, zqual_b1, aband_b1, ID_b1, time_b1 = import_data('B1M33P')
z_c1, error_c1, zqual_c1, aband_c1, ID_c1, time_c1 = import_data('C1M33P')
z_d1, error_d1, zqual_d1, aband_d1, ID_d1, time_d1 = import_data('D1M33P')
z_e1, error_e1, zqual_e1, aband_e1, ID_e1, time_e1 = import_data('E1M33P')

#data from target list (RA, Dec, and photometry)
ref_data = np.genfromtxt('/Users/amandaquirk/Documents/M33/Data/all_target_list.in', dtype=None, names='ID, ras, decs, frame, magnitude1, filter1, magnitude2, filter2, priority, list_assignment, selection_flag, HST_F110W, HST_F160W, HST_F275W, HST_F336W')

ref_ID = ref_data['ID']
ref_ID = [a.decode("utf-8") for a in ref_ID]
ra = ref_data['ras']
ra = [a.decode("utf-8") for a in ra]
dec = ref_data['decs']
dec = [a.decode("utf-8") for a in dec]
mag1 = ref_data['magnitude1']
mag2 = ref_data['magnitude2']
mag3 = ref_data['HST_F110W']
mag4 = ref_data['HST_F160W']
mag5 = ref_data['HST_F275W']
mag6 = ref_data['HST_F336W']

#combining all of the zspec files together
zspec_IDs = list(ID_a1) + list(ID_b1) + list(ID_c1) + list(ID_d1) + list(ID_e1)
zspec_errors = list(error_a1) + list(error_b1) + list(error_c1) + list(error_d1) + list(error_e1)
zspec_zquals = list(zqual_a1) + list(zqual_b1) + list(zqual_c1) + list(zqual_d1) + list(zqual_e1)
zspec_aband = list(aband_a1) + list(aband_b1) + list(aband_c1) + list(aband_d1) + list(aband_e1)
zspec_zs = list(z_a1) + list(z_b1) + list(z_c1) + list(z_d1) + list(z_e1)
zspec_times = list(time_a1) + list(time_b1) + list(time_c1) + list(time_d1) + list(time_e1)

#matching the photometry to the spectroscopy -- ONLY FOR HST SELECTED STARS RIGHT NOW
F475W = []
F814W = []
F110W = []
F160W = []
F275W = []
F336W = []
redshift = []
zquals = []
errors = []
abands = []
RA = []
Dec = []
IDs = []
times = []
for i in range(len(zspec_IDs)):
	if zspec_IDs[i].startswith('serendip') == False: #later go in and add the CFHT data and the serendips
		N = ref_ID.index(zspec_IDs[i])
		F475W.append(mag2[N])
		F814W.append(mag1[N])
		F110W.append(mag3[N])
		F160W.append(mag4[N])
		F275W.append(mag5[N])
		F336W.append(mag6[N])
		RA.append(ra[N])
		Dec.append(dec[N])
		IDs.append(ref_ID[N])
		zquals.append(zspec_zquals[i])
		errors.append(zspec_errors[i])
		abands.append(zspec_aband[i])
		redshift.append(zspec_zs[i])
		times.append(zspec_times[i])
		#print(zspec_IDs[i], ref_ID[N])

#some plotting to test match===================================================================================================================
# import matplotlib.pyplot as plt 
# color = np.array(F475W) - np.array(F814W)
# zquals = np.array(zquals)
# F814W = np.array(F814W)
# good_data = (zquals == 1) | (zquals > 2)
# bad_data = (zquals == 2) | (zquals < 0)
# plt.scatter(color[good_data], F814W[good_data], c='orange', edgecolors='peru')
# plt.scatter(color[bad_data], F814W[bad_data], c='deeppink')#c=zquals, cmap=cmap, norm=norm, edgecolor='none', vmin=1, s=25, alpha=.9)
# plt.xlim(-2 ,8)
# plt.ylim(24.6, 13.5)
# plt.xlabel('F475W-F814W')
# plt.ylabel('F814W')
# plt.savefig('/Users/amandaquirk/Desktop/test_CMD.png')
# plt.close()

# from astropy import units as u 
# from astropy.coordinates import SkyCoord
# c = SkyCoord(RA, Dec, unit=(u.hourangle, u.hourangle))
# ra = c.ra.value 
# dec = c.dec.value
# plt.scatter(ra, dec)
# plt.savefig('/Users/amandaquirk/Desktop/test_map.png')
# plt.close()
#===========================================================================================================================================

'''
calculate the velocity from the redshift of the star, given by zspec
most of the code below comes from Karrie's Jupyter Notebook. she calculates the aband and heliocentric correction
for the repeated stars weighted average for their velocities
'''
# keck = EarthLocation.of_site('Keck')  # the easiest way... but requires internet
keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)

# make an M33 sky coordinate object
m33coord = SkyCoord.from_name('M33')
# and set systemic velocity of m33
m33_sys = -180. # pm 1 km/s; vanderMarel et al. 2008

# produce heliocentric correction
sc = SkyCoord(ra=RA, dec=Dec, unit=(u.hourangle, u.deg))

# radial_velocity_correction returns: The correction with a positive sign. 
# I.e., add this to an observed radial velocity 
# to get the barycentric (or heliocentric) velocity.
heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(times, format='mjd'), location=keck)  
heliocorr_km_s = heliocorr.to(u.km/u.s) 

vraw = redshift * const.c.to(u.km/u.s)
verr = errors * const.c.to(u.km/u.s)
vhelio = vraw + heliocorr_km_s 
vcorr = vraw + heliocorr_km_s - abands * const.c.to(u.km/u.s)

#save to file
#np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_phot_spec.txt', np.c_[IDs, RA, Dec, F275W, F336W, F475W, F814W, F110W, F160W, redshift, vcorr, verr, zquals, abands, times], fmt='%s', delimiter='\t', header='ID, RA, Dec, F275W, F336W, F475W, F814W, F110W, F160W, redshift, corrected vel (km/s), velocity error (km/s), zquality, A band, MJD') 

#file for Laurent to give me HI data
np.savetxt('/Users/amandaquirk/Desktop/M33_2018b_phot_spec.txt', np.c_[IDs, RA, Dec, vcorr], fmt='%s', delimiter='\t', header='ID, RA, Dec, vel') 

















