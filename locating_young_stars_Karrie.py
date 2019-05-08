import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as const
from astropy.table import vstack
from matplotlib import rc 

'''
first, let's read in all of the data: the IDs of the carbon and weak CN stars in the 16B and 18B masks and the zspec files. we care about the redshift, verr, zquality, maskname, objectname, ra, dec,  aband, and slitname
'''
#I have removed the serendips from this file by hand -- see yong_CMD.py for when I had python ignore these stars
young_star_mask_ref, young_stars_ID_ref = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/carbon_weak_CN_M33_all.txt', usecols=(0,1,), dtype='str', unpack=True)

def import_data(maskname):
	hdu = fits.open('/Users/amandaquirk/Documents/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	data = Table(hdu[1].data)
	ra = data['RA']
	dec = data['DEC']
	z = data['Z']
	error = data['Z_ERR']
	zqual = data['ZQUALITY']
	slit = data['SLITNAME']
	aband = data['ABAND']
	ID = data['OBJNAME']
	mask = data['MASKNAME']
	time = data['DATE']

	return ra, dec, z, error, zqual, slit, aband, ID, mask, time  

data1 = import_data('M33D2A')
data2 = import_data('M33D2B')
data3 = import_data('M33D3A')
data4 = import_data('M33D3B')
data5 = import_data('M33D3D')
data6 = import_data('M33D4A')
data7 = import_data('M33D4B')
data8 = import_data('A1M33P')
data9 = import_data('B1M33P')

'''
below grabs the data for the carbon and weak cn stars 
sanity check that ID is matched with correct star: check the masks
cut out any stars that don't have a zquality of 3 or 4
'''
#combine the data arrays
all_ra = list(data1[0]) + list(data2[0]) + list(data3[0]) + list(data4[0]) + list(data5[0]) + list(data6[0]) + list(data7[0]) + list(data8[0]) + list(data9[0])
all_dec = list(data1[1]) + list(data2[1]) + list(data3[1]) + list(data4[1]) + list(data5[1]) + list(data6[1]) + list(data7[1]) + list(data8[1]) + list(data9[1])
all_z = list(data1[2]) + list(data2[2]) + list(data3[2]) + list(data4[2]) + list(data5[2]) + list(data6[2]) + list(data7[2]) + list(data8[2]) + list(data9[2])
all_error = list(data1[3]) + list(data2[3]) + list(data3[3]) + list(data4[3]) + list(data5[3]) + list(data6[3]) + list(data7[3]) + list(data8[3]) + list(data9[3])
all_zqual = list(data1[4]) + list(data2[4]) + list(data3[4]) + list(data4[4]) + list(data5[4]) + list(data6[4]) + list(data7[4]) + list(data8[4]) + list(data9[4])
all_slit = list(data1[5]) + list(data2[5]) + list(data3[5]) + list(data4[5]) + list(data5[5]) + list(data6[5]) + list(data7[5]) + list(data8[5]) + list(data9[5])
all_aband = list(data1[6]) + list(data2[6]) + list(data3[6]) + list(data4[6]) + list(data5[6]) + list(data6[6]) + list(data7[6]) + list(data8[6]) + list(data9[6])
all_ID = list(data1[7]) + list(data2[7]) + list(data3[7]) + list(data4[7]) + list(data5[7]) + list(data6[7]) + list(data7[7]) + list(data8[7]) + list(data9[7])
all_mask = list(data1[8]) + list(data2[8]) + list(data3[8]) + list(data4[8]) + list(data5[8]) + list(data6[8]) + list(data7[8]) + list(data8[8]) + list(data9[8])
all_time = list(data1[9]) + list(data2[9]) + list(data3[9]) + list(data4[9]) + list(data5[9]) + list(data6[9]) + list(data7[9]) + list(data8[9]) + list(data9[9])

#get the data for the young stars in above list
#remove extra spaces from IDs in all_ID
all_ID = [a.strip() for a in all_ID]

young_ra = np.zeros_like(young_stars_ID_ref) #strings
young_dec = np.zeros_like(young_stars_ID_ref)
young_z = np.zeros(len(young_stars_ID_ref)) #not strongs
young_error = np.zeros(len(young_stars_ID_ref))
young_zqual = np.zeros(len(young_stars_ID_ref))
young_slit = np.zeros(len(young_stars_ID_ref))
young_aband = np.zeros(len(young_stars_ID_ref))
young_ID = np.zeros_like(young_stars_ID_ref)
young_mask = np.zeros_like(young_stars_ID_ref)
young_time = np.zeros_like(young_stars_ID_ref) #strings
for i in range(len(young_stars_ID_ref)):
	N = all_ID.index(young_stars_ID_ref[i])
	young_ra[i] = all_ra[N]
	young_dec[i] = all_dec[N]
	young_z[i] = all_z[N]
	young_error[i] = all_error[N]
	young_zqual[i] = all_zqual[N]
	young_slit[i] = all_slit[N]
	young_aband[i] = all_aband[N]
	young_ID[i] = all_ID[N]
	young_mask[i] = all_mask[N]
	young_time[i] = all_time[N]
	#print(young_stars_ID_ref[i], all_ID[N])
	#print(young_star_mask_ref[i], all_mask[N])
	#if young_star_mask_ref[i] != all_mask[N]:
	#	print(young_stars_ID_ref[i], all_ID[N], young_star_mask_ref[i], all_mask[N])

#only want high zquality
good_quality = (young_zqual == 3) | (young_zqual == 4)
young_ra = young_ra[good_quality] 
young_dec = young_dec[good_quality]
young_z = young_z[good_quality]
young_error = young_error[good_quality]
young_zqual = young_zqual[good_quality]
young_slit = young_slit[good_quality]
young_aband = young_aband[good_quality]
young_ID = young_ID[good_quality]
young_mask = young_mask[good_quality]
young_time = young_time[good_quality]

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
sc = SkyCoord(ra=young_ra, dec=young_dec, unit=(u.hourangle, u.deg))

# radial_velocity_correction returns: The correction with a positive sign. 
# I.e., add this to an observed radial velocity 
# to get the barycentric (or heliocentric) velocity.
heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(young_time), location=keck)  
heliocorr_km_s = heliocorr.to(u.km/u.s) 

vraw = young_z * const.c.to(u.km/u.s)
verr = young_error * const.c.to(u.km/u.s)
vhelio = vraw + heliocorr_km_s 
vcorr = vraw + heliocorr_km_s - young_aband * const.c.to(u.km/u.s)

#replace duplicates with the weighted average of their corrected velocity
duplicates = [x for x in list(young_ID) if list(young_ID).count(x) > 1] #ID of repeated stars
duplicates2 = [i for i, x in enumerate(list(young_ID)) if list(young_ID).count(x) > 1] #indices of repeated stars

ra_unique = []
dec_unique = []
error_unique = []
slit_unique = []
mask_unique = []
ID_unique = []
vcorr_unique = []
#non duplicates
for i in range(len(young_ID)):
	if i in duplicates2:
		print(i)
	else:
		ra_unique.append(sc.ra[i].value)
		dec_unique.append(sc.dec[i].value)
		error_unique.append(verr[i].value)
		slit_unique.append(young_slit[i])
		mask_unique.append(young_mask[i])
		ID_unique.append(young_ID[i])
		vcorr_unique.append(vcorr[i].value)

#duplicates
duplicated_stars = list(set(duplicates)) #each duplicated star only listed once
for i in range(len(duplicated_stars)):
	N = np.where(duplicated_stars[i] == young_ID) #N will be an array of the indices of everywhere a duplicate star appears
	#calculate the weighted mean of the velocity and the standard error
	errors = verr[N]
	velocities = vcorr[N]
	weighted_vel_sum = sum(velocities / errors**2)
	weight_sum = sum(1 / errors**2)
	weighted_mean_vel = weighted_vel_sum / weight_sum
	standard_error = np.sqrt(1 / weight_sum)
	#print(velocities)
	#print(weighted_mean_vel)
	ra_unique.append(sc.ra[N[0][0]].value) #just looking at the first coordinate value since will be the same for all entries
	dec_unique.append(sc.dec[N[0][0]].value)
	error_unique.append(standard_error.value)
	slit_unique.append(young_slit[N[0][0]]) #will just have the slit of the first mask its on
	mask_unique.append(young_mask[N[0][0]]) #will just have listed the first mask its on
	ID_unique.append(duplicated_stars[i])
	vcorr_unique.append(weighted_mean_vel.value)

'''
outputs: data file and velocity histogram
'''
#spatial map for sanity check
# plt.scatter(ra_unique, dec_unique, c =vcorr_unique, vmin=-620, vmax=620)
# plt.xlim(23.7, 23.1)
# plt.colorbar()
# plt.show()

#velocity histogram
def single_plot():
	rc('font', family = 'serif')
	fig, ax=plt.subplots(1)
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(2)
	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	plt.tick_params(which='both', width=2)
	plt.tick_params(which='major', length=7)
	plt.tick_params(which='minor', length=4)
	plt.tick_params(labelsize=12) 
	plt.minorticks_on()

single_plot()
plt.hist(vcorr_unique, bins=20)
plt.plot([m33_sys, m33_sys], [0,50], color = 'darkgrey', linestyle = '--')
plt.xlabel('Velocity (km/s)')
plt.ylim(0, 40)
plt.savefig('/Users/amandaquirk/Desktop/weakCN_C_velocities.png', bbox_inches='tight')
plt.close()

#datafile
np.savetxt('/Users/amandaquirk/Desktop/M33_weakCN_C_stars.txt', np.c_[ID_unique, ra_unique, dec_unique, vcorr_unique, mask_unique, slit_unique], fmt='%s', delimiter=' ', header='ID, RA (deg), Dec (deg), corrected LOS vel (km/s), Mask, Slit')

#fits file
c1 = fits.Column(name='ID', array=ID_unique, format='16A')
c2 = fits.Column(name='RA_deg', array=ra_unique, format='D')
c3 = fits.Column(name='Dec_deg', array=dec_unique, format='D')
c4 = fits.Column(name='V_km_s', array=vcorr_unique, format='D')
c5 = fits.Column(name='Mask', array=mask_unique, format='6A')
c6 = fits.Column(name='Slit', array=slit_unique, format='K')
data_table = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6])
data_table.writeto('/Users/amandaquirk/Desktop/M33_weakCN_C_stars.fits')



