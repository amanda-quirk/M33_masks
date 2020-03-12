'''
makes plot for 2020B keck proposal - histogram of velocities of minor axis masks
'''

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt 
from matplotlib import rc 
import numpy as np 

#import zspec data ==================================================================================================================================
def import_data(maskname): #import spectroscopy data from zspec files
	hdu = fits.open('/Volumes/Titan/M33/Data/zspecs/zspec.{}.fits'.format(maskname))
	data = hdu[1].data
	z = data['Z'] #redshift
	slit = data['slitname']
	zqual = data['ZQUALITY'] #eventually only want to use 1, 3, or 4
	aband = data['ABAND']
	time = data['MJD']
	ra = data['ra']
	dec = data['dec']
	ID = data['OBJNAME']
	return z, slit, ra, dec, time, aband, zqual, ID

#mask details
minor_axis_masks = ['M33MA1', 'M33MA2']
MW_stars_A1 = [0, 7, 25, 41, 88, 107, 122, 144, 145, 197, 198, 202, 204, 205] #slit numbers that Raja flagged as MWFG stars; from Raja's notes
MW_stars_A2 = [24, 57, 219, 1, 182, 200, 203, 219, 158, 61, 37, 58, 36, 63, 70, 0, 7, 25, 41, 88, 107, 122, 144, 145, 197, 198, 202, 204, 205, 13, 44, 52, 140, 150, 196, 197, 198, 199, 201, 202, 205, 29, 62, 63, 81, 103, 142, 189, 73, 48, 114, 46, 78, 151, 200, 0, 39, 151, 160, 8, 17, 32, 100, 135, 142, 64, 2, 132]

#read in mask data
A1_data = import_data(minor_axis_masks[0])
A2_data = import_data(minor_axis_masks[1])

#eliminate MWFG and serendips
A1_z = []
A1_ra = []
A1_dec = []
A1_time = []
A1_aband = []
A1_zqual = []
MW_A1_z = []
MW_A1_ra = []
MW_A1_dec = []
MW_A1_time = []
MW_A1_aband = []
MW_A1_zqual = []
for i in range(len(A1_data[0])):
	if int(A1_data[1][i]) not in MW_stars_A1 and A1_data[7][i] != 'serendip1' and A1_data[7][i] != 'serendip2' and A1_data[7][i] != 'serendip3':
		A1_z.append(A1_data[0][i])
		A1_ra.append(A1_data[2][i])
		A1_dec.append(A1_data[3][i])
		A1_time.append(A1_data[4][i])
		A1_aband.append(A1_data[5][i])
		A1_zqual.append(A1_data[6][i])
	elif int(A1_data[1][i]) in MW_stars_A1 and A1_data[7][i] != 'serendip1' and A1_data[7][i] != 'serendip2' and A1_data[7][i] != 'serendip3':
		MW_A1_z.append(A1_data[0][i])
		MW_A1_ra.append(A1_data[2][i])
		MW_A1_dec.append(A1_data[3][i])
		MW_A1_time.append(A1_data[4][i])
		MW_A1_aband.append(A1_data[5][i])
		MW_A1_zqual.append(A1_data[6][i])

A2_z = []
A2_ra = []
A2_dec = []
A2_time = []
A2_aband = []
A2_zqual = []
MW_A2_z = []
MW_A2_ra = []
MW_A2_dec = []
MW_A2_time = []
MW_A2_aband = []
MW_A2_zqual = []
for i in range(len(A2_data[0])):
	if A2_data[1][i] not in MW_stars_A2 and A2_data[7][i] != 'serendip1' and A2_data[7][i] != 'serendip2' and A2_data[7][i] != 'serendip3':
		A2_z.append(A2_data[0][i])
		A2_ra.append(A2_data[2][i])
		A2_dec.append(A2_data[3][i])
		A2_time.append(A2_data[4][i])
		A2_aband.append(A2_data[5][i])
		A2_zqual.append(A2_data[6][i])
	elif int(A2_data[1][i]) in MW_stars_A2 and A2_data[7][i] != 'serendip1' and A2_data[7][i] != 'serendip2' and A2_data[7][i] != 'serendip3':
		MW_A2_z.append(A2_data[0][i])
		MW_A2_ra.append(A2_data[2][i])
		MW_A2_dec.append(A2_data[3][i])
		MW_A2_time.append(A2_data[4][i])
		MW_A2_aband.append(A2_data[5][i])
		MW_A2_zqual.append(A2_data[6][i])

#calculate corrected LOS velocity 
def correct_vel(ra, dec, time, redshift, aband):
	# from astropy.utils.iers import conf
	# conf.auto_max_age = None #astropy told me to do this 
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
	keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)
	heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(time, format='mjd'), location=keck) 
	heliocorr_km_s = heliocorr.to(u.km/u.s) 
	vraw = redshift * const.c.to(u.km/u.s)
	vcorr = vraw + heliocorr_km_s - aband * const.c.to(u.km/u.s)

	return vcorr.value #km/s

#M33
good_data_A1 = (np.array(A1_zqual) == 1) | (np.array(A1_zqual) == 3) | (np.array(A1_zqual) == 4)
good_data_A2 = (np.array(A2_zqual) == 1) | (np.array(A2_zqual) == 3) | (np.array(A2_zqual) == 4)
print('calculating the velocities')
A1_vel = correct_vel(A1_ra, A1_dec, A1_time, A1_z, A1_aband)
A2_vel = correct_vel(A2_ra, A2_dec, A2_time, A2_z, A2_aband)
all_vels = list(A1_vel[good_data_A1]) + list(A2_vel[good_data_A2])

#MW
good_data_MW_A1 = (np.array(MW_A1_zqual) == 1) | (np.array(MW_A1_zqual) == 3) | (np.array(MW_A1_zqual) == 4)
#good_data_MW_A2 = (np.array(MW_A2_zqual) == 1) | (np.array(MW_A2_zqual) == 3) | (np.array(MW_A2_zqual) == 4)
print('calculating the velocities')
MW_A1_vel = correct_vel(MW_A1_ra, MW_A1_dec, MW_A1_time, MW_A1_z, MW_A1_aband)
MW_A1_vel = MW_A1_vel[good_data_MW_A1]
#MW_A2_vel = correct_vel(MW_A2_ra, MW_A2_dec, MW_A2_time, MW_A2_z, MW_A2_aband)
#MW_all_vels = list(MW_A1_vel[good_data_MW_A1]) + list(MW_A2_vel[good_data_MW_A2])

#plot a histogram of these values
def single_plot():
	rc('font', family = 'serif')
	fig, ax=plt.subplots(1)
	for axis in ['top','bottom','left','right']:
	        ax.spines[axis].set_linewidth(1)
	ax.tick_params(axis='x',which='both',bottom='on',top='on', direction='out')
	ax.tick_params(axis='x',which='both',top='on', direction='in')
	ax.tick_params(axis='y',which='both',left='on',top='off', direction='out')
	ax.tick_params(axis='y',which='both',right='on', direction='in')
	plt.tick_params(which='both', width=1)
	plt.tick_params(which='major', length=7)
	plt.tick_params(which='minor', length=4)
	plt.tick_params(labelsize=12) 
	plt.minorticks_on()

single_plot()
plt.hist(all_vels, bins=range(-350, 0, 10), histtype='step', stacked=True, fill=True, color='darkred', alpha=0.6)
plt.hist(MW_A1_vel, bins=range(-350, 0, 10), histtype='step', stacked=True, fill=False, color='purple', alpha=0.5)
plt.ylabel(r'$ \rm N$', fontsize = 13)
plt.xlabel(r'$ \rm Heliocentric\ Velocity\ (km\ s^{-1})$', fontsize = 13)
plt.title(r'$ \rm LOS\ Velocity\ on\ Minor\ Axis$', fontsize = 14)
plt.savefig('/Users/amandaquirk/Desktop/mam_vel_hist_MW.png')
plt.close()




