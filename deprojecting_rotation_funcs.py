import numpy as np 
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import constants as const
from astropy import units as u
from astropy.time import Time

def correct_vel(ra, dec, time, redshift, aband):
	# from astropy.utils.iers import conf
	# conf.auto_max_age = None #astropy told me to do this 
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
	keck = EarthLocation.from_geodetic(lat=19.8283*u.deg, lon=-155.4783*u.deg, height=4160*u.m)
	heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(time, format='mjd'), location=keck) 
	heliocorr_km_s = heliocorr.to(u.km/u.s) 
	vraw = redshift * const.c.to(u.km/u.s)
	vcorr_aband = vraw + heliocorr_km_s - aband * const.c.to(u.km/u.s)
	vcorr = vraw + heliocorr_km_s
	return vcorr.value, vcorr_aband.value #km/s

def deprojection_geo(ra, dec, galaxy, unit='hourangle'):
	#get all the galaxy info
	if galaxy == 'M33':
		galaxy_ra = 23.4583 #deg
		galaxy_dec = 30.6602 #deg
		galaxy_pa = 22.5 #deg
		galaxy_inc = 52 #deg
		deg_2_kpc = 14.12
		dmod = 24.67 # what Anil assumed  # Kam et al 2017 uses 840 kpc; but tilted ring model in arcmin, so this is not crucial
		galaxy_sys_v = -180 #km/s
	elif galaxy == 'M31':
		galaxy_ra = 10.6847083 #deg4
		galaxy_dec = 41.26875 #deg
		galaxy_pa = 37 #deg
		galaxy_inc = 77 #deg
		deg_2_kpc = 13.67
		dmod = 24.38 # Riess 2012
		galaxy_sys_v = -300 #km/s
	else:
		print('You do not study that galaxy!!')

	#find xi and eta
	gal = SkyCoord(galaxy_ra * u.deg, galaxy_dec * u.deg)

	if unit == 'hourangle':
		sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
	else:
		sc = SkyCoord(ra=ra, dec=dec, unit=(u.deg,u.deg))
	c_in_gal = sc.transform_to(gal.skyoffset_frame())
	xi, eta = c_in_gal.lon, c_in_gal.lat
	xi = xi.degree
	eta = eta.degree
	
	#first, find alpha and beta with global parameters in order to calculate the distance from center ========================
	sine0 = np.sin(galaxy_pa * np.pi / 180) #sine of galaxy's PA
	cosine0 = np.cos(galaxy_pa * np.pi / 180) #cosine of galaxy's PA
	beta0 = xi * cosine0 - eta * sine0 #deg, y major axis beta
	alpha0 = eta * cosine0 + xi * sine0 #deg, x minor axis alpha
	inclination_factor0 = np.cos(galaxy_inc * np.pi / 180) #inclination of galaxy to deproject
	ang_dist0 = np.sqrt(alpha0**2 +  (beta0 / inclination_factor0)**2)
	dist0 = ang_dist0 * deg_2_kpc
	
	#recalculate alpha and beta with tilted ring PA and i values =======================
	#read in the titled ring parameters
	if galaxy == 'M31':
		r_kpc, PA_tr, i_tr = np.loadtxt('data/HI_PA_i_vrot.txt', usecols=(0, 1, 2,), unpack=True)
	elif galaxy == 'M33':
		r_kpc, i_tr, PA_tr = np.loadtxt('../Data/M33_HI_tilted_ring.txt', usecols=(1, 4, 5,), unpack=True)
	#find the ring that each star fits in
	assigned_i = np.zeros(len(dist0)) #deg
	assigned_PA = np.zeros(len(dist0)) #deg
	for i in range(len(dist0)):
		difference = r_kpc - dist0[i] #difference between the star radii and all the rings
		pos_diff = np.array([a for a in difference if a > 0]) #eliminate rings inner to the star
		if len(pos_diff) == 0:
			assigned_i[i] = i_tr[-1]
			assigned_PA[i] = PA_tr[-1]
		else:
			smallest_diff = min(pos_diff) #find the smallest positive difference
			ring = list(difference).index(smallest_diff) #find where in the original list the right ring is
			assigned_i[i] = i_tr[ring]
			assigned_PA[i] = PA_tr[ring] #Kam seems to define PA as my definition + 180 but this doesn't matter for tangent
	
	sine = np.sin(assigned_PA * np.pi / 180) #sine of TR PA
	cosine = np.cos(assigned_PA * np.pi / 180) #cosine of TR PA
	beta = xi * cosine - eta * sine #deg, y major axis beta
	alpha = eta * cosine + xi * sine #deg, x minor axis alpha
	inclination_factor = np.cos(assigned_i * np.pi / 180) #inclination of galaxy to deproject
	ang_dist = np.sqrt(alpha**2 +  (beta / inclination_factor)**2)
	dist = ang_dist * deg_2_kpc #kpc
	
	#calculate PA and theta ========================================================
	PA = np.arctan2(beta, alpha) 
	theta = np.arctan2(beta, alpha * np.cos(np.deg2rad(assigned_i)))	

	return xi, eta, alpha, beta, dist, PA, theta, assigned_PA, assigned_i #kpc, kpc, deg, deg, kpc, rad, rad, deg, deg

def vrot_0(vel, ra, dec, galaxy):
	#get galaxy info
	if galaxy == 'M33':
		galaxy_sys_v = -180 #km/s
	elif galaxy == 'M31':
		galaxy_sys_v = -300 #km/s
	else:
		print('You do not study that galaxy!!')

	vsys = galaxy_sys_v #km/s
	pos_data = deprojection_geo(ra, dec, galaxy)
	vrot = abs((vel - vsys) / (np.cos(pos_data[6]) * np.sin(np.deg2rad(pos_data[-1]))))
	return vrot

def vrot_tr(vel, ra, dec, galaxy):
	'''
	Things to figure out: 
		-for tan(PA) term: should it be tan(PA_star - PA_TR)? That would be np.tan(pos_data[5] - np.deg2rad(pos_data[-2])). It gives different vrot values than vrot_0, so I am guessing not to do that and that it worked with a different definition of PA?
	'''
	#get galaxy info
	if galaxy == 'M33':
		galaxy_sys_v = -180 #km/s
	elif galaxy == 'M31':
		galaxy_sys_v = -300 #km/s
	else:
		print('You do not study that galaxy!!')

	vsys = galaxy_sys_v #km/s
	pos_data = deprojection_geo(ra, dec, galaxy)
	first_term = abs(vel - vsys) / np.sin(np.deg2rad(pos_data[-1])) 
	sqrt_term = np.sqrt(1 + (np.tan(pos_data[5])**2 / np.cos(np.deg2rad(pos_data[-1]))**2)) 
	vrot = first_term * sqrt_term
	return vrot 


