import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib import rc 

def deprojection_geo(ra, dec):
	#find xi and eta
	m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
	c_inm33 = sc.transform_to(m33.skyoffset_frame())
	xi, eta = c_inm33.lon, c_inm33.lat
	xi = xi.degree
	eta = eta.degree
	
	#first, find alpha and beta with global parameters in order to calculate the distance from center ========================
	M33_i = 54 #deg
	M33_PA = 22 #deg
	sine = np.sin(M33_PA * np.pi / 180) #sine of M33's PA
	cosine = np.cos(M33_PA * np.pi / 180) #cosine of M33's PA
	beta = xi * cosine - eta * sine #deg, y major axis beta
	alpha = eta * cosine + xi * sine #deg, x minor axis alpha
	inclination_factor = np.cos(M33_i * np.pi / 180)**2 #inclination of M33 to deproject
	ang_dist = np.sqrt(alpha**2 +  (beta / inclination_factor)**2)
	dist = ang_dist * 14.12
	
	#recalculate alpha and beta with tilted ring PA and i values =======================
	#read in the titled ring parameters
	r_kpc, i_tr, PA_tr = np.loadtxt('../Data/M33_HI_tilted_ring.txt', usecols=(1, 4, 5,), unpack=True)
	#find the ring that each star fits in
	assigned_i = np.zeros(len(dist)) #deg
	assigned_PA = np.zeros(len(dist)) #deg
	for i in range(len(dist)):
		difference = r_kpc - dist[i] #difference between the star radii and all the rings
		pos_diff = np.array([a for a in difference if a > 0]) #eliminate rings inner to the star
		smallest_diff = min(pos_diff) #find the smallest positive difference
		ring = list(difference).index(smallest_diff) #find where in the original list the right ring is
		assigned_i[i] = i_tr[ring]
		assigned_PA[i] = PA_tr[ring] #Kam seems to define PA as my definition + 180 but this doesn't matter for tangent
	
	sine = np.sin(assigned_PA * np.pi / 180) #sine of TR PA
	cosine = np.cos(assigned_PA * np.pi / 180) #cosine of TR PA
	beta = xi * cosine - eta * sine #deg, y major axis beta
	alpha = eta * cosine + xi * sine #deg, x minor axis alpha
	inclination_factor = np.cos(assigned_i * np.pi / 180)**2 #inclination of M33 to deproject
	ang_dist = np.sqrt(alpha**2 +  (beta / inclination_factor)**2)
	dist = ang_dist * 14.12 #kpc
	
	#calculate PA and theta ========================================================
	PA = np.arctan2(beta, alpha) 
	theta = np.arctan2(beta, alpha * np.cos(np.deg2rad(M33_i)))	

	return xi, eta, alpha, beta, dist, PA, theta, assigned_PA, assigned_i #kpc, kpc, deg, deg, kpc, rad, rad, deg, deg

def vrot_0(vel, ra, dec):
	vsys = -180 #km/s
	pos_data = deprojection_geo(ra, dec)
	vrot = abs((vel - vsys) / (np.cos(pos_data[6]) * np.sin(np.deg2rad(pos_data[-1])))) #use 54 or the TR inc?
	return vrot

def vrot_tr(vel, ra, dec):
	'''
	Things to figure out: 
		-for tan(PA) term: should it be tan(PA_star - PA_TR)? That would be np.tan(pos_data[5] - np.deg2rad(pos_data[-2])). It gives different vrot values than vrot_0, so I am guessing not to do that and that it worked with a different definition of PA?
	'''
	vsys = -180 #km/s
	pos_data = deprojection_geo(ra, dec)
	first_term = abs(vel - vsys) / np.sin(np.deg2rad(pos_data[-1])) 
	sqrt_term = np.sqrt(1 + (np.tan(pos_data[5])**2 / np.cos(np.deg2rad(pos_data[-1]))**2)) 
	vrot = first_term * sqrt_term
	return vrot 


