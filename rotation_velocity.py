import numpy as np 
from astropy.coordinates import SkyCoord
from astropy import units as u

'''
file for Karrie; also used for calculating projected and deprojected radii of potential masks for a proposal 
takes the corrected LOS velocity and ra/dec of a data point and returns the rotation velocity from the tilted ring model, as described by Kam et al. 2017

to use, just call vrot_tr(velocity (km/s), ra (ha), dec (deg)); it will read in the earlier functions
'''

#calculate distance from star to center =======================================================================================
def deprojected_r(ra, dec): #ha and deg 
	#convert to xi and eta centered on M33
	m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
	c_inm33 = sc.transform_to(m33.skyoffset_frame())
	xi, eta = c_inm33.lon, c_inm33.lat
	xi = xi.degree
	eta = eta.degree

	#first deproject xi and eta
	sine = np.sin(22 * np.pi / 180) #sine of M33's PA
	cosine = np.cos(22 * np.pi / 180) #cosine of M33's PA
	x = xi * cosine - eta * sine #deg
	y = eta * cosine + xi * sine #deg

	#throwing in PA here to be used later
	angle = np.arctan(y / x ) * 180 / np.pi #deg
	PA = [] #east of north
	for i in range(len(angle)):
		if y[i] < 0:
			PA.append(angle[i] + 180)
		else:
			PA.append(angle[i])

	#calculate the distance from the star to the center 
	inclination_factor = np.cos(54 * np.pi / 180)**2 #inclination of M33 to deproject
	ang_dist = np.sqrt(x**2 + y**2 / inclination_factor)
	dist = ang_dist * 14.12 #convet from deg to kpc

	return x, y, PA, dist #deg, deg, deg, kpc

#apply the tilted ring model to calculate v_rot ==============================================================================
#first we need to assign an inclination and PA from the tilted ring to each star
def assign_TR_params(dist): #stellar distance to center in kpc
	#read in the titled ring parameters
	r_kpc, i_tr, PA_tr = np.loadtxt('M33_HI_tilted_ring.txt', usecols=(1, 4, 5,), unpack=True)

	#find the ring that each star fits in
	assigned_i = np.zeros(len(dist))
	assigned_PA = np.zeros(len(dist))
	for i in range(len(dist)):
		difference = r_kpc - dist[i] #difference between the star radii and all the rings
		pos_diff = np.array([a for a in difference if a > 0]) #eliminate rings inner to the star
		smallest_diff = min(pos_diff) #find the smallest positive difference
		ring = list(difference).index(smallest_diff) #find where in the original list the right ring is
		assigned_i[i] = i_tr[ring]
		assigned_PA[i] = PA_tr[ring] #Kam seems to define PA as my definition + 180 but this doesn't matter for tangent

	return assigned_i, assigned_PA

#now we do the tilted ring model and also calculate AD
def vrot_tr(vel, ra, dec): #km/s, ha, deg
	#first use the previous functions to get PA and radius of the star
	deproj_data = deprojected_r(ra, dec)
	PA_star = np.array(deproj_data[2]) #deg
	dist = np.array(deproj_data[3]) #kpc

	#grab the TR values
	i_TR, PA_TR = assign_TR_params(dist) #deg, deg

	#use the formula separately for star vels and gases vels
	vsys = -180 #km/s
	deg_2_rad = np.pi / 180 #to be used for trig

	#tilted ring model --------
	sqrt_term = np.sqrt(1 + (np.tan((PA_TR - PA_star) * deg_2_rad)**2 / np.cos(i_TR * deg_2_rad)**2))
	vsys_term = (vel - vsys) / np.sin(i_TR * deg_2_rad)
	vrot = abs(vsys_term) * sqrt_term

	return vrot

#=========================================================================================================================
#calculate radii

def deprojected_radii(ra, dec): #deg and deg 
	#convert to xi and eta centered on M33
	m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.deg,u.deg))
	c_inm33 = sc.transform_to(m33.skyoffset_frame())
	xi, eta = c_inm33.lon, c_inm33.lat
	xi = xi.degree
	eta = eta.degree

	#first deproject xi and eta
	sine = np.sin(22 * np.pi / 180) #sine of M33's PA
	cosine = np.cos(22 * np.pi / 180) #cosine of M33's PA
	x = xi * cosine - eta * sine #deg
	y = eta * cosine + xi * sine #deg

	#calculate the distance from the star to the center 
	inclination_factor = np.cos(54 * np.pi / 180)**2 #inclination of M33 to deproject
	ang_dist = np.sqrt(x**2 + y**2 / inclination_factor)
	dist = ang_dist * 14.12 #convet from deg to kpc

	return dist #kpc

