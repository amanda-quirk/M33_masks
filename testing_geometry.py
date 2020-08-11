import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u

'''
trying to get the vrot tr and karrie's vrot equation to make
'''

#starting with tan(PA) = tan(theta) * cos(i)
# alpha = np.array([2, 3, 4, 5, 7]) #dummy variables to test that at least the math works
# beta = np.array([0, 5, 6, 9, 2]) #check -- the math works
# xi = np.array([2, -1, 0.5, 1, 1.75]) / 14.12 #deg #dummy variables to test that at least the math works
# eta = np.array([-2.75, 3.01, .45, -1.1, 0.01]) / 14.12
#======================================================================================
MS_data = np.genfromtxt('../Data/M33_2018b_kinematics_MS.txt', dtype=None, names='ra, dec, vel, HI, CO, Ha, ID, zqual')
MS_ra = MS_data['ra']
MS_ra = np.array([a.decode("utf-8") for a in MS_ra])
MS_dec = MS_data['dec']
MS_dec = np.array([a.decode("utf-8") for a in MS_dec])
MS_vel = MS_data['vel']
MS_HI = MS_data['HI']
MS_CO = MS_data['CO']
MS_Ha = MS_data['Ha']
MS_ID = MS_data['ID']
MS_ID = np.array([a.decode("utf-8") for a in MS_ID])

#find xi and eta
m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
sc = SkyCoord(ra=MS_ra, dec=MS_dec, unit=(u.hourangle,u.deg))
c_inm33 = sc.transform_to(m33.skyoffset_frame())
xi, eta = c_inm33.lon, c_inm33.lat
xi = xi.degree
eta = eta.degree

#first, find alpha and beta with global parameters in order to calculate the distance from center ======
M33_i = 54 #deg
M33_PA = 22 #deg
sine = np.sin(M33_PA * np.pi / 180) #sine of M33's PA
cosine = np.cos(M33_PA * np.pi / 180) #cosine of M33's PA
beta = xi * cosine - eta * sine #deg, y major axis beta
alpha = eta * cosine + xi * sine #deg, x minor axis alpha
inclination_factor = np.cos(M33_i * np.pi / 180)**2 #inclination of M33 to deproject
ang_dist = np.sqrt(alpha**2 +  (beta / inclination_factor)**2)
dist = ang_dist * 14.12

#recalculate alpha and beta with tilted ring PA and i values ===========
#read in the titled ring parameters
r_kpc, i_tr, PA_tr = np.loadtxt('../Data/M33_HI_tilted_ring.txt', usecols=(1, 4, 5,), unpack=True)
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

beta = xi * cosine - eta * sine #deg, y major axis beta
alpha = eta * cosine + xi * sine #deg, x minor axis alpha
inclination_factor = np.cos(M33_i * np.pi / 180)**2 #inclination of M33 to deproject
ang_dist = np.sqrt(alpha**2 +  (beta / inclination_factor)**2)
dist = ang_dist * 14.12

#calculate PA and theta
PA = np.arctan2(beta, alpha)
theta = np.arctan2(beta, alpha * np.cos(np.deg2rad(M33_i)))

PA_part = np.tan(PA)
theta_part = np.tan(theta) * np.cos(np.deg2rad(M33_i))
print(PA_part - theta_part)
# plt.hist(PA_part, color='k', label='tan(PA)')
# plt.hist(theta_part, color='b', alpha=.5, label='tan(theta) * cos(i)')
# plt.legend()
# plt.show()
# plt.close()
