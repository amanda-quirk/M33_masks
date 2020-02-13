import numpy as np 
import matplotlib.pyplot as plt 
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib import rc 

'''
takes the output of smoothing_pos_vel and returns rotation curves
-uses the tilted ring
-also creates AD hists
'''
#read in the smoothed data ====================================================================================================
#data has header: ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter
#MS_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_MS.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
HeB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_HeB.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
AGB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_AGB.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
RGB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_RGB.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')

#grab the parameters we care about here
# MS_ra = MS_data['ra']
# MS_ra = np.array([a.decode("utf-8") for a in MS_ra])
# MS_dec = MS_data['dec']
# MS_dec = np.array([a.decode("utf-8") for a in MS_dec])
# MS_vel = MS_data['vel']
# MS_HI = MS_data['HI']
# MS_CO = MS_data['CO']
# MS_Ha = MS_data['Ha']
# MS_ID = MS_data['ID']
# MS_ID = np.array([a.decode("utf-8") for a in MS_ID])

HeB_ra = HeB_data['ra']
HeB_ra = np.array([a.decode("utf-8") for a in HeB_ra])
HeB_dec = HeB_data['dec']
HeB_dec = np.array([a.decode("utf-8") for a in HeB_dec])
HeB_vel = HeB_data['vel']
HeB_HI = HeB_data['HI']
HeB_CO = HeB_data['CO']
HeB_Ha = HeB_data['Ha']
HeB_ID = HeB_data['ID']
HeB_ID = np.array([a.decode("utf-8") for a in HeB_ID])

AGB_ra = AGB_data['ra']
AGB_ra = np.array([a.decode("utf-8") for a in AGB_ra])
AGB_dec = AGB_data['dec']
AGB_dec = np.array([a.decode("utf-8") for a in AGB_dec])
AGB_vel = AGB_data['vel']
AGB_HI = AGB_data['HI']
AGB_CO = AGB_data['CO']
AGB_Ha = AGB_data['Ha']
AGB_ID = AGB_data['ID']
AGB_ID = np.array([a.decode("utf-8") for a in AGB_ID])

RGB_ra = RGB_data['ra']
RGB_ra = np.array([a.decode("utf-8") for a in RGB_ra])
RGB_dec = RGB_data['dec']
RGB_dec = np.array([a.decode("utf-8") for a in RGB_dec])
RGB_vel = RGB_data['vel']
RGB_HI = RGB_data['HI']
RGB_CO = RGB_data['CO']
RGB_Ha = RGB_data['Ha']
RGB_ID = RGB_data['ID']
RGB_ID = np.array([a.decode("utf-8") for a in RGB_ID])

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
	PA = []
	for i in range(len(angle)):
		if y[i] < 0:
			PA.append(angle[i] + 180)
		else:
			PA.append(angle[i])

	#calculate the distance from the star to the center using angular distance
	inclination_factor = np.cos(54 * np.pi / 180)**2 #inclination of M33 to deproject
	ang_dist = np.sqrt(x**2 + y**2 / inclination_factor)
	dist = ang_dist * 14.12 #convet from deg to kpc

	return x, y, PA, dist #deg, deg, deg, kpc

#since I am paranoid, below tests the PA assignment
# import matplotlib as mpl
# HeB_deproj = deprojected_r(HeB_ra, HeB_dec)
# plt.plot([0,0], [-4, 4], c='k')
# plt.plot([-4,4], [0,0], c='k')
# plt.xlim(4, -4)
# bounds = np.array([-90, 0, 90, 180, 270, 360]) #discrete bounds for color map
# cmap = plt.cm.jet  #define the colormap
# cmaplist = [cmap(i) for i in range(cmap.N)] #extract all the colors
# cmaplist[0] = (.5, .5, .5, 1.0) #set first entry to be grey
# cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N) #create new colormap
# norm = mpl.colors.BoundaryNorm(bounds, cmap.N) 
# plt.scatter(HeB_deproj[0] * 14.12, HeB_deproj[1] * 14.12, c = HeB_deproj[2], norm=norm)
# plt.colorbar()
# plt.show()

#apply the tilted ring model to calculate v_rot ==============================================================================
#first we need to assign an inclination and PA from the tilted ring to each star
def assign_TR_params(dist): #stellar distance to center; kpc
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
		assigned_PA = PA_tr[ring] #should i subtract 180 from these?

	return assigned_i, assigned_PA

#now we do the tilted ring model and also calculate AD
def vrot_tr_AD(vstar, vHI, vCO, vHa, ra, dec): #km/s, ha, deg
	#first use the previous functions to get PA and radius of the star
	deproj_data = deprojected_r(ra, dec)
	PA_star = deproj_data[2] #deg
	dist = deproj_data[3] #kpc

	#grab the TR values
	i_TR, PA_TR = assign_TR_params(dist)

	#use the formula separately for star vels and gases vels
	vsys = -180 #km/s
	deg_2_rad = np.pi / 180 #to be used for trig

	sqrt_term = np.sqrt(1 + (np.tan((PA_TR - PA_star) * deg_2_rad)**2 / np.cos(i_TR * deg_2_rad)**2))
	star_vel_term = (vstar - vsys) / np.sin(i_TR * deg_2_rad)
	HI_vel_term = (vHI - vsys) / np.sin(i_TR * deg_2_rad)
	CO_vel_term = (vCO - vsys) / np.sin(i_TR * deg_2_rad)
	Ha_vel_term = (vHa - vsys) / np.sin(i_TR * deg_2_rad)

	star_vrot = abs(star_vel_term) * sqrt_term
	HI_vrot = abs(HI_vel_term) * sqrt_term
	CO_vrot = abs(CO_vel_term) * sqrt_term
	Ha_vrot = abs(Ha_vel_term) * sqrt_term

	AD_HI = HI_vrot - star_vrot
	AD_CO = CO_vrot - star_vrot
	AD_Ha = Ha_vrot - star_vrot

	return dist, star_vrot, HI_vrot, CO_vrot, Ha_vrot, AD_HI, AD_CO, AD_Ha

#plotting and saving data ====================================================================================================
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

def outputs(vstar, vHI, vCO, vHa, ra, dec, age, IDs): #will return a file of vrots and AD values, 1 RC, and 1 AD hist 
	#use the titled ring function to get dist, vrots, and ADs
	dist, star_vrot, HI_vrot, CO_vrot, Ha_vrot, AD_HI, AD_CO, AD_Ha = vrot_tr_AD(vstar, vHI, vCO, vHa, ra, dec)

	#save the data here
	np.savetxt('/Users/amandaquirk/Desktop/{}_vrot_ad_data.txt'.format(age), np.c_[IDs, dist, star_vrot, HI_vrot, CO_vrot, Ha_vrot, AD_HI, AD_CO, AD_Ha], fmt='%s', delimiter='\t', header='ID, dist (kpc), star vrot (km/s), HI vrot, CO vrot, Ha vrot, AD wrt HI, AD wrt CO, AD wrt Ha')

	#rotation curve
	single_plot()
	plt.scatter(dist, star_vrot, c='r', label='{}'.format(age))
	plt.scatter(dist, HI_vrot, c='darkgrey', label='HI')
	plt.scatter(dist, CO_vrot, c='teal', label='CO')
	plt.scatter(dist, Ha_vrot, c='b', label='Ha')
	plt.xlabel('Deprojected R [kpc]')
	plt.ylabel('V rot titled ring [kpc]')
	plt.legend()
	plt.savefig('/Users/amandaquirk/Desktop/{}_rc.png'.format(age))
	plt.close()

	#ad hist
	#get rid of nans in the AD list
	AD_HI = AD_HI[~np.isnan(AD_HI)]
	AD_CO = AD_CO[~np.isnan(AD_CO)]
	AD_Ha = AD_Ha[~np.isnan(AD_Ha)]

	single_plot()
	plt.hist(AD_HI, bins=range(-300, 350, 10), label='HI' + r'$={}$'.format(round(np.median(AD_HI),2)) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, linestyle='--', stacked=True, fill=False, color='darkgrey')
	plt.hist(AD_CO, bins=range(-300, 350, 10), label='CO' + r'$={}$'.format(round(np.median(AD_CO),2)) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, linestyle='--', stacked=True, fill=True, hatch='//', color='teal')
	plt.hist(AD_Ha, bins=range(-300, 350, 10), label='Ha' + r'$={}$'.format(round(np.median(AD_Ha),2)) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, stacked=True, fill=False, color='blue')
	plt.legend()
	plt.savefig('/Users/amandaquirk/Desktop/{}_ad.png'.format(age))
	plt.close()

outputs(HeB_vel, HeB_HI, HeB_CO, HeB_Ha, HeB_ra, HeB_dec, 'HeB', HeB_ID)





