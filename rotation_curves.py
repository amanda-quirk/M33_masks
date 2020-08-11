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
# MS_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_MS.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
# HeB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_HeB_all.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
# AGB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_AGB.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
# RGB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_RGB.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')

#NOT smoothed data
MS_data = np.genfromtxt('../Data/M33_2018b_kinematics_MS.txt', dtype=None, names='ra, dec, vel, HI, CO, Ha, ID, zqual')
HeB_data = np.genfromtxt('../Data/M33_2018b_kinematics_HeB_all.txt', dtype=None, names='ra, dec, vel, HI, CO, Ha, ID, zqual')
AGB_data = np.genfromtxt('../Data/M33_2018b_kinematics_AGB.txt', dtype=None, names='ra, dec, vel, HI, CO, Ha, ID, zqual')
RGB_data = np.genfromtxt('../Data/M33_2018b_kinematics_RGB.txt', dtype=None, names='ra, dec, vel, HI, CO, Ha, ID, zqual')

#grab the parameters we care about here
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
	sine = np.sin(201.3 * np.pi / 180) #sine of M33's PA
	cosine = np.cos(201.3 * np.pi / 180) #cosine of M33's PA
	beta = xi * cosine - eta * sine #deg, y major axis beta
	alpha = eta * cosine + xi * sine #deg, x minor axis alpha

	#throwing in PA here to be used later
	angle = np.degrees(np.arctan2(beta, alpha)) #deg
	# PA = [] #east of north
	# for i in range(len(angle)):
	# 	if beta[i] < 0:
	# 		PA.append(angle[i] + 180)
	# 	else:
	# 		PA.append(angle[i])

	#calculate the distance from the star to the center 
	inclination_factor = np.cos(54 * np.pi / 180)**2 #inclination of M33 to deproject
	ang_dist = np.sqrt(alpha**2 +  (beta / inclination_factor)**2)
	dist = ang_dist * 14.12 #convet from deg to kpc

	# #calculating azimuthal angle
	# theta = np.degrees(np.arctan2(beta / inclination_factor, alpha))

	return alpha, beta, angle, dist #deg, deg, deg, kpc

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
# plt.scatter(HeB_deproj[0] * 14.12, HeB_deproj[1] * 14.12, c = HeB_deproj[4], norm=norm)
# plt.colorbar()
# plt.show()

#apply the tilted ring model to calculate v_rot ==============================================================================
#first we need to assign an inclination and PA from the tilted ring to each star
def assign_TR_params(dist): #stellar distance to center in kpc
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

	return assigned_i, assigned_PA - 180

#now we do the tilted ring model and also calculate AD
def vrot_tr_AD(vstar, vHI, vCO, vHa, ra, dec): #km/s, ha, deg
	#first use the previous functions to get PA and radius of the star
	deproj_data = deprojected_r(ra, dec)
	PA_star = np.array(deproj_data[2]) #deg
	dist = np.array(deproj_data[3]) #kpc

	#grab the TR values
	i_TR, PA_TR = assign_TR_params(dist)

	#use the formula separately for star vels and gases vels
	vsys = -180 #km/s
	deg_2_rad = np.pi / 180 #to be used for trig

	#tilted ring model --------
	sqrt_term = np.sqrt(1 + (np.tan((PA_TR) * deg_2_rad)**2 / np.cos(i_TR * deg_2_rad)**2))
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

def theta_vrot_tr(vel, ra, dec):#km/s, ha, deg
	#tilted ring model, as Karrie is using it -- Karrie first uses M33's overall inc and PA to calculate alpha and beta and r, then she assigns it to a ring, then she recalculates alpha beta and r using the PA and inc of the ring i/o M33 and then she calculates theta and vrot
	m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
	c_inm33 = sc.transform_to(m33.skyoffset_frame())
	xi, eta = c_inm33.lon, c_inm33.lat
	xi = xi.degree
	eta = eta.degree

	#first use the previous functions to get PA and radius of the star
	deproj_data = deprojected_r(ra, dec)
	dist = np.array(deproj_data[3]) #kpc

	#grab the TR values
	i_TR, PA_TR = assign_TR_params(dist) #deg, deg

	#first deproject xi and eta
	sine = np.sin(PA_TR * np.pi / 180) #sine of M33's PA
	cosine = np.cos(PA_TR * np.pi / 180) #cosine of M33's PA
	beta = xi * cosine - eta * sine #deg, y major axis beta
	alpha = eta * cosine + xi * sine #deg, x minor axis alpha

	#calculate the distance from the star to the center 
	inclination_factor = np.cos(i_TR * np.pi / 180)**2 #inclination of M33 to deproject

	#calculating azimuthal angle
	theta = np.degrees(np.arctan2(beta / inclination_factor, alpha))

	#use the formula separately for star vels and gases vels
	vsys = -180 #km/s
	deg_2_rad = np.pi / 180 #to be used for trig
	vrot = (vel - vsys) / (np.cos(theta * deg_2_rad) * np.sin(i_TR * deg_2_rad))
	return vrot

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

def outputs(vstar, vHI, vCO, vHa, ra, dec, age, IDs, smoothed=True): #will return a file of vrots and AD values, 1 RC, and 1 AD hist 
	#use the titled ring function to get dist, vrots, and ADs
	dist, star_vrot, HI_vrot, CO_vrot, Ha_vrot, AD_HI, AD_CO, AD_Ha = vrot_tr_AD(vstar, vHI, vCO, vHa, ra, dec)

	if smoothed == False:
		tag = 'all'
	else:
		tag = 'smoothed'

	#save the data here
	np.savetxt('/Volumes/Titan/M33/Data/{}_vrot_ad_data_{}.txt'.format(age, tag), np.c_[IDs, dist, star_vrot, HI_vrot, CO_vrot, Ha_vrot, AD_HI, AD_CO, AD_Ha], fmt='%s', delimiter='\t', header='ID, dist (kpc), star vrot (km/s), HI vrot, CO vrot, Ha vrot, AD wrt HI, AD wrt CO, AD wrt Ha')

	#rotation curve
	single_plot()
	plt.scatter(dist, star_vrot, c='r', label='{}'.format(age))
	plt.scatter(dist, HI_vrot, c='darkgrey', label='HI')
	plt.scatter(dist, CO_vrot, c='teal', label='CO')
	plt.scatter(dist, Ha_vrot, c='b', label='Ha')
	plt.xlabel('Deprojected R [kpc]')
	plt.ylabel('V rot titled ring [kpc]')
	plt.ylim(0, 300)
	plt.legend()
	plt.savefig('/Volumes/Titan/M33/Plots/{}_rc_{}.png'.format(age, tag))
	plt.close()

	#ad hist
	#get rid of nans in the AD list
	AD_HI = AD_HI[~np.isnan(AD_HI)]
	AD_CO = AD_CO[~np.isnan(AD_CO)]
	AD_Ha = AD_Ha[~np.isnan(AD_Ha)]

	single_plot()
	plt.hist(AD_HI, bins=range(-100, 200, 10), label='HI' + r'$={}$'.format(round(np.median(AD_HI),2)) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, linestyle='--', stacked=True, fill=False, color='darkgrey')
	plt.hist(AD_CO, bins=range(-100, 200, 10), label='CO' + r'$={}$'.format(round(np.median(AD_CO),2)) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, linestyle='--', stacked=True, fill=True, hatch='//', color='teal')
	plt.hist(AD_Ha, bins=range(-100, 200, 10), label='Ha' + r'$={}$'.format(round(np.median(AD_Ha),2)) + r'$\rm \ km \ s^{-1}$', normed=1, histtype='step', linewidth=1.6, stacked=True, fill=False, color='blue')
	plt.legend()
	plt.savefig('/Volumes/Titan/M33/Plots/{}_ad_{}.png'.format(age, tag))
	plt.close()

# outputs(MS_vel, MS_HI, MS_CO, MS_Ha, MS_ra, MS_dec, 'MS', MS_ID)
# outputs(HeB_vel, HeB_HI, HeB_CO, HeB_Ha, HeB_ra, HeB_dec, 'HeB', HeB_ID)
# outputs(AGB_vel, AGB_HI, AGB_CO, AGB_Ha, AGB_ra, AGB_dec, 'AGB', AGB_ID)
# outputs(RGB_vel, RGB_HI, RGB_CO, RGB_Ha, RGB_ra, RGB_dec, 'RGB', RGB_ID)

# outputs(MS_vel, MS_HI, MS_CO, MS_Ha, MS_ra, MS_dec, 'MS', MS_ID, smoothed=False)
# outputs(HeB_vel, HeB_HI, HeB_CO, HeB_Ha, HeB_ra, HeB_dec, 'HeB', HeB_ID, smoothed=False)
# outputs(AGB_vel, AGB_HI, AGB_CO, AGB_Ha, AGB_ra, AGB_dec, 'AGB', AGB_ID, smoothed=False)
# outputs(RGB_vel, RGB_HI, RGB_CO, RGB_Ha, RGB_ra, RGB_dec, 'RGB', RGB_ID, smoothed=False)


#compares the tilted ring model I have been using to the one Karrie is using (except does not keep her iteration) ================
def rc_comp_plot(ra, dec, vel, vHI, vCO, vHa, age):
	r = deprojected_r(ra, dec)[3] #alpha, beta, PA, dist

	vrot_karrie = theta_vrot_tr(vel, ra, dec)
	vrot_0 = vrot_tr_AD(vel, vHI, vCO, vHa, ra, dec)[1] #dist, star_vrot, HI_vrot, CO_vrot, Ha_vrot, AD_HI, AD_CO, AD_Ha

	single_plot()
	plt.scatter(r, vrot_karrie, c='black', label='K')
	plt.scatter(r, vrot_0, c='blue', label='o')
	plt.xlabel('Deprojected R [kpc]')
	plt.ylabel('V_rot [kpc]')
	plt.ylim(0, 300)
	plt.legend()
	plt.savefig('/Users/amandaquirk/Desktop/{}_vrot_test.png'.format(age))
	plt.close()

# rc_comp_plot(MS_ra, MS_dec, MS_vel, MS_HI, MS_CO, MS_Ha, 'MS')
# rc_comp_plot(AGB_ra, AGB_dec, AGB_vel, AGB_HI, AGB_CO, AGB_Ha, 'AGB')
# rc_comp_plot(RGB_ra, RGB_dec, RGB_vel, RGB_HI, RGB_CO, RGB_Ha, 'RGB')
# rc_comp_plot(HeB_ra, HeB_dec, HeB_vel, HeB_HI, HeB_CO, HeB_Ha, 'HeB')

def theta(ra, dec):
	m33 = SkyCoord(23.4583 * u.deg, 30.6602 * u.deg)
	sc = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle,u.deg))
	c_inm33 = sc.transform_to(m33.skyoffset_frame())
	xi, eta = c_inm33.lon, c_inm33.lat
	xi = xi.degree
	eta = eta.degree

	#grab the TR values
	deproj_data = deprojected_r(ra, dec)
	dist = np.array(deproj_data[3]) 
	i_TR, PA_TR = assign_TR_params(dist) #deg, deg

	#first deproject xi and eta
	sine = np.sin(PA_TR * np.pi / 180) #sine of M33's PA
	cosine = np.cos(PA_TR * np.pi / 180) #cosine of M33's PA
	beta = xi * cosine - eta * sine #deg, y major axis beta
	alpha = eta * cosine + xi * sine #deg, x minor axis alpha

	#calculate the distance from the star to the center 
	inclination_factor = np.cos(54 * np.pi / 180)**2 #inclination of M33 to deproject

	#calculating azimuthal angle
	theta = np.degrees(np.arctan2(beta, alpha * inclination_factor))
	return theta

theta_MS = theta(MS_ra, MS_dec)
deprojected_MS = deprojected_r(MS_ra, MS_dec)
MS_i, MS_PAs = assign_TR_params(deprojected_MS[3])
MS_PA = np.array(deprojected_MS[2])

#plt.scatter(np.tan(np.deg2rad(MS_PA)), np.cos(np.deg2rad(MS_i)) * np.tan(np.deg2rad(theta_MS)))
PA_part = np.tan(MS_PA * np.pi / 180)
theta_part = np.cos(54 * np.pi / 180) * np.tan(theta_MS * np.pi / 180)
print(PA_part - theta_part)
# plt.hist(PA_part, color='k', label='tan(PA)')
# plt.legend()
# plt.show()
# plt.close()
# plt.hist(theta_part, color='b', alpha=.5, label='tan(theta) * cos(i)')
# plt.legend()
# plt.show()
# plt.close()


