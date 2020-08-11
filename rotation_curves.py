import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rc 
from deprojecting_rotation_funcs import deprojection_geo, vrot_0, vrot_tr

'''
takes the output of smoothing_pos_vel and returns rotation curves and AD hists
'''

#read in the smoothed data ====================================================================================================
#data has header: ra_goodcenter (ha), dec_goodcenter (deg), smoothed_v (km/s), smoothed_err (km/s), dispersion, HI_goodcenter, CO_goodcenter, Ha_goodcenter, ID_goodcenter, zqual_goodcenter
#MS_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_MS.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
HeB_data = np.genfromtxt('../Data/M33_2018b_smoothed_kinematics_HeB_all.txt', dtype=None, names='ra, dec, vel, err, disp, HI, CO, Ha, ID, zqual')
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
	#use the titled ring function to get dist, vrots
	xi, eta, alpha, beta, dist, PA, theta, assigned_PA, assigned_i = deprojection_geo(ra, dec)
	vrot0_star = vrot_0(vstar, ra, dec)
	vrot_star = vrot_tr(vstar, ra, dec)
	#vrot0_HI = vrot_0(vHI, ra, dec)
	vrot_HI = vrot_tr(vHI, ra, dec)
	#vrot0_Ha = vrot_0(vHa, ra, dec)
	vrot_Ha = vrot_tr(vHa, ra, dec)
	#vrot0_CO = vrot_0(vCO, ra, dec)
	vrot_CO = vrot_tr(vCO, ra, dec)

	#calculate AD
	AD_HI = vrot_HI - vrot_star
	AD_CO = vrot_CO - vrot_star
	AD_Ha = vrot_Ha - vrot_star

	#save the data here
	np.savetxt('/Volumes/Titan/M33/Data/{}_vrot_ad_data.txt'.format(age), np.c_[IDs, dist, vrot_star, vrot_HI, vrot_CO, vrot_Ha, AD_HI, AD_CO, AD_Ha], fmt='%s', delimiter='\t', header='ID, dist (kpc), star vrot (km/s), HI vrot, CO vrot, Ha vrot, AD wrt HI, AD wrt CO, AD wrt Ha')

	#rotation curve
	single_plot()
	plt.scatter(dist, vrot_star, c='r', label='{}'.format(age))
	plt.scatter(dist, vrot0_star, c='b') #comp for sanity check of the rotation velocities
	# plt.scatter(dist, vrot_HI, c='darkgrey', label='HI')
	# plt.scatter(dist, vrot_CO, c='teal', label='CO')
	# plt.scatter(dist, vrot_Ha, c='b', label='Ha')
	plt.xlabel('Deprojected R [kpc]')
	plt.ylabel('V rot titled ring [kpc]')
	plt.ylim(0, 300)
	plt.legend()
	plt.savefig('/Volumes/Titan/M33/Plots/{}_rc.png'.format(age))
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
	plt.savefig('/Volumes/Titan/M33/Plots/{}_ad.png'.format(age))
	plt.close()

# outputs(MS_vel, MS_HI, MS_CO, MS_Ha, MS_ra, MS_dec, 'MS', MS_ID)
outputs(HeB_vel, HeB_HI, HeB_CO, HeB_Ha, HeB_ra, HeB_dec, 'HeB', HeB_ID)
outputs(AGB_vel, AGB_HI, AGB_CO, AGB_Ha, AGB_ra, AGB_dec, 'AGB', AGB_ID)
outputs(RGB_vel, RGB_HI, RGB_CO, RGB_Ha, RGB_ra, RGB_dec, 'RGB', RGB_ID)


