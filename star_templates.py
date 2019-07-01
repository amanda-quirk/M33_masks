from astropy.io import fits
import matplotlib.pyplot as plt 
import numpy as np 
from astropy import constants as const
from astropy import units as u

hdu = fits.open('../Data/deimos_templates.fits')
data = hdu[0].data

'''
template names:
NAME0   = 'HD066665 B1V'       /                                                
NAME1   = 'HD087737 A0Ib'      /                                                
NAME2   = 'HR6161 A0III'       /                                                
NAME3   = 'HD070569 A9V'       /                                                
NAME4   = 'HD107513 A9V'       /                                                
NAME5   = 'HD089449 F6IV'      /                                                
NAME6   = 'HD237846 F8'        /                                                
NAME7   = 'HD105546 G2III'     /                                                
NAME8   = 'HD065934 G8III'     /                                                
NAME9   = 'HD085990 K0III'     /                                                
NAME10  = 'HD052071 K2III'     /                                                
NAME11  = 'HD088230 K7V'       /                                                
NAME12  = 'HD061913 M2II-III'  /                                                
NAME13  = 'HD113285 M8III'     /                                                
NAME14  = 'HD41117 B2Iaev'     /                                                
NAME15  = 'HD164353 B5I'       /                                                
NAME16  = 'HD65900 A1V'        /                                                
NAME17  = 'HD164136 F2II'      /                                                
NAME18  = 'HD197572 F7Ib'      /                                                
NAME19  = 'HD172365 F8Ib-II'   /                                                
NAME20  = 'HD74395 G1Ib'       /                                                
NAME21  = 'HD188727 G5Ibv'     /                                                
NAME22  = 'HD164349 K0.5IIb'   /                                                
NAME23  = 'HD52005 K4Iab'      /                                                
NAME24  = 'HD151217 K5III'     /                                                
NAME25  = 'HD60522 M0III'      /                                                
NAME26  = 'SAO98230 C5'        /                                                
NAME27  = 'SAO152427 C6'       /                                                
NAME28  = 'SAO162465 C6'       /                                                
NAME29  = 'SAO96548 C7'        /                                                
NAME30  = 'SAO69636 C8'  
'''

#will contain the number of the templates Karrie is using based on the dictionary above
RGB_temps = [24, 13, 0, 28, 12, 7, 30, 17, 23, 11, 20, 26, 14, 8, 25, 9, 6, 18, 22, 15, 10, 16, 29, 3, 21, 27]

#arbitrary wavelength units -- doesn't matter since just looking at the zeropoint
x_values = range(len(data[0]))

#first get the wavelength
coe0 = 3.6442 
coe1 = 3.6189e-05
def wavelength(pixnum):  
	exp = coe0 + pixnum * coe1 
	return 10**exp 

#array of pixel numbers
indices0 = np.linspace(0, 9381, 9382)
wavelengths = wavelength(indices0)

#function to plot all of the templates on top of each other
def plot_temp(temp_num):
	plt.plot(wavelengths, data[temp_num], alpha = .2)
	return 

def plot_single_temp(temp_num):
	y_vals = data[temp_num]
	#plt.plot([0, 9400], [1, 1], c='b')
	plt.plot(wavelengths, y_vals, alpha = .5, c='k')

	# zero_point = y_vals == 0
	# print(sum(zero_point))
	# plt.scatter(x_values[zero_point], y_vals[zero_point], s=2, c='k')

	plt.xlabel('Wavelength')
	#plt.ylim(.95, 1.05)
	plt.savefig('/Users/amandaquirk/Desktop/template_{}.png'.format(temp_num))
	plt.close()
	return 

# for num in RGB_temps:
# 	plot_single_temp(num)

#plt.plot([0, 9400], [1, 1], c='k')
# plt.xlabel('Wavelength')
# #plt.ylim(.8, 1.2)
# plt.savefig('/Users/amandaquirk/Desktop/RGB_Templates.png')
# plt.close()

#cross correlating the templates to see if the spectral features line up ==============================================
shifts = np.linspace(-2, 2, 201)
def ccf(temp_num1, temp_num2):
	plt.clf()
	temp1 = data[temp_num1]
	temp2 = data[temp_num2]
	#temp2_2d = np.vstack((wavelengths, temp2))
	#temp2_2d = temp2_2d.T

	ccs = []
	#want to shift one of the templates to evaluate if the templates are aligned; the .2 here is the step size of shifts
	for i in range(len(shifts)):
		if shifts[i] < 0:
			shitfted_temp = list(temp1[abs(int(round(shifts[i] /.2))): len(temp1)]) 
			for j in range(int(round(shifts[i] / .2) - 1), len(temp1)): 
				if len(shitfted_temp) < len(temp1):
					shitfted_temp.append(1)
		elif shifts[i] == 0:
			shitfted_temp = temp1
		else:
			shitfted_temp = np.ones_like(temp1)
			indices = np.linspace(0, len(temp1), len(temp1) + 1) 
			shitfted_indices = indices - round(shifts[i] /.2) 
			for j in range(len(shitfted_indices) - 1):
				if shitfted_indices[j] >= 0:
					shitfted_temp[j] = temp1[int(shitfted_indices[j])]
		ccs.append(np.correlate(temp2, shitfted_temp)) #1D ignoring wavelength since all on same scale
		#temp1_2d = np.vstack((wavelengths, shitfted_temp))
		#temp1_2d = temp1_2d.T
		#ccs.append(signal.correlate2d(temp2_2d, temp1_2d))#2D not ignoring wavelength so can get the velocity shift

	n = np.argmax(ccs)
	#calculate velocity shift
	shift_need = shifts[n]
	v = const.c.to(u.km/u.s) * (10**(shift_need * coe1) - 1)
	velocity_shift = const.c.to(u.km/u.s) * (10**(shifts * coe1) - 1)
	plt.plot(velocity_shift.value, ccs, label='Max CCF at shift {}'.format(v))
	plt.xlabel('Velocity Shift (km/s)')
	plt.ylabel('CCF')
	plt.legend(frameon=False, loc=2)
	plt.savefig('/Users/amandaquirk/Desktop/CCFs/{}_{}.png'.format(temp_num1, temp_num2))
	plt.close()

	if abs(v.value) > 10:
		print(temp_num1, temp_num2, v)
	return

#want to cross match all of the templates against each other (and themselves for a sanity check)
for i in range(len(RGB_temps)):
 	for j in range(len(RGB_temps)):
 		if j >= i:
 			ccf(RGB_temps[i], RGB_temps[j])


