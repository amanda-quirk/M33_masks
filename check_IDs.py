import numpy as np 

#read in the data
data = np.genfromtxt('/Users/amandaquirk/Desktop/K1M33P.out', dtype=None, names='ID, ras, decs, frame, mag, filter, priority, list, selection')

ID = data['ID']
RA = data['ras'] #does this work or nah? break up into ra and dec?
Dec = data['decs']
coordinate_reference_frame = data['frame']
magnitude = data['mag']
filter_tag = data['filter']
priority = data['priority']
list_assignment = data['list']
selection_flag = data['selection']

total = []
RGB = []
MS = []
AGB = []
RHeb = []
BHeb = []
CFHT = []
Heb = []
Pne = []
KG = []
for i in range(len(ID)):
	if selection_flag[i] == 1:
		total.append(i)

		if ID[i].startswith(b'RGB') == True:
			RGB.append(i)
		elif ID[i].startswith(b'MS') == True:
			MS.append(i)
		elif ID[i].startswith(b'AGB') == True:
			AGB.append(i)
		elif ID[i].startswith(b'RHeB') == True:
			RHeb.append(i)
		elif ID[i].startswith(b'BHeB') == True:
			BHeb.append(i)
		elif ID[i].startswith(b'HeB') == True:
			Heb.append(i)
		elif ID[i].startswith(b'CF') == True:
			CFHT.append(i)
		elif ID[i].startswith(b'PN') == True:
			Pne.append(i)
		elif ID[i].startswith(b'KG') == True:
			KG.append(i)

RGB_per = len(RGB) / len(total) * 100
MS_per = len(MS) / len(total) * 100
AGB_per = len(AGB) / len(total) * 100
RHeb_per = len(RHeb) / len(total) * 100
BHeb_per = len(BHeb) / len(total) * 100
CFHT_per = len(CFHT) / len(total) * 100
Heb_per = len(Heb) / len(total) * 100
Pne_per = len(Pne) / len(total) * 100
KG_per = len(KG) / len(total) * 100

print('total = {}'.format(len(total)))
print('RGB = {}, {} %'.format(len(RGB), RGB_per))
print('MS = {}, {} %'.format(len(MS), MS_per))
print('AGB = {}, {} %'.format(len(AGB), AGB_per))
print('RHeB = {}, {} %'.format(len(RHeb), RHeb_per))
print('BHeB = {}, {} %'.format(len(BHeb), BHeb_per))
print('CFHT = {}, {} %'.format(len(CFHT), CFHT_per))
print('HeB = {}, {} %'.format(len(Heb), Heb_per))
print('PNe = {}, {} %'.format(len(Pne), Pne_per))
print('KG = {}, {} %'.format(len(KG), KG_per))





