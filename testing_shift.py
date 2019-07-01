import numpy as np 

shifts = np.array((-.2, 0, .2, .4, -.4, -.6, .6))
temp1 = [1, 2, 3, 4, 5]
for i in range(len(shifts)):
	#print(i)
	if shifts[i] < 0:
		shitfted_temp = list(temp1[abs(int(round(shifts[i] /.2))): len(temp1)]) #added the /.2 here and removed the -1 and added the len(temp1)
		for j in range(int(round(shifts[i] / .2) - 1), len(temp1)): #added the /.2 here
			if len(shitfted_temp) < len(temp1):
				shitfted_temp.append(1)
	elif shifts[i] == 0:
		shitfted_temp = temp1
	else:
		shitfted_temp = np.ones_like(temp1)
		indices = np.linspace(0, len(temp1), len(temp1) + 1) 
		shitfted_indices = indices - round(shifts[i] /.2) #added the /.2
		for j in range(len(shitfted_indices) - 1):
			if shitfted_indices[j] >= 0:
				shitfted_temp[j] = temp1[int(shitfted_indices[j])]
	print(shifts[i], shitfted_temp)