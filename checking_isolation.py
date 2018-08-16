import numpy as np 

python_isolation = np.loadtxt('/Users/amandaquirk/Documents/M33/Data/CFHT_target_criteria.txt', usecols=(4), unpack=True) #0 means not isolated
julia_isolation = np.loadtxt('/Users/amandaquirk/Desktop/isolation/CFHT_julia_isolated.txt') #1 means not isolated 

added_together = python_isolation + julia_isolation
print(sum(added_together))
print(len(python_isolation))

