'''
scrapes energy of final relaxed geometries.
'''


import sys
import os
import numpy as np
import re

energies = []

for i in range(25):
	# file names scan* example: scan1
        iteration = 'scan'+str(i)
        print(iteration)    
        os.chdir(iteration)
        
        relax_out = open('relax.out', 'r')
        lines = relax_out.readlines()
        relax_out.close()

        line = [line for line in lines if 'Final energy' in line]
        final_energy = float(re.search('[-+]?\d{1,}\.\d{1,}', line[0]).group()) 

        energies.append(final_energy)
        
        os.chdir('..')


scan_test_data = np.column_stack((np.arange(25), energies))
np.savetxt('scan_long', scan_test_data)

print('scan = ', range(25))
print('Energies = ', energies)
