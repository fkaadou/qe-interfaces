'''
Script to resubmit stopped or expired calculations using latest atomic positions
in geometry optimization calculation to schedular on Compute Canada cluster.
'''

import re
import os
import numpy as np

# item in scans corresponds to a folder 'scan*' example: scan1
scans = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
#scans = [1,3,4,5,11]
#scans = 1 + np.arange(11)

for i in scans:

    di = 'scan' + str(i)
    print(di)
    os.chdir(di)
    
    runs = [int(dr[-3:]) for dr in os.listdir() if 'relax.out_' in dr]
    if not runs:
    	runs=[-1]

    print(f'Run number {max(runs)+1}')

    os.system(f'cp relax.in relax.in_{str(max(runs)+1).zfill(3)}')
    os.system(f'cp relax.out relax.out_{str(max(runs)+1).zfill(3)}')
    os.system('rm core*')
    print('core files removed.')

    # get number of atoms
    relax_in = open('relax.in', 'r')
    nat_line = [line for line in relax_in if 'nat' in line]
    nat = int(re.search('\d{1,}', str(nat_line)).group())
    relax_in.close()

    # read out latest coords from relax.out
    relax_out = open('relax.out', 'r')  
    lines = relax_out.readlines()

    index_list = np.zeros([])

    for i, line in enumerate(lines):
        if 'ATOMIC_POSITION' in line:
            index_list = np.append(index_list, i)

    max_index = int(np.amax(index_list))

    new_coords = lines[max_index:max_index+nat+1]

    relax_out.close()

    if max_index != 0:
        # update relax.in with latest coords
        relax_in = open('relax.in', 'r+')
        lines = relax_in.readlines()

        for i, line in enumerate(lines):
            if 'ATOMIC_POSITION' in line:
                upto = i
                break
                
        relax_in.seek(0)
        relax_in.truncate()
        relax_in.writelines(lines[:upto])
        relax_in.writelines(new_coords)
        relax_in.writelines(lines[upto+nat+1:])

        relax_in.close() 
    
    else: print('***no new coords. relax.in not updated.***')
    
    os.system('sbatch relax.sh')
    print('\n')
    os.chdir('..')
