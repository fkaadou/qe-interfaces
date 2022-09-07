import sys
import os
import numpy as np

for i in range(25):
    iteration = 'scan' + str(i)
    os.chdir(iteration)
    print(iteration)
    os.system('rm core*')
    os.system('sbatch relax.sh')
    os.chdir('..')

