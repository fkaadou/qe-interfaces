'''
creates and submits geometry relaxation calculations with top material
translated along the long diagonal of a hexagonal cell relative to 
bottom material.
'''

import sys
import os
import numpy as np

sys.path.insert(0, "/home/scratch/makers")
from kpts_makers import make_auto_kpts
from param_makers import make_relax_param
from bash_makers import make_bash_graham
from struct_makers import make_Mo_Ca2N_var_z


k1 = np.array([12.70964067,  3.371e-05, 0.000000000])
k2 = np.array([-6.35479125,       11.00685249,  0.000000000])

q = -1*k1 + k2

make_auto_kpts('KPTS', grid = [2,2,1,0,0,0])

for i in range(13):
    iteration = 'scan' + str(i)
    os.system('mkdir ' + iteration)
    make_relax_param('PARAMS', pseudo_dir = '../../../../../PP/', ecutwfc = 120, ecutrho = 800, nat = 153, ntyp = 3, nbnd = 1210, conv_thr = 1e-6, mixing_beta=0.4)
    # tranlates top layer by (xd,xd). see maker in struct_maker.
    make_Mo_Ca2N_var_z('STRUCT', xd = (i/12)*q[0], yd = (i/12)*q[1]) 
    os.system('cat PARAMS STRUCT KPTS >> relax.in')
    os.system('dos2unix relax.in')
    os.system('mv relax.in  ' + iteration)

    make_bash_graham('relax.sh', qe_inputfile = 'relax.in', qe_outputfile = 'relax.out', nodes=3, job_name = 'Mo_Ca2N_rlx_S_' + iteration, time = '1-00:00:00')
    os.system('mv relax.sh ' + iteration)

    os.chdir(iteration)
    os.system('sbatch relax.sh')
    os.chdir('..')

os.system('rm -f PARAMS KPTS STRUCT')


