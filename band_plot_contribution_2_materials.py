#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Plots the bandstructure of a heterostructure constructed from two materials with 
the contribution from each material to each state given using a colorbar.
Initially developed for metal-MoS2 heterostructures but can be adapted to others.

*REQUIRED FILES:
    projwfc.x output file containged the projected atmomic states (ex: pdos.out)
    pw.x input file (ex: scf.in)
    pw.x output file (ex: scf.out)
    bands.x output file containing plotable bandstrucure (ex: bands.dat.gnu)
    
calculation specific parameters to be modified are labelled with ***
''' 

import re
import os
import numpy as np
import matplotlib.pyplot as plt
import time

tic = time.perf_counter()

path = os.path.dirname(__file__)
metal = 'Ag'

# location of all files
path_data = '../' + metal +'_MoS2/pdos'


################################################################################

def contribution(path_data):
    path_band    = os.path.relpath(path_data + '/bands.dat.gnu', path)
    path_pdos    = os.path.relpath(path_data + '/pdos.out', path)
    path_scf_in  = os.path.relpath(path_data + '/scf.in', path)
    path_scf_out = os.path.relpath(path_data + '/scf.out', path)
    
    
    bands_file = np.loadtxt(fname=path_band)
    
    # get number of bands    
    scf_in = open(path_scf_in, 'r') 
    nbnd = [line for line in scf_in if 'nbnd' in line]
    nbnd = int(re.search('\d{1,}', str(nbnd)).group())
    scf_in.close()
    # get Fermi level
    scf_out = open(path_scf_out, 'r')
    fermi = [line for line in scf_out if 'Fermi energy' in line]
    Ef = float(re.search('[-+]?\d{1,}\.\d{1,}', str(fermi)).group())
    scf_out.close()
    
    pdos_file = open(path_pdos, 'r')
    
    pdos = pdos_file.readlines()
    
    nbnd_plot = nbnd
    nkpt = 211      # number of kpnts (x-axis) depends on your choice during band calculation ***
    
    E = [line for line in pdos if "e = " in line]
    E = E[1:nbnd+1]
    
    # number of paths in k-space (calculation file dependent) ***
    npaths = 3      
    # number of points in k-space per segment (calculaiton file dependent) ***
    pnts = 70
    a = npaths*pnts+1
    k = bands_file[0:a,0]
    
    tic = time.perf_counter()
    # extract the energy bands and sort by band number
    bands = np.zeros((nbnd, a))
    for i in range(nbnd):
        bands[i,:] = bands_file[(a*i):(a*(i+1)), 1]
    toc = time.perf_counter()
    print(f"Elapsed time get bands{toc-tic:0.4f}")
    
    
    # extract wavefunction contribution of each state
    tic = time.perf_counter()
    wfc = ''
    all_wfc = []
    
    for line in pdos:
        if 'psi = ' in line:
            wfc = line
        elif '|psi|' in line:
            all_wfc.append(wfc)
        else:
            wfc = wfc + line   
    
       
    toc = time.perf_counter()
    print(f"Elapsed time get wfcs:{toc-tic:0.4f}")
    
    tic = time.perf_counter()
    # nbnds x kpnts      
    metal_wfc_cont = np.empty((nbnd, nkpt))
            
    for i, psi in enumerate(all_wfc):
        psi_split = psi.split('+')
        comp_sum = 0
        for x in psi_split:
            if re.search('\d{1,}\]', x):
                comp_state = int(re.search('\d{1,}\]', x).group()[:-1])
                comp = 1e-3*float(re.search('\d{3,}',x).group())
                
                # ***
                # state contributions are composed from individual atomic states
                # find up to which numbered atomic state corresponds to the first material
                # in order to properly sum up contribution from said material
                if comp_state < 241:
                    comp_sum = comp_sum + comp
        
        metal_wfc_cont[(i)%nbnd, (i)//nbnd] = comp_sum
     
        
     
    toc = time.perf_counter()
    print(f"Elapsed time get component:{toc-tic:0.4f}")
        
    pdos_file.close()
    
    return k, bands, metal_wfc_cont, nbnd_plot, Ef

#########################################################################


k, bands, metal_wfc_cont, nbnd_plot, Ef = contribution(path_data)

fs = 30
fig, axs = plt.subplots(figsize=(12,8))

for i in range(nbnd_plot):
    plt.scatter(k, bands[i,:]-Ef, cmap='jet', c=metal_wfc_cont[i,:])
    plt.clim(0, 1)

cbar = plt.colorbar()
cbar.set_ticks([0,1])
cbar.set_ticklabels(['$MoS_2$', metal])
cbar.ax.tick_params(labelsize=26)

# Estimated CBE - Ef. Not always reliable. Chech on plot
print('CBE - Ef = ' + str(round((bands[:,70][bands[:,70]> Ef] - Ef)[0]*1000,3)) + ' meV')
print('^^^^^^^^^^  check  ^^^^^^^^^^^^^')
lw = 2

# plot Fermi level (0)
plt.hlines(0, k[0], k[len(k)-1], linestyles='dashed', linewidth=lw, color='k')

# plot vertical lines to separate kspace segments ***
plt.vlines(k[0],   -5, 5, color='k')
plt.vlines(k[70],  -5, 5, color='k', linewidth=lw)
plt.vlines(k[140], -5, 5, color='k', linewidth=lw)
plt.vlines(k[210], -5, 5, color='k')

plt.text(0.05,0.05,'$E_f$',fontsize=fs)

# plotting parameters
ymin = -2.5
ymax =  2.0
l = 0.25

# labels symmetry points on x-axis ***
plt.text(k[0]-0.02,        ymin-l,'$K$',fontsize=fs)
plt.text(k[70]-0.02,      ymin-l,'$\Gamma$',     fontsize=fs)
plt.text(k[140]-0.035,       ymin-l,'$M$',     fontsize=fs)
plt.text(k[len(k)-1]-0.02, ymin-l,'$K$',fontsize=fs)

plt.xlim(k[0], k[len(k)-1])
plt.ylim(ymin, ymax)
plt.xlabel('k',fontsize=fs)
plt.ylabel('$E-E_f$ (eV)', fontsize=fs)

axs.get_xaxis().set_visible(False)
plt.xticks(fontsize=fs-2)
plt.yticks(fontsize=fs-2)
axs.set_yticks([-2,-1,0,1,2])
axs.set_yticklabels([-2,-1,0,1,2])

plt.show()

toc = time.perf_counter()
print(f"Elapsed time:{toc-tic:0.4f}")

# plt.savefig(fname = path_data + '/pdos_bands.png', format='png')

