#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Plots the bandstructure of a heterostructure constructed from three materials with 
the contribution from each material to each state given using a colorbar.
Initially developed for metal-Ca2N-MoS2 heterostructures but can be adapted to others.

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

path = os.path.dirname(__file__)
metal = 'Au'
# 
path_data = '../' + metal + '_Ca2N_MoS2/pdos'

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
nkpt = 121   # number of kpnts (x-axis) depends on your choice during band calculation ***

k_vector = [line for line in pdos if "k = " in line]
E = [line for line in pdos if "e = " in line]
E = E[1:nbnd+1]

# number of paths in k-space (calculation file dependent) ***
npaths = 3    
# number of points in k-space per segment (calculaiton file dependent) ***
pnts = 40
a = npaths*pnts+1
k = bands_file[0:a,0]

# extract the energy bands and sort by band number
bands = np.zeros((nbnd, a))
for i in range(nbnd):
    bands[i,:] = bands_file[(a*i):(a*(i+1)), 1]

# extract wavefunction contribution of each state
wfc = ''
all_wfc = []
all_norm = []

for line in pdos:
    if 'psi = ' in line:
        wfc = line
    elif '|psi|' in line:
        all_wfc.append(wfc)
        
        wfc_norm = 1e-3*float(re.search('\d{3,}', line).group())
        
        if wfc_norm > 0:
            all_norm.append(wfc_norm)
        else:
            all_norm.append(1.0)
    else:
        wfc = wfc + line   
        
# nbnds x kpnts   
wfc_cont = np.empty((nbnd, nkpt, 3))

for i, psi in enumerate(all_wfc):
    psi_split = psi.split('+')
    comp_metal=0
    comp_ca2n=0
    comp_mos2=0
    for x in psi_split:
        if re.search('\d{1,}\]', x):
            comp_state = int(re.search('\d{1,}\]', x).group()[:-1])
            comp = 1e-3*float(re.search('\d{3,}',x).group())
            
            # ***
            # state contributions are composed from individual atomic states
            # find up to which numbered atomic state corresponds to the first material
            # in order to properly sum up contribution from said material
            if comp_state < 865:
                comp_metal = comp + comp_metal
            elif comp_state > 865 and comp_state < 1045:
                comp_ca2n = comp + comp_ca2n 
            else:
                comp_mos2 = comp + comp_mos2 
    
    # atomic contribution of 3 materials is given as an RGB value 
    wfc_cont[(i)%nbnd, (i)//nbnd,:] = [comp_metal, comp_ca2n, comp_mos2]/np.linalg.norm([comp_metal, comp_ca2n, comp_mos2])
#   wfc_cont[(i)%nbnd, (i)//nbnd,:] = [comp_metal, comp_ca2n, comp_mos2]        

print('contribution calculations done.')
    
pdos_file.close()

fs=28

fig, axs = plt.subplots(figsize=(12,8))

for i in range(nbnd_plot):
    plt.scatter(k, bands[i,:]-Ef, s=26, c=wfc_cont[i,:])

# plot Fermi level (0)
plt.hlines(0, k[0], k[len(k)-1], color='k', linestyles='dashed', linewidth=3)

# plot vertical lines to separate kspace segments ***
plt.vlines(k[0],   -5, 5)
plt.vlines(k[40],  -5, 5)
plt.vlines(k[80], -5, 5)
plt.vlines(k[120], -5, 5)

plt.text(0.05,0.05,'$E_f$',fontsize=fs)

# ploting parameters
ymin = -2
ymax =  1
l = 0.2

# labels symmetry points on x-axis ***
plt.text(k[0]-0.02,        ymin-l,'$K$',fontsize=fs)
plt.text(k[40]-0.035,      ymin-l,'$\Gamma$',     fontsize=fs)
plt.text(k[80]-0.05,       ymin-l,'$M$',     fontsize=fs)
plt.text(k[len(k)-1]-0.05, ymin-l,'$K$',fontsize=fs)

plt.xlim(k[0], k[len(k)-1])
plt.ylim(ymin, ymax)
plt.xlabel('k',fontsize=fs)
plt.ylabel('$E-E_f$ (eV)', fontsize=fs)

axs.get_xaxis().set_visible(False)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.show()

#plt.savefig(fname = '../outline/figures/bands_au_ca2n_mos2.png', format='png')
