'''
Plots the z-axis (perpendicular to material) average of a given cube file. 
Useful for investigating potentials accross a heterostructure.
Initially developed for metal-Ca2N-MoS2 heterostructures but can be adapted to others.

*REQUIRED FILES:
    cube file containg V(x,y,x) across cell. Can be calculated using pp.x (ex: ionic_hartree_potential.cube)
    pw.x output file (ex: scf.out)
    
calculation specific parameters to be modified are labelled with ***
''' 

import re
import os
import numpy as np
import matplotlib.pyplot as plt

path_current = os.path.dirname(__file__)

# location of scf pw.x output file
path_scf_out = os.path.relpath('../input_file_examples/scf.out', path_current)
# figure save location
path_save    = os.path.relpath('./')
# location of one or multiple cube files to be averaged and plotted
sources = ['../input_file_examples/ionic_hartree_potential.cube']

# labels of one or multiple averaged cube files
labels = ['ionic + hartree potential']


# get Fermi energy from pw.x output
scf_out = open(path_scf_out, 'r')
fermi = [line for line in scf_out if 'Fermi energy' in line]
Ef = float(re.search('[-+]?\d{1,}\.\d{1,}', str(fermi)).group())
scf_out.close()

# get parameters from cube flie
file = open(os.path.relpath(sources[0], path_current))
cube00 = file.readlines()
# num of atoms in cell
natom = int(cube00[2].split(' ')[3])
# size of V(x,y,z) file (Nx,Ny,Nz)
Nx = int(cube00[3].split(' ')[2])
Ny = int(cube00[4].split(' ')[2])
Nz = int(cube00[5].split(' ')[2])
# spatial spacing (dx,dy,dz)
dx = float(cube00[3].split(' ')[6])
dy = dx
dz = float(cube00[5].split(' ')[-1][:-1])

# get one or multiple V(x,y,z) from cube file(s)
for i, src in enumerate(sources):
    path = os.path.relpath(src, path_current)
    file = open(path,'r')
    exec('cube%s = file.readlines()[(5+natom+1):-1]' % (i))  
    exec('data%s = np.zeros((Nx,Ny,Nz))' % (i))
    exec('data%s_z = np.zeros((Nz,1))' % (i))

r = Nz//6
m = Nz%6

prog='''
if m != 0:
    temp1 = np.loadtxt(cube%s[k*r+k:(k+1)*r+k]) 
    b=cube%s[(k+1)*r+k][:-1].split(' ')[1:]
    a=np.array([float(x) for x in b if x != ''])
    temp2=np.append(a,np.zeros(6-m)).reshape(1,6)
    
    temp=np.concatenate((temp1,temp2),axis=0).reshape(Nz + 6 - m,1)
else:
    temp = np.loadtxt(cube%s[k*r:(k+1)*r]).reshape(Nz,1)
                   
for z in range(0,Nz): 
    data%s[x,y,z] = temp[z]''' 
    
# read data in xy-plane along z     
for  k in range(Nx*Ny-1):
    
    x = k//Nx
    y = k%Ny
    
    for i in range(len(sources)):
       exec( prog % (i,i,i,i))

# average data in z-diretion
for i in range(len(sources)):
    for k in range(Nz):
        exec('data%s_z[k] = np.mean(data%s[:,:,k])' % (i,i))
          
              
b_to_A = 0.529177               # bohr to angstrom
Ry_to_eV = 13.60566             # Ry to eV
z = b_to_A*dz*np.arange(Nz)     #z-axis

plt.figure(figsize=(8,6))
# plot all averaged potentials relative to Fermi level Ef
for i in range(len(sources)):
    exec('plt.plot(z, Ry_to_eV*data%s_z - Ef, label=labels[%s])' % (i,i))

plt.xlim(0,15)

plt.hlines(0, 0, 45, linestyles='dashed')

# veritical lines and labels to indicate atomic layer ***
plt.vlines(3.75,  -20, -5, linestyles='dashed')
plt.vlines(2.24, -20, -5, linestyles='dashed')
plt.vlines(5.27, -20, -5, linestyles='dashed')
plt.text(3.75-0.5,   -21, 'Mo',fontsize=10)
plt.text(2.24-0.5,  -21, 'S',fontsize=10)
plt.text(5.27-0.5,  -21, 'S',fontsize=10)
plt.vlines(10.1, -20, -5, linestyles='dashed')
plt.vlines(7.67,  -20, -5, linestyles='dashed')
plt.vlines(8.94,  -20, -5, linestyles='dashed')
plt.text(10.10-0.5, -21, 'Ca',fontsize=10)
plt.text(7.67-0.5,  -21, 'Ca',fontsize=10)
plt.text(8.94-0.5,  -21, 'N',fontsize=10)


plt.xlabel('z (A)',fontsize=14)
plt.ylabel('$V - Ef$ (eV)', fontsize=14)

plt.legend()
plt.show()

#plt.savefig(fname = path_save + '/potential_ionic_hartree.png', format='png')
