""" make_scf_param(): Creates a Quantum Espresso input file for a self-consistent calculation

REQUIRED INPUTS
filename :       Name of SCF input file to be created
prefix :         Name given to charge density file output by SCF calculation

OPTIONAL INPUTS
pseudo_dir : Directory in which required pseudopotential files are located
outdir     : Directory into which output files from SCF calculation will be written
conv_thr   : Threshold dictating acceptable convergence in SCF cycle
ecutwfc    : Cutoff energy for wavefunctions
ecutrho    : Cutoff energy for charge density
nbnd       : Number of bands to include in SCF calculation

"""
import os
import sys

def make_scf_param(filename, calculation = 'scf', prefix='GS', restart_mode = 'from_scratch', wf_collect = '.false.', pseudo_dir = '../PP/', outdir = './', verbosity ='high', tprnfor = '.true.', tstress= '.true.', diagonalization = 'david', mixing_beta = 0.7, conv_thr = 1e-8, ecutwfc = 110, ecutrho = 1320, occupations = 'fixed', noncolin = '.true.', lspinorb = '.true.', london = '.true.', nbnd = 20, nat = 2, ntyp = 2, ibrav = 0):


    f = open(filename, 'w+')
    f.write(' &control\n')
    f.write("   calculation     = '" + calculation + "'\n" )
    f.write("   prefix          = '" + prefix + "'\n" )
    f.write("   restart_mode    = '" + restart_mode + "'\n" )
    f.write('   wf_collect      = ' + wf_collect + "\n" )
    f.write("   pseudo_dir      = '" + pseudo_dir + "'\n" )
    f.write("   outdir          = '" + outdir + "'\n" )
    f.write("   verbosity       = '" + verbosity + "'\n" )
    f.write('   tprnfor         = ' + tprnfor + "\n" )
    f.write('   tstress         = ' + tstress + "\n" )
    f.write(' /\n')

    f.write(' &system\n')
    f.write('   ecutwfc     =' +str(ecutwfc) +'\n')
    f.write('   ecutrho     =' +str(ecutrho) +'\n')
    f.write('   nbnd        =' +str(nbnd) +'\n')
    f.write('   nat         =' +str(nat) +'\n')
    f.write('   ntyp        =' +str(ntyp) +'\n')
    f.write('   ibrav       =' +str(ibrav) +'\n')
    f.write('   noncolin    = ' + noncolin + "\n" )
    f.write('   lspinorb    = ' + lspinorb + "\n" )
    f.write('   london      = ' + london + "\n" )
    f.write("   occupations     = '" + occupations +"'\n" )
    f.write(' /\n')
    
    f.write(' &electrons\n')
    f.write("   diagonalization = '" + diagonalization + "'\n")
    f.write('   mixing_beta     =' +str(mixing_beta) +'\n')
    f.write('   conv_thr        =' +str(conv_thr) +'\n')
    f.write(' /\n')
    f.close()

def make_scf_param_metal(filename, calculation = 'scf', prefix='GS', restart_mode = 'from_scratch', wf_collect = '.false.', pseudo_dir = '../PP/', outdir = './', verbosity ='high', tprnfor = '.true.', tstress= '.true.', diagonalization = 'david', mixing_mode = 'plain', mixing_beta = 0.7, conv_thr = 1e-8, ecutwfc = 110, ecutrho = 1320, occupations = 'smearing', noncolin = '.false.', lspinorb = '.false.', vdw_corr = 'xdm', xdm_a1 = 0.3275, xdm_a2 = 2.7673, degauss = 1.0e-3, smearing = 'gauss', nbnd = 20, nat = 2, ntyp = 2, ibrav = 0):


    f = open(filename, 'w+')
    f.write(' &control\n')
    f.write("   calculation     = '" + calculation + "'\n" )
    f.write("   prefix          = '" + prefix + "'\n" )
    f.write("   restart_mode    = '" + restart_mode + "'\n" )
    f.write('   wf_collect      = ' + wf_collect + "\n" )
    f.write("   pseudo_dir      = '" + pseudo_dir + "'\n" )
    f.write("   outdir          = '" + outdir + "'\n" )
    f.write("   verbosity       = '" + verbosity + "'\n" )
    f.write('   tprnfor         = ' + tprnfor + "\n" )
    f.write('   tstress         = ' + tstress + "\n" )
    f.write(' /\n')

    f.write(' &system\n')
    f.write('   ecutwfc     =' +str(ecutwfc) +'\n')
    f.write('   ecutrho     =' +str(ecutrho) +'\n')
    f.write('   nbnd        =' +str(nbnd) +'\n')
    f.write('   nat         =' +str(nat) +'\n')
    f.write('   ntyp        =' +str(ntyp) +'\n')
    f.write('   ibrav       =' +str(ibrav) +'\n')
    f.write('   noncolin    = ' + noncolin + "\n" )
    f.write('   lspinorb    = ' + lspinorb + "\n" )
    f.write("   occupations = '" + occupations +"'\n" )
    f.write('   degauss     = ' + str( degauss) + '\n' )
    f.write("   smearing    = '" + smearing +"'\n")
    f.write("   vdw_corr    = '" + vdw_corr +"'\n")
    f.write('   xdm_a1      =' +str(xdm_a1) + '\n ')
    f.write('   xdm_a2      =' +str(xdm_a2) + '\n ')
    f.write(' /\n')
    
    f.write(' &electrons\n')
    f.write("   diagonalization = '" + diagonalization + "'\n")
    f.write('   mixing_beta     =' +str(mixing_beta) +'\n')
    f.write('   conv_thr        =' +str(conv_thr) +'\n')
    f.write("   mixing_mode     = '" + mixing_mode + "'\n")
    f.write(' /\n')
    f.close()


def make_vc_relax_param(filename, calculation = 'vc-relax', prefix='GS', restart_mode = 'from_scratch', wf_collect = '.false.', pseudo_dir = '../PP/', outdir = './', verbosity ='high', tprnfor = '.true.', tstress= '.true.', etot_conv_thr = 1e-5, forc_conv_thr = 1e-4, nstep = 100, electron_maxstep = 100, mixing_mode = 'plain', diagonalization = 'david', mixing_beta = 0.7, conv_thr = 1e-8, ecutwfc = 120, ecutrho = 800, occupations = 'smearing', degauss = 1.0e-3 , smearing = 'gauss', vdw_corr = 'xdm', xdm_a1 = 0.3275, xdm_a2 = 2.7673, noncolin = '.false.', lspinorb = '.false.', london = '.false.', nbnd = 82, nat = 5, ntyp = 2, ibrav = 0,ion_dynamics = 'bfgs', pot_extrapolation = 'atomic', cell_dynamics = 'bfgs', cell_dofree = 'all', cell_factor = 2.0, magnet = False, mag = 0, nspin=2, mixing_ndim = 8):

    f = open(filename, 'w+')
    f.write(' &control\n')
    f.write("   calculation     = '" + calculation +"'\n" )
    f.write("   prefix          = '" + prefix + "'\n" )
    f.write("   restart_mode    = '" + restart_mode + "'\n" )
    f.write('   wf_collect      = ' + wf_collect + "\n" )
    f.write("   pseudo_dir      = '" + pseudo_dir + "'\n" )
    f.write("   outdir          = '" + outdir + "'\n" )
    f.write("   verbosity       = '" + verbosity + "'\n" )
    f.write('   tprnfor         = ' + tprnfor + "\n" )
    f.write('   tstress         = ' + tstress + "\n" )
    f.write('   etot_conv_thr   = ' + str(etot_conv_thr) + '\n')
    f.write('   forc_conv_thr   = ' + str(forc_conv_thr) + '\n')
    f.write('   nstep           = ' + str(nstep) + '\n')
    f.write(' /\n')

    f.write(' &system\n')
    f.write('   ecutwfc     =' +str(ecutwfc) +'\n')
    f.write('   ecutrho     =' +str(ecutrho) +'\n')
    f.write('   nbnd        =' +str(nbnd) +'\n')
    f.write('   nat         =' +str(nat) +'\n')
    f.write('   ntyp        =' +str(ntyp) +'\n')
    f.write('   ibrav       =' +str(ibrav) +'\n')
    f.write('   noncolin    = ' + noncolin + "\n" )
    f.write('   lspinorb    = ' + lspinorb + "\n" )
    f.write('   london      = ' + london + "\n" )
    f.write("   occupations = '" + occupations +"'\n" )
    f.write('   degauss     = ' + str( degauss) + '\n' )
    f.write("   smearing    = '" + smearing +"'\n")
    f.write("   vdw_corr    = '" + vdw_corr +"'\n")
    f.write('   xdm_a1      =' +str(xdm_a1) + '\n ')
    f.write('   xdm_a2      =' +str(xdm_a2) + '\n ') 
    if magnet == True:
        f.write('   nspin      =' +str(nspin) + '\n ')
        f.write('   starting_magnetization(1)      =' +str(mag) + '\n ') 
    f.write(' /\n')
    
    f.write(' &electrons\n')
    f.write("   diagonalization = '" + diagonalization + "'\n")
    f.write('   mixing_beta     =' +str(mixing_beta) +'\n')
    f.write('   conv_thr        =' +str(conv_thr) +'\n')
    f.write('   electron_maxstep =' +str(electron_maxstep) + '\n')
    f.write("   mixing_mode       ='" + mixing_mode + "'\n")
    f.write('   mixing_ndim       =' + str(mixing_ndim) + '\n')
    f.write(' /\n')
    
    f.write(' &IONS\n')
    f.write("   ion_dynamics = '" + ion_dynamics +"'\n")
    f.write("   pot_extrapolation = '" + pot_extrapolation +"'\n")
    f.write(' /\n')

    f.write(' &CELL\n')
    f.write("   cell_dynamics = '" + cell_dynamics +"'\n")
    f.write("   cell_dofree         = '" + cell_dofree + "'\n")
    f.write('   cell_factor   = ' + str(cell_factor) + '\n')
    f.write(' /\n')
    
    f.close()


def make_relax_param(filename, calculation = 'relax', prefix='GS', restart_mode = 'from_scratch', wf_collect = '.false.', pseudo_dir = '../PP/', outdir = './', verbosity ='high', tprnfor = '.true.', tstress= '.true.', etot_conv_thr = 1e-4, forc_conv_thr = 1e-3, nstep = 100, electron_maxstep = 100, mixing_mode = 'plain', diagonalization = 'david', mixing_beta = 0.7, conv_thr = 1e-8, ecutwfc = 120, ecutrho = 800, occupations = 'smearing', degauss = 1.0e-3 , smearing = 'gauss', vdw_corr = 'xdm', xdm_a1 = 0.3275, xdm_a2 = 2.7673, nosym='.false.' , noncolin = '.false.', lspinorb = '.false.', london = '.false.', nbnd = 82, nat = 5, ntyp = 2, ibrav = 0,ion_dynamics = 'bfgs', pot_extrapolation = 'atomic', cell_dynamics = 'bfgs', cell_factor = 2.0):

    f = open(filename, 'w+')
    f.write(' &control\n')
    f.write("   calculation     = '" + calculation +"'\n" )
    f.write("   prefix          = '" + prefix + "'\n" )
    f.write("   restart_mode    = '" + restart_mode + "'\n" )
    f.write('   wf_collect      = ' + wf_collect + "\n" )
    f.write("   pseudo_dir      = '" + pseudo_dir + "'\n" )
    f.write("   outdir          = '" + outdir + "'\n" )
    f.write("   verbosity       = '" + verbosity + "'\n" )
    f.write('   tprnfor         = ' + tprnfor + "\n" )
    f.write('   tstress         = ' + tstress + "\n" )
    f.write('   etot_conv_thr   = ' + str(etot_conv_thr) + '\n')
    f.write('   forc_conv_thr   = ' + str(forc_conv_thr) + '\n')
    f.write('   nstep           = ' + str(nstep) + '\n')
    f.write(' /\n')

    f.write(' &system\n')
    f.write('   ecutwfc     =' +str(ecutwfc) +'\n')
    f.write('   ecutrho     =' +str(ecutrho) +'\n')
    f.write('   nbnd        =' +str(nbnd) +'\n')
    f.write('   nat         =' +str(nat) +'\n')
    f.write('   ntyp        =' +str(ntyp) +'\n')
    f.write('   ibrav       =' +str(ibrav) +'\n')
    f.write('   noncolin    = ' + noncolin + "\n" )
    f.write('   lspinorb    = ' + lspinorb + "\n" )
    f.write('   london      = ' + london + "\n" )
    f.write("   occupations = '" + occupations +"'\n" )
    f.write('   degauss     = ' + str( degauss) + '\n' )
    f.write("   smearing    = '" + smearing +"'\n")
    if vdw_corr == 'xdm':
        f.write("   vdw_corr    = '" + vdw_corr +"'\n")
        f.write('   xdm_a1      =' +str(xdm_a1) + '\n ')
        f.write('   xdm_a2      =' +str(xdm_a2) + '\n ') 
    f.write('   nosym       =' + nosym + '\n')
    f.write(' /\n')
    
    f.write(' &electrons\n')
    f.write("   diagonalization = '" + diagonalization + "'\n")
    f.write('   mixing_beta     =' +str(mixing_beta) +'\n')
    f.write('   conv_thr        =' +str(conv_thr) +'\n')
    f.write('   electron_maxstep =' +str(electron_maxstep) + '\n')
    f.write("   mixing_mode       ='" + mixing_mode + "'\n")
    f.write(' /\n')
    
    f.write(' &IONS\n')
    f.write("   ion_dynamics = '" + ion_dynamics +"'\n")
    f.write("   pot_extrapolation = '" + pot_extrapolation +"'\n")
    f.write(' /\n')

    f.write(' &CELL\n')
    f.write("   cell_dynamics = '" + cell_dynamics +"'\n")
    f.write('   cell_factor   = ' + str(cell_factor) + '\n')
    f.write(' /\n')
    
    f.close()



    #os.system('dos2unix ' + filename)


if __name__ == "__main__":
    filename = sys.argv[1]
    prefix = sys.argv[2]
    make_scf_param(filename)
