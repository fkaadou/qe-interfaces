import numpy as np
import sys 


def make_MoS2_mono(struc_filename, vac_layer = 15):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.960 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   3.1603999138' + '  0.0000000000  ' + '  0.0000000000\n')
    f.write('  -1.5801999569' + '  2.7369866115  ' + '  0.0000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Mo     0.000000000  1.824657798  3.073750019\n')
    f.write('S      1.580199838  0.912328839  1.586054325\n')
    f.write('S      1.580199838  0.912328839  4.561444759\n')
    f.close()

# atomic positions at +1% change
def make_MoS2_mono_super_1(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.960 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
     
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.301802312)   + '   ' + str(scale*0.000016714)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.150886734)  + '   ' + str(scale*5.457511372)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Mo     -0.008293297   1.825093408   3.727508092\n')
        f.write('S      1.581661617   0.904632653   2.168790429\n')
        f.write('S      1.581689979   0.904649392   5.286203407\n')
       
        f.write('Mo     3.170984191   1.829648382   3.727476782\n')
        f.write('S      4.761258251   0.911512233   2.168724253\n')
        f.write('S      4.761229519   0.911528299   5.286268822\n')
        
        f.write('Mo    -1.594046472   4.580724405   3.727508237\n')
        f.write('S      -0.002118552   3.661632434   2.168861212\n')
        f.write('S      -0.002117879   3.661631685   5.286175910\n')

        f.write('Mo     1.585283446   4.585256381   3.727508331\n')
        f.write('S      3.177399477   3.668546998   2.168790253\n')
        f.write('S      3.177400410   3.668513712   5.286203940\n')
    
    f.close()



def make_MoS2_mono_super_2(struc_filename, scale=1.0):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.960 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
     
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.301024880)   + '   ' + str(scale*0.000009649)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.150504087)  + '   ' + str(scale*5.456837176)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Mo     0.010713542   1.833657208   3.504155378\n')
    f.write('S      1.584231226   0.923095684   1.932318784\n')
    f.write('S      1.587517387   0.925138382   5.074023157\n')
   
    f.write('Mo     3.161019863   1.833914980   3.500960395\n')
    f.write('S      4.738196401   0.923710348   1.932317244\n')
    f.write('S      4.734762240   0.925548230   5.074012248\n')
    
    f.write('Mo    -1.565048438   4.562221205   3.504149569\n')
    f.write('S      0.010438499   3.652773507   1.933638763\n')
    f.write('S      0.010448942   3.652801112   5.074375878\n')

    f.write('Mo     1.585767313   4.562571095   3.504152268\n')
    f.write('S      3.160609330   3.654902803   1.932310782\n')
    f.write('S      3.160749993   3.650999588   5.074015310\n')
    
    f.close()



def make_MoS2_mono_super_3(struc_filename, scale=1.0):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.960 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
     
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.304282910)   + '   ' + str(scale*-0.000067833)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.152200199)  + '   ' + str(scale*5.459632941)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    41.000000000\n')
 
    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Mo     0.010485393   1.834407217   3.498412111\n')
    f.write('S      1.586381082   0.924414095   1.928046099\n')
    f.write('S      1.586782344   0.924686585   5.068378376\n')
   
    f.write('Mo     3.162503038   1.834595863   3.497932783\n')
    f.write('S      4.738825683   0.924692916   1.928050416\n')
    f.write('S      4.738381199   0.924906794   5.068375527\n')
    
    f.write('Mo    -1.565973807   4.564418514   3.498409434\n')
    f.write('S      0.010305947   3.654457449   1.928883028\n')
    f.write('S      0.010307925   3.654462859   5.068164257\n')

    f.write('Mo     1.586469697   4.564593643   3.498411175\n')
    f.write('S      3.162264080   3.654670157   1.928046369\n')
    f.write('S      3.162302982   3.654174204   5.068375996\n')
    
    f.close()




def make_MoS2_mono_super_4(struc_filename, scale=1.0):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.960 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
     
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.300927955)   + '   ' + str(scale*-0.000029454)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.150489485)  + '   ' + str(scale*5.456731895)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    45.000000000\n')
 
    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Mo     0.010586774   1.833161605   3.493328237\n')
    f.write('S      1.584164058   0.922626624   1.922143916\n')
    f.write('S      1.587471030   0.924620324   5.061975254\n')
   
    f.write('Mo     3.160726955   1.833501016   3.489972164\n')
    f.write('S      4.737827436   0.923475900   1.922143732\n')
    f.write('S      4.734436449   0.925348511   5.061968716\n')
    
    f.write('Mo    -1.565402723   4.561672118   3.493326032\n')
    f.write('S      0.010155498   3.652290858   1.924630185\n')
    f.write('S      0.010163920   3.652306519   5.062816144\n')

    f.write('Mo     1.585427824   4.562215162   3.493328040\n')
    f.write('S      3.160125150   3.654331602   1.922140435\n')
    f.write('S      3.160202494   3.650452244   5.061970785\n')
    
    f.close()


# atomic positions at +3% in xy and +5% in z
def make_Ca2N_bulk_hex(struc_filename, scale_xy=1.0, scale_z=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale_xy*1.728879606)  + '   ' + str(scale_xy*-2.994506898)  + '    0.000000000\n')
    f.write('   ' +  str(scale_xy*1.728879606)  + '   ' + str(scale_xy*2.994506898)   + '    0.000000000\n')
    f.write('     0.000000000      0.000000000      ' + str(scale_z*17.617661378)+ '\n')
 
    f.write('\n')

    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (crystal)\n')
        f.write('Ca     0.000000000  -0.000000000   0.265525260\n')
        f.write('Ca     0.000000000  -0.000000000   0.734474740\n')
        f.write('Ca     0.333333000   0.666667000   0.598861184\n')
        f.write('Ca     0.333333000   0.666667000   0.067807223\n')
        f.write('Ca    -0.333333000   0.333333000   0.932192777\n')
        f.write('Ca    -0.333333000   0.333333000   0.401137816\n')

        f.write('N      0.333333000   0.666667000   0.333331346\n')
        f.write('N     -0.333333000   0.333333000   0.666668654\n')
        f.write('N     -0.000000000  -0.000000000   0.000000000\n')
    
    f.close()




# atomic postions at +3% stretch
def make_Ca2N_bulk_prm(struc_filename, scale=1.0):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*-1.728802504) + '   ' + str(scale*0.998124591)  + '   ' + str(scale*5.878346101)+ '\n')
    f.write('   ' +  str(scale*1.728802736)  + '   ' + str(scale*0.998124588)  + '   ' + str(scale*5.878347052)+ '\n')
    f.write('   ' +  str(scale*0.000000115)  + '   ' + str(scale*-1.996249382) + '   ' + str(scale*5.878347052)+ '\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (crystal)\n')
    f.write('Ca     0.264693794   0.264693810   0.264693851\n')
    f.write('Ca     0.735306190   0.735306195   0.735306156\n')
    f.write('N      0.000000000  -0.000000000   0.000000000\n')
    f.close()

# atomic positions at +4% stretch
def make_Ca2N_mono_prm(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*3.465408795)   + '   ' + str(scale*0.00000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.732704398)  + '   ' + str(scale*3.001132051)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    20.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca    -0.024312497  -0.035864378   5.092951895\n')
        f.write('Ca    -0.024311186   2.036619297   7.555050535\n')
        f.write('N      1.781328080   1.000377131   6.323999047\n')
    
    f.close()

# atomic positions at +3% stretch
def make_Ca2N_bi_prm(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*3.472084049)   + '   ' + str(scale*0.00000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.736042025)  + '   ' + str(scale*3.006912990)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    25.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca     -0.020470766  -0.028078143   5.387247874\n')
        f.write('Ca     -0.016847474   2.030267730   7.852979741\n')
        f.write('N       1.771758703   1.005645372   6.508504989\n')

        f.write('Ca     -0.014080524  -0.031765977  13.584747850\n')
        f.write('Ca      1.766688841   1.000543954  11.119022906\n')
        f.write('N      -0.014964501   2.037212925  12.463498359\n')
    
    f.close()

# atomic postitions at +3% stretch
def make_Ca2N_tri_prm(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*3.467623579)   + '   ' + str(scale*0.00000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.733811790)  + '   ' + str(scale*3.003050110)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    33.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca     -0.017558272  -0.029951001   5.657801123\n')
        f.write('Ca     -0.017531676   2.032296737   8.122189289\n')
        f.write('N       1.768016771   1.001095905   6.794698157\n')

        f.write('Ca     -0.016878386  -0.030433282  13.871526054\n')
        f.write('Ca     1.768601967   1.000414852  11.424474995\n')
        f.write('N      -0.017144474   2.031722854  12.648002778\n')
 
        f.write('Ca     -0.017438702   2.032345245  17.173811953\n')
        f.write('Ca     1.768534760   1.001244220  19.638194531\n')
        f.write('N      -0.017166617  -0.029585196  18.501303563\n')    
   
    f.close()

#atomic positions at +3% stretch
def make_Ca2N_quad_prm(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*3.465088733)   + '   ' + str(scale*0.00000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.732544367)  + '   ' + str(scale*3.000854870)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    37.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca     -0.017507463  -0.029875885   5.916825565\n')
        f.write('Ca      -0.017405539   2.030231071   8.382278718\n')
        f.write('N       1.767006945   1.000363222   7.056699600\n')

        f.write('Ca     1.767325186   1.000010557  11.681418999\n')
        f.write('Ca     -0.017191132  -0.030210172  14.133129827\n')
        f.write('N      -0.017105269   2.030177391  12.922383342\n')
 
        f.write('Ca     -0.017206193   2.030751159  17.486839321\n')
        f.write('Ca     1.767312998   1.000538313  19.938578834\n')
        f.write('N      -0.017118096  -0.029626874  18.697580648\n')    
  
        f.write('Ca     -0.017415176  -0.029638331  23.237754886\n')
        f.write('Ca     -0.017517315   2.030469181  25.703196051\n')
        f.write('N      1.766998525   1.000229851  24.563318454\n')    
   
  
    f.close()



# atomic positions at +3%
def make_Ca2N_mono_super(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*5.999630868)   + '   ' + str(scale*0.000002115)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.999813651)  + '   ' + str(scale*5.195831457)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca     -0.028732759  -0.049762465   7.714853856\n')
        f.write('Ca     4.090751468  -0.052375465  10.190918128\n')
        f.write('Ca    -2.090740230   3.516508383  10.190918269\n')
        f.write('Ca     -0.030514384   3.514561894   7.711640508\n')
        f.write('Ca     1.000424109   1.732789344  10.194791336\n')
        f.write('Ca     3.058951750   1.730856944   7.711640421\n')
        f.write('N     -1.060265722   1.733668918   8.949621131\n')
        f.write('N      2.031529420  -0.051378493   8.949622008\n')
        f.write('N      2.029143595   3.514584480   8.956972891\n')
    
    f.close()


# atomic positions at +3% stretch
def make_Ca2N_bi_super(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.007049839)   + '   ' + str(scale*0.000027010)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.003501532)  + '   ' + str(scale*5.202279440)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')

    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca     -0.028993757  -0.050780750   7.610572400\n')
        f.write('Ca     4.095830194  -0.051799473  10.072708701\n')
        f.write('Ca    -2.092305363   3.520926630  10.072708465\n')
        f.write('Ca     -0.029608133   3.521078596   7.607265277\n')
        f.write('Ca     1.002857716   1.736456829  10.070767021\n')
        f.write('Ca     3.064631009   1.734616526   7.607264873\n')
        f.write('N     -1.061499257   1.736105024   8.733954286\n')
        f.write('N      2.034745157  -0.051515725   8.733950629\n')
        f.write('N      2.032889129   3.520510946   8.738782857\n')
   
        f.write('Ca     -0.029150718  -0.050555064  15.789650952\n')
        f.write('Ca    -1.062668036   1.735115661  13.327565276\n')
        f.write('Ca     2.034057508  -0.052781839  13.327565435\n')
        f.write('Ca     -0.029360099   3.519058286  15.791471645\n')
        f.write('Ca     2.032191444   3.519781539  13.326597643\n')
        f.write('Ca     3.062324141   1.734071587  15.791470962\n')
        f.write('N      1.002188975   1.735781179  14.662948101\n')
        f.write('N     -2.093068654   3.520118299  14.663661448\n')
        f.write('N      4.095097488  -0.052627100  14.663663189\n')
    
    f.close()




# atomic positions at +3% stretch
def make_Ca2N_tri_super(struc_filename, scale=1.0, include_atoms=True):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.000381506)   + '   ' + str(scale*-0.000529853)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.000649618)  + '   ' + str(scale*5.196226883)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    41.000000000\n')
 
    f.write('\n')

    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ca     -0.028500375  -0.049908832   7.570332200\n')
        f.write('Ca     4.091067160  -0.051450375  10.038176330\n')
        f.write('Ca    -2.089631910   3.516984556  10.038176088\n')
        f.write('Ca     -0.029462721   3.516579698   7.564992337\n')
        f.write('Ca     1.002101872   1.735165600  10.037452313\n')
        f.write('Ca     3.060644311   1.732501615   7.564995653\n')
        f.write('N     -1.060795495   1.735053574   8.701239959\n')
        f.write('N      2.033467040  -0.051425596   8.701239002\n')
        f.write('N      2.030028090   3.515575452   8.707964356\n')
       
        f.write('Ca     -0.026689984  -0.046321575  15.820789324\n')
        f.write('Ca    -1.057400728   1.736383589  13.365134394\n')
        f.write('Ca     2.032542909  -0.047598396  13.365134636\n')
        f.write('Ca     -0.027152062   3.520869028  15.820071500\n')
        f.write('Ca     2.032107182   3.519612400  13.364471030\n')
        f.write('Ca     3.062812915   1.736874324  15.820070595\n')
        f.write('N      1.002627573   1.736508011  14.592596834\n')
        f.write('N     -2.089417206   3.521725382  14.592611002\n')
        f.write('N      4.094691609  -0.048678244  14.592612461\n')
       
        f.write('Ca    -1.041992844   1.763382763  21.620340686\n')
        f.write('Ca    -2.072901716   3.546484605  19.147157958\n')
        f.write('Ca     2.048135116  -0.020706681  21.620340502\n')
        f.write('Ca     4.107800879  -0.021951522  19.147157894\n')
        f.write('Ca     1.016061195   1.759861300  19.147871970\n')
        f.write('Ca     2.047147083   3.545760805  21.615018465\n')
        f.write('N      -0.011415065  -0.019779570  20.477410722\n')
        f.write('N      -0.014860038   3.547225343  20.484078842\n')
        f.write('N      3.079418788   1.760739426  20.484078487\n')
   

    f.close()





def make_Ca2N_quad_super(struc_filename, scale=1.0):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*5.997035699)   + '   ' + str(scale*-0.000321006)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.998795846)  + '   ' + str(scale*5.193424578)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    45.000000000\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca     0.001222592   0.002163900   7.466864853\n')
    f.write('Ca     3.999140778   0.001126805   9.979865218\n')
    f.write('Ca    -1.997380102   3.464587115   9.979811213\n')
    f.write('Ca     0.001050719   3.463790636   7.466631797\n')
    f.write('Ca     1.000844811   1.733442742   9.979828651\n')
    f.write('Ca     3.000128610   1.732826300   7.466856838\n')
    f.write('N     -0.998113550   1.733158142   8.649214810\n')
    f.write('N      2.000474143   0.001588937   8.647709570\n')
    f.write('N      1.999996683   3.463967569   8.647744220\n')
   
    f.write('Ca     0.004930647   0.008362996  15.795903536\n')
    f.write('Ca    -0.994796017   1.739260119  13.295589180\n')
    f.write('Ca     2.003810342   0.007875860  13.295714151\n')
    f.write('Ca     0.004711880   3.470513241  15.796069831\n')
    f.write('Ca     2.003486075   3.470131849  13.295715140\n')
    f.write('Ca     3.003368166   1.739300864  15.795969301\n')
    f.write('N      1.004220474   1.739178813  14.548623036\n')
    f.write('N     -1.994396581   3.470379944  14.548568900\n')
    f.write('N      4.002831157   0.007994406  14.548664676\n')
   
    f.write('Ca    -0.982053235   1.761579377  21.679826286\n')
    f.write('Ca    -1.982511148   3.491968662  19.179278785\n')
    f.write('Ca     2.016613158   0.029599321  21.679832520\n')
    f.write('Ca     4.015402512   0.029760135  19.179272233\n')
    f.write('Ca     1.016723684   1.760306974  19.179266006\n')
    f.write('Ca     2.015767476   3.492091710  21.679978129\n')
    f.write('N      0.017248117   0.029857997  20.426446257\n')
    f.write('N      0.017117922   3.492062935  20.424304830\n')
    f.write('N      3.015668622   1.760842470  20.426238189\n')
   
    f.write('Ca     0.024370654   0.042216298  24.996492222\n')
    f.write('Ca     4.022618654   0.041663202  27.509818775\n')
    f.write('Ca    -1.975287527   3.504491837  27.509821664\n')
    f.write('Ca     0.023661293   3.504097769  24.996788178\n')
    f.write('Ca     1.023974804   1.773518810  27.509838472\n')
    f.write('Ca     3.022802302   1.772545208  24.996497269\n')
    f.write('N     -0.975270774   1.773289010  26.329354500\n')
    f.write('N      2.023343057   0.042112721  26.329355522\n')
    f.write('N      2.023206391   3.504225154  26.327478117\n')
   
    f.close()



def make_Ca2N_mono(struc_filename, vac_layer = 20):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   1.728879606' + '  -2.994506898  ' + '  0.000000000\n')
    f.write('   1.728879606' + '   2.994506898  ' + '  0.000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca     0.000000000  0.000000000  0.000000000\n')
    f.write('Ca     0.000000119  2.094107151  2.543133259\n')
    f.write('N      1.813549995  1.047053576  1.271566868\n')
    f.close()


def make_Ca2N_bi(struc_filename, vac_layer = 20):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   3.6270999908' + '    0.0000000000' + '     0.0000000000\n')
    f.write('  -1.8135499954' + '    3.1411607341' + '     0.0000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca       0.000000000000   0.000000000000   5.052434680556\n')
    f.write('Ca       0.000000000000   2.094107156067   7.595566754180\n')
    f.write('N        1.813549995400   1.047053578033   6.324000040700\n')

    f.write('Ca       0.000000000   0.000000000   13.919567108\n')
    f.write('Ca       1.813550234   1.047053695   11.376433372\n')
    f.write('N        0.000000000   2.094106913   12.647999763\n')
    f.close()



def make_Ca2N_tri(struc_filename, vac_layer = 20):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   3.6270999908' + '    0.0000000000' + '     0.0000000000\n')
    f.write('  -1.8135499954' + '    3.1411607341' + '     0.0000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca       0.000000000000   0.000000000000   5.052434680556\n')
    f.write('Ca       0.000000000000   2.094107156067   7.595566754180\n')
    f.write('N        1.813549995400   1.047053578033   6.324000040700\n')

    f.write('Ca       0.000000000   0.000000000   13.919567108\n')
    f.write('Ca       1.813550234   1.047053695   11.376433372\n')
    f.write('N        0.000000000   2.094106913   12.647999763\n')

    f.write('Ca        0.000000000000   2.094107156067   17.7004337310\n')
    f.write('Ca        1.813549995400   1.047053578033   20.2435668711\n')
    f.write('N	      0.000000000000   0.000000000000   18.9720001221\n')
    f.close()



def make_Ca2N_quad(struc_filename, vac_layer = 42):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   3.6270999908' + '    0.0000000000' + '     0.0000000000\n')
    f.write('  -1.8135499954' + '    3.1411607341' + '     0.0000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca       0.000000000000   0.000000000000   5.052434680556\n')
    f.write('Ca       0.000000000000   2.094107156067   7.595566754180\n')
    f.write('N        1.813549995400   1.047053578033   6.324000040700\n')

    f.write('Ca       0.000000000   0.000000000   13.919567108\n')
    f.write('Ca       1.813550234   1.047053695   11.376433372\n')
    f.write('N        0.000000000   2.094106913   12.647999763\n')

    f.write('Ca       0.000000000000   2.094107156067   17.7004337310\n')
    f.write('Ca       1.813549995400   1.047053578033   20.2435668711\n')
    f.write('N	      0.000000000000   0.000000000000   18.9720001221\n')
    
    f.write('Ca       0.000000000000   0.000000000000    24.0244340900\n')
    f.write('Ca       0.000000000000   2.094107156067    26.5675672301\n')
    f.write('N        1.813549995400   1.047053578033    25.2960004811\n')
    f.close()



def make_Ca2N_quin(struc_filename, vac_layer = 48):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   3.6270999908' + '    0.0000000000' + '     0.0000000000\n')
    f.write('  -1.8135499954' + '    3.1411607341' + '     0.0000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca       0.000000000000   0.000000000000   5.052434680556\n')
    f.write('Ca       0.000000000000   2.094107156067   7.595566754180\n')
    f.write('N        1.813549995400   1.047053578033   6.324000040700\n')

    f.write('Ca       0.000000000   0.000000000   13.919567108\n')
    f.write('Ca       1.813550234   1.047053695   11.376433372\n')
    f.write('N        0.000000000   2.094106913   12.647999763\n')

    f.write('Ca       0.000000000000   2.094107156067   17.7004337310\n')
    f.write('Ca       1.813549995400   1.047053578033   20.2435668711\n')
    f.write('N	      0.000000000000   0.000000000000   18.9720001221\n')
    
    f.write('Ca       0.000000000000   0.000000000000    24.0244340900\n')
    f.write('Ca       0.000000000000   2.094107156067    26.5675672301\n')
    f.write('N        1.813549995400   1.047053578033    25.2960004811\n')
 
    f.write('Ca       0.000000000   0.000000000   32.891567589\n')
    f.write('Ca       1.813550234   1.047053695   30.348434449\n')
    f.write('N        0.000000000   2.094106913   31.620000840\n')
    f.close()




def make_Ca2N_hex(struc_filename, vac_layer = 54):
    
    f = open(struc_filename, 'w+')
    
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')
    
    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   3.6270999908' + '    0.0000000000' + '     0.0000000000\n')
    f.write('  -1.8135499954' + '    3.1411607341' + '     0.0000000000\n')
    f.write('   0.0000000000' + '  0.0000000000  ' + '  ' + str(vac_layer)+ '\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca       0.000000000000   0.000000000000   5.052434680556\n')
    f.write('Ca       0.000000000000   2.094107156067   7.595566754180\n')
    f.write('N        1.813549995400   1.047053578033   6.324000040700\n')

    f.write('Ca       0.000000000   0.000000000   13.919567108\n')
    f.write('Ca       1.813550234   1.047053695   11.376433372\n')
    f.write('N        0.000000000   2.094106913   12.647999763\n')

    f.write('Ca       0.000000000000   2.094107156067   17.7004337310\n')
    f.write('Ca       1.813549995400   1.047053578033   20.2435668711\n')
    f.write('N	      0.000000000000   0.000000000000   18.9720001221\n')
    
    f.write('Ca       0.000000000000   0.000000000000    24.0244340900\n')
    f.write('Ca       0.000000000000   2.094107156067    26.5675672301\n')
    f.write('N        1.813549995400   1.047053578033    25.2960004811\n')
 
    f.write('Ca       0.000000000   0.000000000   32.891567589\n')
    f.write('Ca       1.813550234   1.047053695   30.348434449\n')
    f.write('N        0.000000000   2.094106913   31.620000840\n')

    f.write('Ca       0.000000000000   2.094107156067   36.6724348079\n')
    f.write('Ca       1.813549995400   1.047053578033   39.2155679481\n')
    f.write('N	      0.000000000000   0.000000000000   37.9440011991\n')
    f.close()




# atomic positions at +2% stretch
def make_mono_Ca2N_MoS2(struc_filename, scale=1.0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.212547717)   + '   ' + str(scale*0.000089109)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.106196739)  + '   ' + str(scale*5.380172221)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')   
    f.write('Ca      -0.018740752  -0.016834681   7.625755402\n')
    f.write('Ca      4.214708073  -0.010157419  10.093283019\n')
    f.write('Ca     -2.165386108   3.619473106  10.085191840\n')
    f.write('Ca     -0.021746507   3.613378058   7.805509353\n')
    f.write('Ca      1.055562488   1.767863543  10.093119303 \n')
    f.write('Ca      3.150659105   1.780871776   7.608524091\n')
    f.write('N      -1.074814037   1.789246180   8.953015217\n')
    f.write('N       2.091890476  -0.036123957   8.947277008\n')
    f.write('N       2.083443058   3.623032512   8.972809607\n')

    f.write('Mo      -0.033526902   1.784076213   3.751846242\n')
    f.write('S       1.588583378   0.876749335   2.226094733\n')
    f.write('S       1.557858646   0.886447290   5.334867993\n')
    f.write('Mo      3.144908829   1.933132850   3.778510767\n')
    f.write('S       4.688752387   0.860048190   2.231655079\n')
    f.write('S       4.712713726   0.854059453   5.331214954\n')
    f.write('Mo     -1.741185656   4.484108548   3.746639280\n')
    f.write('S       -0.027797849   3.648535798   2.213883791\n')
    f.write('S       -0.029287402   3.656173327   5.195200604\n')
    f.write('Mo       1.696198522   4.465021960   3.750219748\n')
    f.write('S       3.149789005   3.602955612   2.056556556\n')
    f.write('S       3.174004961   3.633669324   5.499823623\n')    
   
    f.close()





def make_bi_Ca2N_MoS2(struc_filename, scale=1.0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.154456204)   + '   ' + str(scale*0.000033056)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.077199478)  + '   ' + str(scale*5.329915163)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    
    f.write('Ca      0.001195365   0.001875137   7.543973046\n')
    f.write('Ca      4.112653787   0.021928441  10.032585089\n')
    f.write('Ca     -2.073636831   3.551672206  10.033032923\n')
    f.write('Ca      0.000274988   3.553324129   7.678389606\n')
    f.write('Ca      1.040884152   1.759098999  10.032886517\n')
    f.write('Ca      3.079129799   1.778713226   7.521743744\n')
    f.write('N      -1.030277849   1.768384797   8.786689272\n')
    f.write('N       2.047299669   0.010321352   8.784709695\n')
    f.write('N       2.063030212   3.554229099   8.784760849\n')
 
    f.write('Ca      0.000023839   0.000033923  15.805988793\n')
    f.write('Ca     -1.026285579   1.775604122  13.333343546\n')
    f.write('Ca      2.050930520   0.001119178  13.333498914\n')
    f.write('Ca      0.000034682   3.553286760  15.805175060\n')
    f.write('Ca      2.052766195   3.553364586  13.333481796\n')
    f.write('Ca      3.077278527   1.776688526  15.805921604\n')
    f.write('N       1.025727784   1.776681068  14.658795487\n')
    f.write('N      -2.051411473   3.553275943  14.658794381\n')
    f.write('N       4.103009107   0.000048826  14.658798839\n')

    f.write('Mo      0.005946779   1.761930581   3.491860323\n')
    f.write('S       1.559002707   0.913472245   1.884218048\n')
    f.write('S       1.537167900   0.907871622   5.155870686\n')
    f.write('Mo      3.086401480   1.790890593   3.516587274\n')
    f.write('S       4.611298214   0.907088905   1.884971700\n')
    f.write('S       4.626134193   0.891952895   5.154869078\n')
    f.write('Mo     -1.552046117   4.475030458   3.491937922\n')
    f.write('S       0.010646015   3.568399838   1.902005080\n')
    f.write('S       0.011089393   3.569033901   5.026202118\n')
    f.write('Mo      1.577978045   4.467138650   3.491788253\n')
    f.write('S       3.091904507   3.552878659   1.884359600\n')
    f.write('S       3.098082770   3.572402202   5.155759694\n')    
   
    f.close()




def make_tri_Ca2N_MoS2(struc_filename, scale=1.0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.122397896)   + '   ' + str(scale*-0.000340373)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.061493719)  + '   ' + str(scale*5.301969885)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    41.000000000\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    
    f.write('Ca      0.001391969   0.002072960   7.541578215\n')
    f.write('Ca      4.090672819   0.020720495  10.036200240\n')
    f.write('Ca      -2.061703228   3.533217854  10.036634802\n')
    f.write('Ca      0.000220424   3.534813600   7.673355897\n')
    f.write('Ca      1.034861210   1.751048505  10.036676155\n')
    f.write('Ca      3.063081047   1.769360182   7.519578557\n')
    f.write('N       -1.024723839   1.759388614   8.791205050\n')
    f.write('N       2.036733262   0.010236590   8.789180341\n')
    f.write('N       2.052054429   3.535631821   8.789227972\n')
 
    f.write('Ca      0.001051581   0.001790911  15.824141829\n')
    f.write('Ca      -1.020609204   1.767354934  13.362238012\n')
    f.write('Ca      2.040947226   0.001634943  13.362356249\n')
    f.write('Ca      0.000542273   3.536668164  15.822894938\n')
    f.write('Ca      2.042424243   3.535707458  13.362507952\n')
    f.write('Ca      3.062640615   1.768799949  15.823372524\n')
    f.write('N       1.020864462   1.768466607  14.590426556\n')
    f.write('N       -2.040228988   3.536102869  14.591101117\n')
    f.write('N       4.082560058   0.000894707  14.591106521\n')

    f.write('Ca      -0.999315645   1.801907986  21.596108307\n')
    f.write('Ca      -2.020483735   3.568068000  19.119006044\n')
    f.write('Ca      2.060160153   0.035520383  21.596112120\n')
    f.write('Ca      4.100289415   0.034287774  19.119010036\n')
    f.write('Ca      1.039682528   1.800736978  19.114391783\n')
    f.write('Ca      2.060566901   3.568997620  21.593624164\n')
    f.write('N       0.020460639   0.035433591  20.439381191\n')
    f.write('N       0.019241964   3.569349104  20.439959655\n')
    f.write('N       3.081534045   1.801337911  20.440139316\n')

    f.write('Mo      0.006005496   1.753189702   3.486699424\n')
    f.write('S       1.551141739   0.908164236   1.873035041\n')
    f.write('S       1.529668318   0.902613511   5.156969203\n')
    f.write('Mo      3.070690039   1.780826767   3.512113513\n')
    f.write('S       4.587400590   0.902297900   1.873133901\n')
    f.write('S       4.602014729   0.887374774   5.156758488\n')
    f.write('Mo      -1.543975200   4.450329477   3.486744847\n')
    f.write('S       0.010179409   3.549382108   1.890456500\n')
    f.write('S       0.010587835   3.550019903   5.026930484\n')
    f.write('Mo      1.568203585   4.443686079   3.486616783\n')
    f.write('S       3.074733494   3.533560454   1.873098564\n')
    f.write('S       3.080854346   3.552764341   5.156928821\n')    
   
    f.close()




def make_quad_Ca2N_MoS2(struc_filename, scale=1.0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*6.102673502)   + '   ' + str(scale*-0.000257905)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-3.051560101)  + '   ' + str(scale*5.284935537)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    45.000000000\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca      0.001509829   0.002254152   7.539948521\n')
    f.write('Ca      4.078192420   0.021400811  10.038182224\n')
    f.write('Ca     -2.055571544   3.522286653  10.038620090\n')
    f.write('Ca      0.000402196   3.523723320   7.671452032\n')
    f.write('Ca      1.031858539   1.744990777  10.038463628\n')
    f.write('Ca      3.053387812   1.763906760   7.518185742\n')
    f.write('N      -1.021549088   1.753685936   8.789733187\n')
    f.write('N       2.030216515   0.010800566   8.787750348\n')
    f.write('N       2.046045047   3.524447409   8.787765778\n')
 
    f.write('Ca      0.001488133   0.002651769  15.823049029\n')
    f.write('Ca     -1.017030240   1.761809172  13.355765310\n')
    f.write('Ca      2.034409575   0.003760404  13.355926339\n')
    f.write('Ca      0.001405662   3.525834223  15.821696981\n')
    f.write('Ca      2.037394197   3.525223267  13.355924929\n')
    f.write('Ca      3.052882854   1.764104861  15.823116665\n')
    f.write('N       1.018200529   1.763822211  14.594066391\n')
    f.write('N      -2.032949017   3.525429679  14.59403120\n')
    f.write('N       4.069698859   0.001886852  14.594011292\n')
   
    f.write('Ca     -0.997327119   1.795899922  21.616351381\n')
    f.write('Ca     -2.015790176   3.556757276  19.152523986\n')
    f.write('Ca      2.053969862   0.034045479  21.616331696\n')
    f.write('Ca      4.088178373   0.033505745  19.152550005\n')
    f.write('Ca      1.036502355   1.794486323  19.152433867\n')
    f.write('Ca      2.053865383   3.557619468  21.616354624\n')
    f.write('N       0.019585362   0.033892129  20.379830830\n')
    f.write('N       0.019436418   3.557203561  20.377532353\n')
    f.write('N       3.070911708   1.795473380  20.379266856\n')
 
    f.write('Ca      0.026273590   0.045525759  24.909502164\n')
    f.write('Ca      4.096763449   0.044676812  27.391068643\n')
    f.write('Ca     -2.009723847   3.570226732  27.391066828\n')
    f.write('Ca      0.025457503   3.569333149  24.908718947\n')
    f.write('Ca      1.043886411   1.808034229  27.390234960\n')
    f.write('Ca      3.078401241   1.806702919  24.908620047\n')
    f.write('N      -0.990701104   1.807506268  26.232433818\n')
    f.write('N       2.060692425   0.045800052  26.232432350\n')
    f.write('N       2.060619051   3.569085524  26.230289827\n')

    f.write('Mo      0.006041213   1.749004624   3.480969487\n')
    f.write('S       1.545701840   0.904694882   1.864506739\n')
    f.write('S       1.524760806   0.899246041   5.153636530\n')
    f.write('Mo      3.060642969   1.775061618   3.506572554\n')
    f.write('S       4.572994494   0.899230512   1.864667184\n')
    f.write('S       4.587277223   0.884711548   5.153289810\n')
    f.write('Mo     -1.537775361   4.434937080   3.480999108\n')
    f.write('S       0.010046190   3.537790947   1.880687488\n')
    f.write('S       0.010445330   3.538411494   5.025478862\n')
    f.write('Mo      1.561670003   4.428544380   3.480834254\n')
    f.write('S       3.064635482   3.522568189   1.864512290\n')
    f.write('S       3.070563383   3.541333079   5.153589333\n')    
   
    f.close()



def make_interface_mono(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  6.212547717   0.000089109   0.000000000\n')
    f.write(' -3.106196739   5.380172221   0.000000000\n')
    f.write('  0.000000000   0.000000000  30.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    
    f.write('Ca      ' + str(0.000662319) + '   '  + str(0.001409588) + '     7.669178666\n')
    f.write('Ca      ' + str(4.149574561) + '   '  + str(0.023595636) + '    10.151824472\n')
    f.write('Ca      ' + str(-2.095865967) + '   ' + str(3.582458548) + '    10.152713382\n')
    f.write('Ca      ' + str(-0.000948011) + '   ' + str(3.585326583) + '     7.819199532\n')
    f.write('Ca      ' + str(1.052553784) + '   '  + str(1.774239387) + '    10.152518190\n')
    f.write('Ca      ' + str(3.107913050) + '   '  + str(1.795444767) + '     7.641559945\n')
    f.write('N       ' + str(-1.034174075) + '   ' + str(1.790720842) + '     8.996233733\n')
    f.write('N       ' + str(2.067949179) + '   '  + str(0.000365830) + '     8.993747279\n')
    f.write('N       ' + str(2.072522604) + '   '  + str(3.589244948) + '     8.994003348\n')
    f.write('Mo      ' + str(0.007660054  + xd) + '   ' + str(1.764397297 + yd) + '     3.707531671\n')
    f.write('S       ' + str(1.580892939  + xd) + '   ' + str(0.925789980 + yd) + '     2.115566773\n')
    f.write('S       ' + str(1.551543247  + xd) + '   ' + str(0.920100419 + yd) + '     5.388381505\n')
    f.write('Mo      ' + str(3.115040170  + xd) + '   ' + str(1.809404121 + yd) + '     3.745889857\n')
    f.write('S       ' + str(4.649483734  + xd) + '   ' + str(0.920261486 + yd) + '     2.117438882\n')
    f.write('S       ' + str(4.667745235  + xd) + '   ' + str(0.899851134 + yd) + '     5.386420938\n')
    f.write('Mo      ' + str(-1.580840488 + xd) + '   ' + str(4.525349524 + yd) + '     3.707411197\n')
    f.write('S       ' + str(0.011408623  + xd) + '   ' + str(3.603316579 + yd) + '     2.142180273\n')
    f.write('S       ' + str(0.011985395  + xd) + '   ' + str(3.604260804 + yd) + '     5.207606058\n')
    f.write('Mo      ' + str(1.608493324  + xd) + '   ' + str(4.516432175 + yd) + '     3.707304303\n')
    f.write('S       ' + str(3.121654354  + xd) + '   ' + str(3.579588878 + yd) + '     2.114913844\n')
    f.write('S       ' + str(3.131333413  + xd) + '   ' + str(3.604168493 + yd) + '     5.389374365\n')

    f.close()



def make_interface_mono_var_z_old(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  6.212547717   0.000089109   0.000000000\n')
    f.write(' -3.106196739   5.380172221   0.000000000\n')
    f.write('  0.000000000   0.000000000  30.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca      ' + str(0.000662319) + '   '  + str(0.001409588) + '     7.669178666  0   0   1\n')
    f.write('Ca      ' + str(4.149574561) + '   '  + str(0.023595636) + '    10.151824472  0   0   1\n')
    f.write('Ca      ' + str(-2.095865967) + '   ' + str(3.582458548) + '    10.152713382  0   0   1\n')
    f.write('Ca      ' + str(-0.000948011) + '   ' + str(3.585326583) + '     7.819199532  0   0   1\n')
    f.write('Ca      ' + str(1.052553784) + '   '  + str(1.774239387) + '    10.152518190  0   0   1\n')
    f.write('Ca      ' + str(3.107913050) + '   '  + str(1.795444767) + '     7.641559945  0   0   1\n')
    f.write('N       ' + str(-1.034174075) + '   ' + str(1.790720842) + '     8.996233733  0   0   1\n')
    f.write('N       ' + str(2.067949179) + '   '  + str(0.000365830) + '     8.993747279  0   0   1\n')
    f.write('N       ' + str(2.072522604) + '   '  + str(3.589244948) + '     8.994003348  0   0   1\n')
    f.write('Mo      ' + str(0.007660054  + xd) + '   ' + str(1.764397297 + yd) + '     3.707531671  0   0   1\n')
    f.write('S       ' + str(1.580892939  + xd) + '   ' + str(0.925789980 + yd) + '     2.115566773  0   0   1\n')
    f.write('S       ' + str(1.551543247  + xd) + '   ' + str(0.920100419 + yd) + '     5.388381505  0   0   1\n')
    f.write('Mo      ' + str(3.115040170  + xd) + '   ' + str(1.809404121 + yd) + '     3.745889857  0   0   1\n')
    f.write('S       ' + str(4.649483734  + xd) + '   ' + str(0.920261486 + yd) + '     2.117438882  0   0   1\n')
    f.write('S       ' + str(4.667745235  + xd) + '   ' + str(0.899851134 + yd) + '     5.386420938  0   0   1\n')
    f.write('Mo      ' + str(-1.580840488 + xd) + '   ' + str(4.525349524 + yd) + '     3.707411197  0   0   1\n')
    f.write('S       ' + str(0.011408623  + xd) + '   ' + str(3.603316579 + yd) + '     2.142180273  0   0   1\n')
    f.write('S       ' + str(0.011985395  + xd) + '   ' + str(3.604260804 + yd) + '     5.207606058  0   0   1\n')
    f.write('Mo      ' + str(1.608493324  + xd) + '   ' + str(4.516432175 + yd) + '     3.707304303  0   0   1\n')
    f.write('S       ' + str(3.121654354  + xd) + '   ' + str(3.579588878 + yd) + '     2.114913844  0   0   1\n')
    f.write('S       ' + str(3.131333413  + xd) + '   ' + str(3.604168493 + yd) + '     5.389374365  0   0   1\n')

    f.close()



def make_interface_mono_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  6.354820335   0.000016855     0.000000000\n')
    f.write(' -3.177395624   5.503426247     0.000000000\n')
    f.write('  0.000000000   0.000000000  30.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Ca      ' + str(-0.600793983) + '   '  + str(-0.104373504) + '     7.669537504  0   0   1\n')
    f.write('Ca      ' + str(3.606777368) + '   '  + str(-0.047541328) + '    10.060486917  0   0   1\n')
    f.write('Ca      ' + str(-2.748390914) + '   ' + str(3.622367208) + '    10.109649163  0   0   1\n')
    f.write('Ca      ' + str(-0.695300055) + '   ' + str(3.623767791) + '     7.671373196  0   0   1\n')
    f.write('Ca      ' + str(0.429409415) + '   '  + str(1.786891756) + '    10.058644407  0   0   1\n')
    f.write('Ca      ' + str(2.583402540) + '   '  + str(1.842518960) + '     7.671569092  0   0   1\n')
    f.write('N       ' + str(-1.688554740) + '   ' + str(1.787454229) + '     8.950560834  0   0   1\n')
    f.write('N       ' + str(1.488011763) + '   '  + str(-0.046870451) + '     8.948861835  0   0   1\n')
    f.write('N       ' + str(1.488267502) + '   '  + str(3.621125102) + '     8.949298965  0   0   1\n')
    f.write('Mo      ' + str(0.424533022  + xd) + '   ' + str(1.797975037 + yd) + '     3.746738378  0   0   1\n')
    f.write('S       ' + str(2.032502248  + xd) + '   ' + str(0.832031981 + yd) + '     2.235698858  0   0   1\n')
    f.write('S       ' + str(1.992945365  + xd) + '   ' + str(0.809605352 + yd) + '     5.270608443  0   0   1\n')
    f.write('Mo      ' + str(3.555311573  + xd) + '   ' + str(1.952921329 + yd) + '     3.772201328  0   0   1\n')
    f.write('S       ' + str(5.139279068  + xd) + '   ' + str(0.887687956 + yd) + '     2.236031738  0   0   1\n')
    f.write('S       ' + str(5.139991272  + xd) + '   ' + str(0.932707273 + yd) + '     5.272643784  0   0   1\n')
    f.write('Mo      ' + str(-1.275049258 + xd) + '   ' + str(4.430543190 + yd) + '     3.771931862  0   0   1\n')
    f.write('S       ' + str(0.457607627  + xd) + '   ' + str(3.673329941 + yd) + '     2.235907804  0   0   1\n')
    f.write('S       ' + str(0.496407365  + xd) + '   ' + str(3.651254775 + yd) + '     5.271009533  0   0   1\n')
    f.write('Mo      ' + str(2.172621967  + xd) + '   ' + str(4.512275205 + yd) + '     3.772290014  0   0   1\n')
    f.write('S       ' + str(3.602357996  + xd) + '   ' + str(3.631799015 + yd) + '     2.065051651  0   0   1\n')
    f.write('S       ' + str(3.603292531  + xd) + '   ' + str(3.630220152 + yd) + '     5.560902905  0   0   1\n')

    f.close()



def make_au_mos2_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  5.501797669   0.00000000   0.000000000\n')
    f.write(' -2.750898834    4.764696549   0.000000000\n')
    f.write('  0.000000000   0.000000000  35.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Au      ' + str(0.103429615) + '   '  + str(0.180567724) + '  6.974424176   0   0   1\n')
    f.write('Au      ' + str(1.460581256) + '   '  + str(0.939492512) + '  9.552173402  0   0   1\n')
    f.write('Au      ' + str(0.063298774) + '   '  + str(1.693489462) + '  12.159996736  0   0   1\n')
    f.write('Au      ' + str(0.046060407) + '   '  + str(0.081415799) + '    14.779085962  0   0   1\n')
    
    f.write('Au      ' + str(2.854995107) + '   '  + str(0.186193126) + '     6.981024464  0   0   1\n')
    f.write('Au      ' + str(4.211944759) + '   '  + str(0.936670863) + '     9.596714887  0   0   1\n')
    f.write('Au      ' + str(2.812070697) + '   '  + str(1.693572054) + '    12.156236638  0   0   1\n')
    f.write('Au      ' + str(2.796506868) + '   '  + str(0.081216106) + '    14.778024244  0   0   1\n')

    f.write('Au      ' + str(-1.265442627) + '   '  + str(2.563719348) + '     6.967749688  0   0   1\n')
    f.write('Au      ' + str(0.088930708)  + '   '  + str(3.316476895) + '     9.541845134  0   0   1\n')
    f.write('Au      ' + str(-1.313749034) + '   '  + str(4.074557010) + '    12.133707817  0   0   1\n')
    f.write('Au      ' + str(-1.329739034) + '   '  + str(2.464214067) + '    14.780257999  0   0   1\n')
    
    f.write('Au      ' + str(1.481790054) + '   '  + str(2.565391402) + '     6.798986632  0   0   1\n')
    f.write('Au      ' + str(2.832197220) + '   '  + str(3.315917231) + '     9.530992672  0   0   1\n')
    f.write('Au      ' + str(1.437586724) + '   '  + str(4.073753810) + '    12.150432452  0   0   1\n')
    f.write('Au      ' + str(1.421269348) + '   '  + str(2.463825656) + '    14.791857708  0   0   1\n')


    f.write('Mo      ' + str(0.582858358  + xd) + '   ' + str(7.363455835 + yd) + '     2.501816532  0   0   1\n')
    f.write('S       ' + str(1.497764805  + xd) + '   ' + str(5.775020235 + yd) + '     0.938149936  0   0   1\n')
    f.write('S       ' + str(1.498321535  + xd) + '   ' + str(5.775475567 + yd) + '     4.059599421  0   0   1\n')

    f.write('Mo      ' + str(3.330860480  + xd) + '   '  + str(5.773395247 + yd) + '     2.501769734  0   0   1\n')
    f.write('S       ' + str(4.248966523  + xd) + '   '  + str(4.186893156 + yd) + '     0.938221819  0   0   1\n')
    f.write('S       ' + str(4.248617773  + xd) + '   '  + str(4.186986368 + yd) + '     4.060357295  0   0   1\n')
    
    f.write('Mo      ' + str(6.081757663  + xd) + '   ' + str(4.188508287 + yd)  + '     2.502113578  0   0   1\n')
    f.write('S       ' + str(4.248762461  + xd) + '   ' + str(7.363557560 + yd)  + '     0.940667778  0   0   1\n')
    f.write('S       ' + str(4.248476578  + xd) + '   ' + str(7.362628429 + yd)  + '     4.057430256  0   0   1\n')
    
    f.close()



def make_au_mos2_var_z_old(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  5.584610679   0.000078444   0.000000000\n')
    f.write(' -2.792286805   4.837072617   0.000000000\n')
    f.write('  0.000000000   0.000000000  35.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Au      ' + str(0.047731266) + '   '  + str(0.069552751) + '     7.011231851  0   0   1\n')
    f.write('Au      ' + str(1.420500919) + '   '  + str(0.871127421) + '     9.542297987  0   0   1\n')
    f.write('Au      ' + str(0.037914242) + '   '  + str(1.660387123) + '    12.068940762  0   0   1\n')
    f.write('Au      ' + str(0.032867198) + '   '  + str(0.067525773) + '    14.637223111  0   0   1\n')
    
    f.write('Au      ' + str(2.833155207) + '   '  + str(0.070720916) + '     6.966639275  0   0   1\n')
    f.write('Au      ' + str(4.209848189) + '   '  + str(0.867600749) + '     9.543267226  0   0   1\n')
    f.write('Au      ' + str(2.830820074) + '   '  + str(1.660397076) + '    12.077491850  0   0   1\n')
    f.write('Au      ' + str(2.825113155) + '   '  + str(0.067878547) + '    14.642628395  0   0   1\n')

    f.write('Au      ' + str(-1.344276267) + '   '  + str(2.482155111) + '     6.963362607  0   0   1\n')
    f.write('Au      ' + str(0.024119304)  + '   '  + str(3.286404519) + '     9.521606556  0   0   1\n')
    f.write('Au      ' + str(-1.356802403) + '   '  + str(4.078544363) + '    12.065783297  0   0   1\n')
    f.write('Au      ' + str(-1.362865104) + '   '  + str(2.486664572) + '    14.641983066  0   0   1\n')
    
    f.write('Au      ' + str(1.448613852) + '   '  + str(2.494861941) + '     6.950931895  0   0   1\n')
    f.write('Au      ' + str(2.817699798) + '   '  + str(3.285254482) + '     9.538846239  0   0   1\n')
    f.write('Au      ' + str(1.434626523) + '   '  + str(4.078417300) + '    12.067450718  0   0   1\n')
    f.write('Au      ' + str(1.429310731) + '   '  + str(2.487484922) + '    14.643565445  0   0   1\n')


    f.write('Mo      ' + str(-1.861947027  + xd) + '   ' + str(3.226586574 + yd) + '     2.589148406  0   0   1\n')
    f.write('S       ' + str(-0.931573977  + xd) + '   ' + str(1.614283657 + yd) + '     1.042163755  0   0   1\n')
    f.write('S       ' + str(-0.931492800  + xd) + '   ' + str(1.617088554 + yd) + '     4.132263063  0   0   1\n')

    f.write('Mo      ' + str(0.930428027  + xd) + '   '  + str(1.611871148 + yd) + '     2.588468574  0   0   1\n')
    f.write('S       ' + str(1.863368850  + xd) + '   '  + str(0.000370285 + yd) + '     1.042078201  0   0   1\n')
    f.write('S       ' + str(1.865485867  + xd) + '   '  + str(0.001729644 + yd) + '     4.132241928  0   0   1\n')
    
    f.write('Mo      ' + str(3.724573288  + xd) + '   ' + str(0.000865878 + yd)  + '     2.589231313  0   0   1\n')
    f.write('S       ' + str(1.861158508  + xd) + '   ' + str(3.224658164 + yd)  + '     1.042364152  0   0   1\n')
    f.write('S       ' + str(1.858876042  + xd) + '   ' + str(3.220321789 + yd)  + '     4.132427285  0   0   1\n')
    
    f.close()


#this struct is scan11 of first scan along short diagonal, vc-relaxed
def make_au_mos2_var_z_new_short(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  5.583627657   0.000042450   0.000000000\n')
    f.write(' -2.791826456   4.835634908   0.000000000\n')
    f.write('  0.000000000   0.000000000  35.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Au      ' + str(0.049011116) + '   '  + str(0.071460529) + '     6.969973562  0   0   1\n')
    f.write('Au      ' + str(1.423352439) + '   '  + str(0.867157034) + '     9.517131289  0   0   1\n')
    f.write('Au      ' + str(0.037757780) + '   '  + str(1.660580787) + '    12.065853122  0   0   1\n')
    f.write('Au      ' + str(0.033363472) + '   '  + str(0.067466223) + '    14.649486900  0   0   1\n')
    
    f.write('Au      ' + str(2.837729521) + '   '  + str(0.067598517) + '     6.902150564  0   0   1\n')
    f.write('Au      ' + str(4.209548662) + '   '  + str(0.864323497) + '     9.527483007  0   0   1\n')
    f.write('Au      ' + str(2.829118108) + '   '  + str(1.660598996) + '    12.074003172  0   0   1\n')
    f.write('Au      ' + str(2.824866803) + '   '  + str(0.067135841) + '    14.653605098  0   0   1\n')

    f.write('Au      ' + str(-1.350794990) + '   '  + str(2.484123888) + '     7.005368761  0   0   1\n')
    f.write('Au      ' + str(0.027257343)  + '   '  + str(3.287256597) + '     9.528780180  0   0   1\n')
    f.write('Au      ' + str(-1.357482821) + '   '  + str(4.079434520) + '    12.072473430  0   0   1\n')
    f.write('Au      ' + str(-1.362505710) + '   '  + str(2.485082836) + '    14.651685388  0   0   1\n')
    
    f.write('Au      ' + str(1.444414019) + '   '  + str(2.491085798) + '     6.968495450  0   0   1\n')
    f.write('Au      ' + str(2.817234744) + '   '  + str(3.287236707) + '     9.547966466  0   0   1\n')
    f.write('Au      ' + str(1.432942969) + '   '  + str(4.079144665) + '    12.075461695  0   0   1\n')
    f.write('Au      ' + str(1.429068232) + '   '  + str(2.484981104) + '    14.650667819  0   0   1\n')


    f.write('Mo      ' + str(-9.538463362  + xd) + '   ' + str(7.656783201 + yd) + '     2.594812308  0   0   1\n')
    f.write('S       ' + str(-8.609320583  + xd) + '   ' + str(6.046474740 + yd) + '     1.044096769  0   0   1\n')
    f.write('S       ' + str(-8.603399538  + xd) + '   ' + str(6.035550467 + yd) + '     4.135421100  0   0   1\n')

    f.write('Mo      ' + str(-6.747351593  + xd) + '   '  + str(6.046002788 + yd) + '     2.589678606  0   0   1\n')
    f.write('S       ' + str(-5.816530886  + xd) + '   '  + str(4.432166694 + yd) + '     1.045050045  0   0   1\n')
    f.write('S       ' + str(-5.810545725  + xd) + '   '  + str(4.444499131 + yd) + '     4.134989528  0   0   1\n')
    
    f.write('Mo      ' + str(-3.954251367  + xd) + '   ' + str(4.433132777 + yd)  + '     2.590076857  0   0   1\n')
    f.write('S       ' + str(-5.814340040  + xd) + '   ' + str(7.656932790 + yd)  + '     1.044038328  0   0   1\n')
    f.write('S       ' + str(-5.825733501  + xd) + '   ' + str(7.656613406 + yd)  + '     4.134887514  0   0   1\n')
    
    f.close()




#this struct is scan21 of first scan along long diagonal, vc-relaxed
def make_au_mos2_var_z_new_long(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  5.583670411   0.000034064   0.000000000\n')
    f.write(' -2.791855099   4.836009119   0.000000000\n')
    f.write('  0.000000000   0.000000000  35.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Au      ' + str(0.046424570) + '   '  + str(0.068910279) + '     6.988818205  0   0   1\n')
    f.write('Au      ' + str(1.419347385) + '   '  + str(0.873568590) + '     9.523457115  0   0   1\n')
    f.write('Au      ' + str(0.038943978) + '   '  + str(1.660609445) + '    12.075381096  0   0   1\n')
    f.write('Au      ' + str(0.033211307) + '   '  + str(0.067619813) + '    14.650320575  0   0   1\n')
    
    f.write('Au      ' + str(2.834244450) + '   '  + str(0.071991793) + '     6.991911435  0   0   1\n')
    f.write('Au      ' + str(4.211529501) + '   '  + str(0.866719069) + '     9.555123098  0   0   1\n')
    f.write('Au      ' + str(2.828982797) + '   '  + str(1.660899286) + '    12.076087215  0   0   1\n')
    f.write('Au      ' + str(2.824527022) + '   '  + str(0.067839294) + '    14.652108997  0   0   1\n')

    f.write('Au      ' + str(-1.343565153) + '   '  + str(2.482741994) + '     6.989767750  0   0   1\n')
    f.write('Au      ' + str(0.028047918)  + '   '  + str(3.283677917) + '     9.520836231  0   0   1\n')
    f.write('Au      ' + str(-1.357287785) + '   '  + str(4.078258686) + '    12.063675032  0   0   1\n')
    f.write('Au      ' + str(-1.362706883) + '   '  + str(2.485888412) + '    14.651918247  0   0   1\n')
    
    f.write('Au      ' + str(1.446683343) + '   '  + str(2.492349091) + '     6.874420839  0   0   1\n')
    f.write('Au      ' + str(2.813529472) + '   '  + str(3.284085795) + '     9.521235596  0   0   1\n')
    f.write('Au      ' + str(1.434303401) + '   '  + str(4.076893350) + '    12.074012190  0   0   1\n')
    f.write('Au      ' + str(1.429065746) + '   '  + str(2.486138435) + '    14.655306824  0   0   1\n')


    f.write('Mo      ' + str(0.582003636  + xd) + '   ' + str(7.457205650 + yd) + '     2.591296379  0   0   1\n')
    f.write('S       ' + str(1.511548948  + xd) + '   ' + str(5.845143917 + yd) + '     1.043767055  0   0   1\n')
    f.write('S       ' + str(1.513082125  + xd) + '   ' + str(5.846299416 + yd) + '     4.133538506  0   0   1\n')

    f.write('Mo      ' + str(3.372980160  + xd) + '   '  + str(5.842841573 + yd) + '     2.590738942  0   0   1\n')
    f.write('S       ' + str(4.305546409  + xd) + '   '  + str(4.231930256 + yd) + '     1.043880345  0   0   1\n')
    f.write('S       ' + str(4.305458020  + xd) + '   '  + str(4.233620142 + yd) + '     4.134305526  0   0   1\n')
    
    f.write('Mo      ' + str(6.166097468  + xd) + '   ' + str(4.232812528 + yd)  + '     2.591434394  0   0   1\n')
    f.write('S       ' + str(4.303834697  + xd) + '   ' + str(7.455735427 + yd)  + '     1.045357421  0   0   1\n')
    f.write('S       ' + str(4.302284488  + xd) + '   ' + str(7.452613591 + yd)  + '     4.134937947  0   0   1\n')
    
    f.close()

def make_Au_MoS2(struc_filename, scale=1.0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*5.583670411)   + '   ' + str(scale*0.000034064)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.791855099)  + '   ' + str(scale*4.836009119)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write('Au     0.046424570   0.068910279   6.988818205\n')
    f.write('Au      1.419347385   0.873568590   9.523457115\n')
    f.write('Au      0.038943978   1.660609445  12.075381096\n')
    f.write('Au      0.033211307   0.067619813  14.650320575\n')
   
    f.write('Au      2.834244450   0.071991793   6.991911435\n')
    f.write('Au      4.211529501   0.866719069   9.555123098\n')
    f.write('Au      2.828982797   1.660899286  12.076087215\n')
    f.write('Au      2.824527022   0.067839294  14.652108997\n')
   
    f.write('Au      -1.343565153   2.482741994   6.989767750\n')
    f.write('Au      0.028047918   3.283677917   9.520836231\n')
    f.write('Au      -1.357287785   4.078258686  12.063675032\n')
    f.write('Au      -1.362706883   2.485888412  14.651918247\n')
   
    f.write('Au      1.446683343   2.492349091   6.874420839\n')
    f.write('Au      2.813529472   3.284085795   9.521235596\n')
    f.write('Au      1.434303401   4.076893350  12.074012190\n')
    f.write('Au      1.429065746   2.486138435  14.655306824\n')
   

    f.write('Mo      0.582003636   7.457205650   2.591296379\n')
    f.write('S       1.511548948   5.845143917   1.043767055\n')
    f.write('S       1.513082125   5.846299416   4.133538506\n')

    f.write('Mo      3.372980160   5.842841573   2.590738942\n')
    f.write('S       4.305546409   4.231930256   1.043880345\n')
    f.write('S       4.305458020   4.233620142   4.134305526\n')

    f.write('Mo      6.166097468   4.232812528   2.591434394\n')
    f.write('S       4.303834697   7.455735427   1.045357421\n')
    f.write('S       4.302284488   7.452613591   4.134937947\n')
 
    f.close()

#atomic positions at +3% stretch
def make_Au(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*5.644955091)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.822477545)  + '   ' + str(scale*4.888674511)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Au      -0.031313305  -0.044414948   6.959612717\n')
        f.write('Au      1.416302879   0.778986280   9.359990950\n')
        f.write('Au      -0.033534931   1.616042242  11.703810630\n')
        f.write('Au      -0.022802495  -0.049334617  14.104060274\n')
   
        f.write('Au      2.875260823  -0.046225306   6.959916339\n')
        f.write('Au      4.324789241   0.777821984   9.357985966\n')
        f.write('Au      2.877560840   1.614065842  11.706171045\n')
        f.write('Au      2.886781609  -0.053502649  14.103405513\n')
   
        f.write('Au      -1.489714877   2.473272076   6.960263085\n')
        f.write('Au      -0.040962361   3.299081562   9.357636183\n')
        f.write('Au      -1.488793436   4.134291097  11.705821573\n')
        f.write('Au      -1.477660751   2.466924259  14.103744491\n')
   
        f.write('Au      1.416674397   2.471546541   6.960146847\n')
        f.write('Au      2.867422392   3.298393709   9.354905633\n')
        f.write('Au      1.422771927   4.132458925  11.708898548\n')
        f.write('Au      1.432083320   2.462640078  14.103519852\n')
   
    f.close()

#atomic positions at +3% stretch
def make_Au_prm(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*5.644955091/2)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.822477545/2)  + '   ' + str(scale*4.888674511/2)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    27.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Au      -0.031313305  -0.044414948   6.959612717\n')
        f.write('Au      1.416302879   0.778986280   9.359990950\n')
        f.write('Au      -0.033534931   1.616042242  11.703810630\n')
        f.write('Au      -0.022802495  -0.049334617  14.104060274\n')
  
    f.close()


#atomic positions at +1% stretch
def make_MoS2_super_5(struc_filename, scale=1.0, include_atoms=True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*5.457225415)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.728612707)  + '   ' + str(scale*4.726095845)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')

    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
    
        f.write('Mo      -1.845213025   3.166172529   3.073749946\n')
        f.write('S       -0.927908818   1.575144195   1.514194651\n')
        f.write('S       -0.927908966   1.575144078   4.633304334\n')

        f.write('Mo      0.909229494   1.574831615   3.073750091\n')
        f.write('S       1.828069295  -0.016020479   1.514194651\n')
        f.write('S       1.828069268  -0.016020666   4.633304334\n')

        f.write('Mo      3.664592316  -0.014915065   3.073749946\n')
        f.write('S       1.828454038   3.166975547   1.515889228\n')
        f.write('S       1.828454047   3.166975564   4.631610129\n')
 
    f.close()


def make_ca2n_mos2_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  11.006871693   0.000029194   0.000000000\n')
    f.write(' -5.503410657   9.532213876   0.000000000\n')
    f.write('  0.000000000   0.000000000  40.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write(' Ca      -1.843385883   2.874324619   7.577162505  0   0   1\n')
    f.write(' Ca       0.012087505  -0.334471068   7.600697942  0   0   1\n')
    f.write(' Ca       1.870605782   2.732879593   7.730626628  0   0   1\n')
    f.write(' Ca       3.666308064  -0.298367384   7.550708180  0   0   1\n')
    f.write(' Ca       7.394080080  -0.459606105   7.716683118  0   0   1\n')
    f.write(' Ca       0.005665063   6.019838445   7.595232002  0   0   1\n')
    f.write(' Ca       3.672882331   6.049921353   7.581787077  0   0   1\n')
    f.write(' Ca       5.519984195   2.834434101   7.551349228  0   0   1\n')
    f.write(' Ca       7.390535666   5.896638324   7.706654414  0   0   1\n')
    f.write(' Ca       0.013724175   3.879743922  10.028021211  0   0   1\n')
    f.write(' Ca       1.864213284   0.642333082  10.035163579  0   0   1\n')
    f.write(' Ca       5.535378364   0.712536111  10.008925513  0   0   1\n')
    f.write(' Ca      -1.810475227   7.080754054  10.008820883  0   0   1\n')
    f.write(' Ca       3.719054161   3.882301055  10.008834031  0   0   1\n')
    f.write(' Ca       1.848807738   7.010734654  10.047578941  0   0   1\n')
    f.write(' Ca       9.190757875   0.723794167  10.025072013  0   0   1\n')
    f.write(' Ca       5.544106781   7.053516826  10.008997167  0   0   1\n')
    f.write(' Ca       7.344884454   3.815030132  10.038536514  0   0   1\n')
    f.write(' N       -1.809369148   4.929726459   9.012698296  0   0   1\n')
    f.write(' N        0.032196163   1.754452160   9.002293798  0   0   1\n')
    f.write(' N        3.697265396   1.755312590   8.988702889  0   0   1\n')
    f.write(' N       -3.645918504   8.095704777   8.905584345  0   0   1\n')
    f.write(' N        0.036152332   8.105776329   8.982064163  0   0   1\n')
    f.write(' N        1.856459747   4.914680002   8.963635961  0   0   1\n')
    f.write(' N        3.695516295   8.105059232   9.034064655  0   0   1\n')
    f.write(' N        5.534227025   4.925700625   8.958085731  0   0   1\n')
    f.write(' N        7.360609670   1.742977950   8.899992651  0   0   1\n')
    f.write(' Mo       '+str(-2.028033252  +xd)   + '   '  +str(  3.195127933  +yd) +'    3.806149302  0   0   1\n')
    f.write(' Mo       '+str(-4.664513504  +xd)   + '   '  +str(  8.293947504  +yd) +'    3.781907120  0   0   1\n')
    f.write(' Mo       '+str(-1.927639006  +xd)   + '   '  +str(  6.457919968  +yd) +'    3.769630392  0   0   1\n')
    f.write(' Mo        '+str(6.432700260  +xd)   + '   '  +str(  4.778043188  +yd) +'    3.767745739  0   0   1\n')
    f.write(' Mo        '+str(0.859122792  +xd)   + '   '  +str(  1.617001873  +yd) +'    3.774534825  0   0   1\n')
    f.write(' Mo        '+str(3.408091305  +xd)   + '   '  +str(  0.040683028  +yd) +'    3.797784048  0   0   1\n')
    f.write(' Mo        '+str(0.837283493  +xd)   + '   '  +str(  5.085078970  +yd) +'    3.778922564  0   0   1\n')
    f.write(' Mo        '+str(3.573668429  +xd)   + '   '  +str(  3.269930036  +yd) +'    3.761574372  0   0   1\n')
    f.write(' Mo        '+str(6.374299121  +xd)   + '   '  +str(  1.910360671  +yd) +'    3.774979310  0   0   1\n')
    f.write(' Mo        '+str(9.072180387  +xd)   + '   '  +str(  0.085868149  +yd) +'    3.772200265  0   0   1\n')
    f.write(' Mo        '+str(0.967321983  +xd)   + '   '  +str(  7.934268765  +yd) +'    3.775840681  0   0   1\n')
    f.write(' Mo        '+str(3.400814799  +xd)   + '   '  +str(  6.450156543  +yd) +'    3.799787348  0   0   1\n')
    f.write(' S        '+str(-1.008557108  +xd)   + '   '  +str(  1.670529374  +yd) +'    2.217156141  0   0   1\n')
    f.write(' S        '+str(-0.989356652  +xd)   + '   '  +str(  1.594599415  +yd) +'    5.381239102  0   0   1\n')
    f.write(' S        '+str(-3.809507744  +xd)   + '   '  +str(  6.534585298  +yd) +'    2.252335311  0   0   1\n')
    f.write(' S        '+str(-3.820778287  +xd)   + '   '  +str(  6.571035042  +yd) +'    5.211898288  0   0   1\n')
    f.write(' S        '+str(-1.027092222  +xd)   + '   '  +str(  4.821363777  +yd) +'    2.260568725  0   0   1\n')
    f.write(' S        '+str(-1.066316631  +xd)   + '   '  +str(  4.825522029  +yd) +'    5.371409181  0   0   1\n')
    f.write(' S         '+str(1.743096222  +xd)   + '   '  +str(  0.164817775  +yd) +'    2.077952433  0   0   1\n')
    f.write(' S         '+str(1.747472379  +xd)   + '   '  +str(  0.242313568  +yd) +'    5.569849697  0   0   1\n')
    f.write(' S         '+str(1.691524711  +xd)   + '   '  +str(  3.346864441  +yd) +'    2.244579615  0   0   1\n')
    f.write(' S         '+str(1.694491870  +xd)   + '   '  +str(  3.380157651  +yd) +'    5.220309934  0   0   1\n')
    f.write(' S         '+str(4.480165022  +xd)   + '   '  +str(  1.630501890  +yd) +'    2.268252211  0   0   1\n')
    f.write(' S         '+str(4.446190336  +xd)   + '   '  +str(  1.633463718  +yd) +'    5.337789716  0   0   1\n')
    f.write(' S         '+str(7.202744634  +xd)   + '   '  +str(  0.161016181  +yd) +'    2.240469138  0   0   1\n')
    f.write(' S         '+str(7.192678234  +xd)   + '   '  +str(  0.209892321  +yd) +'    5.229780689  0   0   1\n')
    f.write(' S        '+str(-0.966336043  +xd)   + '   '  +str(  8.064715719  +yd) +'    2.260199519  0   0   1\n')
    f.write(' S        '+str(-0.943626933  +xd)   + '   '  +str(  7.984910473  +yd) +'    5.344512922  0   0   1\n')
    f.write(' S         '+str(1.757807197  +xd)   + '   '  +str(  6.473785962  +yd) +'    2.057235678  0   0   1\n')
    f.write(' S         '+str(1.770981220  +xd)   + '   '  +str(  6.555577683  +yd) +'    5.599031035  0   0   1\n')
    f.write(' S         '+str(4.516501031  +xd)   + '   '  +str(  4.896849901  +yd) +'    2.260135063  0   0   1\n')
    f.write(' S         '+str(4.526344613  +xd)   + '   '  +str(  4.830955654  +yd) +'    5.325232178  0   0   1\n')
    f.write(' S         '+str(4.478613904  +xd)   + '   '  +str(  8.002041467  +yd) +'    2.250454865  0   0   1\n')
    f.write(' S         '+str(4.446727104  +xd)   + '   '  +str(  8.011388640  +yd) +'    5.388158198  0   0   1\n')
    f.write(' S         '+str(7.302334577  +xd)   + '   '  +str(  3.316736431  +yd) +'    2.085355305  0   0   1\n')
    f.write(' S         '+str(7.295953344  +xd)   + '   '  +str(  3.385033873  +yd) +'    5.548055660  0   0   1\n')

    f.close()








def make_au_ca2n_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Au 196.67 Au.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  11.006871693   0.000029194   0.000000000\n')
    f.write(' -5.503410657   9.532213876   0.000000000\n')
    f.write('  0.000000000   0.000000000  40.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    
    f.write(' Ca      '+str(-1.980394696 +xd)   + '   '  +str(3.254716857 +yd) +'   7.720562931  0   0   1\n')
    f.write(' Ca      '+str(-0.141074525 +xd)   + '   '    +str(0.031095580 +yd) +'  7.475331540  0   0   1\n')
    f.write(' Ca       '+str(1.697749504 +xd)   + '   '    +str(3.256912509 +yd) +'  7.687694493  0   0   1\n')
    f.write(' Ca       '+str(3.612632786 +xd)   + '   '    +str(0.006578461 +yd) +'  7.744874670  0   0   1\n')
    f.write(' Ca       '+str(7.154922011 +xd)   + '   '   +str(-0.006220291 +yd) +'  7.737352144  0   0   1\n')
    f.write(' Ca      '+str(-0.112830884 +xd)   + '   '    +str(6.377264464 +yd) +'  7.842326594  0   0   1\n')
    f.write(' Ca       '+str(3.488066618 +xd)   + '   '    +str(6.332779357 +yd) +'  7.703390501  0   0   1\n')
    f.write(' Ca       '+str(5.377144214 +xd)   + '   '    +str(3.178790854 +yd) +'  7.735175137  0   0   1\n')
    f.write(' Ca       '+str(7.261347371 +xd)   + '   '    +str(6.348064607 +yd) +'  7.682074227  0   0   1\n')
    f.write(' Ca      '+str(-0.126569111 +xd)   + '   '    +str(4.272186130 +yd) +' 10.251605044  0   0   1\n')
    f.write(' Ca       '+str(1.647520294 +xd)   + '   '    +str(1.069773416 +yd) +'  9.864335432  0   0   1\n')
    f.write(' Ca       '+str(5.429912498 +xd)   + '   '    +str(1.045912275 +yd) +' 10.254667895  0   0   1\n')
    f.write(' Ca      '+str(-1.930623221 +xd)   + '   '    +str(7.398783420 +yd) +' 10.214729740  0   0   1\n')
    f.write(' Ca       '+str(3.506006032 +xd)   + '   '    +str(4.198985099 +yd) +' 10.155974588  0   0   1\n')
    f.write(' Ca       '+str(1.712617058 +xd)   + '   '    +str(7.411227324 +yd) +' 10.275585732  0   0   1\n')
    f.write(' Ca       '+str(9.092423933 +xd)   + '   '    +str(1.055009279 +yd) +'  9.903238477  0   0   1\n')
    f.write(' Ca       '+str(5.367736725 +xd)   + '   '    +str(7.485620211 +yd) +'  9.842794380  0   0   1\n')
    f.write(' Ca       '+str(7.193573920 +xd)   + '   '    +str(4.293674471 +yd) +' 10.164635296  0   0   1\n')
    f.write(' N       '+str(-1.956067844 +xd)   + '   '    +str(5.320175826 +yd) +'  8.826946943  0   0   1\n')
    f.write(' N       '+str(-0.137962391 +xd)   + '   '    +str(2.125251046 +yd) +'  8.561760039  0   0   1\n')
    f.write(' N        '+str(3.522410620 +xd)   + '   '    +str(2.153500620 +yd) +'  8.702140151  0   0   1\n')
    f.write(' N       '+str(-3.824379405 +xd)   + '   '    +str(8.516854075 +yd) +'  8.512149587  0   0   1\n')
    f.write(' N       '+str(-0.110256218 +xd)   + '   '    +str(8.483560473 +yd) +'  8.892710012  0   0   1\n')
    f.write(' N        '+str(1.705619818 +xd)   + '   '    +str(5.306895449 +yd) +'  8.833110801  0   0   1\n')
    f.write(' N        '+str(3.557156825 +xd)   + '   '    +str(8.502828176 +yd) +'  8.562214710  0   0   1\n')
    f.write(' N        '+str(5.380127405 +xd)   + '   '    +str(5.316967579 +yd) +'  8.660697812  0   0   1\n')
    f.write(' N        '+str(7.210078993 +xd)   + '   '    +str(2.132601003 +yd) +'  8.740716652  0   0   1\n')
    f.write(' Au      -0.030983916   0.012852927  12.244727564  0   0   1\n')
    f.write(' Au       1.564775861   0.862898886  14.691027781  0   0   1\n')
    f.write(' Au       0.059517782   1.739404428  19.055253919  0   0   1\n')
    f.write(' Au       0.191459587  -0.062862269  21.589828916  0   0   1\n')
    f.write(' Au       2.808241085  -0.077147494  12.189179794  0   0   1\n')
    f.write(' Au       4.265839931   0.867303365  15.452629197  0   0   1\n')
    f.write(' Au       2.458517811   1.554660178  17.519105836  0   0   1\n')
    f.write(' Au       3.018065791  -0.077393352  21.598804035  0   0   1\n')
    f.write(' Au      -1.375336215   2.462031697  12.368293952  0   0   1\n')
    f.write(' Au       0.198879655   2.745521485  16.235057613  0   0   1\n')
    f.write(' Au      -1.250198039   4.118361806  18.188416301  0   0   1\n')
    f.write(' Au      -1.265063538   2.357719012  21.655103106  0   0   1\n')
    f.write(' Au       1.380261599   2.560168989  12.380580202  0   0   1\n')
    f.write(' Au       2.687586431   3.378593579  14.982838156  0   0   1\n')
    f.write(' Au       1.653468638   3.958176894  18.730138714  0   0   1\n')
    f.write(' Au       1.634593365   2.299304946  21.495821599  0   0   1\n')
    f.write(' Au       5.361751279  -0.441245299  13.114438781  0   0   1\n')
    f.write(' Au       7.024322161   0.610883154  15.276934436  0   0   1\n')
    f.write(' Au       5.134973458   1.832396002  18.352648218  0   0   1\n')
    f.write(' Au       5.785512380  -0.103870365  20.874054508  0   0   1\n')
    f.write(' Au       7.946754678   0.072067245  12.367247251  0   0   1\n')
    f.write(' Au       9.777564208   0.809028724  14.762505604  0   0   1\n')
    f.write(' Au       8.778626720   1.549209827  17.454615244  0   0   1\n')
    f.write(' Au       8.395017467  -0.009486725  21.743807548  0   0   1\n')
    f.write(' Au       4.172848955   2.259092724  12.644796671  0   0   1\n')
    f.write(' Au       5.580005938   3.340296084  15.203905774  0   0   1\n')
    f.write(' Au       4.233757226   4.277289165  17.525274420  0   0   1\n')
    f.write(' Au       4.461794380   2.257334531  21.355537089  0   0   1\n')
    f.write(' Au       6.868616714   2.648486019  12.658895864  0   0   1\n')
    f.write(' Au       8.392224663   3.168396368  15.121612277  0   0   1\n')
    f.write(' Au       6.969496798   3.833243487  17.901501762  0   0   1\n')
    f.write(' Au       7.247402987   2.220380699  20.227890250  0   0   1\n')
    f.write(' Au      -2.385909505   4.933079831  12.847299645  0   0   1\n')
    f.write(' Au      -0.839000841   5.248131600  15.372983741  0   0   1\n')
    f.write(' Au      -2.863080301   6.356163421  18.470578269  0   0   1\n')
    f.write(' Au      -2.484023048   4.663949655  20.677222672  0   0   1\n')
    f.write(' Au      -0.162301708   6.396240261  12.476855059  0   0   1\n')
    f.write(' Au       1.659557698   5.605336126  16.400327168  0   0   1\n')
    f.write(' Au       0.158638171   6.495409201  18.571438650  0   0   1\n')
    f.write(' Au       0.227834668   4.669863046  21.514138550  0   0   1\n')
    f.write(' Au      -4.140981198   6.966713321  12.184287804  0   0   1\n')
    f.write(' Au      -3.160199973   8.420007588  16.507421204  0   0   1\n')
    f.write(' Au      -4.096238156   8.834279712  19.086177978  0   0   1\n')
    f.write(' Au      -3.872418965   7.011966834  21.576915180  0   0   1\n')
    f.write(' Au      -2.203462513   7.450195448  14.119747614  0   0   1\n')
    f.write(' Au      -0.233548064   7.872258882  16.103265664  0   0   1\n')
    f.write(' Au      -1.362241850   8.848872178  18.639526773  0   0   1\n')
    f.write(' Au      -1.135496628   7.054795479  21.130837811  0   0   1\n')
    f.write(' Au       2.119605732   5.100053558  12.739976071  0   0   1\n')
    f.write(' Au       4.055678800   5.816630257  15.030395716  0   0   1\n')
    f.write(' Au       3.019926003   6.388778687  18.838409020  0   0   1\n')
    f.write(' Au       3.062711549   4.637852906  21.187855286  0   0   1\n')
    f.write(' Au       5.174061943   4.799768745  12.481871474  0   0   1\n')
    f.write(' Au       7.249185291   5.639589930  15.608346520  0   0   1\n')
    f.write(' Au       5.535233939   6.865249974  17.560365994  0   0   1\n')
    f.write(' Au       5.846515853   4.586609616  21.404068621  0   0   1\n')
    f.write(' Au       1.829752952   7.461775499  14.251519309  0   0   1\n')
    f.write(' Au       3.415793854   8.282393118  16.417111580  0   0   1\n')
    f.write(' Au       1.416764665   9.159413551  18.241809480  0   0   1\n')
    f.write(' Au       1.633119210   7.057151484  21.526940516  0   0   1\n')
    f.write(' Au       3.919971573   7.226736007  12.368306324  0   0   1\n')
    f.write(' Au       5.771180204   7.912449854  14.729912766  0   0   1\n')
    f.write(' Au       4.112022731   8.883408524  19.090786286  0   0   1\n')
    f.write(' Au       4.401211645   7.037735915  21.432581929  0   0   1\n')

    f.close()



def make_Ag_MoS2_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.96 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('S  32.065 S.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Ag 107.868 Ag.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  5.50179766915   0.0   0.000000000\n')
    f.write(' -2.750898834070916   4.764696549409532    0.000000000\n')
    f.write('  0.000000000   0.000000000  40.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write(' Ag       0.056541381   0.097561818   7.817935490  0   0   1\n')
    f.write(' Ag       2.788990426   0.097475899   7.788519829  0   0   1\n')
    f.write(' Ag      -1.310139906   2.463171768   7.788675596  0   0   1\n')
    f.write(' Ag       1.441400910   2.495274593   7.780377715  0   0   1\n')
    f.write(' Ag       1.433357005   0.896338908  10.234724986  0   0   1\n')
    f.write(' Ag       4.181448083   0.892725522  10.238337212  0   0   1\n')
    f.write(' Ag       0.057611889   3.275391285  10.223567313  0   0   1\n')
    f.write(' Ag       2.809391053   3.274319854  10.234590762  0   0   1\n')
    f.write(' Ag       0.057996373   1.689673357  12.660073489  0   0   1\n')
    f.write(' Ag       2.808555438   1.689348019  12.665429200  0   0   1\n')
    f.write(' Ag      -1.317371440   4.071306569  12.657867971  0   0   1\n')
    f.write(' Ag       1.432829642   4.071288227  12.659954802  0   0   1\n')
    f.write(' Ag       0.058516976   0.101782037  15.079023887  0   0   1\n')
    f.write(' Ag       2.810407692   0.101558984  15.080302408  0   0   1\n')
    f.write(' Ag      -1.317502855   2.485271685  15.080548177  0   0   1\n')
    f.write(' Ag       1.433196559   2.483283961  15.081161794  0   0   1\n')
    f.write(' Ag       1.434344325   0.896847499  17.518430008  0   0   0\n')
    f.write(' Ag       4.185125940   0.895763353  17.518918312  0  0   0\n')
    f.write(' Ag       0.059186783   3.278625416  17.518429763  0  0   0\n')
    f.write(' Ag       2.810692775   3.278570515  17.518989438  0  0   0\n')
    f.write(' Ag       0.059601943   1.690757871  19.987355324  0  0   0\n')
    f.write(' Ag       2.810425129   1.690712270  19.988463049  0  0   0\n')
    f.write(' Ag      -1.315884967   4.073032903  19.986518830  0  0   0\n')
    f.write(' Ag       1.434936038   4.073068515  19.987263756  0  0   0\n')

    f.write(' Mo        '+str(-1.839085155  +xd)   + '   '  +str(   3.169859232  +yd) +'    3.613283627  0   0   1\n')
    f.write(' S         '+str(-0.922806320  +xd)   + '   '  +str(   1.581746492  +yd) +'    2.050117024  0   0   1\n')
    f.write(' S         '+str(-0.921149976  +xd)   + '   '  +str(   1.583960146  +yd) +'    5.180792050  0   0   1\n')
    f.write(' Mo        '+str( 0.911709684  +xd)   + '   '  +str(   1.578837015  +yd) +'    3.611601229  0   0   1\n')
    f.write(' S         '+str( 1.831294445  +xd)   + '   '  +str(  -0.008622610  +yd) +'    2.050115660  0   0   1\n')
    f.write(' S         '+str( 1.832303560  +xd)   + '   '  +str(  -0.005799479  +yd) +'    5.180830371  0   0   1\n')
    f.write(' Mo        '+str( 3.664734958  +xd)   + '   '  +str(  -0.007625999  +yd) +'    3.613321707  0   0   1\n')
    f.write(' S         '+str( 1.828791227  +xd)   + '   '  +str(   3.167793215  +yd) +'    2.050182926  0   0   1\n')
    f.write(' S         '+str( 1.826333194  +xd)   + '   '  +str(   3.163128439  +yd) +'    5.177893353  0   0   1\n')
    f.close()



def make_Ag_Ca2N_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Ag 107.868 Ag.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  11.003595338   0.000000000   0.000000000\n')
    f.write(' -5.501797668   9.529393098    0.000000000\n')
    f.write('  0.000000000   0.000000000  40.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write(' Ca            '+str(-0.5840608020  +xd)   + '   '  +str(0.2571531679   +yd) +'        8.0053598572  0   0   1\n')
    f.write(' Ca            '+str(-0.6066733150  +xd)   + '   '  +str(4.4899396935   +yd) +'       10.4316935545  0   0   1\n')
    f.write(' Ca            '+str(1.2676247638   +xd)   + '   '  +str(3.4710520672   +yd) +'        8.0488478544  0   0   1\n')
    f.write(' Ca            '+str(1.2489833001   +xd)   + '   '  +str(1.3184323216   +yd) +'       10.4018702120  0   0   1\n')
    f.write(' Ca            '+str(3.0905831040   +xd)   + '   '  +str(0.2737478282   +yd) +'        8.0076300674  0   0   1\n')
    f.write(' Ca            '+str(4.9284353956   +xd)   + '   '  +str(1.3315196133   +yd) +'       10.3866320396  0   0   1\n')
    f.write(' Ca            '+str(6.7581870149   +xd)   + '   '  +str(0.2452735834   +yd) +'        8.0034313135  0   0   1\n')
    f.write(' Ca            '+str(-2.4224860797  +xd)   + '   '  +str(7.6817422543   +yd) +'       10.3931448514  0   0   1\n')
    f.write(' Ca            '+str(-0.5519899275  +xd)   + '   '  +str(6.6179308375   +yd) +'        8.0419543986  0   0   1\n')
    f.write(' Ca            '+str(3.0994924470   +xd)   + '   '  +str(4.4860386662   +yd) +'       10.5180981035  0   0   1\n')
    f.write(' Ca            '+str(1.2383343231   +xd)   + '   '  +str(7.7022189024   +yd) +'       10.4996399019  0   0   1\n')
    f.write(' Ca            '+str(3.0758291175   +xd)   + '   '  +str(6.6130280609   +yd) +'        8.0859047168  0   0   1\n')
    f.write(' Ca            '+str(4.8995777121   +xd)   + '   '  +str(3.4633434990   +yd) +'        8.0112076588  0   0   1\n')
    f.write(' Ca            '+str(8.5861009563   +xd)   + '   '  +str(1.3465621429   +yd) +'       10.3624313992  0   0   1\n')
    f.write(' Ca            '+str(4.9388782617   +xd)   + '   '  +str(7.6920823758   +yd) +'       10.4520255322  0   0   1\n')
    f.write(' Ca            '+str(6.7412085380   +xd)   + '   '  +str(6.6252078786   +yd) +'        8.0055016403  0   0   1\n')
    f.write(' Ca            '+str(6.7719703060   +xd)   + '   '  +str(4.4952755722   +yd) +'       10.3652942307  0   0   1\n')
    f.write(' Ca            '+str(8.6017851007   +xd)   + '   '  +str(3.4528123301   +yd) +'        7.9639319128  0   0   1\n')
    f.write(' N             '+str(4.9215931923   +xd)   + '   '  +str(5.5649849455   +yd) +'        9.1376350592  0   0   1\n')
    f.write(' N             '+str(6.7585460474   +xd)   + '   '  +str(2.3893675073   +yd) +'        9.0290041702  0   0   1\n')
    f.write(' N             '+str(-0.5775451281  +xd)   + '   '  +str(8.7385722375   +yd) +'        9.0984812908  0   0   1\n')
    f.write(' N             '+str(1.2552624002   +xd)   + '   '  +str(5.5635857236   +yd) +'        9.2066313619  0   0   1\n')
    f.write(' N             '+str(3.0865166224   +xd)   + '   '  +str(8.7397523300   +yd) +'        9.1317898293  0   0   1\n')
    f.write(' N             '+str(-2.4119912413  +xd)   + '   '  +str(5.5622350412   +yd) +'        9.0674333537  0   0   1\n')
    f.write(' N             '+str(-0.5805156191  +xd)   + '   '  +str(2.3894502494   +yd) +'        9.0709208481  0   0   1\n')
    f.write(' N             '+str(3.0882215607   +xd)   + '   '  +str(2.3900724042   +yd) +'        9.1089714078  0   0   1\n')
    f.write(' N             '+str(-4.2489445050  +xd)   + '   '  +str(8.7401076778   +yd) +'        9.0717644331  0   0   1\n')
    f.write(' Ag            0.0378794661        0.0790957734       13.0059455036  0   0   1\n')
    f.write(' Ag            5.5900329179        0.0702325847       13.0498490797  0   0   1\n')
    f.write(' Ag           -2.6379644055        4.8871076934       12.9631311534  0   0   1\n')
    f.write(' Ag            2.8400727175        4.8345294276       13.4471872293  0   0   1\n')
    f.write(' Ag            2.8251211902        0.0900566532       12.9207740156  0   0   1\n')
    f.write(' Ag            8.3057984011        0.0472770458       13.0214547610  0   0   1\n')
    f.write(' Ag            0.1291445318        4.8623811797       13.2389925312  0   0   1\n')
    f.write(' Ag            5.5502204341        4.8717727647       13.0328496252  0   0   1\n')
    f.write(' Ag           -1.2949778224        2.4904469905       12.9675214469  0   0   1\n')
    f.write(' Ag            4.1807462352        2.4859705623       13.0562397544  0   0   1\n')
    f.write(' Ag           -4.0581615333        7.2425064292       12.9949187589  0   0   1\n')
    f.write(' Ag            1.4796215231        7.2074057351       13.3961008910  0   0   1\n')
    f.write(' Ag            1.4473031445        2.4857618982       13.1406586086  0   0   1\n')
    f.write(' Ag            6.9433402374        2.4843612632       12.8737451529  0   0   1\n')
    f.write(' Ag           -1.2385097729        7.2320036849       13.1032213663  0   0   1\n')
    f.write(' Ag            4.1839053441        7.2080046640       13.2872766556  0   0   1\n')
    f.write(' Ag            1.4381793633        0.9123866315       15.5658990234  0   0   1\n')
    f.write(' Ag            6.9396010823        0.8495656144       15.5137199633  0   0   1\n')
    f.write(' Ag           -1.3167216986        5.6668923484       15.6083742924  0   0   1\n')
    f.write(' Ag            4.2301502529        5.6404904945       15.7964369262  0   0   1\n')
    f.write(' Ag            4.1898751355        0.9012633631       15.5305830144  0   0   1\n')
    f.write(' Ag            9.6872771208        0.8639442851       15.5132444368  0   0   1\n')
    f.write(' Ag            1.4350479883        5.6266024111       15.8469611759  0   0   1\n')
    f.write(' Ag            6.9584179329        5.6778961125       15.5028022106  0   0   1\n')
    f.write(' Ag            0.0779529451        3.2649049729       15.6290368858  0   0   1\n')
    f.write(' Ag            5.6065456105        3.2695995543       15.4964985495  0   0   1\n')
    f.write(' Ag           -2.7111394558        8.0740875716       15.5127237704  0   0   1\n')
    f.write(' Ag            2.8002444128        8.0024797395       15.7632021949  0   0   1\n')
    f.write(' Ag            2.8569367815        3.2634334778       15.7842525022  0   0   1\n')
    f.write(' Ag            8.3436270528        3.2599103506       15.4567245155  0   0   1\n')
    f.write(' Ag            0.0260928270        8.0281512710       15.7060022153  0   0   1\n')
    f.write(' Ag            5.5573029518        8.0457583917       15.6456122606  0   0   1\n')
    f.write(' Ag            0.0895980004        1.6952686763       18.0129549246  0   0   1\n')
    f.write(' Ag            5.5584012191        1.6867280888       17.9401472198  0   0   1\n')
    f.write(' Ag           -2.7242151424        6.4519255765       17.9768976876  0   0   1\n')
    f.write(' Ag            2.7972931493        6.4728461672       18.2599299697  0   0   1\n')
    f.write(' Ag            2.8266080130        1.6988232487       18.1176010973  0   0   1\n')
    f.write(' Ag            8.3298490555        1.6701062196       17.9207824633  0   0   1\n')
    f.write(' Ag            0.0084560435        6.4576572166       18.1832530522  0   0   1\n')
    f.write(' Ag            5.5614511800        6.4625207730       18.1304020252  0   0   1\n')
    f.write(' Ag           -1.3247822809        4.0595452590       17.9930568876  0   0   1\n')
    f.write(' Ag            4.1957861878        4.0796495390       18.1886598102  0   0   1\n')
    f.write(' Ag           -4.0547349853        8.8487135733       18.0182608867  0   0   1\n')
    f.write(' Ag            1.4184346008        8.8447854244       18.1120580785  0   0   1\n')
    f.write(' Ag            1.4290122478        4.0840285160       18.2281338470  0   0   1\n')
    f.write(' Ag            6.9245499906        4.0618879146       17.9091881092  0   0   1\n')
    f.write(' Ag           -1.3395182880        8.8433708415       18.0486049358  0   0   1\n')
    f.write(' Ag            4.2028484484        8.8465043281       18.0918686009  0   0   1\n')
    f.write(' Ag            0.0610505352        0.1154357202       20.5395774895  0   0   1\n')
    f.write(' Ag            5.5511631705        0.1203603921       20.5401954393  0   0   1\n')
    f.write(' Ag           -2.7025674088        4.8752755349       20.4407893210  0   0   1\n')
    f.write(' Ag            2.8025258625        4.8792829419       20.7317337630  0   0   1\n')
    f.write(' Ag            2.8056241886        0.1072582161       20.5754640064  0   0   1\n')
    f.write(' Ag            8.3074602349        0.1265139036       20.5500477676  0   0   1\n')
    f.write(' Ag            0.0390265896        4.8757248873       20.6655819029  0   0   1\n')
    f.write(' Ag            5.5604885593        4.8759061577       20.6175488724  0   0   1\n')
    f.write(' Ag           -1.3237120714        2.4965669762       20.4591927741  0   0   1\n')
    f.write(' Ag            4.1823381674        2.4917983828       20.6197477064  0   0   1\n')
    f.write(' Ag           -4.0674367728        7.2593853945       20.5500968005  0   0   1\n')
    f.write(' Ag            1.4217542042        7.2728703405       20.6993415385  0   0   1\n')
    f.write(' Ag            1.4258901542        2.4890603334       20.6485335519  0   0   1\n')
    f.write(' Ag            6.9240751294        2.4976694346       20.3931298875  0   0   1\n')
    f.write(' Ag           -1.3304377430        7.2601961086       20.6023779786  0   0   1\n')
    f.write(' Ag            4.1872391620        7.2707310282       20.6788949364  0   0   1\n')

    f.close()




def make_Cu_Ca2N_var_z(struc_filename, xd=0, yd=0):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ca 40.0780 Ca.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('N  14.0067 N.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('Cu 63.546 Cu.pbe-dn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('  12.70964067       0.00003371   0.000000000\n')
    f.write(' -6.35479125      11.00685249     0.000000000\n')
    f.write('  0.000000000   0.000000000  40.000000000\n')
    f.write('\n')

    f.write('ATOMIC_POSITIONS (angstrom)\n')
    f.write(' Ca            '+str(-0.0242227927  +xd)   + '   '  +str(-0.0419806875 +yd) +'        7.9591918407  0   0   1\n')
    f.write(' Ca            '+str(6.3262966146   +xd)   + '   '  +str(-0.0401433467 +yd) +'        7.9605984593  0   0   1\n')
    f.write(' Ca            '+str(-3.1979029537  +xd)   + '   '  +str(5.4586392262  +yd) +'        7.9606003774  0   0   1\n')
    f.write(' Ca            '+str(3.1555862227   +xd)   + '   '  +str(5.4655962302  +yd) +'        7.9579392924  0   0   1\n')
    f.write(' Ca            '+str(4.2108880229   +xd)   + '   '  +str(-0.0444527359 +yd) +'       10.3854093892  0   0   1\n')
    f.write(' Ca            '+str(10.5693872783  +xd)   + '   '  +str(-0.0389872752 +yd) +'       10.3543289449  0   0   1\n')
    f.write(' Ca            '+str(1.0363731959   +xd)   + '   '  +str(5.4648797927  +yd) +'       10.3543129320  0   0   1\n')
    f.write(' Ca            '+str(7.3894578538   +xd)   + '   '  +str(5.4609647406  +yd) +'       10.3516829505  0   0   1\n')
    f.write(' Ca            '+str(-2.1439300932  +xd)   + '   '  +str(3.6244769612  +yd) +'       10.3853832805  0   0   1\n')
    f.write(' Ca            '+str(4.2145816248   +xd)   + '   '  +str(3.6299432420  +yd) +'       10.3543199759  0   0   1\n')
    f.write(' Ca            '+str(-5.3184423542  +xd)   + '   '  +str(9.1338242841  +yd) +'       10.3543388838  0   0   1\n')
    f.write(' Ca            '+str(1.0346382019   +xd)   + '   '  +str(9.1299233810  +yd) +'       10.3516887646  0   0   1\n')
    f.write(' Ca            '+str(-0.0285103849  +xd)   + '   '  +str(3.6287913760  +yd) +'        7.9605885957  0   0   1\n')
    f.write(' Ca            '+str(6.3306030001   +xd)   + '   '  +str(3.6269895713  +yd) +'        7.9591684434  0   0   1\n')
    f.write(' Ca            '+str(-3.1992233772  +xd)   + '   '  +str(9.1345487190  +yd) +'        7.9579625479  0   0   1\n')
    f.write(' Ca            '+str(3.1569282452   +xd)   + '   '  +str(9.1276179778  +yd) +'        7.9605870782  0   0   1\n')
    f.write(' Ca            '+str(1.0346400045   +xd)   + '   '  +str(1.7920066906  +yd) +'       10.3517072822  0   0   1\n')
    f.write(' Ca            '+str(7.3911709174   +xd)   + '   '  +str(1.7959517661  +yd) +'       10.3543296997  0   0   1\n')
    f.write(' Ca            '+str(-2.1402290265  +xd)   + '   '  +str(7.2988919307  +yd) +'       10.3543356080  0   0   1\n')
    f.write(' Ca            '+str(4.2109061017   +xd)   + '   '  +str(7.2934545304  +yd) +'       10.3853706205  0   0   1\n')
    f.write(' Ca            '+str(3.1569106603   +xd)   + '   '  +str(1.7896938983  +yd) +'        7.9605904974  0   0   1\n')
    f.write(' Ca            '+str(9.5104009903   +xd)   + '   '  +str(1.7966626615  +yd) +'        7.9579493121  0   0   1\n')
    f.write(' Ca            '+str(-0.0242129351  +xd)   + '   '  +str(7.2959337449  +yd) +'        7.9591766749  0   0   1\n')
    f.write(' Ca            '+str(6.3263296344   +xd)   + '   '  +str(7.2977714560  +yd) +'        7.9605899863  0   0   1\n')
    f.write(' N             '+str(-1.0841031195  +xd)   + '   '  +str(1.7949950332  +yd) +'        9.0084568169  0   0   1\n')
    f.write(' N             '+str(5.2709670486   +xd)   + '   '  +str(1.7916546822  +yd) +'        9.0137268322  0   0   1\n')
    f.write(' N             '+str(-4.2601935023  +xd)   + '   '  +str(7.2968936032  +yd) +'        9.0045939606  0   0   1\n')
    f.write(' N             '+str(2.0966016797   +xd)   + '   '  +str(7.2965424471  +yd) +'        9.0084392057  0   0   1\n')
    f.write(' N             '+str(2.0965838119   +xd)   + '   '  +str(-0.0413722436 +yd) +'        9.0084661260  0   0   1\n')
    f.write(' N             '+str(8.4494324175   +xd)   + '   '  +str(-0.0409785102 +yd) +'        9.0045799704  0   0   1\n')
    f.write(' N             '+str(-1.0838488721  +xd)   + '   '  +str(5.4605909756  +yd) +'        9.0137305769  0   0   1\n')
    f.write(' N             '+str(5.2707245601   +xd)   + '   '  +str(5.4639672908  +yd) +'        9.0084396908  0   0   1\n')
    f.write(' N             '+str(2.0946187253   +xd)   + '   '  +str(3.6279574755  +yd) +'        9.0045740852  0   0   1\n')
    f.write(' N             '+str(8.4514107540   +xd)   + '   '  +str(3.6275990042  +yd) +'        9.0084498903  0   0   1\n')
    f.write(' N             '+str(-1.0840889570  +xd)   + '   '  +str(9.1329021060  +yd) +'        9.0084585893  0   0   1\n')
    f.write(' N             '+str(5.2709797715   +xd)   + '   '  +str(9.1295553923  +yd) +'        9.0137259723  0   0   1\n')

    f.write(' Cu            0.0936795788        0.1622833032       12.8227677378  0   0   1\n')
    f.write(' Cu            4.3366751610        0.1734556686       13.0753613881  0   0   1\n')
    f.write(' Cu            8.5692979725        0.1667148637       12.8128846172  0   0   1\n')
    f.write(' Cu           -2.0181371436        3.8423769801       13.0753483364  0   0   1\n')
    f.write(' Cu            2.2144993539        3.8356114716       12.8128866863  0   0   1\n')
    f.write(' Cu            6.4484908335        3.8311798152       12.8227232538  0   0   1\n')
    f.write(' Cu           -4.1402872985        7.5045740569       12.8128908736  0   0   1\n')
    f.write(' Cu            0.0936749724        7.5001544158       12.8227391206  0   0   1\n')
    f.write(' Cu            4.3366792958        7.5113441753       13.0753275138  0   0   1\n')
    f.write(' Cu            0.0842376328        2.6309254885       12.9123967509  0   0   1\n')
    f.write(' Cu            4.3336784175        2.5958003049       12.9576170771  0   0   1\n')
    f.write(' Cu            8.5849114487        2.6198234219       12.8840786610  0   0   1\n')
    f.write(' Cu           -2.0211214948        6.2647409035       12.9576353387  0   0   1\n')
    f.write(' Cu            2.2300886868        6.2887422258       12.8840715276  0   0   1\n')
    f.write(' Cu            6.4390659838        6.2998637072       12.9123543461  0   0   1\n')
    f.write(' Cu           -4.1246904340        9.9577380383       12.8841151733  0   0   1\n')
    f.write(' Cu            0.0842394268        9.9688408137       12.9123784725  0   0   1\n')
    f.write(' Cu            4.3336997358        9.9337579121       12.9576448770  0   0   1\n')
    f.write(' Cu            2.2363271502        1.3884283644       12.9124066286  0   0   1\n')
    f.write(' Cu            6.4359975947        1.3820404004       12.9576154686  0   0   1\n')
    f.write(' Cu           10.6859968133        1.4067969358       12.8840977320  0   0   1\n')
    f.write(' Cu            0.0811952203        5.0509701079       12.9576259790  0   0   1\n')
    f.write(' Cu            4.3311788980        5.0756704621       12.8840617563  0   0   1\n')
    f.write(' Cu            8.5911584318        5.0573578476       12.9123599956  0   0   1\n')
    f.write(' Cu           -2.0236146588        8.7446615711       12.8841063913  0   0   1\n')
    f.write(' Cu            2.2363233025        8.7263295827       12.9123809535  0   0   1\n')
    f.write(' Cu            6.4360388060        8.7199685226       12.9576388905  0   0   1\n')
    f.write(' Cu            1.5074747818        0.1676572812       15.0066508492  0   0   1\n')
    f.write(' Cu            5.7544828896        0.1672570377       15.0962059196  0   0   1\n')
    f.write(' Cu            9.9776071618        0.1727923486       15.0173276111  0   0   1\n')
    f.write(' Cu           -0.6003264245        3.8362146434       15.0961819565  0   0   1\n')
    f.write(' Cu            3.6227713260        3.8417000216       15.0172882725  0   0   1\n')
    f.write(' Cu            7.8622941008        3.8365801370       15.0065866373  0   0   1\n')
    f.write(' Cu           -2.7319959789        7.5106595344       15.0173390495  0   0   1\n')
    f.write(' Cu            1.5074601620        7.5055314620       15.0065994023  0   0   1\n')
    f.write(' Cu            5.7544835430        7.5051794875       15.0961778021  0   0   1\n')
    f.write(' Cu           -0.6085631582        1.3893458300       15.0066316264  0   0   1\n')
    f.write(' Cu            3.6224188449        1.3982069731       15.0961989805  0   0   1\n')
    f.write(' Cu            7.8704405968        1.3893370790       15.0173136686  0   0   1\n')
    f.write(' Cu           -2.7323780945        5.0671455946       15.0961871673  0   0   1\n')
    f.write(' Cu            1.5156087211        5.0582806958       15.0172893099  0   0   1\n')
    f.write(' Cu            5.7462408621        5.0582627737       15.0065715320  0   0   1\n')
    f.write(' Cu           -4.8391726329        8.7272453701       15.0173409916  0   0   1\n')
    f.write(' Cu           -0.6085484684        8.7272271256       15.0066335114  0   0   1\n')
    f.write(' Cu            3.6224205965        8.7361205737       15.0961842966  0   0   1\n')
    f.write(' Cu            1.5123435093        2.6194775737       15.0109030106  0   0   1\n')
    f.write(' Cu            5.7488136048        2.6193263219       15.0402948898  0   0   1\n')
    f.write(' Cu            9.9776300870        2.6059421280       15.0746886147  0   0   1\n')
    f.write(' Cu           -0.6059841315        6.2882906250       15.0403161896  0   0   1\n')
    f.write(' Cu            3.6227740343        6.2748601228       15.0746533214  0   0   1\n')
    f.write(' Cu            7.8671820918        6.2884313409       15.0108861835  0   0   1\n')
    f.write(' Cu           -2.7319890852        9.9438274012       15.0747277586  0   0   1\n')
    f.write(' Cu            1.5123566697        9.9573725464       15.0108978337  0   0   1\n')
    f.write(' Cu            5.7488447291        9.9572946593       15.0403345216  0   0   1\n')
    f.write(' Cu            2.9185304998        0.1675122934       17.1486979702  0   0   1\n')
    f.write(' Cu            7.1646073279        0.1691122351       17.1358589585  0   0   1\n')
    f.write(' Cu           11.3929762133        0.1706547095       17.1158511487  0   0   1\n')
    f.write(' Cu            0.8097914970        3.8380814025       17.1358405912  0   0   1\n')
    f.write(' Cu            5.0381508026        3.8395586262       17.1158136402  0   0   1\n')
    f.write(' Cu            9.2733758644        3.8364700338       17.1486704990  0   0   1\n')
    f.write(' Cu           -1.3166283156        7.5085180422       17.1158638941  0   0   1\n')
    f.write(' Cu            2.9185131894        7.5054170692       17.1486535022  0   0   1\n')
    f.write(' Cu            7.1646359670        7.5070453200       17.1358541870  0   0   1\n')
    f.write(' Cu            0.8041118259        1.3927956094       17.1048215442  0   0   1\n')
    f.write(' Cu            5.0417466309        1.3946584058       17.1616688515  0   0   1\n')
    f.write(' Cu            9.2746510621        1.3883568062       17.1319547400  0   0   1\n')
    f.write(' Cu           -1.3130429737        5.0636190831       17.1616651700  0   0   1\n')
    f.write(' Cu            2.9198072202        5.0572983831       17.1319212385  0   0   1\n')
    f.write(' Cu            7.1589496846        5.0617314307       17.1047917346  0   0   1\n')
    f.write(' Cu           -3.4349554589        8.7262459687       17.1319886527  0   0   1\n')
    f.write(' Cu            0.8041186292        8.7306829114       17.1048123430  0   0   1\n')
    f.write(' Cu            5.0417559138        8.7325940346       17.1616631064  0   0   1\n')
    f.write(' Cu           -1.3141983596        2.6112806478       17.1486732007  0   0   1\n')
    f.write(' Cu            2.9189482720        2.6203477647       17.1358459114  0   0   1\n')
    f.write(' Cu            7.1608978487        2.6139977638       17.1158329358  0   0   1\n')
    f.write(' Cu           -3.4358213171        6.2892845798       17.1358629512  0   0   1\n')
    f.write(' Cu            0.8060729936        6.2829655005       17.1158288030  0   0   1\n')
    f.write(' Cu            5.0406001707        6.2802328216       17.1486473035  0   0   1\n')
    f.write(' Cu           -5.5487192252        9.9519206178       17.1158645775  0   0   1\n')
    f.write(' Cu           -1.3141833594        9.9491738013       17.1486997543  0   0   1\n')
    f.write(' Cu            2.9189677569        9.9582602992       17.1358494470  0   0   1\n')
    f.write(' Cu            0.0978869892        0.1695928843       19.2203722089  0   0   1\n')
    f.write(' Cu            4.3336173873        0.1681296986       19.2598157061  0   0   1\n')
    f.write(' Cu            8.5707540114        0.1691781407       19.2468865036  0   0   1\n')
    f.write(' Cu           -2.0211914606        3.8370744099       19.2598052438  0   0   1\n')
    f.write(' Cu            2.2159087761        3.8381225962       19.2468766586  0   0   1\n')
    f.write(' Cu            6.4527251946        3.8385258146       19.2203532556  0   0   1\n')
    f.write(' Cu           -4.1388489187        7.5070796206       19.2469113833  0   0   1\n')
    f.write(' Cu            0.0979107892        7.5074902544       19.2203416344  0   0   1\n')
    f.write(' Cu            4.3335989159        7.5060671314       19.2597759514  0   0   1\n')
    f.write(' Cu            0.0996663862        2.6133041318       19.2390816705  0   0   1\n')
    f.write(' Cu            4.3337179914        2.6177838868       19.2470625107  0   0   1\n')
    f.write(' Cu            8.5690385861        2.6136798067       19.2420752723  0   0   1\n')
    f.write(' Cu           -2.0210585146        6.2867494061       19.2470719868  0   0   1\n')
    f.write(' Cu            2.2141768922        6.2826457952       19.2420374773  0   0   1\n')
    f.write(' Cu            6.4545056630        6.2822800184       19.2390797893  0   0   1\n')
    f.write(' Cu           -4.1405875530        9.9515978951       19.2421006652  0   0   1\n')
    f.write(' Cu            0.0997067705        9.9512024297       19.2390820717  0   0   1\n')
    f.write(' Cu            4.3337374412        9.9557338769       19.2470563067  0   0   1\n')
    f.write(' Cu            2.2133257334        1.3929842032       19.2390949903  0   0   1\n')
    f.write(' Cu            6.4550321469        1.3930519024       19.2470658727  0   0   1\n')
    f.write(' Cu           10.6886378933        1.3899484276       19.2420729289  0   0   1\n')
    f.write(' Cu            0.1001986010        5.0620234325       19.2470500966  0   0   1\n')
    f.write(' Cu            4.3337944776        5.0588808225       19.2420483773  0   0   1\n')
    f.write(' Cu            8.5681883634        5.0619488162       19.2390958823  0   0   1\n')
    f.write(' Cu           -2.0209628154        8.7278316628       19.2420885123  0   0   1\n')
    f.write(' Cu            2.2133383526        8.7309014575       19.2390587862  0   0   1\n')
    f.write(' Cu            6.4550393038        8.7310173669       19.2470716581  0   0   1\n')
    f.close()






def make_In_prm_slab(struc_filename, scale_c=1.0, scale_a=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('In 114.818 In.pbe-dn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale_c*4.346597952)   + '   ' + str(scale_c*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale_a*0.0000000000)   + '   ' + str(scale_a*3.290719868)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('In       0.000000000   0.000000000  -0.156260207\n')
        f.write('In       2.173298975   1.645359934   1.511548741\n')
        f.write('In       0.000000000   0.000000000   3.191645046\n')
        f.write('In       2.173298975   1.645359934   4.840186918\n')
        f.write('In       0.000000000   0.000000000   6.508376684\n') 
        f.write('In       2.173298975   1.645359934   8.176603130\n')
        f.write('In       0.000000000   0.000000000   9.825147067\n')
        f.write('In       2.173298975   1.645359934  11.505294219\n')
        f.write('In       0.000000000   0.000000000  13.173058074\n') 
    f.close()


def make_Mo_prm_slab(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Mo 95.960 Mo.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*4.337749386)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-2.168874765)   + '   ' + str(scale*3.756601698)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    35.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Mo       0.000007083   0.000029965   2.935190559\n')
        f.write('Mo       2.168881741   1.252230441   3.686842375\n')
        f.write('Mo       0.000007008   2.504430943   4.447395227\n')
        f.write('Mo       0.000006945   0.000029737   5.449411497\n')
        f.write('Mo       2.168881554   1.252230188   6.328666764\n')
        f.write('Mo       0.000006768   2.504430731   7.254855679\n')
        f.write('Mo       0.000006767   0.000029648   8.176922876\n')
        f.write('Mo       2.168881410   1.252230276   9.098990129\n')
        f.write('Mo       0.000006671   2.504431024  10.025178622\n')
        f.write('Mo       0.000006730   0.000030164  10.904431649\n')
        f.write('Mo       2.168881376   1.252230952  11.906451019\n')
        f.write('Mo       0.000006672   2.504431579  12.667004938\n')
        f.write('Mo       0.000006756   0.000030640  13.418658666\n')
    
    f.close()


def make_Cu_prm_slab(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Cu 63.546 Cu.pbe-dn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*2.490606683)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.245303336)   + '   ' + str(scale*2.156928659)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Cu      -0.000000000  -0.000000000   0.053132603\n')
        f.write('Cu       1.245303337   0.718976223   2.103156959\n')
        f.write('Cu       0.000000000   1.437952436   4.166845023\n')
        f.write('Cu      -0.000000000  -0.000000000   6.216913264\n') 
    f.close()


def make_Ni_prm_slab(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ni 58.693 Ni.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*2.489015870)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.24450794)   + '   ' + str(scale*2.15555097)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ni       0.00000000       0.00000000       0.00000000\n')
        f.write('Ni       1.24450794       0.71851699       2.03227295\n')
        f.write('Ni       0.00000000       1.43703398       4.06454590\n')
        f.write('Ni       0.00000000       0.00000000       6.09681884\n') 
    f.close()

def make_Ni_Hex_bulk(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ni 58.693 Ni.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*2.489015870)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.24450794)   + '   ' + str(scale*2.15555097)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    4.06454590\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ni       0.00000000       0.00000000       0.00000000\n')
        f.write('Ni       1.24450794       0.71851699       2.03227295\n')
        f.write('Ni       0.00000000       1.43703398       4.06454590\n')
 
    f.close()




def make_Ti_prm_slab(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ti 47.867 Ti.pbe-spn-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*2.815611995)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.407805997)   + '   ' + str(scale*2.438391515)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ti       0.000000000   1.625594430   1.382905069\n')
        f.write('Ti       1.407805884   0.812797101   7.989088448\n')
        f.write('Ti       1.407805884   0.812797101   3.533642751\n')
        f.write('Ti       0.000000000   1.625594430   5.838363215\n') 
    f.close()



def make_Ag_prm_slab(struc_filename, scale=1.0, include_atoms = True):

    f = open(struc_filename, 'w+')
    f.write('ATOMIC_SPECIES\n')
    f.write('Ag 107.868 Ag.pbe-n-kjpaw_psl.1.0.0.UPF\n')
    f.write('\n')

    f.write('CELL_PARAMETERS (angstrom)\n')
    f.write('   ' +  str(scale*2.854038817)   + '   ' + str(scale*0.000000000)  + '    0.000000000\n')
    f.write('   ' +  str(scale*-1.427019408)   + '   ' + str(scale*2.471670112)  + '    0.000000000\n')
    f.write('    0.000000000    0.000000000    30.000000000\n')
 
    f.write('\n')
    
    if include_atoms == True:
        f.write('ATOMIC_POSITIONS (angstrom)\n')
        f.write('Ag      -0.000000000   0.000000000  -0.007753730\n')
        f.write('Ag       1.427019408   0.823890034   2.356338912\n')
        f.write('Ag       0.000000000   1.647780078   4.727675141\n')
        f.write('Ag      -0.000000000   0.000000000   7.091915276\n') 
    f.close()


if __name__ == "__main__":
    filename = sys.argv[1]
    prefix = sys.argv[2]
