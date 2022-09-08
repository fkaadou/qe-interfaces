# qe-interfaces

This repo contains python tools that are useful in understanding interfaces in heterostructures consisting of two or more materials within a density-functional theory framework. These tools are meant to be used in conjunction with Quantum ESPRESSO input and output files. I developed them during my MSc when I was studying Ca2N-MoS2 and Au-Ca2N-MoS2 heterostrucures, but they can be modified for use with other supercells. 

Check out the resulting paper: [Improved Charge Transfer and Barrier Lowering across a Au–MoS2 Interface through Insertion of a Layered Ca2N Electride](https://doi-org.ezproxy.library.dal.ca/10.1021/acs.jpcc.1c02142). It will come in handy in explaining what these tools are useful for.


## Bandstructure Plots with Atomic Contributions

Bandstructure plots can tell us a lot about materials; however, they can be very difficult to understand when the cell is consists of multiple materials interfaced togethor like in a heterostructure supercell. `band_plot_contribution_2_materials.py` and `band_plot_contribution_3_materials.py` get around the spaghetti mess of these kind of bandstructures by assigning a colored weight to each state in order to distinguish its composition as a function of contribution from the indivudual materials. They are used for 2 and 3 material heterostructures, respectively.

`band_plot_contribution_2_materials.py` and `band_plot_contribution_3_materials.py` require `pw.x` input and output files for `scf` calculations. As well, you will need to provide the file containing the projected bandstructure states. You can use Quantum ESPRESSO's `projwfc.x` to get the states projected onto localized atomic orbitals (`pdos.out` in the input_file_examples folder). 

## Potential Energy Curves

Constructing a supercell with two different materials is not a trivial matter. When interfacing two surfaces, it is important that the resulting contact corresponds to a minimum on the potential energy or sliding enegy curve. In order to insure that is the case, one must calculate this potential energy curve by constructing many cells where the position of the top material is translated relative to the bottom material in incremental amounts and then finding the optimized geometry or 'relaxing' the cell. Once all the different cell configuration are optimized, their total energy is extracted and plotted as a function of the translation distance to produce the potential energy curve.

The potential_energy_curve folder contains a set of useful scripts to be run on a compute canada cluster for automatically creating and submitting geometry relaxation calculations to the scheduler given a starting configuration. Due to symmetry considerations, it is only necessary to perform these translations along the diagonals of the cell. Since the cells I was dealing with were hexagonal, I had to create potential energy curves along both the long diagonal (`pec_long_diagonal.py`) and short diagonal (`pec_short_diagonal.py`). These can easily be modified to work with any other cell. Note: these scripts make use of input file makers (in the makers folder) from the [qecc](https://github.com/edmontoneuler/qecc) repo by [edmontoneuler](https://github.com/edmontoneuler). Some makers are unchanged while other are slightly modified to fit the task better.

The folder also contains scripts to automate resubmition of failed or expired jobs on the cluster (`resub.py`) and a scraper to extract the final energies (`scrape.py`).

## Potential Across an Interface

Another very useful measurement when characterizing the type of contact between two materials is the potential energy relative to the Fermi level at the the interface. Once you've calculated the optimized geometry of your cell, Quantum ESPRESSO allows you to easily calculate several different potentials across the cell with its post-processing tool `pp.x`. However, this returns 3D potentials V(x,y,x) which can be cumbersome to investigate. However, averageing V(x,y,z) in the direction perpendicular to the atomic layers (z-direction) and plotting it along said direction provides useful insight to the nature of the interface contact. This is precisely the purpose of `cube_average_plot.py`.
