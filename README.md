# qe-interfaces

This repo contains python tools that are useful in understanding interfaces in heterostructures consisting of two or more materials within a density-functional theory framework. These tools are meant to be used in conjunction with Quantum ESPRESSO input and output files. I developed them during my MSc when I was studying Ca2N-MoS2 and Au-Ca2N-MoS2 heterostrucures, but they can be modified for use with other supercells. 

Check out the resulting paper: [Improved Charge Transfer and Barrier Lowering across a Au–MoS2 Interface through Insertion of a Layered Ca2N Electride](https://doi-org.ezproxy.library.dal.ca/10.1021/acs.jpcc.1c02142). It will come in handy in explaining what these tools are useful for.


## Bandstructure Plots with Atomic Contributions

Bandstructure plots can tell us a lot about materials; however, they can be very difficult to understand when the cell is consists of multiple materials interfaced togethor like in a heterostructure supercell. `band_plot_contribution_2_materials.py` and `band_plot_contribution_3_materials.py` get around the spaghetti mess of these kind of bandstructures by assigning a colored weight to each state in order to distinguish its composition as a function of contribution from the indivudual materials. They are used for 2 and 3 material heterostructures, respectively.

`band_plot_contribution_2_materials.py` and `band_plot_contribution_3_materials.py` require `pw.x` input and output files for `scf` calculations. As well, you will need to provide the file containing the projected bandstructure states. You can use Quantum ESPRESSO's `projwfc.x` to get the states projected onto localized atomic orbitals (`pdos.out` in the input_file_examples folder). 

## Potential Energy Curves

## Other
