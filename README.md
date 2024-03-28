# DomainWalls

This repo includes scripts and programs I wrote for the simulation of **ferroelectric domain walls** (DW) in lithium niobate (LN) and lithium tantalate (LT). The work was part of my Master's thesis at [Justus-Lieibig-University Gießen](https://www.uni-giessen.de/de/fbz/fb07/fachgebiete/physik/institute/theorie/agsanna/people). 

&nbsp;
>The work consists of a multi-scale approach and is devided in three major parts: 
1. Phenomenological description according to Ginzburg-Landau-Devonshire-Theory (GLD), following [Scrymgeour et al](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.71.184110)
   
   *followed by*
   
2. The density functional theory (DFT) modeling employing large supercells and the *Vienna ab initio simulation package* [(VASP)](https://www.vasp.at)

   *which gives input to*
   
2. A derivate of the Ising model, to simulate the formation of domains in thin films on larger scales

&nbsp;

The folder `gld` contains a jupyter notebook to numerically solve the equations given by GLD formalism (1.).

&nbsp;

The folder `dft_scripts` contains 2 important python scripts for the DFT calculations (2.): The first one (`create_supercell.py`) builds the supercell containing (DW).
The second script (`get_pol.py`) analyses the DFT-relaxed cell and calculates the ionic contribution to spontaneous polarization in an electrostatic definition.

&nbsp;

The folder `mc_new` contains C programs, implementing the Ising model for the simulaiton of thin films (3.). It contains a Makefile (that worked for me at least), as well as a sample input file (`input.dat`). There is also a shell script (`run_multiple.sh`), that allows to run multiple instances of the simulation for different temperatures, in order to have kind of an embarrising *parallelization*.

&nbsp;

