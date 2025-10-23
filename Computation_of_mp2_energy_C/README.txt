## COMPUTATION OF THE MP2 ENERGY
This project computes the Hatree-Fock energy and MP2 energy for a closed shell system using input data from a trexio file. 
The program reads molecular orbitals, orbital_energies and other parameters to perform the calculations. 

## Directory structure
- INSTALL_MP2.pdf: Provides instructions on how to compile and run the program
- Makefile: handles the compilation process for the source files
- data/: contains the input files for the program for methane, water and benzene
- src/: Source code for the simulation 
     - mp2_energy.c: implements the core dynamics
     - utils.c: contains utility functions for memory allocation, reading inputs, defining functions, etc.
     - utils.h: header file for utility functions declarations
- tests/: contains the output files from the test runs
