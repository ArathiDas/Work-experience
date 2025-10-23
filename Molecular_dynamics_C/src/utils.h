
#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

// Function to allocate the memory
double** malloc_2d(size_t m, size_t n);

// Functions to free the allocated memory
void free_2d(double** a);

// Function to read the number of atoms from the input file (inp.txt)
size_t read_Natoms(FILE* input_file);

// Functions to read the coordinate and mass of the molecule from the inpiut file (inp.txt)
int read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass, char** symbols);

// Function to compute the internuclear distance between pairs of atoms
void compute_distances(size_t Natoms, double** coord, double** distance);

// Function to compute the potential energy 
double potential_energy(double epsilon, double sigma, size_t Natoms, double** distance);

// Function to compute the kinetic energy
double kinetic_energy(size_t Natoms, double** velocity, double* mass);

//Fucntion to compute the total energy of the system by summing T and V
double Total_energy( double V, double T);

// Function compute to the acceleration of the atoms
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration, double sigma, double epsilon);

// The verlet algorithm 
void verlet_update(size_t Natoms, double dt, double** coord, double** velocity, double** acceleration, double** distance, double* mass, double sigma, double epsilon);

// The file wrting function
void write_trajectory(FILE* trajectory_file, size_t Natoms, double** coord, char** symbols, double kinetic_energy, double potential_energy, double total_energy, size_t step);
#endif

