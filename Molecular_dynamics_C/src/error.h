#ifndef ERROR_H
#define ERROR_H

#include <stdio.h>
#include <stdlib.h>

// Function to check the error in the input file opening.
void error_file_open(FILE* input_file);

// Function to check the error in the memory alliocation
void error_memory_allocation(const char* var_name);

// Function to check the error in reading the number of atom from the input file(input_file)
void error_read_atoms();

// Function to check the error in the reading coordinates and mass from the input file(input_file)
void error_read_molecule();

#endif

