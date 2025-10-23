#include <stdio.h>
#include <stdlib.h>

// File opening error
void error_file_open(const char* filename)
{
    printf("Error: Could not open the file: inp.txt \n");
    exit(EXIT_FAILURE);  // Exit the program
}

// Memory allocation error
void error_memory_allocation(const char* var_name)
{
    printf("Error: Memory allocation failed for %s\n", var_name);
    exit(EXIT_FAILURE);  // Exit the program
}

// Atom reading error
void error_read_atoms()
{
    printf("Error: Failed to read the number of atoms\n");
    exit(EXIT_FAILURE);  // Exit the program
}

// Molecule data reading error
void error_read_molecule()
{
    printf("Error: Failed to read data for atoms \n");
    exit(EXIT_FAILURE);  // Exit the program
}

