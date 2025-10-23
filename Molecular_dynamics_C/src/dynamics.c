#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "utils.h"
#include "error.h"

// --------------------------------------------------------------------------------------------- //
// ********************************** THE MAIN PROGRAM ***************************************** //
// --------------------------------------------------------------------------------------------- //

// Constants for Lennard-Jones potential
const double epsilon = 0.0661;    // Lennard-Jones epsilon in J/mol
const double sigma   = 0.3345;    // Lennard-Jones sigma in nm
const int WRITE_FREQUENCY= 1;     // Frequency of writing trajectory to file

int main()
{
    // Open input file and read the number of atoms
    const char* input_filename = "data/CH4.txt";
    FILE* input_file = fopen(input_filename, "r");
    if (input_file == NULL) error_file_open(input_file);
    
    // Dynamically create the output file name
    char output_file[100]; 						// Assuming max file name length is 100
    strncpy(output_file, input_filename, sizeof(output_file));
    char* dot = strrchr(output_file, '.'); 				// Find the last dot in the file name
    if (dot != NULL) 
	    *dot = '\0'; 					// Remove the extension
    strcat(output_file, ".xyz"); 					// Append ".xyz"
				     

    size_t Natoms = read_Natoms(input_file);
    if (Natoms == 0) error_read_atoms();

    // Allocate memory for coordinates, masses, symbols, and other arrays
    double** coord = malloc_2d(Natoms, 3);
    double* mass = malloc(Natoms * sizeof(double));
    char** symbols = malloc(Natoms * sizeof(char*));
    double** distance = malloc_2d(Natoms, Natoms);
    double** velocity = malloc_2d(Natoms, 3);
    double** acceleration = malloc_2d(Natoms, 3);

    if (coord == NULL || mass == NULL || distance == NULL || velocity == NULL || acceleration == NULL) {
        error_memory_allocation("Allocation error");
        return EXIT_FAILURE;
    }

    // Allocate memory for symbols and read molecule data
    for (size_t i = 0; i < Natoms; i++) 
    {
        symbols[i] = malloc(10 * sizeof(char)); // Assume max symbol length is 10
        if (symbols[i] == NULL) 
	{
            error_memory_allocation("symbols");
            return EXIT_FAILURE;
        }
    }

    if (!read_molecule(input_file, Natoms, coord, mass, symbols)) 
	    error_read_atoms();
    fclose(input_file);

    // Initialize velocities to zero
    for (size_t i = 0; i < Natoms; i++)
        for (size_t j = 0; j < 3; j++)
            velocity[i][j] = 0.0;

    // Compute initial distances and accelerations
    compute_distances(Natoms, coord, distance);
    compute_acc(Natoms, coord, mass, distance, acceleration, sigma, epsilon);

    // Open trajectory file
    FILE* trajectory_file = fopen(output_file, "w");
    if (trajectory_file == NULL) 
    {
        printf("Error opening trajectory file");
        return EXIT_FAILURE;
    }

    // Simulation parameters
    double dt = 0.2;               // Time step
    size_t total_steps = 1000;    // Total number of simulation step
    size_t progress_interval = 50;
    printf("Starting molecular dynamics simulation.........\n"); // Start message

    // Molecular dynamics simulation loop
    for (size_t step = 0; step < total_steps; step++) 
    {
        // Compute kinetic, potential, and total energies
        double kinetic   = kinetic_energy(Natoms, velocity, mass);
        double potential = potential_energy(epsilon, sigma, Natoms, distance);
        double total     = Total_energy(potential, kinetic);

        // Write trajectory every WRITE_FREQUENCY steps
        if (step % WRITE_FREQUENCY == 0) 
	{
            write_trajectory(trajectory_file, Natoms, coord, symbols, kinetic, potential, total, step);
            
	    // Update progress bar
        	if ((step + 1) % progress_interval == 0) 
		{
            		printf("#"); // Print one `#` for every 50 steps
            		fflush(stdout); // Ensure output is printed immediately
        	}
        }

        // Update positions, velocities, and accelerations using Verlet algorithm
        verlet_update(Natoms, dt, coord, velocity, acceleration, distance, mass, sigma, epsilon);
    }
    printf("\n\n");
    printf("Molecular dynamics simulation completed successfully.\n"); // End message

    // Close trajectory file and free allocated memory
    fclose(trajectory_file);
    free_2d(coord); 
    free_2d(distance); 
    free_2d(velocity); 
    free_2d(acceleration); 
    free(mass);
    for (size_t i = 0; i < Natoms; i++) free(symbols[i]);
    free(symbols);

    return 0;
}

