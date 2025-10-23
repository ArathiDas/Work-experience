#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ---------------------------------------------------------------------------------------------//
//                      TO ALLOCATE AND FREE THE MEMORY FOR 2-D ARRAY AND MASS                  //
// ---------------------------------------------------------------------------------------------//

double** malloc_2d(size_t m, size_t n)
{
    // Allocate an array of double pointers with size m
    double** a = malloc(m * sizeof(double*));
    if (a == NULL)
    {
        return NULL;
    }

    // Allocate a contiguous block of memory for the 2D array elements
    a[0] = malloc(n * m * sizeof(double));
    if (a[0] == NULL)
    {
        free(a);
        return NULL;
    }

    // Set the pointers in the array of double pointers to point to the correct locations
    for (size_t i = 1; i < m; i++)
    {
        a[i] = a[i - 1] + n;
    }

    return a;
}

void free_2d(double** a)
{
    free(a[0]);
    a[0] = NULL;
    free(a);
}



// ---------------------------------------------------------------------------------------------//
// 			TO READ THE NUMBER OF ATOMS FROM THE INP.TXT FILE                       // 
// ---------------------------------------------------------------------------------------------//

size_t read_Natoms(FILE* input_file)
{
    size_t Natoms;
	
    // Read the number of atoms from the first line of the file
    if(fscanf(input_file, "%zu", &Natoms) != 1)
    {
	    return 0;
    }
    
    return Natoms;
}



// ---------------------------------------------------------------------------------------------//
//               TO READ THE SYMBOL, COORDINATES AND MASS OF ATOMS FROM THE INP.TXT FILE        //
// ---------------------------------------------------------------------------------------------//

int read_molecule(FILE* input_file, size_t Natoms, double** coord, double* mass, char** symbols)
{
    
    // Read the atomic coordinates and masses for each atom
    for (size_t i = 0; i < Natoms; i++)
    {
        // Read x, y, z coordinates and mass
        if(fscanf(input_file, "%s %lf %lf %lf %lf", symbols[i], &coord[i][0], &coord[i][1], &coord[i][2], &mass[i]) != 5)
	{
		return 0;
	}
    }

    return 1;
}

// ---------------------------------------------------------------------------------------------//
//               		TO CALCULATE THE DISTANCE BETWEEN ATOMS                		//		
// ---------------------------------------------------------------------------------------------//

void compute_distances(size_t Natoms, double** coord, double** distance) 
{
    for (size_t i = 0; i < Natoms; i++) 
    {
        for (size_t j = 0; j < Natoms; j++) 
	{
            if (i == j) 
	    {
                distance[i][j] = 0.0;  // Distance from atom to itself
            } 
	    else 
	    {
                // The distances are multiplied by 0.1 unit to make the "nm" scale
		double dx = (coord[i][0] - coord[j][0]);
                double dy = (coord[i][1] - coord[j][1]);
                double dz = (coord[i][2] - coord[j][2]);
                
                // Euclidean distance
                distance[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------//
//                              TO CALCULATE THE POTENTIAL ENERGY                               //
// ---------------------------------------------------------------------------------------------//

double potential_energy(double epsilon, double sigma, size_t Natoms, double** distance)
{
	double V_LJ = 0; 				
	double V_total = 0; 				

	// Loop over all unique pairs of atoms
	for (int i=0; i<Natoms; i++)
	{
		for (int j=i+1; j<Natoms; j++) 			// Only consider j > i to avoid double counting
		{
			double r = distance[i][j]; 		// Distance between atom i and atom j
			double sigma_over_r_value;  		// Ratio of sigma to r
			double power_12_term_value; 		// Term for the 1/r^12 contribution
			double power_6_term_value;  		// Term for the 1/r^6 contribution
			
			if (r > 0) // Ensure distance is positive to avoid division by zero
			{
				sigma_over_r_value  = sigma/r;  				
				power_6_term_value  = pow(sigma_over_r_value,6);  		// (sigma/r)^6 term
				power_12_term_value = power_6_term_value * power_6_term_value; 	// (sigma/r)^12 term

				// Calculate Lennard-Jones potential for this pair
				V_LJ = 4*epsilon*(power_12_term_value - power_6_term_value);
			}

			V_total += V_LJ;			 // Accumulate the potential energy for this pair 
		}
	}
	return V_total; 					 
}


// ---------------------------------------------------------------------------------------------//
//                              TO CALCULATE THE KINETIC ENERGY                                 //
// ---------------------------------------------------------------------------------------------//

double kinetic_energy(size_t Natoms, double** velocity, double* mass)        // function to calculate the total kinetic energy  
{
	double T_total = 0;

	for (size_t i = 0; i < Natoms; i++)                                  // iterating over each atom
	{
		double vx = velocity[i][0];                                  // extracting the x, y and z components from the array
		double vy = velocity[i][1];
		double vz = velocity[i][2];

		T_total += 0.5 *mass[i]*(vx*vx + vy*vy + vz*vz);      	     // calculating the kinetic energy
	}
	return T_total;
}


// ---------------------------------------------------------------------------------------------//
//                              TO CALCULATE THE TOTAL ENERGY                                   //
// ---------------------------------------------------------------------------------------------//

double Total_energy(double V, double T)
{
	return T+V;
}


// ---------------------------------------------------------------------------------------------//
//                              COMPUTING THE ACCELERATION                                      //
// ---------------------------------------------------------------------------------------------//
void compute_acc(size_t Natoms, double** coord, double* mass, double** distance, double** acceleration, double sigma, double epsilon) 
{
    // Reset accelerations to zero
    for (size_t i = 0; i < Natoms; i++) 
    {
        acceleration[i][0] = 0.0;
        acceleration[i][1] = 0.0;
        acceleration[i][2] = 0.0;
    }

    double r_min = 0.1; // Minimum allowed distance to avoid division by zero

    // Compute accelerations based on pairwise interactions
    for (size_t i = 0; i < Natoms; i++) 
    {
        for (size_t j = 0; j < Natoms; j++) 
	{
	   if(i != j)
	   { 
            double r = distance[i][j];

            // Apply minimum distance threshold
            if (r < r_min) 
	    {
                r = r_min;
            }

            // Compute Lennard-Jones force magnitude
            double sigma_over_r = sigma / r;
            double t6 = pow(sigma_over_r, 6);  // (sigma / r)^6
            double t12 = t6 * t6;              // (sigma / r)^12
            double force = (24.0 * epsilon / r) * (t6 - 2.0 * t12);

            // Compute force components
            double fx = force * (coord[i][0] - coord[j][0]) / r;
            double fy = force * (coord[i][1] - coord[j][1]) / r;
            double fz = force * (coord[i][2] - coord[j][2]) / r;

            // Update acceleration for atom i
            acceleration[i][0] += (-1.0 / mass[i]) * fx;
            acceleration[i][1] += (-1.0 / mass[i]) * fy;
            acceleration[i][2] += (-1.0 / mass[i]) * fz;
	   }
       }
    }
}

// ---------------------------------------------------------------------------------------------//
//                                  VERLET AlGORITHM                                            //
// ---------------------------------------------------------------------------------------------//


void verlet_update(size_t Natoms, double dt, double** coord, double** velocity, double** acceleration, double** distance, double* mass, double sigma, double epsilon)
{
    // Updating the positions of the atoms
    for (size_t i = 0; i < Natoms; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            coord[i][j] += velocity[i][j] * dt + 0.5 * acceleration[i][j] * dt * dt;
        }
    }

    // Computing the distances between the atoms
    compute_distances(Natoms, coord, distance);

    // Computing  accelerations
    double** new_acceleration = malloc_2d(Natoms, 3); // Allocate temporary array for new accelerations
    if (!new_acceleration)
    {
        perror("Error allocating memory for new_acceleration");
        exit(EXIT_FAILURE);
    }
    compute_acc(Natoms, coord, mass, distance, new_acceleration, sigma, epsilon);

    // Updating the velocity vectors
    for (size_t i = 0; i < Natoms; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            velocity[i][j] += 0.5 * (acceleration[i][j] + new_acceleration[i][j]) * dt;
            acceleration[i][j] = new_acceleration[i][j]; // Update old acceleration to the new one
        }
    }

    free_2d(new_acceleration); // Free temporary acceleration array
}


// ---------------------------------------------------------------------------------------------//
//   				THE FILE WRITING FUNCTION                                       //
// ---------------------------------------------------------------------------------------------//
void write_trajectory(FILE* trajectory_file, size_t Natoms, double** coord, char** symbols, double kinetic_energy, double potential_energy, double total_energy, size_t step)
{	
    // The output details
 
    // Coordinates in XYZ format
    fprintf(trajectory_file, "%zu\n", Natoms); // Number of atoms
    fprintf(trajectory_file, "Step: %zu --- Kinetic Energy: %.8f J/mol --- Potential Energy: %.8f J/mol --- Total Energy: %.8f J/mol\n", step, kinetic_energy, potential_energy, total_energy);
    for (size_t i = 0; i < Natoms; i++)
    {
        fprintf(trajectory_file, "%-2s %10.5f %10.5f %10.5f\n", symbols[i], coord[i][0], coord[i][1], coord[i][2]);
    }
}

