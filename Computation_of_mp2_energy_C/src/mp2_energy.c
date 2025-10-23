#include <stdio.h>
#include <trexio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

//--------------------------------------------------------------------------------//
//                                 MAIN PROGRAM                                   //
//--------------------------------------------------------------------------------//

int main() 
{
    trexio_exit_code rc;
    double energy;
    int32_t n_up;
    int32_t mo_num;
    int64_t n_integrals;

    //--------------------------------------------------------------------------------//
    //                                 OPENING FILE                                   //
    //--------------------------------------------------------------------------------//
    
    const char* input_filename = "data/hcn.h5";
    // Check if the filename has the .h5 extension
    const char *extension = strrchr(input_filename, '.');
    if (extension == NULL || strcmp(extension, ".h5") != 0) 
    {
        printf("Error: Input file must have a .h5 extension.\n");
        return -1;
    }

    trexio_t* file = trexio_open(input_filename, 'r', TREXIO_AUTO, &rc);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error: %s\n", trexio_string_of_error(rc));
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                        READING NUCLEAR REPULSION ENERGY                        //
    //--------------------------------------------------------------------------------//
    printf("READING THE DATA in progess . . .");
    rc = trexio_read_nucleus_repulsion(file, &energy);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error reading nuclear repulsion energy:\n%s\n", trexio_string_of_error(rc));
        trexio_close(file);
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                        OBTAINING THE NUMBER OF OCCUPIED ORBITALS               //
    //--------------------------------------------------------------------------------//

    rc = trexio_read_electron_up_num(file, &n_up);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error reading number of up-spin electrons:\n%s\n", trexio_string_of_error(rc));
        trexio_close(file);
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                          READING THE NUMBER OF MOLECULAR ORBITALS              //
    //--------------------------------------------------------------------------------//

    rc = trexio_read_mo_num(file, &mo_num);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error reading number of molecular orbitals:\n%s\n", trexio_string_of_error(rc));
        trexio_close(file);
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                          READING ONE-ELECTRON INTEGRALS                        //
    //--------------------------------------------------------------------------------//

    double* data = malloc(mo_num * mo_num * sizeof(double));
    if (data == NULL) 
    {
        printf("Memory allocation failed for one-electron integrals!\n");
        trexio_close(file);
        exit(1);
    }

    rc = trexio_read_mo_1e_int_core_hamiltonian(file, data);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error reading one-electron integrals:\n%s\n", trexio_string_of_error(rc));
        free(data);
        trexio_close(file);
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                          READING TWO-ELECTRON INTEGRALS SIZE                   //
    //--------------------------------------------------------------------------------//

    rc = trexio_read_mo_2e_int_eri_size(file, &n_integrals);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error reading two-electron integrals size:\n%s\n", trexio_string_of_error(rc));
        free(data);
        trexio_close(file);
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                          READING TWO-ELECTRON INTEGRALS                        //
    //--------------------------------------------------------------------------------//

    int32_t* const index = malloc(4 * n_integrals * sizeof(int32_t));
    double* const value = malloc(n_integrals * sizeof(double));
    if (index == NULL || value == NULL) 
    {
        printf("Memory allocation failed for two-electron integrals!\n");
        free(data);
        if (index) free(index);
        if (value) free(value);
        trexio_close(file);
        exit(1);
    }

    rc = trexio_read_mo_2e_int_eri(file, 0, &n_integrals, index, value);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error reading two-electron integrals:\n%s\n", trexio_string_of_error(rc));
        free(data);
        free(index);
        free(value);
        trexio_close(file);
        exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                         READING THE INTEGRAL VALUES                            //
    //--------------------------------------------------------------------------------//

    double ****integral_array = malloc_4d(mo_num, mo_num, mo_num, mo_num);

    for (int n = 0; n < n_integrals; n++) 
    {
        int i = index[4 * n + 0];
        int j = index[4 * n + 1];
        int k = index[4 * n + 2];
        int l = index[4 * n + 3];

        // The integral values considering the 8-fold symmetry
        integral_array[i][j][k][l] = value[n];
        integral_array[k][l][i][j] = value[n];
        integral_array[k][j][i][l] = value[n];
        integral_array[i][l][k][j] = value[n];
        integral_array[j][i][l][k] = value[n];
        integral_array[l][k][j][i] = value[n];
        integral_array[j][k][l][i] = value[n];
        integral_array[l][i][j][k] = value[n];
    }

    //--------------------------------------------------------------------------------//
    //                         READING THE ORBITAL ENERGIES                           //
    //--------------------------------------------------------------------------------//

    double* mo_energy = malloc(mo_num * sizeof(double));
    if (mo_energy == NULL)
    {
        printf("Memory allocation failed for orbital energies!\n");
        free(data);
        free(index);
        free(value);
        trexio_close(file);
	return -1;
    }

    rc = trexio_read_mo_energy(file, mo_energy);
    if (rc != TREXIO_SUCCESS)
    {
        printf("TREXIO Error reading molecular orbital energies:\n%s\n",trexio_string_of_error(rc));
	exit(1);
    }

    //--------------------------------------------------------------------------------//
    //                          CALCULATING THE HARTREE FOCK ENERGY                   //
    //--------------------------------------------------------------------------------//

    double hf_energy = HF_energy(energy, data, integral_array, mo_num, n_up);

    //--------------------------------------------------------------------------------//
    //                          CALCULATING THE MP2 ENERGY                            //
    //--------------------------------------------------------------------------------//

    double mp2_energy = calculate_MP2_energy(integral_array, mo_energy, n_up, mo_num);
 
    //--------------------------------------------------------------------------------//
    //                          WRITING THE OUTPUT FILE                               //
    //--------------------------------------------------------------------------------//
    
    // Dynamically create the output file name
    char output_file[100];                                              // Assuming max file name length is 100
    strncpy(output_file, input_filename, sizeof(output_file));
    char* dot = strrchr(output_file, '.');                              // Find the last dot in the file name
    if (dot != NULL)
            *dot = '\0';                                                // Remove the extension
    strcat(output_file, ".txt");                                        // Append ".txt"
									
    // Create output file
    create_output_file(output_file, input_filename, energy, n_up, mo_num, hf_energy, mp2_energy);
    printf("\nCalculation completed \nResults written to %s\n", output_file);

    
    //--------------------------------------------------------------------------------//
    //                                 CLEANING UP                                    //
    //--------------------------------------------------------------------------------//

    free(data);
    free(index);
    free(value);
    free_4d(integral_array, mo_num, mo_num, mo_num);

    rc = trexio_close(file);
    if (rc != TREXIO_SUCCESS) 
    {
        printf("TREXIO Error closing file: %s\n", trexio_string_of_error(rc));
        return -1;
    }

    return 0;
}

