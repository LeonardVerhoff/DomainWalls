#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ran.h"
#include "monte_carlo.h"

/*
general, mathematical functions
*/

int pos_mod(int a, int b)
{
    // return "positive" mod. If negative: what needs to get to next multiple of number!
    int res;
    res = a % b;
    if (res >= 0)
        return res;
    else
        return b + res;
}

int sign(int x)
{
    // return sign of an integer
    return (x > 0) - (x < 0);
}

int metropolis(double delta, double temperature)
{
    // decide wether to accept step. 1: success, 0: dont accept

    // if energy is lowered, we always accept
    if (delta <= 0)
        return 1;

    // otherwise:
    double ran = dran();
    double alpha = exp(-delta / temperature);

    if (ran < alpha)
        return 1;

    return 0;
}

/*
Functions for I/O
*/

void print_lattice(Par *par, int *lattice)
{
    // print lattice in terminal
    printf("\n");

    for (int i = 0; i < par->L; i++)
    {
        for (int j = 0; j < par->L; j++)
        {
            printf("%d ", lattice[i * par->L + j]);
        }
        printf("\n");
    }
}

void read_parameters(Par *par)
{
    // read input parameters from file "input.dat"
    FILE *file = fopen("input.dat", "r");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }

    char variable_name[100]; // Adjust the size as needed
    char equals_sign[2];     // To read the '=' character

    while (fscanf(file, "%s %1s", variable_name, equals_sign) == 2)
    {
        if (variable_name[0] == '#')
            continue;
        else if (equals_sign[0] == '=')
        {
            if (strcmp(variable_name, "L") == 0)
            {
                fscanf(file, "%d", &(par->L));
            }
            else if (strcmp(variable_name, "t") == 0)
            {
                fscanf(file, "%lf", &(par->t));
            }
            else if (strcmp(variable_name, "num_print_lattice") == 0)
            {
                fscanf(file, "%d", &(par->num_print_lattice));
            }
            else if (strcmp(variable_name, "print_measurement") == 0)
            {
                fscanf(file, "%d", &(par->print_measurement));
            }
            else if (strcmp(variable_name, "nsteps") == 0)
            {
                fscanf(file, "%ld", &(par->nsteps));
            }
            else if (strcmp(variable_name, "U_I") == 0)
            {
                fscanf(file, "%lf", &(par->U_I));
            }
            else if (strcmp(variable_name, "U_II") == 0)
            {
                fscanf(file, "%lf", &(par->U_II));
            }
            else if (strcmp(variable_name, "bc_width") == 0)
            {
                fscanf(file, "%d", &(par->bc_width));
            }
        }
    }
    fclose(file);
}

void read_cmd_parameters(Par *par, int argc, char *argv[])
{
    // read optional input from command line. Read: temperature (float) and ctr!
    double temperature;
    int ctr;

    if (argc != 3)
    {
        printf("Input error! Either provide no additional input or temperature and ctr! Using values from input.dat.");
        return;
    }

    temperature = atof(argv[1]);
    ctr = atoi(argv[2]);

    par->t = temperature;
    par->ctr = ctr;

    return;
}

void write_lattice(Par *par, int *lattice, int step)
{
    // printf("I was called with %d and %d \n", step, par->print_between);
    //  Define the size of your array
    int array_size = par->L; // Adjust this to your actual array size

    // Specify the output file name
    // const char *filename = "lattice.out";
    char filename[80];
    sprintf(filename, "outputs/lattice_%d.out", step);

    // Open the file for writing in binary mode
    FILE *file = fopen(filename, "wb");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }

    // write simulation parameters to header
    fprintf(file, "# nsteps=%ld istep=%d L=%d t=%f U_I=%lf U_II=%lf bc_width=%d \n", par->nsteps, step, par->L, par->t, par->U_I, par->U_II, par->bc_width);
    // Write the array to the file with L columns and L rows
    for (int i = 0; i < array_size; i++)
    {
        for (int j = 0; j < array_size; j++)
        {
            fprintf(file, "%d", sign(lattice[i * array_size + j])); // use sign, since we have +-2 for boundary conditions
            if (j < array_size - 1)
            {
                fprintf(file, "\t"); // Use a tab (or any other delimiter) between columns
            }
        }
        fprintf(file, "\n"); // Start a new row
    }
    // Close the file
    fclose(file);
}

void write_measurement(Par *par, Measurement *measurement, long step)
{
    char filename[80];
    sprintf(filename, "outputs/measure_%.2f_.out", par->t);

    // if first iteration, delete file if existent:
    if (step == 0)
    {
        char command[100];
        sprintf(command, "rm -f %s", filename);
        // system( "rm -f outputs/lattice*" );
        system(command);
    }

    FILE *file = fopen(filename, "a");
    if (file == NULL)
    {
        perror("Error opening file");
        return;
    }
    // printf("%ld \n", step);
    //  if first step, add file header
    if (step == 0)
        fprintf(file, "# Step   Energy   Polarization   Triangularity \n");

    // print measurements
    fprintf(file, "%ld      %lf      %lf      %lf \n", step, measurement->energy, measurement->polarization, measurement->triangularity);

    // Close the file
    fclose(file);
}

/*
Functions for the lattice
*/

void init_lattice(Par *par, int *lattice)
{
    int L = par->L;
    // initialize the lattice with random domains
    double ran;
    for (int i = 0; i < L * L; i++)
    {
        ran = dran_sign();
        if (ran > 0)
            lattice[i] = 1;
        else
            lattice[i] = 1;
    }

    // in case we apply fixed bundary conditions. Then we have +-2 to indicate it is fixed!
    if (par->bc_width > 0)
    {
        // Set 1st rows and 1st columns to +2 -> stripe
        for (int i = 0; i < L; i++)
        {
            lattice[i] = 2;
            lattice[i * L] = 2;
            lattice[i * L + 1] = 2;
        }
        // set impurity in the center:
        int center = (L / 2) * L + L / 2 - 1;
        int width = par->bc_width;

        lattice[center] = -2;
        for (int r1 = 0; r1 < width; r1++)
        {

            for (int r2 = 0; r2 < width * 2; r2++)
            {
                lattice[center + r2 + L * r1] = -2;
            }
        }
    }
}

/*
Functions for energy calculation and dynamics
*/

double energy_stance(Par *par, int *lattice, int i, int j)
{
    // calculate energy of a single stance
    int L = par->L;
    double energy = 0;

    int spin_ij = sign(lattice[i * L + j]);
    int spin_k;
    int neighbors[3];

    /* even newer idea: put it on triangular lattice.
    2 Sublattices of triangulars pointing up and down.
    Up triangles (even numbers): I interaction for up-down and II for down-up
    down triangles (odd numbers): II interaction for up-down and I for down-up
    */
    int side; // 1: triangle pointing downwards; 0: triangle poining to upwards
    neighbors[0] = i * L + pos_mod(j + 1, L);
    neighbors[1] = i * L + pos_mod(j - 1, L);
    if ((i * L + j) % 2) // case of triangle pointing downwards -> 3rd interaction is with upper triangle
    {
        neighbors[2] = pos_mod(i + 1, L) * L + j - 1;
        side = 1;
    }
    else
    { // case of triangle pointing upwards -> 3rd interaction is with lower triangle
        neighbors[2] = pos_mod(i - 1, L) * L + j + 1;
        side = 0;
    }

    for (int k = 0; k < 3; k++)
    {
        // give energy penalty only if we have to opposite spins, otherwise energy is 0
        spin_k = sign(lattice[neighbors[k]]);

        if (spin_ij != spin_k)
        {                      // in case of different domains:
            if (spin_ij == +1) // case of up->down: I-type interaction for upwards; II-type interaction for downwards triangle
                energy += ((par->U_I) * (1 - side) + (par->U_II) * side) * spin_ij * spin_k;

            if (spin_ij == -1) // case of down->up: II-type interaction for upwards triangle; I-type interaction for downwards triangle
                energy += ((par->U_I) * (side) + (par->U_II) * (1 - side)) * spin_ij * spin_k;
        }
    }
    return energy;
}

void measure(Par *par, Measurement *measurement, int *lattice)
{
    // total energy as sum over all stances
    double tot_en = 0, p_up = 0, p_down = 0;
    // int edge_ctr = 0;
    int L = par->L;
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            tot_en += energy_stance(par, lattice, i, j);

            // if j even, we have an up triangle, if j odd its down
            if (j % 2 == 0)
                p_up += sign(lattice[i * L + j]);

            // count all edges of
            // implement edge counter???

            else
                p_down += sign(lattice[i * L + j]);
        }
    }

    measurement->energy = tot_en;
    measurement->polarization = (p_up + p_down) / (L * L);
    measurement->triangularity = (p_up - p_down) / (L * L);
}

void random_single_flip(Par *par, int *lattice)
{
    // choose random domain and flip. Accept with prop. exp(-dE/t).
    int ran_i, ran_j = 0;
    double old_en, new_en, delta_en;

    ran_i = par->L * dran();
    ran_j = par->L * dran();

    // in case we pick a cell that shall be fixed:
    while (abs(lattice[ran_i * par->L + ran_j]) > 1)
    {
        ran_i = par->L * dran();
        ran_j = par->L * dran();
    }

    // printf("%d ", ran_pos);
    old_en = energy_stance(par, lattice, ran_i, ran_j);
    lattice[ran_i * par->L + ran_j] *= -1;
    new_en = energy_stance(par, lattice, ran_i, ran_j);

    delta_en = new_en - old_en;

    int do_flip = metropolis(delta_en, par->t);

    // undo spin flip if we dont accept
    if (!do_flip)
    {
        lattice[ran_i * par->L + ran_j] *= -1;
    }

    return;
}