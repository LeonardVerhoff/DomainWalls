#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ran.h"
#include "monte_carlo.h"

void mc_sweep(Par *par, int *lattice)
{
    // do LxL random flips -> one sweep
    int L = par->L;
    for (int sweep_step = 0; sweep_step < L; sweep_step++)
    {
        random_single_flip(par, lattice);
    }
}

void mc_run(Par *par, int *lattice)
{

    Measurement measurement;

    // measure intilized array
    measure(par, &measurement, lattice);
    write_measurement(par, &measurement, 0);

    for (long istep = 1; istep < par->nsteps + 1; istep++)
    {

        mc_sweep(par, lattice);

        // measure and write values after specific number of sweeps, i.e. LxL trials of spin flips
        if (istep % par->print_measurement == 0)
        {
            measure(par, &measurement, lattice);
            write_measurement(par, &measurement, istep);
        }

        long lattice_ctr = par->num_print_lattice * istep;
        if (lattice_ctr % par->nsteps == 0 && par->num_print_lattice > 0)
        {
            write_lattice(par, lattice, lattice_ctr / par->nsteps);
            printf("Wrote lattice of step %ld \n", istep);
        }
    }
    return;
}

int main(int argc, char *argv[])
{

    Par par;
    Measurement measurement;
    int *lattice;

    // itinialize with time:
    init_ran(0);

    par.L = 50;
    par.t = 1.0;
    par.nsteps = (int)1.0e3;
    par.U_I = -0.5;
    par.U_II = -0.5;
    par.num_print_lattice = 1;
    par.print_measurement = 10; // by default, write every 10th step
    par.bc_width = 0;
    par.ctr = 0;

    // read parameters from input.dat
    read_parameters(&par);
    // read command line arguments:
    if (argc > 1)
        read_cmd_parameters(&par, argc, argv);

    printf("Start. \n-------------------------------------------------------------------------------------------\n");
    printf("Parameters: \n      nsteps=%ld print_measurement=%d littices printed=%d \n      L=%d t=%f U_I=%lf U_II=%lf bc_width=%d\n", par.nsteps, par.print_measurement, par.num_print_lattice, par.L, par.t, par.U_I, par.U_II, par.bc_width);
    // maybe write all impportant info to log-file?

    // remove old outputs
    system("mkdir -p outputs");
    system("rm -f outputs/lattice*");

    lattice = malloc(par.L * par.L * sizeof(int));
    init_lattice(&par, lattice);
    measure(&par, &measurement, lattice);

    mc_run(&par, lattice);

    printf("-------------------------------------------------------------------------------------------\n Done.\n");
    return 0;
}
