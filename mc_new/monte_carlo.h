#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "ran.h"

// Structure for simulation parameters
typedef struct Par {
  int nblock, nsamp, ntherm, seed;
  int L; // lentth of lattice in 1 direction
  double t; // temperature
  double U_I, U_II; // interaction of walls of type I and type II
  long nsteps; // number of steps per sweep
  int num_print_lattice; // number of configurations we want to print
  int print_measurement; // only write measurements after this many sweeps
  int bc_width; // mode for implementing boundary conditions. 0: no bc; >0: diameter of central, fixed hexagon
  int ctr; // ctr, for the case we run multiple instances with different temperatures
} Par;

// Structure for measured quantities
typedef struct Measurement {
  double energy;
  double polarization;
  double triangularity;
  double trig2; // For getting triangularity via edge cpuntint?
} Measurement;

// functions defined in common.c:
extern int pos_mod(int a, int b);
int metropolis(double delta, double temperature);
extern void print_lattice(Par *par, int *lattice);
void read_parameters(Par *par);
void read_cmd_parameters(Par *par, int argc, char *argv[]);
void write_lattice(Par *par, int *lattice, int step);
void write_measurement(Par *par, Measurement *measurement, long step);
extern void init_lattice(Par *par, int *lattice);
extern double energy_stance(Par *par, int *lattice, int i, int j);
extern void measure(Par *par, Measurement *measure, int *lattice);
void random_single_flip(Par *par, int *lattice);