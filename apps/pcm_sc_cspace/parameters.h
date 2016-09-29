#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h> 

#include <clue_utils.h>
#include <clue_logs.h>

#include "sd_metropolis_parameters.h"
#include "sd_metropolis_statistics.h"
#include "largeN_QFT_parameters.h"
#include "lattice_stack.h"
#include "actions.h"

extern double meff_sq; //Effective mass term - for resummations 
//Some switches
extern int    save_sampling_hist; //Whether to save the histogram of sampling in sectors of different n and order
extern int    save_correlators;   //Whether to save correlators <gx gy>
extern int    resummation;        //Whether to re-sum the kinetic "random walk" term
//Useful calculable parameters
extern double beta;           // == 1/\lambda, strong-coupling expansion parameter
extern double sigma;          // Self-energy, is trivial without resummation
extern double mass2;          //Squared mass determining the pole of the bare propagator
extern double coord_factor;   //Counts how many neighbors of a given site coincide

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();

#endif
