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

#define NMODELS (3)
//Model numbering
#define FREE_GAUSS  (0)
#define LM_INVERSE  (1)
#define LM_POSITIVE (2)

static const char * const model_names[NMODELS] = {"Free Gaussian", "LM with inverses", "LM with positive powers"};

extern double epsilon;
extern int    model;
//Calculable parameters
extern double   create_amplitudes[NMODELS];
extern double increase_amplitudes[NMODELS];
extern double     join_amplitudes[NMODELS];


int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void free_parameters();
void print_parameters();

#endif
