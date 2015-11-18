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

extern double mean_link;
extern double alpha;
extern int    max_alpha_order;

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void free_parameters();
void print_parameters();

#endif
