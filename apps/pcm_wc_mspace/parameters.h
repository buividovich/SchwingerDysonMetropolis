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

//Some switches
extern int    save_sampling_hist;   //Whether to save the histogram of sampling in sectors of different n and order
extern int    save_correlators;     //Whether to save correlators <gx gy>
extern int    save_metropolis_stat; //Whether to save the data on the performance of Metropolis algorithm

//Useful calculable parameters
extern double stereo_alpha;     //-\lambda/8, expansion parameter 

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();

#endif
