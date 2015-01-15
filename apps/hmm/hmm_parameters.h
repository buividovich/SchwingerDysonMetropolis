#ifndef _HMM_PARAMETERS_H_
#define _HMM_PARAMETERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h> 

#include <clue_utils.h>
#include <clue_logs.h>

#include "sd_metropolis_parameters.h"
#include "largeN_QFT_parameters.h"


int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();

#endif
