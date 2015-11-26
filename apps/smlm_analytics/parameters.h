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
#include "sd_metropolis_parameters.h"
#include "free_hmm_genus_expansion.h"

extern int    mmax;
extern int    s1;            //sigma parameters controlling the thimble structure
extern int    s2s3;          //ratio of s2 to s3
extern double mean_link;    
extern double m2;            //Mass squared 

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();

#endif
