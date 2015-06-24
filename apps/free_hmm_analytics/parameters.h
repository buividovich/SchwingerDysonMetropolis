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
#include "free_hmm_genus_expansion.h"

extern int   nphi;
extern char* nphi_str;
extern int   ntr;
extern int*  ngs;
extern int   npairs;

int parse_nphi_str(char* s, int** ns);
int parse_command_line_options(int argc, char **argv);

#endif
