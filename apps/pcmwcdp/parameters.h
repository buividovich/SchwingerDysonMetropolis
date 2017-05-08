#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h> 

#include <clue_logs.h>

#include "largeN_QFT_parameters.h"

extern char   raw_data_dir[512];
extern char exact_data_dir[512];

extern char** cluster_data_suffixes; 
extern int    calculate_correlators;
extern int    number_of_clusters;

extern char data_suffix[512];   //Whether to save the histogram of sampling in sectors of different n and order
extern char scan_suffix[512];
extern char scan_label[512];
extern int  save_summary;

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();

#endif
