#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h> 

#include <clue_utils.h>
#include <clue_logs.h>

extern char data_suffix[512];   //Whether to save the histogram of sampling in sectors of different n and order
extern char scan_suffix[512];
extern char scan_label[512];
extern int  save_summary;

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();

#endif
