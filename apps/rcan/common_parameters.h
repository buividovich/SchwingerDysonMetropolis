#ifndef _COMMON_PARAMETERS_H_
#define _COMMON_PARAMETERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h> 

#include <clue_utils.h>
#include <clue_logs.h>

#include "sd_metropolis_parameters.h"
#include "sd_metropolis_statistics.h"
#include "largeN_QFT_parameters.h"
#include "sd_metropolis_parameters.h"

#include "recursion.h"

#define  FIX2        //This parameter makes LS=2 a hard-coded choice
#define   PACKED      //This switches on an economical saving of momenta sequences

#define  MAX_M       (16)             //Maximal orders - to define storage of pre-computed values as static array
#define  MAX_MP1     (17)             //MAX_M+1
#define  MAX_MP2     (18)             //MAX_M+2
#define  MAX_MX2     (32)             //Maximal orders - to define storage of pre-computed values as static array
#define  MAX_LS      (16)             //Maximal LS  - for static arrays
#define  MAX_LS2     (256)            //Maximal LS2 - for static arrays containing bare propagators etc.

extern int    mmax;          //Max. order
extern int    mmin_prc;
extern int    consistency_check;

extern int    n_threads;

extern char*  out_file;
extern int    append_mode;
extern int    auto_naming;

int  parse_common_command_line_options(int argc, char **argv);
void init_common_parameters();
void print_common_parameters();

#endif
