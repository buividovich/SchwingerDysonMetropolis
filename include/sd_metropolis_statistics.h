#ifndef _SD_METROPOLIS_STATISTICS_H_
#define _SD_METROPOLIS_STATISTICS_H_

#include <stdlib.h>
#include <stdio.h>
#include <clue_utils.h>
#include <clue_logs.h>

#include "sd_metropolis.h"

extern double                   maxnA; //Max. value of nA
extern double                     anA; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
extern double                     dnA; //Expectation value of its square - necessary to estimate the error!
extern double                   msign; //Expectation value of the sign - also necessary for recovering the correct normalization
extern int                        aac;
extern double                     ans;
extern int                        nmc; //Counter of calls to gather_mc_stat()
extern int*            action_counter;

//If control_max_ampl_sum=1 and max_ampl_sum is set to some nonzero value, 
//an error message is generated everytime nA exceeds max_ampl_sum, 
//this is useful for debugging
extern int       control_max_ampl_sum; 
extern double            max_ampl_sum;
extern double        max_ampl_sum_tol;


//These are the variables which are only set by process_mc_stat
extern double acceptance_rate;
extern double mean_recursion_depth;
extern double mean_nA;
extern double err_nA;
extern double mean_sign;

void init_metropolis_statistics();
void gather_mc_stat();
void process_mc_stat(const char* prefix);

void print_max_amplitudes();

#endif
