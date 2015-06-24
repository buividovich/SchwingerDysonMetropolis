#ifndef _SD_METROPOLIS_STATISTICS_H_
#define _SD_METROPOLIS_STATISTICS_H_

#include <stdlib.h>
#include <stdio.h>
#include <clue_utils.h>
#include <clue_logs.h>
#include <clue_io.h>

#include "sd_metropolis.h"

extern double                   maxnA; //Max. value of nA
extern double                     anA; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
extern double                     dnA; //Expectation value of its square - necessary to estimate the error!
extern double                   msign; //Expectation value of the sign - also necessary for recovering the correct normalization
extern int                        aac;
extern double                     ans;
extern int                        nmc; //Counter of calls to gather_mc_stat()
extern int*            action_counter;
extern int*                ns_history;       //MC history of sequence lengths

//If control_max_ampl_sum=1 and max_ampl_sum is set to some nonzero value, 
//an error message is generated everytime nA exceeds max_ampl_sum, 
//this is useful for debugging
extern int       control_max_ampl_sum; 
extern double            max_ampl_sum;
extern double        max_ampl_sum_tol;

//Variables for estimating return times
extern int     prev_return_time;
extern int     n_returns;
extern double  mean_rt;
extern double  mean_rt2;
extern double  mean_rt4;
extern int      max_rt;


//These are the variables which are only set by process_mc_stat
extern double acceptance_rate;
extern double mean_recursion_depth;
extern double mean_nA;
extern double err_nA;
extern double mean_sign;

void init_metropolis_statistics();
void free_metropolis_statistics();
void gather_mc_stat();
void process_mc_stat(const char* prefix, int save_to_files);

double     f_max_ampl_sum();
void print_max_amplitudes();

void init_pplus_tuning();
void collect_pplus_tuning_data();
void tune_pplus();
void free_pplus_tuning();

#endif
