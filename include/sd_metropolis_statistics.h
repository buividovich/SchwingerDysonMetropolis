#ifndef _SD_METROPOLIS_STATISTICS_H_
#define _SD_METROPOLIS_STATISTICS_H_

#include <stdlib.h>
#include <stdio.h>
#include <clue_utils.h>
#include <clue_logs.h>
#include <clue_io.h>

#include "sd_metropolis.h"
#include "largeN_QFT_parameters.h"

extern double                   maxnA; //Max. value of nA
extern double                     anA; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
extern double                     dnA; //Expectation value of its square - necessary to estimate the error!
extern double                   msign; //Expectation value of the sign - also necessary for recovering the correct normalization
extern int                        aac;
extern double                     ans;
extern int                     max_ns;
extern double                   astop;
extern int                   max_stop;
extern int                        nmc; //Counter of calls to gather_mc_stat()
extern int*            action_counter;

//Variables for estimating return times
extern int     n_returns;
extern double  mean_return_time;

//These are the variables which are only set by process_mc_stat
extern double acceptance_rate;
extern double mean_recursion_depth;
extern double mean_stack_top;
extern double mean_nA;
extern double mean_sign;

void init_metropolis_statistics();
void free_metropolis_statistics();
void gather_mc_stat();
void process_mc_stat();
void save_mc_stat();

void init_pplus_tuning();
void collect_pplus_tuning_data();
void tune_pplus();
void free_pplus_tuning();

#endif
