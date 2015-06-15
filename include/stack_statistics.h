#ifndef _STACK_STATISTICS_H_
#define _STACK_STATISTICS_H_

/* This file implements the tuning of cc and NN coefficients */

#include <clue_io.h>

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>
#include <lattice_stack.h>

extern double* seq_len_hist;
extern double* num_seq_hist;
extern int     nstat;

void init_stack_statistics(int max_hist_size);
void free_stack_statistics();

void gather_stack_statistics(t_lat_stack* X);

void print_stack_statistics(FILE* f);

double mean_nA_prediction(double NN_new, double cc_new, double source_norm);

#endif
