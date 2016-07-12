#ifndef _STACK_STATISTICS_H_
#define _STACK_STATISTICS_H_

/* This file implements the tuning of cc and NN coefficients */

#include <clue_io.h>

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>
#include <lattice_stack.h>

typedef struct 
{
 int* seq_len_hist;
 int* num_seq_hist;
 int  nstat;
 int  my_max_hist_size;
} t_stack_stat;

t_stack_stat* init_stack_statistics(int max_hist_size);
void free_stack_statistics(t_stack_stat* stat);

void gather_stack_statistics(t_stack_stat* stat, t_lat_stack* X);

void print_stack_statistics(t_stack_stat* stat);
void print_stack_histogram(t_stack_stat* stat, FILE* f);

double mean_nA_prediction(t_stack_stat* stat, double NN_new, double cc_new, double source_norm);

#endif
