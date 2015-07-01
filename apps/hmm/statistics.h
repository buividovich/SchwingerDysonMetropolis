#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>

#include "parameters.h"
#include "actions.h"
#include "stack_statistics.h"

#include "free_hmm_genus_expansion.h"

extern double    astop;
extern int**     G_hist;
extern int*      genus_hist;
extern int       actual_max_genus;

void init_observable_stat();
void gather_observable_stat();
void process_observable_stat();
void free_observable_stat();

#endif
