#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>

#include "parameters.h"
#include "actions.h"
#include "stack_statistics.h"

#include "free_hmm_genus_expansion.h"

#define  GMN_HIST_INDEX(_g, _top, _nel)\
  ((_g)*(gmn_hist_max_traces*gmn_hist_max_nel) + (_top)*(gmn_hist_max_nel) + (_nel))

extern double    astop;
extern int**     G_hist;
extern int*      gmn_hist;

extern int       actual_max_genus;
extern int       actual_max_top;
extern int       actual_max_nel;

void init_observable_stat();
void gather_observable_stat();
void process_observable_stat();
void free_observable_stat();

#endif
