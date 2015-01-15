#ifndef _HMM_STATISTICS_H_
#define _HMM_STATISTICS_H_

#include <linalg.h>

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>

#include "hmm_parameters.h"
#include "hmm_actions.h"


extern double    astop;
extern int*      G_hist[2];

void init_observable_stat();
void gather_observable_stat();
void process_observable_stat();
void free_observable_stat();

#endif
