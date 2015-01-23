#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>

#include "parameters.h"
#include "actions.h"


extern double    astop;
extern int*      G_hist[2];

void init_observable_stat();
void gather_observable_stat();
void process_observable_stat();
void free_observable_stat();

#endif
