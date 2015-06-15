#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>

#include <sd_metropolis_parameters.h>
#include <largeN_QFT_parameters.h>

#include "clue_io.h"

#include "parameters.h"
#include "actions.h"
#include "stack_statistics.h"

extern double    astop;
extern int*      G2_hist[2];
extern int       max_G2_order;

void init_observable_stat();
void gather_observable_stat();
void process_observable_stat();
void free_observable_stat();

#endif
