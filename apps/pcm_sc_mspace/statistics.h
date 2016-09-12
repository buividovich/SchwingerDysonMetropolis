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

typedef struct 
{
 double    astop;
 double*   G2[2];
 int*      sampling_hist;
 int       actual_max_beta_order;
 int       max_history_nel;
 int       max_stack_nel;
 int       nstat;
 int       nstat_useless;
} t_observable_stat;

t_observable_stat* init_observable_stat();
void gather_observable_stat(t_observable_stat* stat);
void process_observable_stat(t_observable_stat* stat);
void free_observable_stat(t_observable_stat* stat);

#endif
