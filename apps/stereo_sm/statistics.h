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

typedef struct 
{
 double    astop;
 double*   Gx[2];
 double*   Gxy[2];
 double    Gtotal;
 int*      sampling_hist;
 int       actual_max_alpha_order;
 int       nstat;
 int       nstat_useless;
} t_observable_stat;

t_observable_stat* init_observable_stat();
void gather_observable_stat(t_observable_stat* stat);
void process_observable_stat(t_observable_stat* stat);
void free_observable_stat(t_observable_stat* stat);

#endif
