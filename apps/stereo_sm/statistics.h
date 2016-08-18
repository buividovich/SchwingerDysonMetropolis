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
 double  meanA[2];
 double meanA2[2];
 int        nA[2];       
} t_obs;

typedef struct
{
 double meanA;
 double errA;
 double spA;       
} t_obs_res;

void init_obs(t_obs* obs);
void add_to_obs(t_obs* obs, double A, int s);
void process_obs(t_obs* obs, int n, t_obs_res* res);

typedef struct 
{
 double    astop;
 t_obs* Gx;
 t_obs* Gxy;
 int*      G2_hist[2];
 int*      G4_hist[2]; //This is for debugging
 int       nG4;
 int       actual_max_alpha_order;
 int*      G_hist[2];
 int       nstat;
 int       nstat_useless;
 int       max_useful_Xtop;
} t_observable_stat;

t_observable_stat* init_observable_stat();
void gather_observable_stat(t_observable_stat* stat);
void process_observable_stat(t_observable_stat* stat);
void free_observable_stat(t_observable_stat* stat);

#endif
