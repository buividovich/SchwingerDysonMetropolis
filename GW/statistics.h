#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <stdio.h>
#include <stdlib.h>

#include <linalg.h>

#include "my_stack.h"
#include "metropolis.h"
#include "parameters.h"

extern double   maxnA;
extern double     anA;
extern double     dnA;
extern double   msign;
extern int        aac;
extern double     ans;
extern double   astop;
extern int    *G_hist[2];
extern double *ns_hist;

void init_mc_stat();
void gather_mc_stat();
void process_mc_stat();
void free_mc_stat();

void init_observable_stat();
void gather_observable_stat();
void process_observable_stat_sc();
void process_observable_stat_wc();
void free_observable_stat();

#endif
