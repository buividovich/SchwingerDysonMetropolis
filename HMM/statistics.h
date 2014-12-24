#ifndef _STATISTICS_H_
#define _STATISTICS_H_

#include <stdio.h>
#include <stdlib.h>

#include <linalg.h>

#include "my_stack.h"
#include "parameters.h"
#include "metropolis.h"

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
void process_observable_stat();
void free_observable_stat();

#endif
