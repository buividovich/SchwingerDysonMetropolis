#ifndef _METROPOLIS_H_
#define _METROPOLIS_H_

#include <stdlib.h>
#include <stdio.h>

#include <rand_num_generators.h>

#include "parameters.h"
#include "my_stack.h"

extern int overflow_count;

int mc_step();

double A_create();
double A_evolve_line();
double A_evolve_vertex();
double A_join();

#endif
