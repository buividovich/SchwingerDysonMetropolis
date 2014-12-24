#ifndef _METROPOLIS_H_
#define _METROPOLIS_H_

#include <stdlib.h>
#include <stdio.h>

#include <rand_num_generators.h>

#include "parameters.h"
#include "my_stack.h"
#include "statistics.h"

extern int overflow_count;

int mc_step_sc();
int mc_step_wc();
int tune_parameters_sc(double accuracy);

#endif
