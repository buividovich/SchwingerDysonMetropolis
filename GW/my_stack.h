#ifndef _MY_STACK_H_
#define _MY_STACK_H_

#include <stdlib.h>
#include <stdio.h>

#include <clue_errors.h>

#include "parameters.h"

extern int    **stack; //Storage for the correlators themselves
extern int      *stop; //Sequence length
extern double     *nA; //Here sums of probabilities will be saved - in order to avoid recalculations
extern int         ns; //The number of the topmost element in a sequence
extern int     *asign; //Sign of the configuration

void init_stack();
void free_stack();

#endif
