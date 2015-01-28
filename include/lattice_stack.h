#ifndef _LATTICE_STACK_H_
#define _LATTICE_STACK_H_

#include <stdlib.h>
#include <stdio.h>
#include <clue_logs.h>
#include <limits.h>

//Stack itself
typedef struct 
{
 int    dim;         //Dimensionality of stack elements
 int    max_nel;     //Max number of elements in stack
 int**  stack;       //sint[i] contains the pointer to the array of DIM sint's
 int*   start;       //seq_start[i]  is the first index of i'th sequence in the stack
 int*   len;         //seq_length[i] is the length of i'th sequence in the stack
 int    top;         //Number of sequences in the stack
 int    nel;         //Total number of elements in the stack
 
} t_lat_stack;

void init_lat_stack(t_lat_stack* lat_stack, int dim, int max_nel);
void free_lat_stack(t_lat_stack* lat_stack);

int   check_stack_consistency(t_lat_stack* lat_stack, const char* stack_name);

#endif
