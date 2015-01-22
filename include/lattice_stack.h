#ifndef _LATTICE_STACK_H_
#define _LATTICE_STACK_H_

#include <stdlib.h>
#include <stdio.h>
#include <clue_logs.h>
#include <limits.h>

typedef short    int sint;

//Stack itself
typedef struct 
{
 sint** stack;       //sint[i] contains the pointer to the array of DIM sint's
 int*   start;       //seq_start[i]  is the first index of i'th sequence in the stack
 int*   len;         //seq_length[i] is the length of i'th sequence in the stack
 int    top;         //Number of sequences in the stack
 int    nel;         //Total number of elements in the stack
} t_lat_stack;


extern int DIM;                     //Our stack will contain sequences of lattice points which are the integer coordinates of a DIM-dimensional lattice
extern int lat_stack_max_nel;       //Maximal total number of lattice points in all the sequences in the stack

void  print_lattice_stack_parameters();
void  init_lattice_stack(t_lat_stack* lat_stack);
void  free_lattice_stack(t_lat_stack* lat_stack);
int   check_stack_consistency(t_lat_stack* lat_stack, const char* stack_name);

#define LATTICE_STACK_LONG_OPTIONS                                                  \
 {                     "DIM",  required_argument,                       NULL, 'Z'}, \
 {       "lat-stack-max-nel",  required_argument,                       NULL, 'Y'}

#define PARSE_LATTICE_STACK_OPTIONS                                      \
   case 'Z':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%i", DIM);                                \
    ASSERT(DIM<0);                                                       \
   break;                                                                \
   case 'Y':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%i", lat_stack_max_nel);                  \
   break;                                                                

static const char lattice_stack_short_option_list[] = "Z:Y:";

#endif
