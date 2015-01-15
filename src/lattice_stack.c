#include "lattice_stack.h"

//Parameters of our implementation of the stack
int DIM                     =   0;            //Our stack will contain sequences of lattice points which are the integer coordinates of a DIM-dimensional lattice
int lat_stack_max_nel       =   10000;        //Maximal total number of lattice points in all the sequences in the stack

void print_lattice_stack_parameters()
{
 logs_Write(0, "\tPARAMETERS OF A STACK OF SEQUENCES OF LATTICE POINTS");
 logs_WriteParameter(            "Lattice dimensionality", "%i",       DIM);
 logs_WriteParameter(           "Max. number of elements", "%i",       lat_stack_max_nel);
 logs_WriteParameter(   "Max. range of coordinate values", "%i...%i",  SHRT_MIN, SHRT_MAX);
}

void init_lattice_stack(t_lat_stack* lat_stack)
{
 int i;
 SAFE_MALLOC( lat_stack->stack,      sint*, lat_stack_max_nel);
 SAFE_MALLOC( lat_stack->seq_start,    int, lat_stack_max_nel);
 SAFE_MALLOC( lat_stack->seq_length,   int, lat_stack_max_nel);
 for(i=0; i<lat_stack_max_nel; i++)
  SAFE_MALLOC(lat_stack->stack[i], sint, DIM);
 lat_stack->stack_top = 0;
 lat_stack->stack_nel = 0; 
}

void free_lattice_stack(t_lat_stack* lat_stack)
{
 int i;
 for(i=0; i<lat_stack_max_nel; i++)
  SAFE_FREE(lat_stack->stack[i]);
 SAFE_FREE(lat_stack->stack);
 SAFE_FREE(lat_stack->seq_start);
 SAFE_FREE(lat_stack->seq_length);
}
