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
 SAFE_MALLOC( lat_stack->stack,    sint*, lat_stack_max_nel);
 SAFE_MALLOC( lat_stack->start,      int, lat_stack_max_nel);
 SAFE_MALLOC( lat_stack->len,        int, lat_stack_max_nel);
 for(i=0; i<lat_stack_max_nel; i++)
  SAFE_MALLOC(lat_stack->stack[i], sint, DIM);
 lat_stack->top = 0;
 lat_stack->nel = 0;
}

void free_lattice_stack(t_lat_stack* lat_stack)
{
 int i;
 for(i=0; i<lat_stack_max_nel; i++)
  SAFE_FREE(lat_stack->stack[i]);
 SAFE_FREE(lat_stack->stack);
 SAFE_FREE(lat_stack->start);
 SAFE_FREE(lat_stack->len);
}

int   check_stack_consistency(t_lat_stack* lat_stack, const char* stack_name)
{
 int res = 0, my_nel = 0, i;
 
 if(lat_stack->top>0)
 {
  if(lat_stack->start[0]!=0)
  {
   res = -1;
   logs_WriteError("%s->start[0]!=0", stack_name);
  };
  my_nel += lat_stack->len[0];
 }; 
 
 for(i=1; i<lat_stack->top; i++)
 {
  if(lat_stack->start[i] != (lat_stack->start[i-1] + lat_stack->len[i-1]))
  {
   res = -2;
   logs_WriteError("%s->start[%i] = %i", stack_name,   i, lat_stack->start[i]);
   logs_WriteError("%s->start[%i] = %i", stack_name, i-1, lat_stack->start[i-1]);
   logs_WriteError("%s->len[%i]   = %i", stack_name, i-1, lat_stack->len[i-1]);
  };
  my_nel += lat_stack->len[i];
 }; 
 
 if(my_nel != lat_stack->nel)
 {
  logs_WriteError("my_nel = %i, %s->nel = %i", my_nel, stack_name, lat_stack->nel);
  system("PAUSE");
  res = -3;
 }; 
   
 return res;     
}
