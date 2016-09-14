#include "lattice_stack.h"

void init_lat_momentum_stack(t_lat_momentum_stack* lat_stack, int dim, int max_nel)
{
 lat_stack->dim     = dim;
 lat_stack->max_nel = max_nel;
 SAFE_MALLOC( lat_stack->stack,     int*, max_nel);
 SAFE_MALLOC( lat_stack->start,      int, max_nel);
 SAFE_MALLOC( lat_stack->len,        int, max_nel);
 for(int i=0; i<max_nel; i++)
  SAFE_MALLOC(lat_stack->stack[i],   int, dim);
 lat_stack->top = 0;
 lat_stack->nel = 0;
}

void free_lat_momentum_stack(t_lat_momentum_stack* lat_stack)
{
 for(int i=0; i<lat_stack->max_nel; i++)
  SAFE_FREE(lat_stack->stack[i]);
 SAFE_FREE(lat_stack->stack);
 SAFE_FREE(lat_stack->start);
 SAFE_FREE(lat_stack->len);
}

int   check_stack_consistency(t_lat_momentum_stack* lat_stack, const char* stack_name)
{
 int res = 0, my_nel = 0;
 
 if(lat_stack->top>0)
 {
  if(lat_stack->start[0]!=0)
  {
   res = -1;
   logs_WriteError("STACK CHECK:\t %s->start[0]!=0", stack_name);
  };
  my_nel += lat_stack->len[0];
 }; 
 
 for(int i=1; i<lat_stack->top; i++)
 {
  if(lat_stack->start[i] != (lat_stack->start[i-1] + lat_stack->len[i-1]))
  {
   res = -2;
   logs_WriteError("STACK CHECK:\t %s->start[%i] = %i", stack_name,   i, lat_stack->start[i]);
   logs_WriteError("STACK CHECK:\t %s->start[%i] = %i", stack_name, i-1, lat_stack->start[i-1]);
   logs_WriteError("STACK CHECK:\t %s->len[%i]   = %i", stack_name, i-1, lat_stack->len[i-1]);
  };
  my_nel += lat_stack->len[i];
 }; 
 
 if(my_nel != lat_stack->nel)
 {
  logs_WriteError("my_nel = %i, %s->nel = %i", my_nel, stack_name, lat_stack->nel);
  res = -3;
 }; 
   
 return res;     
}

void init_lat_coordinate_stack(t_lat_coordinate_stack* lat_stack, int max_nel)
{
 lat_stack->max_nel = max_nel;
 SAFE_MALLOC( lat_stack->stack,      int, max_nel);
 SAFE_MALLOC( lat_stack->start,      int, max_nel);
 SAFE_MALLOC( lat_stack->len,        int, max_nel);
 lat_stack->top = 0;
 lat_stack->nel = 0;
}

void free_lat_coordinate_stack(t_lat_coordinate_stack* lat_stack)
{
 SAFE_FREE(lat_stack->stack);
 SAFE_FREE(lat_stack->start);
 SAFE_FREE(lat_stack->len);
}
