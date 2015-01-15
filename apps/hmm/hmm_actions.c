#include "hmm_actions.h"

//Create new factorized-out line
double action_create_amplitude(void* data_in)
{
 return 1.0/NN;
}

int action_create_do(void** data_out, void* data_in)
{
 return 0;
}

int action_create_undo(void* data_out, void* data_in)
{
 return 0;
}

//Create new factorized-in line
double action_evolve_line_amplitude(void* data_in)
{
 return 2.0/cc;
}

int action_evolve_line_do(void** data_out, void* data_in)
{
 return 0;
}

int action_evolve_line_undo(void* data_out, void* data_in)
{
 return 0;
}

//Create new vertex
double action_evolve_vertex_amplitude(void* data_in)
{
 if(stack[ns][stop[ns]]>1)
  return fabs(lambda)*cc;
 return 0.0;
}

int   action_evolve_vertex_do(void** data_out, void* data_in)
{
 return 0;
}

int   action_evolve_vertex_undo(void* data_out, void* data_in)
{
 return 0;
}

//Join two sets of lines
double action_join_amplitude(void* data_in)
{
 if(stop[ns]>0)
  return NN/SQR(cc);
 return 0.0;
}

int   action_join_do(void** data_out, void* data_in)
{
 return 0;
}

int   action_join_undo(void* data_out, void* data_in)
{
 return 0;
}
