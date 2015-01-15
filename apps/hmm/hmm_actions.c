#include "hmm_actions.h"

t_lat_stack lat_stack;
t_lat_stack history_data;

void init_hmm_actions()
{
 init_lattice_stack(&lat_stack);
 init_lattice_stack(&history_data);    
}

void free_hmm_actions()
{
 free_lattice_stack(&lat_stack);
 free_lattice_stack(&history_data);
}

//Create new factorized-out line
DECLARE_ACTION_AMPLITUDE(create)
{
 return 1.0/NN;
}

DECLARE_ACTION_DO(create)
{
 return 0;
}

DECLARE_ACTION_UNDO(create)
{
 return 0;
}

//Create new factorized-in line
DECLARE_ACTION_AMPLITUDE(evolve_line)
{
 return 2.0/cc;
}

DECLARE_ACTION_DO(evolve_line)
{
 return 0;
}

DECLARE_ACTION_UNDO(evolve_line)
{
 return 0;
}

//Create new vertex
DECLARE_ACTION_AMPLITUDE(evolve_vertex)
{
 //if(stack[ns][stop[ns]]>1)
 // return fabs(lambda)*cc;
 return 0.0;
}

DECLARE_ACTION_DO(evolve_vertex)
{
 return 0;
}

DECLARE_ACTION_UNDO(evolve_vertex)
{
 return 0;
}

//Join two sets of lines
DECLARE_ACTION_AMPLITUDE(join)
{
 //if(stop[ns]>0)
 // return NN/SQR(cc);
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 return 0;
}

DECLARE_ACTION_UNDO(join)
{
 return 0;
}
