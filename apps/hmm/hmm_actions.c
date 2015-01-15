#include "hmm_actions.h"

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
