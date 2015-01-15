#ifndef _HMM_ACTIONS_H_
#define _HMM_ACTIONS_H_

#include "hmm_parameters.h"
#include "sd_metropolis.h"

//Create new factorized-out line
DECLARE_ACTION_DO(create);
DECLARE_ACTION_UNDO(create);

//Create new factorized-in line
DECLARE_ACTION_DO(evolve_line);
DECLARE_ACTION_UNDO(evolve_line);

//Create new vertex
DECLARE_ACTION_DO(evolve_vertex);
DECLARE_ACTION_UNDO(evolve_vertex);

//Join two sets of lines
DECLARE_ACTION_DO(join);
DECLARE_ACTION_UNDO(join);

#endif
