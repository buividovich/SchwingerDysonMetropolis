#ifndef _ACTIONS_H_
#define _ACTIONS_H_

#include <sd_metropolis.h>
#include <lattice_stack.h>
#include <largeN_QFT_parameters.h>

#include "parameters.h"

extern t_lat_stack X; //This stack is the current state of the system
extern t_lat_stack H; //This stack will contain the data related to the sequence of actions

//These two functions set up a configuration space on which the actions will be performed
void init_actions();
void free_actions();

//Create new factorized-out line
DECLARE_ACTION_AMPLITUDE( create);
DECLARE_ACTION_DO(        create);
DECLARE_ACTION_UNDO(      create);

//Create new factorized-in line
DECLARE_ACTION_AMPLITUDE( evolve_line);
DECLARE_ACTION_DO(        evolve_line);
DECLARE_ACTION_UNDO(      evolve_line);

//Create new vertex
DECLARE_ACTION_AMPLITUDE( evolve_vertex);
DECLARE_ACTION_DO(        evolve_vertex);
DECLARE_ACTION_UNDO(      evolve_vertex);

//Join two sets of lines
DECLARE_ACTION_AMPLITUDE( join);
DECLARE_ACTION_DO(        join);
DECLARE_ACTION_UNDO(      join);

//Action fetcher
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length, int step_number);

#endif
