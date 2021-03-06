#ifndef _ACTIONS_H_
#define _ACTIONS_H_

#include <sd_metropolis.h>
#include <lattice_stack.h>
#include <largeN_QFT_parameters.h>
#include <lattice_propagator.h>

#include "parameters.h"

extern t_lat_stack      X; //This stack is the current state of the system
extern t_lat_stack      H; //This stack will contain the data related to the sequence of actions

extern int beta_order;

extern t_lat_propagator P;

//These two functions set up a configuration space on which the actions will be performed
void init_actions();
void free_actions();

//Custom rand_momentum which works faster in the case of pure SC expansion
void my_rand_momentum(int* m);

//Create new factorized-out line
DECLARE_ACTION_AMPLITUDE( create);
DECLARE_ACTION_DO(        create);
DECLARE_ACTION_UNDO(      create);

//Create new factorized-in line
DECLARE_ACTION_AMPLITUDE( add_line);
DECLARE_ACTION_DO(        add_line);
DECLARE_ACTION_UNDO(      add_line);

//Join two sets of lines
DECLARE_ACTION_AMPLITUDE( join);
DECLARE_ACTION_DO(        join);
DECLARE_ACTION_UNDO(      join);

//Flip momenta on the two sets of lines
DECLARE_ACTION_AMPLITUDE( exchange_momenta);
DECLARE_ACTION_DO(        exchange_momenta);
DECLARE_ACTION_UNDO(      exchange_momenta);

//Create new vertex
DECLARE_ACTION_AMPLITUDE( vertex);
DECLARE_ACTION_DO(        vertex);
DECLARE_ACTION_UNDO(      vertex);

//The trivial random walk term
DECLARE_ACTION_AMPLITUDE( random_walk);
DECLARE_ACTION_DO(        random_walk);
DECLARE_ACTION_UNDO(      random_walk);


//Action fetcher
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length);

#endif
