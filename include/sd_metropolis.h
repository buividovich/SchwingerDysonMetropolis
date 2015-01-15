#ifndef _SD_METROPOLIS_H_
#define _SD_METROPOLIS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <clue_utils.h>
#include <clue_errors.h>
#include <clue_logs.h>
#include <rand_num_generators.h>

#include "sd_metropolis_statistics.h"
#include "sd_metropolis_parameters.h"

//Functional type for an elementary action on the configuration space
typedef int (*t_action_do)(void** data_out, void* data_in);
typedef int (*t_action_undo)(void* data_out, void* data_in);

#define  DECLARE_ACTION_AMPLITUDE(_action_name) static inline double action##_action_name##_amplitude(void*  data_in)
#define         DECLARE_ACTION_DO(_action_name) int                  action##_action_name##_do(       void** data_out, void* data_in)
#define       DECLARE_ACTION_UNDO(_action_name) int                  action##_action_name##_undo(     void*  data_out, void* data_in)

//Structure which contains all the information related to some single action
typedef struct 
{
 int            action_id;
 double         action_amplitude;
 t_action_do    action_do;
 t_action_undo  action_undo;
 void*          action_data_in;
 void*          action_data_out;
} t_action;

//Functional type for a generic function which returns all possible actions and their probabilities
typedef int (*t_action_fetcher)(t_action** action_list); //Should return the number of possible actions

//Functional type for a generic function which initializes the configuration 
typedef int (*t_state_initializer)(); //Should return the "sign" of the initial state

extern int     overflow_count; //Counts exceptions when the sequence of elements becomes larger than allocated 
extern int     ns;             //Counter of the sequence depth  
extern double* nA;             //Contains the probability normalization factor on every step
extern int*    asign;          //Reweighting signs of every configuration

int   init_metropolis(t_state_initializer state_initializer);
void  free_metropolis();
int   metropolis_step(t_action_fetcher action_fetcher, t_state_initializer state_initializer);

#endif
