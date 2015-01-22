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

//Possible improvements in the code
//Important: Fixed arrays for pointers to actions + do not save amplitudes
//Important: in case of overflow, reset the stack and the history
//Important: Function for generating random momenta, struct t_propagator
//Important: calculation of sigma
//Important: auto-tuning of actions, more automatic
//Important: LS, LT, DIM in large-N QFT parameters
//Important: make dataclean
//Important: make mcstat_plots - arguments for gnuplot?
//Important: update git repository
//Can wait:  Numerical values for characters in getopt
//Can wait:  Coordinates as integers in stack, not as lists - can give some speedup...
//Can wait:  Timing? Compare, how much time is taken by stack manipulations? How fast is ranlux?

//Functional type for an elementary action on the configuration space
typedef int (*t_action)(  int* data_in); //In future, probably more data will be necessary to characterize the actions. For the time being, however, single integer seems enough

#define  DECLARE_ACTION_AMPLITUDE(_action_name) double  action_##_action_name##_amplitude(int* data_in)
#define         DECLARE_ACTION_DO(_action_name) int     action_##_action_name##_do(       int* data_in)
#define       DECLARE_ACTION_UNDO(_action_name) int     action_##_action_name##_undo(     int* data_in)

//Structure which contains all the information related to some single action
typedef struct 
{
 int            action_id;          //Id of this action - included for further extensions
 int            action_data_in;     //Parameter which should be passed to the action_do function
} t_action_data;

#define FETCH_ACTION(_action_name, _action_id, _action_list, _amplitude_list, _list_len, _acounter, _adata, _amp_temp, _step_num)                           \
{                                                                                                                                                           \
 (_amp_temp) = action_##_action_name##_amplitude(&(_adata));                                                                                                \
 if(fabs(_amp_temp)>0.0)                                                                                                                                    \
 {                                                                                                                                                          \
  if((_acounter)>=(_list_len))                                                                                                                              \
  {                                                                                                                                                         \
   SAFE_REALLOC(   (_action_list), t_action_data, ((_acounter)+1));                                                                                         \
   SAFE_REALLOC((_amplitude_list),        double, ((_acounter)+1));                                                                                         \
   logs_Write(0, "Step %08i:", (_step_num));                                                                                                                \
   logs_Write(0, "\t Action Fetcher:\t Reallocating action_list and amplitude_list to hold %i elements instead of %i", ((_acounter)+1), (_list_len));       \
  };                                                                                                                                                        \
  (_amplitude_list)[(_acounter)]                  = (_amp_temp);                                                                                            \
     (_action_list)[(_acounter)].action_id        = (_action_id);                                                                                           \
     (_action_list)[(_acounter)].action_data_in   = (_adata);                                                                                               \
  (_acounter)++;                                                                                                                                            \
 }                                                                                                                                                          \
} 

#define ADD_TO_ACTION_COLLECTION(_action_name, _action_id)                     \
{                                                                              \
 action_collection_do[(_action_id)]   = &(action_##_action_name##_do);         \
 action_collection_undo[(_action_id)] = &(action_##_action_name##_undo);       \
}

//Functional type for a generic function which returns all possible actions and their probabilities
typedef int (*t_action_fetcher)(t_action_data** action_list, double** amplitude_list, int list_length, int step_number); //Should return the number of possible actions, if there are more actions than the length of the currently allocated action_list, should perform realloc

//These pointers to functions should be set in the user code before calling init_metropolis
extern t_action*          action_collection_do;
extern t_action*          action_collection_undo;
extern int                action_collection_size;
extern t_action           state_initializer;
extern t_action_fetcher   action_fetcher;

//Variables characterizing the state of the random process which should be visible outside
extern int            ns;               //Counter of the sequence depth  
extern double*        nA;               //Contains the probability normalization factor on every step
extern int*           asign;            //Reweighting signs of every configuration
extern t_action_data* action_history;   //History of actions, elements up to ns are referenced

int   init_metropolis();
void  free_metropolis();
int   metropolis_step(int step_number);

#endif
