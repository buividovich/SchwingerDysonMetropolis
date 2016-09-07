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
typedef int (*t_action)(  int* data_in); //In future, probably more data will be necessary to characterize the actions. For the time being, however, single integer seems enough
//On success, t_action should return 0 
typedef double (*t_action_amplitude)(int* data_in); //Generic functional type for action amplitudes

#define  DECLARE_ACTION_AMPLITUDE(_action_name) double  action_##_action_name##_amplitude(int* data_in)
#define         DECLARE_ACTION_DO(_action_name) int     action_##_action_name##_do(       int* data_in)
#define       DECLARE_ACTION_UNDO(_action_name) int     action_##_action_name##_undo(     int* data_in)

//Error codes returned by action_do and action_undo
#define  ACTION_SUCCESS             (1)
#define  ERR_WRONG_STATE            (-2)
#define  ERR_HISTORY_OVERFLOW       (-3)
#define  ERR_STACK_OVERFLOW         (-4)
#define  ERR_WRONG_DATA             (-5)
#define  ERR_OTHER                  (-10) 

//Structure which contains all the information related to some single action
typedef struct 
{
 int            action_id;          //Id of this action - included for further extensions
 int            action_data_in;     //Parameter which should be passed to the action_do function
} t_action_data;

#define FETCH_ACTION(_action_name, _action_id, _action_list, _amplitude_list, _list_len, _acounter, _adata, _amp_temp)                           \
{                                                                                                                                                           \
 (_amp_temp) = action_##_action_name##_amplitude(&(_adata));                                                                                                \
 if(fabs(_amp_temp)>0.0)                                                                                                                                    \
 {                                                                                                                                                          \
  if((_acounter)>=(_list_len))                                                                                                                              \
  {                                                                                                                                                         \
   SAFE_REALLOC(   (_action_list), t_action_data, ((_acounter)+1));                                                                                         \
   SAFE_REALLOC((_amplitude_list),        double, ((_acounter)+1));                                                                                         \
   logs_Write(0, "Step %08i:", step_number);                                                                                                                \
   logs_Write(0, "\t Action Fetcher:\t Reallocating action_list and amplitude_list to hold %i elements instead of %i", ((_acounter)+1), (_list_len));       \
  };                                                                                                                                                        \
  (_amplitude_list)[(_acounter)]                  = (_amp_temp);                                                                                            \
     (_action_list)[(_acounter)].action_id        = (_action_id);                                                                                           \
     (_action_list)[(_acounter)].action_data_in   = (_adata);                                                                                               \
  (_acounter)++;                                                                                                                                            \
 }                                                                                                                                                          \
} 

#define ADD_TO_ACTION_COLLECTION(_action_name, _action_id)                               \
{                                                                                        \
 action_collection_do[(_action_id)]        = &(action_##_action_name##_do);              \
 action_collection_undo[(_action_id)]      = &(action_##_action_name##_undo);            \
 action_collection_amplitude[(_action_id)] = &(action_##_action_name##_amplitude);       \
 SAFE_MALLOC(action_collection_name[(_action_id)], char, strlen(#_action_name)+2);       \
 sprintf(action_collection_name[(_action_id)], "%s", #_action_name);                     \
}

//Functional type for a generic function which returns all possible actions and their probabilities
typedef int (*t_action_fetcher)(t_action_data** action_list, double** amplitude_list, int list_length); //Should return the number of possible actions, if there are more actions than the length of the currently allocated action_list, should perform realloc

//These pointers to functions should be set in the user code before calling init_metropolis
extern t_action*            action_collection_do;
extern t_action*            action_collection_undo;
extern t_action_amplitude*  action_collection_amplitude;
extern char**               action_collection_name;
extern int                  action_collection_size;
extern t_action             state_initializer;
extern t_action_fetcher     action_fetcher;

//Variables characterizing the state of the random process which should be visible outside
extern int            step_number;      //Counts the number of calls to metropolis_step() after init_metropolis()
extern int            ns;               //Counter of the sequence depth  
extern double*        nA;               //Contains the probability normalization factor on every step
extern int*           asign;            //Reweighting signs of every configuration
extern t_action_data* action_history;   //History of actions, elements up to ns are referenced

void  init_metropolis();
void  free_metropolis();
int   metropolis_step();

void  call_state_initializer();

void  print_action_history();

#endif
