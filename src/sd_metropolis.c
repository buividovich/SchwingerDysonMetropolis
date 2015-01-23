#include "sd_metropolis.h"

//These pointers to functions should be set in the user code before calling init_metropolis
t_action*          action_collection_do     =   NULL;
t_action*          action_collection_undo   =   NULL;
int                action_collection_size   =   NULL;
t_action           state_initializer        =   NULL;
t_action_fetcher   action_fetcher           =   NULL;

//Variables characterizing the state of the random process which should be visible outside
int            step_number        = 0;    //Counts the number of calls to metropolis_step() after init_metropolis()
int            ns                 = 0;    //Counter of the sequence depth  
double*        nA                 = NULL; //Contains the probability normalization factor on every step
int*           asign              = NULL; //Reweighting signs of every configuration
t_action_data* action_history     = NULL; //History of actions, elements up to ns are referenced

//These are internal variables of sd_metropolis
t_action_data* action_list        = NULL; //List of currently possible actions at every MC step
double*        amplitude_list     = NULL; //List of the amplitudes of currently possible actions
double*        probability_list   = NULL; //LIst of the probabilities of currently possible actions
int            action_list_length = 0;


void   init_metropolis()
{
 step_number    = 0;
 ns             = 0;
 SAFE_MALLOC(             nA,        double, max_recursion_depth);
 SAFE_MALLOC(          asign,           int, max_recursion_depth);
 SAFE_MALLOC( action_history, t_action_data, max_recursion_depth);
 
 action_list        = NULL;
 amplitude_list     = NULL;
 probability_list   = NULL;
 action_list_length = 0;
 
 init_metropolis_statistics();
 
 ASSERT(   action_collection_do == NULL);
 ASSERT( action_collection_undo == NULL);
 ASSERT(      state_initializer == NULL);
 ASSERT( action_collection_size <= 0);
 
 asign[ns] = state_initializer(NULL);      
}

void  free_metropolis()
{
 SAFE_FREE(nA);
 SAFE_FREE(asign);
 SAFE_FREE(action_history);
 SAFE_FREE(action_list);
 SAFE_FREE(amplitude_list);
 SAFE_FREE(probability_list);
}

int   metropolis_step()
{
 int n_actions, iaction, todo, accepted = 0, res;
 double a, alpha;
 
 if(ns>max_recursion_depth-1)
 {
  logs_WriteError("Overflow detected - size of action_history is larger than the allocated memory!!! Resetting the state of the system!!!");
  ns = 0;
  state_initializer(NULL);                           
 };
 
 n_actions = action_fetcher(&action_list, &amplitude_list, action_list_length); //Action fetcher should also perform reallocation of the action_list if necessary
 
 if(n_actions>action_list_length)
 {
  logs_Write(0, "\t     Metropolis:\t Reallocating               probability_list to hold %i elements instead of %i", n_actions, action_list_length);
  action_list_length = n_actions;
  SAFE_REALLOC(probability_list, double, action_list_length);
 };
 
 if(logs_noise_level>=2)
 {
  logs_Write(3, "Fetched the following action amplitudes at step %i: ", step_number);
  for(iaction=0; iaction<n_actions; iaction++)
   logs_Write(3, "Action %i: \t ampl = %2.4E, action_data_in = %i", action_list[iaction].action_id, amplitude_list[iaction], action_list[iaction].action_data_in);
 }; 
 
 nA[ns] = 0.0;
 for(iaction=0; iaction<n_actions; iaction++)
  nA[ns] += fabs(amplitude_list[iaction]);
 
 //Decide whether we move forward or backward in the order of expansion
 a = rand_double(0.0, 1.0);
 if(a <= p_plus)
 {//Attaching new element to the sequence
  //Acceptance probability
  alpha  = MIN(1.0, (1.0 - p_plus)/p_plus*nA[ns]);
  a      = rand_double(0.0, 1.0);      
  if(a<alpha)
  {//Accept and perform the forward step
   accepted = 1;
   asign[ns+1] = asign[ns];
   //Now calculate the normalized probabilities of different steps
   for(iaction=0; iaction<n_actions; iaction++)
    probability_list[iaction] = fabs(amplitude_list[iaction])/nA[ns];
   //And increase the order counter 
   ns++; 
   //Choose at random which action to perform
   todo = rand_choice(probability_list, n_actions);
   //Perform the selected action
   if(action_list[todo].action_id<action_collection_size)
   {
    res = (action_collection_do[action_list[todo].action_id])(&(action_list[todo].action_data_in));
    if(res==-2 || res==-3)
    {
     logs_WriteError("Action %i with action_data_in %i caused %s overflow at mc step %i, resetting the state of the system...", action_list[todo].action_id, action_list[todo].action_data_in, (res==-2? "history" : "state"), step_number);
     ns = 0;
     state_initializer(NULL);
    };
    if(res==-1)
     logs_WriteError("Attempted to apply action %i with action_data_in %i to the improper system state, nothing was really done...", action_list[todo].action_id, action_list[todo].action_data_in);
   } 
   else
   {
    logs_WriteError("action_list[%i].action_id = %i is larger than the total number %i of actions in the collection", todo, action_list[todo].action_id, action_collection_size); 
   }; 
   //Write to the log
   logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t Performed action %i, action_data_in = %i, ns = %i, nA[ns] = %2.4lf", step_number, action_list[todo].action_id, action_list[todo].action_data_in, ns, nA[ns]);
   //And save it to the history
   action_history[ns] = action_list[todo];
   //Modify the reweighting sign/phase
   asign[ns] *= SIGN(amplitude_list[todo]);
  };
 }
 else
 {//Step backward
  if(ns>0)
  {
   alpha  = MIN( 1.0, p_plus/((1.0 - p_plus)*nA[ns-1]) ); //Acceptance probability
   a      = rand_double(0.0, 1.0);
   if(a<alpha)
   {
    accepted = 1;
    //Undo the last action!!!
    logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t Undoing the action %i, action_data_in = %i, ns = %i, nA[ns-1] = %2.4lf", step_number, action_history[ns].action_id, action_history[ns].action_data_in, ns, nA[ns-1]);
    if(action_history[ns].action_id<action_collection_size)
    {
     res = (action_collection_undo[action_history[ns].action_id])(&(action_history[ns].action_data_in));
     if(res==-2 || res==-3)
     {
      logs_WriteError("Undoing action %i with action_data_in %i caused %s overflow at mc step %i, resetting the state of the system...", action_history[ns].action_id, action_history[ns].action_data_in, (res==-2? "history" : "state"), step_number);
      ns = 0;
      state_initializer(NULL);
     };
     if(res==-1)
      logs_WriteError("Attempted to undo action %i with action_data_in %i from the improper system state, nothing was really done...", action_history[ns].action_id, action_history[ns].action_data_in);
    } 
    else
    {
     logs_WriteError("action_history[%i].action_id = %i is larger than the total number %i of actions in the collection", ns, action_history[ns].action_id, action_collection_size);  
    };
    //Finally, decrease the counter
    ns--;
   }; 
  }
  else
  {
   accepted = 1;
   asign[ns] = state_initializer(NULL);
   //Write to the log
   logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t state_initializer has been just called", step_number);
  };
 };
 
 //Updating statistics of the MC process - the mean nA, its dispersion and the mean sign
 gather_mc_stat(accepted);
 
 step_number ++;
 
 return accepted; //Return 1 if some new state is accepted, otherwise 0
}

