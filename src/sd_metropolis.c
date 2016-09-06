#include "sd_metropolis.h"

//These pointers to functions should be set in the user code before calling init_metropolis
t_action*            action_collection_do        =   NULL;
t_action*            action_collection_undo      =   NULL;
t_action_amplitude*  action_collection_amplitude =   NULL;
char**               action_collection_name      =   NULL;
int                  action_collection_size      =   0;
t_action             state_initializer           =   NULL;
t_action_fetcher     action_fetcher              =   NULL;

//Variables characterizing the state of the random process which should be visible outside
int            step_number        = 0;    //Counts the number of calls to metropolis_step() after init_metropolis()
int            ns                 = 0;    //Counter of the sequence depth  
double*        nA                 = NULL; //Contains the probability normalization factor on every step
int*           asign              = NULL; //Reweighting signs of every configuration
t_action_data* action_history     = NULL; //History of actions, elements up to ns are referenced

//These are internal variables of sd_metropolis
t_action_data* action_list        = NULL; //List of currently possible actions at every MC step
double*        amplitude_list     = NULL; //List of the amplitudes of currently possible actions
double*        probability_list   = NULL; //List of the probabilities of currently possible actions
int            action_list_length = 0;

void   init_metropolis()
{
 step_number    = 0;
 ns             = 0;
 SAFE_MALLOC_IF_NULL(             nA,        double, (max_recursion_depth+1));
 SAFE_MALLOC_IF_NULL(          asign,           int, (max_recursion_depth+1));
 SAFE_MALLOC_IF_NULL( action_history, t_action_data, (max_recursion_depth+1));
 
 init_metropolis_statistics();
 
 //Check that the action collection is properly initialized
 ASSERT( action_collection_do        == NULL);
 ASSERT( action_collection_undo      == NULL);
 ASSERT( action_collection_amplitude == NULL);
 ASSERT( action_collection_name      == NULL);
 
 ASSERT( action_collection_size      <= 0);
 for(int iact=0; iact<action_collection_size; iact++)
 {
  ASSERT( action_collection_do[iact]        == NULL);
  ASSERT( action_collection_undo[iact]      == NULL);
  ASSERT( action_collection_amplitude[iact] == NULL);
  ASSERT( action_collection_name[iact]      == NULL);      
 };
  
 ASSERT(      state_initializer      == NULL);
 
 //Add a name for state_initializer
 SAFE_REALLOC(action_collection_name, char*, action_collection_size+1);
 SAFE_MALLOC( action_collection_name[action_collection_size], char, 20);
 sprintf(action_collection_name[action_collection_size], "State initializer");
  
 call_state_initializer();
}

void  free_metropolis()
{
 free_metropolis_statistics();
 SAFE_FREE(nA);
 SAFE_FREE(asign);
 SAFE_FREE(action_history);
 SAFE_FREE(action_list);
 SAFE_FREE(amplitude_list);
 SAFE_FREE(probability_list);
}

void call_state_initializer()
{
 ns = 0;
 action_history[ns].action_id = action_collection_size; //This is the additional action code for state initializer
 action_history[ns].action_data_in = -1;
 asign[ns] = state_initializer(&(action_history[ns].action_data_in));
}

int   metropolis_step()
{
 int n_actions, iaction, todo, accepted = 0, res;
 double a, alpha;
 
 if(ns>max_recursion_depth)
 {
  logs_WriteError("Step %08i:\tOverflow detected - size of action_history (ns = %i) is larger than max_recursion_depth=%i !!! %s", step_number, ns, max_recursion_depth, (exit_upon_overflow? "Stopping the MC process" : "Resetting the state of the system"));
  ns = 0;
  if(exit_upon_overflow)
   exit(EXIT_FAILURE);
  call_state_initializer();
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
   logs_Write(3, "Action %s [id = %i]: \t ampl = %2.4E, action_data_in = %i", action_collection_name[action_list[iaction].action_id], action_list[iaction].action_id, amplitude_list[iaction], action_list[iaction].action_data_in);
 }; 
 
 nA[ns] = 0.0;
 for(iaction=0; iaction<n_actions; iaction++)
  nA[ns] += fabs(amplitude_list[iaction]);
  
 gather_mc_stat(); 
 
 //Controlling the sum of amplitudes
 if(control_max_ampl_sum && (nA[ns] > (max_ampl_sum + max_ampl_sum_tol)))
 {
  logs_WriteError("nA[%i] = %2.4E > max_ampl_sum = %2.4E", ns, nA[ns], max_ampl_sum);
  int mdata = -1;
  for(int i=0; i<n_actions; i++)
   logs_WriteParameter(0, action_collection_name[action_list[i].action_id], "A = %2.4E,\t AMax = %2.4E", amplitude_list[i], (action_collection_amplitude[action_list[i].action_id])(&mdata) );
  logs_Write(0, " ");
 }; 
 
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
   //Now calculate the normalized probabilities of different steps
   for(iaction=0; iaction<n_actions; iaction++)
    probability_list[iaction] = fabs(amplitude_list[iaction])/nA[ns];
   //Choose at random which action to perform
   todo = rand_choice(probability_list, n_actions);
   //Perform the selected action
   logs_Write(2, "Action %i [%s, id = %i] selected by rand_choice...", todo, action_collection_name[action_list[todo].action_id], action_list[todo].action_id);
   if(action_list[todo].action_id<action_collection_size)
   {
    res = (action_collection_do[action_list[todo].action_id])(&(action_list[todo].action_data_in));
    if(abs(res)==ACTION_SUCCESS)
    {
     //Write to the log
     logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t%s Performed action % 10s [id = % 2i], action_data_in = % 3i, ns = % 5i, nA[ns] = %2.4lf", step_number, (ansi_colors? "\x1b[5;40m" : ""), action_collection_name[action_list[todo].action_id], action_list[todo].action_id, action_list[todo].action_data_in, ns, nA[ns]);
     //Increase the order counter 
     ns++;
     //...save it to the history
     action_history[ns] = action_list[todo];
     //...add this choice to the statistics
     action_counter[action_list[todo].action_id]++;
     //Modify the reweighting sign/phase
     asign[ns] = asign[ns-1]*SIGN(amplitude_list[todo]);         
    };
    if(res==ERR_HISTORY_OVERFLOW || res==ERR_STACK_OVERFLOW)
    {
     logs_WriteError("Action %s [id = %i] with action_data_in %i caused %s overflow at mc step %i!!! %s!!!", action_collection_name[action_list[todo].action_id], action_list[todo].action_id, action_list[todo].action_data_in, (res==ERR_HISTORY_OVERFLOW? "history" : "state"), step_number, (exit_upon_overflow? "Stopping the MC process" : "Resetting the state of the system"));
     if(exit_upon_overflow)
      exit(EXIT_FAILURE);
     ns = 0;
     call_state_initializer();
    };
    if(res==ERR_WRONG_STATE)
     logs_WriteError("Attempted to apply action %s [id = %i] with action_data_in %i to the improper system state, nothing was really done...", action_collection_name[action_list[todo].action_id], action_list[todo].action_id, action_list[todo].action_data_in); 
   } 
   else
   {
    logs_WriteError("action_list[%i].action_id = %i is larger than the total number %i of actions in the collection", todo, action_list[todo].action_id, action_collection_size); 
   };
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
    logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t%s Undoing the action %s [id = %i], action_data_in = %i, ns = %i, nA[ns-1] = %2.4lf", step_number, (ansi_colors? "\x1b[5;40m" : ""), action_collection_name[action_history[ns].action_id], action_history[ns].action_id, action_history[ns].action_data_in, ns, nA[ns-1]);
    if(action_history[ns].action_id<action_collection_size)
    {
     res = (action_collection_undo[action_history[ns].action_id])(&(action_history[ns].action_data_in));
     if(res==ERR_HISTORY_OVERFLOW || res==ERR_STACK_OVERFLOW)
     {
      logs_WriteError("Undoing action %s [id = %i] with action_data_in %i caused %s overflow at mc step %i!!! %s!!!", action_collection_name[action_history[ns].action_id], action_history[ns].action_id, action_history[ns].action_data_in, (res==ERR_HISTORY_OVERFLOW? "history" : "state"), step_number, (exit_upon_overflow? "Stopping the MC process" : "Resetting the state of the system"));
      if(exit_upon_overflow)
       exit(EXIT_FAILURE);
      ns = 0;
      call_state_initializer();
     };
     if(res==ERR_WRONG_STATE)
      logs_WriteError("Attempted to undo action %s [id = %i] with action_data_in %i from the improper system state, nothing was really done...", action_collection_name[action_history[ns].action_id], action_history[ns].action_id, action_history[ns].action_data_in);
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
   call_state_initializer();
   //Write to the log
   logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t%s state_initializer has been just called", step_number, (ansi_colors? "\x1b[5;40m" : ""));
  };
 };
 
 //Updating statistics of the MC process - the mean nA, its dispersion and the mean sign

 aac   += accepted;
 step_number ++;
 
 return accepted; //Return 1 if some new state is accepted, otherwise 0
}

void  print_action_history()
{
 logs_Write(0, "\nAction history at step %i:", step_number);
 for(int ins=0; ins<=ns; ins++)
  if(action_history[ins].action_id>=0 && action_history[ins].action_id<=action_collection_size)
   logs_Write(0, "ins=%i, action \"%s\"(id=%i), action_data_in = %i", ins, action_collection_name[action_history[ins].action_id], action_history[ins].action_id, action_history[ins].action_data_in);
  else
   logs_WriteError("ins=%i, action id=%i is out of range...", ins, action_history[ins].action_id);
 logs_Write(0, "");      
}

