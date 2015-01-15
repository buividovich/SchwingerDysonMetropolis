#include "sd_metropolis.h"

int       overflow_count = 0;    //Counts exceptions when the sequence of elements becomes larger than allocated 
int       ns             = 0;    //Counter of the sequence depth  
double*   nA             = NULL; //Contains the probability normalization factor on every step
int*      asign          = NULL; //Reweighting signs of every configuration
t_action* action_history = NULL;

int   init_metropolis(t_state_initializer state_initializer)
{
 int res = 0;
 overflow_count = 0;
 ns             = 0;
 SAFE_MALLOC(             nA,   double, max_recursion_depth);
 SAFE_MALLOC(          asign,      int, max_recursion_depth);
 SAFE_MALLOC( action_history, t_action, max_recursion_depth);
 init_metropolis_statistics();
 asign[ns] = state_initializer();

 return res;      
}

void  free_metropolis()
{
 SAFE_FREE(nA);
 SAFE_FREE(asign);
 SAFE_FREE(action_history);
}

int   metropolis_step(t_action_fetcher action_fetcher, t_state_initializer state_initializer)
{
 int n_actions, iaction, todo, accepted = 0;
 double a, alpha;
 double *p               = NULL;
 
 t_action* action_list = NULL;
 n_actions = action_fetcher(&action_list);
 
 for(iaction=0; iaction<n_actions; iaction++)
  nA[ns] += fabs(action_list[iaction].action_amplitude);
 
 //Decide whether we move forward or backward in the order of expansion
 if(ns<max_recursion_depth-1)
 {
  a = rand_double(0.0, 1.0);
 }
 else
 {
  if(overflow_count<2)
   logs_WriteError("Overflow detected - sequence size is larger than the allocated memory!!! Attachment step will be omitted, statistical averages might be incorrect!!!");
  overflow_count ++;
  a = 1.0;
 };
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
   SAFE_MALLOC(p, double, n_actions);
   for(iaction=0; iaction<n_actions; iaction++)
    p[iaction] = fabs(action_list[iaction].action_amplitude)/nA[ns];
   //And increase the order counter 
   ns++; 
   //Choose at random which action to perform
   todo = rand_choice(p, n_actions);
   //Perform the selected action
   action_list[todo].action_do(&(action_list[todo].action_data_out), action_list[todo].action_data_in);
   //And save it to the history
   action_history[ns] = action_list[todo]; //Should we copy every member of struct?
   //Modify the reweighting sign/phase
   asign[ns] *= SIGN(action_list[todo].action_amplitude);
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
    action_history[ns].action_undo(action_history[ns].action_data_out, action_history[ns].action_data_out);
    //Finally, decrease the counter
    ns--;
   }; 
  }
  else
  {
   accepted = 1;
   asign[ns] = state_initializer();
  };
 };
 
 //Updating statistics of the MC process - the mean nA, its dispersion and the mean sign
 gather_mc_stat(accepted);
 
 return accepted; //Return 1 if some new state is accepted, otherwise 0
}

