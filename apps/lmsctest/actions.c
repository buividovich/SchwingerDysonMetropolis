#include "actions.h"

t_lat_stack      X; //This stack is the current state of the system

void init_actions()
{
 //Initialize the collection of actions to be used in MC process
 action_collection_size = 3;
 SAFE_MALLOC( action_collection_do,        t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_undo,      t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_amplitude, t_action_amplitude, action_collection_size);
 SAFE_MALLOC( action_collection_name,      char*,              action_collection_size);
 
 ADD_TO_ACTION_COLLECTION(            create, 0);
 ADD_TO_ACTION_COLLECTION(          increase, 1); 
 ADD_TO_ACTION_COLLECTION(              join, 2);
 
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
 //Initialize the lattice stack
 init_lat_stack(&X, DIM, max_stack_nel);
}

void free_actions()
{
 free_lat_stack(&X);
 
 SAFE_FREE(action_collection_do);
 SAFE_FREE(action_collection_undo);
 SAFE_FREE(action_collection_amplitude);
 for(int i=0; i<=action_collection_size; i++)
  SAFE_FREE(action_collection_name[i]);
 SAFE_FREE(action_collection_name); 
}

/************* Create new factorized-out line ****************/
DECLARE_ACTION_AMPLITUDE(create)
{
 return create_amplitudes[model];
}

DECLARE_ACTION_DO(create)
{
 if(data_in == NULL)
 {
  X.top  = 0; //If called with NULL, should completely reset the state
  X.nel  = 0;
 };
 RETURN_IF_FALSE(X.nel<X.max_nel-1, ERR_STACK_OVERFLOW);  

 X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
 X.len[   X.top] = 1; //We push a pair of momenta on the top of the stack
 X.top ++;
 X.nel += 1;
  
 return ACTION_SUCCESS; 
}

DECLARE_ACTION_UNDO(create)
{
 RETURN_IF_FALSE(X.top>0, ERR_WRONG_STATE);
 
 X.top --;
 X.nel -= 1;
 
 return ACTION_SUCCESS;
}

/************************ Create new factorized-in line *****************/
DECLARE_ACTION_AMPLITUDE(increase)
{
 return increase_amplitudes[model];
}

DECLARE_ACTION_DO(increase)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-1, ERR_STACK_OVERFLOW);
 
 X.len[X.top-1] += 1; //We are adding two momenta to the topmost sequence
 X.nel          += 1;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(increase)
{
 RETURN_IF_FALSE(               X.len[X.top-1]>1, ERR_WRONG_STATE);

 X.len[X.top-1] -= 1;
 X.nel          -= 1;
 
 return ACTION_SUCCESS;
}

/*************************** Join ******/

DECLARE_ACTION_AMPLITUDE(join)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1)
 {
  return join_amplitudes[model];
 }; 
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 RETURN_IF_FALSE(                X.top>1, ERR_WRONG_STATE);
 
 //Combine the topmost sequences in the stack
 X.len[X.top-2] += X.len[X.top-1];
 if(model==FREE_GAUSS || model==LM_POSITIVE)
 {
  X.len[X.top-2] += 1;
  X.nel          += 1;
 }; 
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
 X.top --;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join)
{
 RETURN_IF_FALSE(             *(data_in)>=1, ERR_WRONG_DATA);
 RETURN_IF_FALSE( X.len[X.top-1]>(*data_in), ERR_WRONG_STATE);
 if(model==FREE_GAUSS || model==LM_POSITIVE)
  RETURN_IF_FALSE( X.len[X.top-1]>(*data_in)+1, ERR_WRONG_STATE);
  
 //And now split the sequences
 if(model==FREE_GAUSS || model==LM_POSITIVE)
 {
  X.len[X.top-1] -= 1;
  X.nel          -= 1;
 }; 
 X.len[X.top-1] -= (*data_in);
 X.start[X.top] = X.start[X.top-1] + X.len[X.top-1];
 X.len[X.top]   = (*data_in);
 
 X.top ++;
 
 return ACTION_SUCCESS;
}

/************** Action fetcher *****************/
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length)
{
 int nact = 0, adata = 0;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = %i, X.nel = %i", step_number, X.top, X.nel);
 
 if(check_stack)
  check_stack_consistency(    &X, "X");
 
 FETCH_ACTION(            create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(          increase, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(              join, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}
