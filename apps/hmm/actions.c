#include "actions.h"

t_lat_stack X; //This stack is the current state of the system
t_lat_stack H; //This stack will contain the data related to the sequence of actions

void init_actions()
{
 //Initialize the collection of actions to be used in MC process
 action_collection_size = 4;
 SAFE_MALLOC(  action_collection_do, t_action, action_collection_size);
 SAFE_MALLOC(action_collection_undo, t_action, action_collection_size);
 
 ADD_TO_ACTION_COLLECTION(        create, 0);
 ADD_TO_ACTION_COLLECTION(   evolve_line, 1);
 ADD_TO_ACTION_COLLECTION( evolve_vertex, 2);
 ADD_TO_ACTION_COLLECTION(          join, 3);
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
 //Initialize the lattice stack
 init_lat_stack(&X, DIM, max_stack_nel);
 init_lat_stack(&H, DIM, max_history_nel);
}

void free_actions()
{
 free_lat_stack(&X);
 free_lat_stack(&H);
 
 SAFE_FREE(action_collection_do);
 SAFE_FREE(action_collection_undo);
}

/************* Create new factorized-out line ****************/
DECLARE_ACTION_AMPLITUDE(create)
{
 return 1.0/NN;
}

DECLARE_ACTION_DO(create)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 if(data_in == NULL)
 {
  X.top = 0; //If called with NULL, should completely reset the state
  X.nel = 0;
  H.top = 0;
  H.nel = 0;
 };  

 X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
 X.len[   X.top] = 2; //We push a pair of momenta on the top of the stack
 X.top ++;
 X.nel += 2;
 
 return ACTION_SUCCESS; 
}

DECLARE_ACTION_UNDO(create)
{
 RETURN_IF_FALSE(X.top>0, ERR_WRONG_STATE);
 
 X.top --;
 X.nel -= 2;
 
 return ACTION_SUCCESS;
}

/************************ Create new factorized-in line *****************/
DECLARE_ACTION_AMPLITUDE(evolve_line)
{
 return 2.0/cc;
}

DECLARE_ACTION_DO(evolve_line)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);

 X.len[X.top-1] += 2;
 X.nel          += 2;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(evolve_line)
{
 RETURN_IF_FALSE(X.len[X.top-1]>2, ERR_WRONG_STATE);

 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 return ACTION_SUCCESS;
}

/*************************** Create new vertex ****************************/
DECLARE_ACTION_AMPLITUDE(evolve_vertex)
{
 if(data_in==NULL || (*data_in)<0 || X.len[X.top-1]>=3)  //Now we think that elements in the stack are really like momenta
  return lambda*cc;
 return 0.0;
}

DECLARE_ACTION_DO(evolve_vertex)
{
 RETURN_IF_FALSE(X.len[X.top-1]>=3, ERR_WRONG_STATE);
 RETURN_IF_FALSE(H.nel<H.max_nel-2, ERR_HISTORY_OVERFLOW);
 //Push the momenta which are being joined into the H(istory)stack
 H.start[ H.top] = (H.top>0? H.start[H.top-1] + H.len[H.top-1] : 0);
 H.len[   H.top] = 2; //We push a pair of momenta on the top of the stack
 H.top ++;
 H.nel += 2;
 //Now modify the stack
 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(evolve_vertex)
{
 RETURN_IF_FALSE(H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE(H.len[H.top-1]>=2, ERR_WRONG_STATE);
 //Pop the momenta incoming to the vertex from the H(istory)stack
 H.top --;
 H.nel -= 2;
 //And add them back to the topmost sequence in the stack
 X.len[X.top-1] += 2;
 X.nel          += 2;
 
 return ACTION_SUCCESS;
}

//Join two sets of lines
DECLARE_ACTION_AMPLITUDE(join)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1) //We can join two sequences if there are more than two elements in the stack
  return NN/SQR(cc);
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 RETURN_IF_FALSE(X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 
 X.len[X.top-2] += (X.len[X.top-1] + 2);
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
 
 X.top --;
 X.nel += 2;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join)
{
 RETURN_IF_FALSE(X.len[X.top-1] > ((*data_in) + 2), ERR_WRONG_STATE);
 
 X.len[X.top-1] -= ((*data_in) + 2);
 X.len[X.top]    = (*data_in);
 X.start[X.top]  = X.start[X.top-1] + X.len[X.top-1];
 X.top ++;
 X.nel -= 2;
 
 return ACTION_SUCCESS;
}

/************** Action fetcher *****************/
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length)
{
 int nact = 0, adata = 0;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = %i, X.nel = %i, H.top = %i, H.nel = %i", step_number, X.top, X.nel, H.top, H.nel);
 
 check_stack_consistency(&X, "X");
 check_stack_consistency(&H, "H");
 
 FETCH_ACTION(         create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(    evolve_line, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(  evolve_vertex, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(           join, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}

//Max. sum of all amplitudes - for the tuning of cc and NN
double f_max_ampl_sum()
{
 int adata = -1;
 double ampl_sum = fabs(action_create_amplitude(&adata))        + 
                   fabs(action_evolve_line_amplitude(&adata))   +
                   fabs(action_evolve_vertex_amplitude(&adata)) +
                   fabs(action_join_amplitude(&adata)); 
                   
 return ampl_sum;
}
