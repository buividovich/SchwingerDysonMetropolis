#include "actions.h"

t_lat_coordinate_stack      X; //This stack is the current state of the system

int  beta_order         = 0; 

void init_actions()
{
 //Initialize the collection of actions to be used in MC process
 action_collection_size = (resummation? 5 : 6);
 SAFE_MALLOC( action_collection_do,        t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_undo,      t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_amplitude, t_action_amplitude, action_collection_size);
 SAFE_MALLOC( action_collection_name,      char*,              action_collection_size);
 
 ADD_TO_ACTION_COLLECTION(            create, 0);
 ADD_TO_ACTION_COLLECTION(          add_line, 1); 
 ADD_TO_ACTION_COLLECTION(              join, 2);
 ADD_TO_ACTION_COLLECTION(      join_contact, 3);
 ADD_TO_ACTION_COLLECTION(            vertex, 4);
 if(!resummation)
  ADD_TO_ACTION_COLLECTION(      random_walk, 5);
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
 //Initialize the lattice stack with sequences
 init_lat_coordinate_stack(&X, 2*(max_order+1));
}

void free_actions()
{
 free_lat_coordinate_stack(&X);

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
 int new_eff_order = (X.nel/2 + beta_order - 1) + 1;
 if(new_eff_order <= max_order)
  return 1.0/(NN*cc);
 return 0.0;
}

DECLARE_ACTION_DO(create)
{
 if((data_in==NULL) || ((*data_in)<0))
 {
  X.top  = 0; //If called with NULL, should completely reset the state
  X.nel  = 0;
  beta_order     = 0;
 };
 
 RETURN_IF_FALSE(X.nel<=X.max_nel-2, ERR_STACK_OVERFLOW);  

 X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
 X.len[   X.top] = 2; //We push a pair of momenta on the top of the stack
 X.top ++;
 X.nel += 2;
 
 STACK_EL(X, 0) = rand_int(lat_vol);
 STACK_EL(X, 1) = STACK_EL(X, 0);
 
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
DECLARE_ACTION_AMPLITUDE(add_line)
{
 int new_eff_order = (X.nel/2 + beta_order - 1) + 1;
 if(new_eff_order <= max_order)
  return 2.0/cc;
 return 0.0; 
}

DECLARE_ACTION_DO(add_line) 
{
 RETURN_IF_FALSE(X.nel<=X.max_nel-2, ERR_STACK_OVERFLOW);
 
 X.len[X.top-1] += 2; //We are adding two momenta to the topmost sequence
 X.nel          += 2;
 
 (*data_in) = rand_int(2); //1 is for prepending to the beginning of the sequence, 0 is to the beginning and to the end

 if((*data_in)==0)
 {
  //Shift the momentum sequence
  for(int i=0; i<X.len[X.top-1]-2; i++)
   STACK_EL(X, i) = STACK_EL(X, i+2);
  STACK_EL(X, X.len[X.top-1]-2) = STACK_EL(X, 0);
 };  

 STACK_EL(X, 0) = rand_int(lat_vol);
 
 if((*data_in)==1)
  STACK_EL(X, 1) = STACK_EL(X, 0);
 
 if((*data_in)==0)
  STACK_EL(X, X.len[X.top-1]-1) = STACK_EL(X, 0); 
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(add_line)
{
 RETURN_IF_FALSE(               X.len[X.top-1]>2, ERR_WRONG_STATE);
 RETURN_IF_FALSE( (*data_in)==0 || (*data_in)==1, ERR_WRONG_DATA);
 
 if((*data_in)==0)
 {
  //Back-shifting the sequence of momenta
  STACK_EL(X,0) = STACK_EL(X, X.len[X.top-1]-2);
  for(int i=X.len[X.top-1]-1; i>=2; i--)
   STACK_EL(X, i) = STACK_EL(X, i-2);
 };
   
 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 return ACTION_SUCCESS;
}

/************************** Join two sets of lines **********************/
DECLARE_ACTION_AMPLITUDE(join)
{
 if(X.top>1) //We can join two sequences if there are more than two elements in the stack
 {
  int new_eff_order = (X.nel/2 + beta_order - 1) + 1;
  if(new_eff_order <= max_order)
   return NN/cc;
 };
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 RETURN_IF_FALSE(                          X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(               X.nel<=X.max_nel-2, ERR_STACK_OVERFLOW);
 
 X.len[X.top-2] += (X.len[X.top-1] + 2);
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
 X.top --;
 X.nel += 2;
 
 //Shift and Exchange the momenta in the ex-topmost sequence
 for(int i=0; i<(*data_in); i++)
  STACK_EL(X, i) = STACK_EL(X,i+2);
 STACK_EL(X, (*data_in)) = STACK_EL(X, 0);
 STACK_EL(X, 0) = rand_int(lat_vol);
 STACK_EL(X, (*data_in)+1) = STACK_EL(X, 0);
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join)
{
 RETURN_IF_FALSE(                        (*data_in) >= 2, ERR_WRONG_DATA);
 RETURN_IF_FALSE(X.len[X.top-1] >= ((*data_in) + 2) +  2, ERR_WRONG_STATE);
 
 //First bringing the topmost sequence into the form which is easy to split = {_, _, x1, y1, ..., p_m, q_m, pt_1, qt_1, ..., pt_n, qt_n}
 STACK_EL(X, 0) = STACK_EL(X, (*data_in)); //now we have {x1, y1, ..., x_m, y_m, _, _, ...}
 
 //Shifting by two elements
 for(int i=(*data_in)+1; i>=2; i--)
  STACK_EL(X, i) = STACK_EL(X,i-2);
 
 X.len[X.top-1] -= ((*data_in) + 2);
 X.len[X.top]    = (*data_in);
 X.start[X.top]  = X.start[X.top-1] + X.len[X.top-1];
 
 X.top ++;
 X.nel -= 2;
 
 return ACTION_SUCCESS;
}

/************************** Join two sets of lines **********************/
DECLARE_ACTION_AMPLITUDE(join_contact)
{
 if(X.top>1 && STACK_EL(X, 0) == STACK_EL_PREV(X, 0) ) //We can join two sequences if there are more than two elements in the stack
  return -1.0*NN;
 return 0.0;
}

DECLARE_ACTION_DO(join_contact)
{
 RETURN_IF_FALSE(                          X.top>1, ERR_WRONG_STATE);
 
 //Combine the topmost sequences in the stack
 X.len[X.top-2] += X.len[X.top-1];
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 //... Diagram orders are also stored
 (*data_in) = X.len[X.top-1];
 
 X.top --;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join_contact)
{
 RETURN_IF_FALSE(                   (*data_in) >= 2, ERR_WRONG_DATA);
 RETURN_IF_FALSE(X.len[X.top-1] >= ((*data_in) + 2), ERR_WRONG_STATE);
  
 //And now split the sequences
 X.len[X.top-1] -= (*data_in);
 X.start[X.top] = X.start[X.top-1] + X.len[X.top-1];
 X.len[X.top]   = (*data_in);
 
 X.top ++;
 
 return ACTION_SUCCESS;
}


/*************************** Create new vertex ****************************/
DECLARE_ACTION_AMPLITUDE(vertex)
{
 if(X.len[X.top-1]>=4 && (STACK_EL(X, 0)==STACK_EL(X, 2)))
 {
  (*data_in) = are_neighbors(STACK_EL(X, 0), STACK_EL(X, 1));
  if((*data_in)>=0) 
   return -cc*alpha*lat_vol;
 }; 
 return 0.0;
}

DECLARE_ACTION_DO(vertex)
{
 RETURN_IF_FALSE(                 X.len[X.top-1]>=4, ERR_WRONG_STATE);
 RETURN_IF_FALSE( (*data_in)>=0 && (*data_in)<2*DIM, ERR_WRONG_DATA);
  
 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 beta_order ++;
  
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(vertex)
{
 RETURN_IF_FALSE( (*data_in)>=0 && (*data_in)<2*DIM, ERR_WRONG_DATA);
 
 X.len[X.top-1] += 2;
 X.nel          += 2;
 
 //Restore the original sequence of coordinates
 STACK_EL(X, 0) = STACK_EL(X, 2);
 STACK_EL(X, 1) = ((*data_in)<DIM? lat_shift_fwd(STACK_EL(X, 0), (*data_in)) : lat_shift_bwd(STACK_EL(X, 0), (*data_in)-DIM) );

 beta_order --;
 
 return ACTION_SUCCESS;
}

/********************* Random walk term, somewhat trivial *********************/
DECLARE_ACTION_AMPLITUDE(random_walk)
{
 int new_eff_order = (X.nel/2 + beta_order - 1) + 1;
 if(new_eff_order<=max_order)
  return alpha*(double)(2*DIM);
 return 0.0; 
}

DECLARE_ACTION_DO(random_walk)
{
 (*data_in) = rand_int(2*DIM);
 
 if((*data_in)<DIM)
  STACK_EL(X, 0) = lat_shift_fwd(STACK_EL(X, 0), (*data_in));
 else
  STACK_EL(X, 0) = lat_shift_bwd(STACK_EL(X, 0), (*data_in)-DIM);
 beta_order ++; 
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(random_walk)
{
 RETURN_IF_FALSE( (*data_in)>=0 && (*data_in)<2*DIM, ERR_WRONG_DATA);
 
 if((*data_in)<DIM)
  STACK_EL(X, 0) = lat_shift_bwd(STACK_EL(X, 0), (*data_in));
 else
  STACK_EL(X, 0) = lat_shift_fwd(STACK_EL(X, 0), (*data_in)-DIM);
 beta_order --;

 return ACTION_SUCCESS;
}

/************** Action fetcher *****************/
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length)
{
 int nact = 0, adata = 0;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = % 5i, X.len[X.top-1] = % 5i,  X.nel = % 5i, beta_order = % 3i, ns = % 3i", step_number, X.top, X.len[X.top-1], X.nel, beta_order, ns);
 
 //if(check_stack)
 // check_stack_consistency(    &X, "X"); TODO: still implement for debugging purposes
 
 FETCH_ACTION(            create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(          add_line, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(              join, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(      join_contact, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(            vertex, 4, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(       random_walk, 5, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}
