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
 ADD_TO_ACTION_COLLECTION(          join, 2);
 ADD_TO_ACTION_COLLECTION( evolve_vertex, 3); 
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
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
 return 1.0/(NN*cc);
}

DECLARE_ACTION_DO(create)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 if(data_in == NULL)
 {
  X.top = 0;
  X.nel = 0;
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
 double amp = 0.0;
 if(data_in==NULL || (*data_in)<0) //Return the maximal possible amplitude of this action
 {
  amp = 2.0*cc*SQR(alpha_wc)/(lambda*(1.0 - cc*alpha_wc))*(4.0 + lambda + 4.0/(1.0 - cc*alpha_wc));
 }
 else
  if(X.len[X.top-1]>(*data_in))  //Now we think that elements in the stack are really like momenta, data_in contains the number of momenta which we need to join - 2, 4, 6, ... in this model
  {
   amp = 2.0*alpha_wc*pow(alpha_wc*cc, 0.5*(double)(*data_in))*((double)(2*(*data_in) + 4) + lambda)/lambda;
   amp *= (((*data_in)/2)%2==0? -1.0 : 1.0); //The sign of the amplitude
  }; 
 return amp;
}

DECLARE_ACTION_DO(evolve_vertex)
{
 RETURN_IF_FALSE(X.len[X.top-1]>(*data_in),  ERR_WRONG_STATE);
 RETURN_IF_FALSE(H.nel<H.max_nel-(*data_in), ERR_HISTORY_OVERFLOW);
 
 //Push the momenta which are being joined into the H(istory)stack
 H.start[ H.top] = (H.top>0? H.start[H.top-1] + H.len[H.top-1] : 0);
 H.len[   H.top] = (*data_in); //We push a pair of momenta on the top of the stack
 H.top ++;
 H.nel += (*data_in);
 //Now modify the stack
 X.len[X.top-1] -= (*data_in);
 X.nel          -= (*data_in);
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(evolve_vertex)
{
 RETURN_IF_FALSE(H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE(H.len[H.top-1]==(*data_in), ERR_WRONG_STATE);
 
 //Pop the momenta incoming to the vertex from the H(istory)stack
 H.top --;
 H.nel -= (*data_in);
 //And add them back to the topmost sequence in the stack
 X.len[X.top-1] += (*data_in);
 X.nel          += (*data_in);
 
 return ACTION_SUCCESS;
}

//Join two sets of lines
DECLARE_ACTION_AMPLITUDE(join)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1) //We can join two sequences if there are more than two elements in the stack
  return NN/cc;
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
 int nact = 0, adata = 0, iact;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = %i, X.nel = %i, H.top = %i, H.nel = %i", step_number, X.top, X.nel, H.top, H.nel);
 
     check_stack_consistency(&X, "X");
     check_stack_consistency(&H, "H");
 check_momentum_conservation(&X, "X");
 
 FETCH_ACTION(           create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(      evolve_line, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(             join, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 for(iact=0; iact<(X.len[X.top-1]-2)/2; iact++)
 {
  adata = 2*(iact + 1);
  FETCH_ACTION(   evolve_vertex, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 };
 
 return nact;
}

double f_max_ampl_sum()
{
 int adata = -1;
 double ampl_sum = fabs(action_create_amplitude(&adata))        + 
                   fabs(action_evolve_line_amplitude(&adata))   +
                   fabs(action_join_amplitude(&adata))          + 
                   fabs(action_evolve_vertex_amplitude(&adata));
 return ampl_sum;
}
