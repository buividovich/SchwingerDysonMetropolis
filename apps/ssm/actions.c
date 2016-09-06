#include "actions.h"

t_lat_stack      X; //This stack is the current state of the system
t_lat_stack      H; //This stack will contain the data related to the sequence of actions

//Data for the formal counting of diagram orders
int* alpha_order_stack   = NULL; //This is necessary to make the full use of factorization
int  alpha_order_stop    = 0;
int* alpha_order_history = NULL;
int  alpha_order_htop    = 0;

int* sign_stack          = NULL;
int  sign_stop           = 0;
int* sign_history        = NULL;
int  sign_htop           = 0;

t_lat_propagator P;

void init_actions()
{
 //Initialize the collection of actions to be used in MC process
 action_collection_size = 5;
 SAFE_MALLOC( action_collection_do,        t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_undo,      t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_amplitude, t_action_amplitude, action_collection_size);
 SAFE_MALLOC( action_collection_name,      char*,              action_collection_size);
 
 ADD_TO_ACTION_COLLECTION(        create, 0);
 ADD_TO_ACTION_COLLECTION(   evolve_line, 1);
 ADD_TO_ACTION_COLLECTION( evolve_vertex, 2);
 ADD_TO_ACTION_COLLECTION(          join, 3);
 ADD_TO_ACTION_COLLECTION(  flip_momenta, 4);
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
 //Initialize the lattice stack
 init_lat_stack(&X, DIM, max_stack_nel);
 init_lat_stack(&H, DIM, max_history_nel);
 
 //Initializing the order counting and sign stack
 SAFE_MALLOC(  alpha_order_stack, int, max_stack_nel);
 SAFE_MALLOC(alpha_order_history, int, max_history_nel);
 alpha_order_stop     = 0;
 alpha_order_htop     = 0;
 
 SAFE_MALLOC(         sign_stack, int, max_stack_nel);
 SAFE_MALLOC(       sign_history, int, max_history_nel);
 sign_stop            = 0;
 sign_htop            = 0;
}

void free_actions()
{
 free_lat_stack(&X);
 free_lat_stack(&H);
 
 SAFE_FREE(alpha_order_stack);
 SAFE_FREE(alpha_order_history);
 SAFE_FREE(sign_stack);
 SAFE_FREE(sign_history);
 
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
 return P.sigma/(NN*cc);
}

DECLARE_ACTION_DO(create)
{
 if((data_in==NULL) || ((*data_in)<0))
 {
  X.top  = 0; //If called with NULL, should completely reset the state
  X.nel  = 0;
  H.top  = 0;
  H.nel  = 0;
  alpha_order_stop     = 0;
  alpha_order_htop     = 0;
         sign_stop     = 0;
         sign_htop     = 0;
 };
 
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);  

 X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
 X.len[   X.top] = 2; //We push a pair of momenta on the top of the stack
 X.top ++;
 X.nel += 2;
 
 rand_momentum(&P, STACK_EL(X, 0) );
 invert_momentum(  STACK_EL(X, 1), STACK_EL(X, 0));
 
 alpha_order_stack[alpha_order_stop] = 0;
 alpha_order_stop ++;
        sign_stack[       sign_stop] = +1;
        sign_stop ++;
 
 return ACTION_SUCCESS; 
}

DECLARE_ACTION_UNDO(create)
{
 RETURN_IF_FALSE(X.top>0, ERR_WRONG_STATE);
 
 X.top --;
 X.nel -= 2;
 
 alpha_order_stop --;
        sign_stop --;        
 
 return ACTION_SUCCESS;
}

/************************ Create new factorized-in line *****************/
DECLARE_ACTION_AMPLITUDE(evolve_line)
{
 return 2.0*P.sigma/cc;
}

DECLARE_ACTION_DO(evolve_line)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 
 X.len[X.top-1] += 2; //We are adding two momenta to the topmost sequence
 X.nel          += 2;
 
 (*data_in) = rand_int(2); //1 is for prepending to the beginning of the sequence, 0 is to the beginning and to the end

 if((*data_in)==0)
 {
  //Shift the momentum sequence
  for(int i=0; i<X.len[X.top-1]-2; i++)
   assign_momentum(STACK_EL(X, i), STACK_EL(X, i+2) );
  assign_momentum(STACK_EL(X, X.len[X.top-1]-2), STACK_EL(X, 0));
 };  
 
 rand_momentum(&P, STACK_EL(X, 0));
 
 if((*data_in)==1)
  invert_momentum(STACK_EL(X, 1), STACK_EL(X, 0));
 
 if((*data_in)==0)
  invert_momentum(STACK_EL(X, X.len[X.top-1]-1), STACK_EL(X, 0)); 
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(evolve_line)
{
 RETURN_IF_FALSE(               X.len[X.top-1]>2, ERR_WRONG_STATE);
 RETURN_IF_FALSE( (*data_in)==0 || (*data_in)==1, ERR_WRONG_DATA);
 
 if((*data_in)==0)
 {
  //Back-shifting the sequence of momenta
  assign_momentum(STACK_EL(X,0), STACK_EL(X, X.len[X.top-1]-2));
  for(int i=X.len[X.top-1]-1; i>=2; i--)
   assign_momentum(STACK_EL(X, i), STACK_EL(X, i-2) );
 };

 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 return ACTION_SUCCESS;
}

//TODO: introduce alpha in the interacting part

/*************************** Create new vertex ****************************/
DECLARE_ACTION_AMPLITUDE(evolve_vertex)
{
 int    Q[4];
 if(data_in==NULL || (*data_in)<0)
 {
  double res = cc*lambda/(meff_sq + lambda);
  if(fabs(4.0*DIM + meff_sq) >= fabs(meff_sq))
   return res*(4.0*DIM + meff_sq);
  else
   return res*meff_sq;
 }
 else
  if(X.len[X.top-1]>=3)
  {
   add3momenta(Q, STACK_EL(X, 0), STACK_EL(X, 1), STACK_EL(X, 2));
   return lambda*cc*lat_propagator(Q, meff_sq + lambda)*(lat_momentum_sq(STACK_EL(X, 1)) + meff_sq);
  };
 return 0.0;
}

DECLARE_ACTION_DO(evolve_vertex)
{
 RETURN_IF_FALSE( X.len[X.top-1]>=3, ERR_WRONG_STATE);
 RETURN_IF_FALSE( H.nel<H.max_nel-2, ERR_HISTORY_OVERFLOW);
 
 //Push the momenta which are being joined into the H(istory)stack
 H.start[ H.top] = (H.top>0? H.start[H.top-1] + H.len[H.top-1] : 0);
 H.len[   H.top] = 2; //We push a pair of momenta on the top of the stack
 H.top ++;
 H.nel += 2;
 
 assign_momentum(STACK_EL(H, 0), STACK_EL(X, 0)); //p_1
 assign_momentum(STACK_EL(H, 1), STACK_EL(X, 1)); //q_1
  
 //Now modify the stack
 addto2momenta(STACK_EL(X, 2), +1, STACK_EL(X, 1), +1, STACK_EL(X, 0)); //p2 -> p2 + p1 + q1
 
 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 O[X.top-1] ++;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(evolve_vertex)
{
 RETURN_IF_FALSE(           H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE( H.len[H.top-1]==2, ERR_WRONG_STATE);
 RETURN_IF_FALSE(      O[X.top-1]>0, ERR_WRONG_STATE);
 
 O[X.top-1] --;
 
 //Increase again the size of the topmost sequence
 X.len[X.top-1] += 2;
 X.nel          += 2;
 
 //Pop the momenta incoming to the vertex from the H(istory)stack
 assign_momentum(STACK_EL(X, 0), STACK_EL(H, 0));
 assign_momentum(STACK_EL(X, 1), STACK_EL(H, 1));
 H.top --;
 H.nel -= 2;

 addto2momenta(STACK_EL(X, 2), -1, STACK_EL(X, 0), -1, STACK_EL(X, 1)); //p2 -> p2 - p1 - q1

 return ACTION_SUCCESS;
}

/************************** Join two sets of lines **********************/
DECLARE_ACTION_AMPLITUDE(join)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1) //We can join two sequences if there are more than two elements in the stack
  return NN*lambda*P.sigma/cc;
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 RETURN_IF_FALSE(                 X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(       X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 RETURN_IF_FALSE(  OH_top<max_history_nel, ERR_OTHER);
 RETURN_IF_FALSE(               OH_top>=0, ERR_OTHER);
 
 O[X.top-2] += O[X.top-1];
 OH_PUSH(O[X.top-1]);
 
 X.len[X.top-2] += (X.len[X.top-1] + 2);
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
 X.top --;
 X.nel += 2;
 
 //Shift and Exchange the momenta in the ex-topmost sequence
 for(int i=0; i<(*data_in); i++)
  assign_momentum(STACK_EL(X, i), STACK_EL(X,i+2));
 assign_momentum(STACK_EL(X, (*data_in)), STACK_EL(X, 0));
 rand_momentum(&P, STACK_EL(X, 0));
 invert_momentum(STACK_EL(X, (*data_in)+1), STACK_EL(X, 0)); 
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join)
{
 RETURN_IF_FALSE(                        (*data_in) >= 2, ERR_WRONG_DATA);
 RETURN_IF_FALSE(X.len[X.top-1] >= ((*data_in) + 2) +  2, ERR_WRONG_STATE);
 RETURN_IF_FALSE(                               OH_top>0, ERR_OTHER);
  
 //First bringing the topmost sequence into the form which is easy to split = {_, _, p1, q1, ..., p_m, q_m, pt_1, qt_1, ..., pt_n, qt_n}
 assign_momentum(STACK_EL(X, 0), STACK_EL(X, (*data_in))); //now we have {p1, q1, ..., p_m, q_m, _, _, ...}
 
 //Shifting by two elements
 for(int i=(*data_in)+1; i>=2; i--)
  assign_momentum(STACK_EL(X, i), STACK_EL(X,i-2));
 
 X.len[X.top-1] -= ((*data_in) + 2);
 X.len[X.top]    = (*data_in);
 X.start[X.top]  = X.start[X.top-1] + X.len[X.top-1];
 
 OH_POP(O[X.top]);
 O[X.top-1] -= O[X.top];
 
 X.top ++;
 X.nel -= 2;
 
 return ACTION_SUCCESS;
}

/************************** Join two sets of lines and flips topmost momentum **********************/
DECLARE_ACTION_AMPLITUDE(flip_momenta)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1) //We can join two sequences if there are more than two elements in the stack
  return -1.0*lambda*NN*P.sigma;
 return 0.0;
}

DECLARE_ACTION_DO(flip_momenta)
{
 RETURN_IF_FALSE(                X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(      H.nel<H.max_nel-2, ERR_HISTORY_OVERFLOW);
 RETURN_IF_FALSE( OH_top<max_history_nel, ERR_OTHER);
 RETURN_IF_FALSE(              OH_top>=0, ERR_OTHER);
 
 //Combine the topmost sequences in the stack
 X.len[X.top-2] += X.len[X.top-1];
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 //... Diagram orders are also stored
 (*data_in) = X.len[X.top-1];
 
 O[X.top-2] += O[X.top-1];
 OH_PUSH(O[X.top-1]);
 
 X.top --;
  
 //Save p_1 to the history stack
 H.start[ H.top] = (H.top>0? H.start[H.top-1] + H.len[H.top-1] : 0);
 H.len[   H.top] = 1; //We push just one momentum on the top of the stack
 H.top ++;
 H.nel += 1;
 assign_momentum(STACK_EL(H, 0), STACK_EL(X, 0));
 
 //Now flip the momenta
 addto_momentum(STACK_EL(X, (*data_in)), +1, STACK_EL(X, 0));
 rand_momentum(&P, STACK_EL(X, 0));
 addto_momentum(STACK_EL(X, (*data_in)), -1, STACK_EL(X, 0));
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(flip_momenta)
{
 RETURN_IF_FALSE(                   H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE(         H.len[H.top-1]==1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(             *(data_in)>=2, ERR_WRONG_DATA);
 RETURN_IF_FALSE( X.len[X.top-1]>(*data_in), ERR_WRONG_STATE);
 RETURN_IF_FALSE(                  OH_top>0, ERR_OTHER);
 
 //First bringing the sequences back to the form p1, q1, ..., p_m, q_m, pt_1, pt_2, ..., qt_1, qt_2
 addto_momentum(STACK_EL(X, (*data_in)), +1, STACK_EL(X, 0));
 
 assign_momentum(STACK_EL(X, 0), STACK_EL(H, 0));
 H.top --;
 H.nel --;

 addto_momentum(STACK_EL(X, (*data_in)), -1, STACK_EL(X, 0));
 
 //And now split the sequences
 X.len[X.top-1] -= (*data_in);
 X.start[X.top] = X.start[X.top-1] + X.len[X.top-1];
 X.len[X.top]   = (*data_in);
 
 OH_POP(O[X.top]);
 O[X.top-1] -= O[X.top];
 
 X.top ++;
 
 return ACTION_SUCCESS;
}

/************** Action fetcher *****************/
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length)
{
 int nact = 0, adata = 0;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = %i, X.nel = %i, H.top = %i, H.nel = %i", step_number, X.top, X.nel, H.top, H.nel);
 
 if(check_stack)
 {
  check_stack_consistency(    &X, "X");
  check_stack_consistency(    &H, "H");
  check_momentum_conservation(&X, "X");
 }; 
 
 FETCH_ACTION(         create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(    evolve_line, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(  evolve_vertex, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(           join, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(   flip_momenta, 4, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}
