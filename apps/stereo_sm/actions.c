#include "actions.h"

t_lat_stack      X; //This stack is the current state of the system
t_lat_stack      H; //This stack will contain the data related to the sequence of actions

int* alpha_order_stack   = NULL;
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
 
 ADD_TO_ACTION_COLLECTION(            create, 0);
 ADD_TO_ACTION_COLLECTION(          add_line, 1); 
 ADD_TO_ACTION_COLLECTION(  exchange_momenta, 2);
 ADD_TO_ACTION_COLLECTION(              join, 3);
 ADD_TO_ACTION_COLLECTION(            vertex, 4);
 
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
 //Initialize the lattice stack
 init_lat_stack(&X, DIM, max_stack_nel);
 init_lat_stack(&H, DIM, max_history_nel);
 
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
DECLARE_ACTION_AMPLITUDE(add_line)
{
 return P.sigma/cc;
}

DECLARE_ACTION_DO(add_line)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 
 X.len[X.top-1] += 2; //We are adding two momenta to the topmost sequence
 X.nel          += 2;
 
 rand_momentum(&P, STACK_EL(X, 0));
 invert_momentum(STACK_EL(X, 1), STACK_EL(X, 0));
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(add_line)
{
 RETURN_IF_FALSE(               X.len[X.top-1]>2, ERR_WRONG_STATE);

 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 return ACTION_SUCCESS;
}

/*************************** Exchange momenta via the interaction vertex ******/

DECLARE_ACTION_AMPLITUDE(exchange_momenta)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1)
  return -lambda*alpha*NN*P.sigma;
 return 0.0;
}

DECLARE_ACTION_DO(exchange_momenta)
{
 RETURN_IF_FALSE(                X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(      H.nel<H.max_nel-1, ERR_HISTORY_OVERFLOW);
 
 //Combine the topmost sequences in the stack
 X.len[X.top-2] += X.len[X.top-1];
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
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
 
 //Now remember the combination of orders and signs
 alpha_order_stack[alpha_order_stop-2] += (alpha_order_stack[alpha_order_stop-1] + 1);
 alpha_order_history[alpha_order_htop]  = alpha_order_stack[alpha_order_stop-1];
 alpha_order_htop ++;
 alpha_order_stop --;
 
 sign_stack[sign_stop-2] *= -sign_stack[sign_stop-1];
 sign_history[sign_htop] = sign_stack[sign_stop-1];
 sign_htop ++;
 sign_stop --;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(exchange_momenta)
{
 RETURN_IF_FALSE(                   H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE(         H.len[H.top-1]==1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(             *(data_in)>=2, ERR_WRONG_DATA);
 RETURN_IF_FALSE( X.len[X.top-1]>(*data_in), ERR_WRONG_STATE);
 
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
 
 X.top ++;
 
 //Now restore the combination of orders
 alpha_order_stack[alpha_order_stop]    = alpha_order_history[alpha_order_htop-1];
 alpha_order_stack[alpha_order_stop-1] -= (alpha_order_history[alpha_order_htop-1] + 1);
 alpha_order_htop --;
 alpha_order_stop ++;
 
 sign_stack[sign_stop]    = sign_history[sign_htop-1];
 sign_stack[sign_stop-1] *= -sign_history[sign_htop-1];
 sign_htop --;
 sign_stop ++;
  
 return ACTION_SUCCESS;
}

/********** Join two sets of lines via an interaction vertex **************/
DECLARE_ACTION_AMPLITUDE(join)
{
 int    Q[4];
 if(data_in==NULL || (*data_in)<0)
 {
  double res = cc*NN*lambda*alpha/P.mass_sq;
  if(fabs(4.0*DIM + meff_sq) >= fabs(meff_sq))
   return res*(4.0*DIM + meff_sq);
  else
   return res*meff_sq;
 };
 
 if(X.top>1 && X.len[X.top-1]>=4)
 {
  add3momenta(Q, STACK_EL(X, 0), STACK_EL(X, X.len[X.top-1]-1), STACK_EL_PREV(X, 0));
  return cc*NN*lambda*alpha*lat_propagator(Q, P.mass_sq)*(lat_momentum_sq(STACK_EL(X, X.len[X.top-1]-1)) + meff_sq);
 };
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 RETURN_IF_FALSE(                X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(      H.nel<H.max_nel-2, ERR_HISTORY_OVERFLOW);
 
 //Save \tilde{q}_A and \tilde{p}_A to the history stack
 H.start[ H.top] = (H.top>0? H.start[H.top-1] + H.len[H.top-1] : 0);
 H.len[   H.top] = 2; //We push two momenta on the top of the stack
 H.top ++;
 H.nel += 2;
 assign_momentum(STACK_EL(H, 0), STACK_EL(     X, X.len[X.top-1]-1)); //\tilde{q}_A
 assign_momentum(STACK_EL(H, 1), STACK_EL_PREV(X,                0)); //\tilde{p}_A
 
 //Add \tilde{q}_A and \tilde{p}_A to \tilde{p}_1 to get p_1
 addto2momenta(STACK_EL(X, 0), 1, STACK_EL(X, X.len[X.top-1]-1), 1, STACK_EL_PREV(X, 0));
 
 //Combine the topmost sequences in the stack
 X.len[X.top-2] += X.len[X.top-1];
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
 X.top --;
  
 //Now we have to shift the elements in the joined sequence by two
 //in order to erase \tilde{q}_A and \tilde{p}_A
 for(int i=(*data_in)-2; i>=0; i--)
  assign_momentum(STACK_EL(X, i+2), STACK_EL(X, i));
 //And decrease the length of the topmost sequence by two
 X.len[X.top-1] -= 2; 
 X.nel          -= 2;
 
 //Now remember the combination of orders + signs
 alpha_order_stack[alpha_order_stop-2] += (alpha_order_stack[alpha_order_stop-1] + 1);
 alpha_order_history[alpha_order_htop]  = alpha_order_stack[alpha_order_stop-1];
 alpha_order_htop ++;
 alpha_order_stop --;
 
 sign_stack[sign_stop-2] *= sign_stack[sign_stop-1];
 sign_history[sign_htop]  = sign_stack[sign_stop-1];
 sign_htop ++;
 sign_stop --;
  
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join)
{
 RETURN_IF_FALSE(                    H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE(          H.len[H.top-1]==2, ERR_WRONG_STATE);
 RETURN_IF_FALSE(              *(data_in)>=2, ERR_WRONG_DATA);
 RETURN_IF_FALSE( X.len[X.top-1]>=(*data_in), ERR_WRONG_STATE);
 
 //First free up the place in X.stack for \tilde{q}_A and \tilde{p}_A
 //Increase X.len[X.top-1]
 X.len[X.top-1] += 2;
 X.nel          += 2;
 for(int i=0; i<(*data_in)-1; i++)
  assign_momentum(STACK_EL(X, i), STACK_EL(X, i+2));
 //Now restore the values of \tilde{q}_A and \tilde{p}_A 
 assign_momentum(STACK_EL(X, (*data_in)-1), STACK_EL(H, 0) );
 assign_momentum(STACK_EL(X, (*data_in)  ), STACK_EL(H, 1) );
 //Finally, restore the value of \tilde{p}_1 by subtracting \tilde{q}_A and \tilde{p}_A
 addto2momenta(STACK_EL(X, 0), -1, STACK_EL(H, 0), -1, STACK_EL(H, 1));
 //Free up the H stack
 H.top --;
 H.nel -= 2;
 
 //And now split the sequence back into two
 X.top ++;
 X.len[X.top-1]   = (*data_in);
 X.len[X.top-2]  -= (*data_in);
 X.start[X.top-1] = X.start[X.top-2] + X.len[X.top-2];
 
 //Now restore the combination of orders
 alpha_order_stack[alpha_order_stop]    = alpha_order_history[alpha_order_htop-1];
 alpha_order_stack[alpha_order_stop-1] -= (alpha_order_history[alpha_order_htop-1] + 1);
 alpha_order_htop --;
 alpha_order_stop ++;
 
 sign_stack[sign_stop]    = sign_history[sign_htop-1];
 sign_stack[sign_stop-1] *= sign_history[sign_htop-1];
 sign_htop --;
 sign_stop ++;
  
 return ACTION_SUCCESS;
}

/*************************** Create new vertex ****************************/
DECLARE_ACTION_AMPLITUDE(vertex)
{
 int    Q[4];
 if(data_in==NULL || (*data_in)<0)
 {
  double res = lambda*alpha*cc/P.mass_sq;
  if(fabs(4.0*DIM + meff_sq) >= fabs(meff_sq))
   return res*(4.0*DIM + meff_sq);
  else
   return res*meff_sq;
 };
 
 if(X.len[X.top-1]>=3)
 {
  add3momenta(Q, STACK_EL(X, 0), STACK_EL(X, 1), STACK_EL(X, 2));
  return lambda*alpha*cc*lat_propagator(Q, P.mass_sq)*(lat_momentum_sq(STACK_EL(X, 1)) + meff_sq);
 };
 
 return 0.0;
}

DECLARE_ACTION_DO(vertex)
{
 RETURN_IF_FALSE( X.len[X.top-1]>=3, ERR_WRONG_STATE);
 RETURN_IF_FALSE( H.nel<H.max_nel-2, ERR_HISTORY_OVERFLOW);
 
 //Push the momenta which are being joined into the H(istory)stack
 H.start[ H.top] = (H.top>0? H.start[H.top-1] + H.len[H.top-1] : 0);
 H.len[   H.top] = 2; //We push a pair of momenta on the top of the stack
 H.top ++;
 H.nel += 2;
 
 assign_momentum(STACK_EL(H, 0), STACK_EL(X, 0)); //\tilde{p}_1
 assign_momentum(STACK_EL(H, 1), STACK_EL(X, 1)); //\tilde{q}_1
  
 //Now modify the stack
 addto2momenta(STACK_EL(X, 2), +1, STACK_EL(X, 1), +1, STACK_EL(X, 0)); //p2 -> \tilde{p}2 + \tilde{p}1 + \tilde{q}1
 
 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 alpha_order_stack[alpha_order_stop-1] ++;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(vertex)
{
 RETURN_IF_FALSE(           H.top>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE( H.len[H.top-1]==2, ERR_WRONG_STATE);
 
 //Increase again the size of the topmost sequence
 X.len[X.top-1] += 2;
 X.nel          += 2;
 
 //Pop the momenta incoming to the vertex from the H(istory)stack
 assign_momentum(STACK_EL(X, 0), STACK_EL(H, 0));
 assign_momentum(STACK_EL(X, 1), STACK_EL(H, 1));
 H.top --;
 H.nel -= 2;

 addto2momenta(STACK_EL(X, 2), -1, STACK_EL(X, 0), -1, STACK_EL(X, 1)); //p2 -> p2 - p1 - q1
 
 alpha_order_stack[alpha_order_stop-1] --;

 return ACTION_SUCCESS;
}

/************** Action fetcher *****************/
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length)
{
 int nact = 0, adata = 0;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = % 5i, X.nel = % 5i, H.top = % 5i, H.nel = % 5i, alpha_order = % 3i", step_number, X.top, X.nel, H.top, H.nel, alpha_order_stack[alpha_order_stop-1]);
 
 if(check_stack)
 {
  check_stack_consistency(    &X, "X");
  check_stack_consistency(    &H, "H");
  check_momentum_conservation(&X, "X");
 }; 
 
 FETCH_ACTION(            create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(          add_line, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(  exchange_momenta, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(              join, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(            vertex, 4, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}
