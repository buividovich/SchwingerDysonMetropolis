#include "actions.h"

t_lat_stack      X; //This stack is the current state of the system
t_lat_stack      H; //This stack will contain the data related to the sequence of actions

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
 action_collection_size = 4;
 SAFE_MALLOC( action_collection_do,        t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_undo,      t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_amplitude, t_action_amplitude, action_collection_size);
 SAFE_MALLOC( action_collection_name,      char*,              action_collection_size);
 
 ADD_TO_ACTION_COLLECTION(            create, 0);
 ADD_TO_ACTION_COLLECTION(          add_line, 1); 
 ADD_TO_ACTION_COLLECTION(              join, 2);
 ADD_TO_ACTION_COLLECTION(            vertex, 3);
 
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
 return 2.0*P.sigma/cc;;
}

DECLARE_ACTION_DO(add_line) 
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 
 X.len[X.top-1] += 2; //We are adding two momenta to the topmost sequence
 X.nel          += 2;
 
 (*data_in) = rand_int(2); //1 is for prepending to the beginning of the sequence, 0 is to the beginning and to the end

 if((*data_in)==0) //Shift the momentum sequence
  for(int i=1; i<=X.len[X.top-1]-2; i++)
   assign_momentum(STACK_EL(X, i), STACK_EL(X, i+1) );
 
 rand_momentum(&P, STACK_EL(X, 0));
 if((*data_in)==1)
  invert_momentum(STACK_EL(X, 1), STACK_EL(X, 0))
 else
  invert_momentum(STACK_EL(X, X.len[X.top-1]-1), STACK_EL(X, 0)); 
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(add_line)
{
 RETURN_IF_FALSE(               X.len[X.top-1]>2, ERR_WRONG_STATE);
 RETURN_IF_FALSE( (*data_in)==0 || (*data_in)==1, ERR_WRONG_DATA);
 
 if((*data_in)==0) //Back-shifting the sequence of momenta
  for(int i=X.len[X.top-1]-1; i>=2; i--)
   assign_momentum(STACK_EL(X, i), STACK_EL(X, i-1) );
   
 X.len[X.top-1] -= 2;
 X.nel          -= 2;
 
 return ACTION_SUCCESS;
}

/************************** Join two sets of lines **********************/
DECLARE_ACTION_AMPLITUDE(join)
{
 if(data_in==NULL || (*data_in)<0 || X.top>1) //We can join two sequences if there are more than two elements in the stack
  return NN*P.sigma/cc;
 return 0.0;
}

DECLARE_ACTION_DO(join)
{
 RETURN_IF_FALSE(                          X.top>1, ERR_WRONG_STATE);
 RETURN_IF_FALSE(                X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 RETURN_IF_FALSE( alpha_order_htop<max_history_nel, ERR_HISTORY_OVERFLOW);
 RETURN_IF_FALSE(              alpha_order_htop>=0, ERR_WRONG_STATE);
 RETURN_IF_FALSE( alpha_order_htop<max_history_nel, ERR_HISTORY_OVERFLOW);
 RETURN_IF_FALSE(              alpha_order_htop>=0, ERR_WRONG_STATE);
 
 alpha_order_stack[X.top-2]           += alpha_order_stack[X.top-1];
 alpha_order_history[alpha_order_htop] = alpha_order_stack[X.top-1];
 alpha_order_htop ++;
 
 sign_stack[X.top-2]     *= sign_stack[X.top-1];
 sign_history[sign_htop] =  sign_stack[X.top-1];
 sign_htop ++;
 
 X.len[X.top-2] += (X.len[X.top-1] + 2);
 //... and now we have to remember what was the length of both sequences in order to perform undo
 //... to this end we store X.seq_length[X.stack_top-1] in action_data_in
 (*data_in) = X.len[X.top-1];
 X.top --;
 X.nel += 2;
 
 //Shift the momenta in the ex-topmost sequence
 for(int i=1; i<=(*data_in); i++)
  assign_momentum(STACK_EL(X, i), STACK_EL(X,i+1));
 //Generate the momenta on the new legs 
 rand_momentum(&P, STACK_EL(X, 0));
 invert_momentum(STACK_EL(X, (*data_in)+1), STACK_EL(X, 0)); 
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(join)
{
 RETURN_IF_FALSE(                        (*data_in) >= 2, ERR_WRONG_DATA);
 RETURN_IF_FALSE(X.len[X.top-1] >= ((*data_in) + 2) +  2, ERR_WRONG_STATE);
 RETURN_IF_FALSE(                     alpha_order_htop>0, ERR_WRONG_STATE);
 RETURN_IF_FALSE(                            sign_htop>0, ERR_WRONG_STATE);
 
 //Shifting the elements back
 for(int i=(*data_in)+1; i>=2; i--)
  assign_momentum(STACK_EL(X, i), STACK_EL(X,i-1));
 
 X.len[X.top-1] -= ((*data_in) + 2);
 X.len[X.top]    = (*data_in);
 X.start[X.top]  = X.start[X.top-1] + X.len[X.top-1];
 
 alpha_order_stack[X.top]    = alpha_order_history[alpha_order_htop-1];
 alpha_order_stack[X.top-1] -= alpha_order_history[alpha_order_htop-1];
        sign_stack[X.top]    =        sign_history[sign_htop-1];
        sign_stack[X.top-1] *=        sign_history[sign_htop-1];
 alpha_order_htop --;
        sign_htop --;
 
 X.top ++;
 X.nel -= 2;
 
 return ACTION_SUCCESS;
}

/*************************** Create new vertex ****************************/
DECLARE_ACTION_AMPLITUDE(vertex)
{
 int    Q[4];
 if(data_in==NULL || (*data_in)<0)
 {
  //TODO: change here
  double res = alpha*cc/P.mass_sq;
  if(fabs(4.0*DIM + meff_sq) >= fabs(meff_sq))
   return res*(4.0*DIM + meff_sq);
  else
   return res*meff_sq;
 };
 
 if(adata>=3)
 {
  //TODO: here we'll need the sum of all momenta
  //and the vertex amplitude
  double res = vertex(&(STACK_EL(X, 0), Q, adata)
  return alpha*cc*lat_propagator(Q, P.mass_sq)*res;
 };
 
 return 0.0;
}

DECLARE_ACTION_DO(vertex)
{
 RETURN_IF_FALSE( X.len[X.top-1]>=4, ERR_WRONG_STATE);
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
 int nact = 0, adata = 0, iact;
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
 FETCH_ACTION(              join, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 for(adata=3; adata<X.len[X.top-1]; adata+=2) //adata contains the number of momenta to be joined
  FETCH_ACTION(           vertex, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}

//TODO [not urgent]: given the max. number of momenta to be contracted, can we calculate all the vertex function at the cost of max.? Will this speed up?
//TODO [not urgent]: can we speed up the calculation by a factor of two using momentum conservation?
/************* Vertex functions **************/
double vertex(int** P, int* Pt, int n)
{
 double res = 0.0;
 for(int m=0; m<n; m++)
 {
  int Qt[4] = {0, 0, 0, 0};
  int s = +1;
  for(int l=0; l<n-m; l++)
  {
   addto_momentum(Qt, +1, P[m+l])
   res += s*lat_momentum_sq(Qt);
   s *= -1;
  };
  if(m==0)
   assign_momentum(Pt, Qt);
 };
 return res;      
}
