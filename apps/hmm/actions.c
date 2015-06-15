#include "actions.h"

t_lat_stack X; //This stack is the current state of the system
t_lat_stack H; //This stack will contain the data related to the sequence of actions

double* A_choice_probs = NULL;
int     max_A_choices  = 0; 
double* a_choice_probs = NULL;
int     max_a_choices  = 0;

void init_actions()
{
 //Initialize the collection of actions to be used in the MC process
 action_collection_size = 5;
 SAFE_MALLOC( action_collection_do,        t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_undo,      t_action,           action_collection_size);
 SAFE_MALLOC( action_collection_amplitude, t_action_amplitude, action_collection_size);
 SAFE_MALLOC( action_collection_name,      char*,              action_collection_size);
 
 ADD_TO_ACTION_COLLECTION(        create, 0);
 ADD_TO_ACTION_COLLECTION(   evolve_line, 1);
 ADD_TO_ACTION_COLLECTION( evolve_vertex, 2);
 ADD_TO_ACTION_COLLECTION(          join, 3);
 ADD_TO_ACTION_COLLECTION(         split, 4);
 
 state_initializer       = &action_create_do;
 action_fetcher          = &my_action_fetcher;
 
 //Initialize the lattice stack
 init_lat_stack(&X, DIM, max_stack_nel);
 init_lat_stack(&H, DIM, max_history_nel);
 
 //Initializing the order counting
 genus = 0;
}

void free_actions()
{
 free_lat_stack(&X);
 free_lat_stack(&H);
 
 SAFE_FREE(action_collection_do);
 SAFE_FREE(action_collection_undo);
 SAFE_FREE(action_collection_amplitude);
 for(int i=0; i<action_collection_size; i++)
  SAFE_FREE(action_collection_name[i]);
 SAFE_FREE(action_collection_name); 
 
 SAFE_FREE(A_choice_probs);
 SAFE_FREE(a_choice_probs);
 
 //SAFE_FREE(genus);
 //SAFE_FREE(genus_history);
}

/************* Create new factorized-out line ****************/
DECLARE_ACTION_AMPLITUDE(create)
{                                
 double aG2  = 1.0/(NN_genus[genus]*SQR(cc_genus[genus])*f_genus[genus]);
 double aG11 = (genus<max_genus? 1.0/(SQR(NN_genus[genus+1])*SQR(cc_genus[genus+1])*f_genus[genus+1]) : 0.0);
 return (aG2 + aG11); 
}

DECLARE_ACTION_DO(create)
{
 if(data_in == NULL)
 {
  X.top = 0; //If called with NULL, should completely reset the state
  X.nel = 0;
  H.top = 0;
  H.nel = 0;
  genus = 0;
 };
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);

 //For higher genera, we can create either tr\phi^2 or \tr\phi\tr\phi, here we calculate the relative probabilities
 double aG2  = 1.0/(NN_genus[genus]*SQR(cc_genus[genus])*f_genus[genus]);
 double aG11 = (genus<max_genus? 1.0/(SQR(NN_genus[genus+1])*SQR(cc_genus[genus+1])*f_genus[genus+1]) : 0.0);
 double p[2] = {aG2/(aG2 + aG11), aG11/(aG2 + aG11)};
 
 int c = rand_choice(p, 2);
 
 if(c==0) //Create tr\phi^2
 {
  X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
  X.len[   X.top] = 2; //We push a pair of momenta on the top of the stack
  X.top ++;
  X.nel += 2;
 }
 else //Create tr\phi tr\phi
 {
  //The first trace...
  X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
  X.len[   X.top] = 1; //We push a pair of momenta on the top of the stack
  X.top ++;
  X.nel += 1;
  //And the second one
  X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
  X.len[   X.top] = 1; //We push a pair of momenta on the top of the stack
  X.top ++;
  X.nel += 1;
  //Finally, increase the genus counter
  genus ++;   
 };
 //We also have to remember which operation was done
 if(data_in!=NULL)
  (*data_in) = c;
 return ACTION_SUCCESS; 
}

DECLARE_ACTION_UNDO(create)
{
 RETURN_IF_FALSE((*data_in)==0 || (*data_in)==1, ERR_WRONG_STATE);
 if((*data_in)==0)
 {
  RETURN_IF_FALSE(          X.top>0, ERR_WRONG_STATE);
  RETURN_IF_FALSE(X.len[X.top-1]==2, ERR_WRONG_STATE);
  X.top --;
  X.nel -= 2;
 }
 else
 {
  RETURN_IF_FALSE(          X.top>1, ERR_WRONG_STATE);
  RETURN_IF_FALSE(X.len[X.top-1]==1, ERR_WRONG_STATE);
  RETURN_IF_FALSE(X.len[X.top-2]==1, ERR_WRONG_STATE);
  X.top -= 2;
  X.nel -= 2;
  genus --;
 };
 
 return ACTION_SUCCESS;
}

/************************ Create new factorized-in line *****************/
DECLARE_ACTION_AMPLITUDE(evolve_line)
{
 return 2.0/(SQR(cc_genus[genus]));
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
  return lambda*SQR(cc_genus[genus]);
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
  return NN_genus[genus]/SQR(cc_genus[genus]);
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

//Split singlet operators and increase genus
DECLARE_ACTION_AMPLITUDE(split)
{
 if(genus<max_genus)
 {
  double W = 0.0;
  if(data_in==NULL || (*data_in)<0) 
  {
   W = n2an_sup(cc_genus[genus]/cc_genus[genus+1]);
  } 
  else
  {
   for(int A=0; A<X.top; A++)
    W += 0.5*(double)((X.len[A]+1)*(X.len[A]+2));
   W *= pow(cc_genus[genus]/cc_genus[genus+1], (double)(X.nel));
   W *= pow(NN_genus[genus]/NN_genus[genus+1], (double)(X.top));
  };
  W *= f_genus[genus]/f_genus[genus+1];
  W *= 1.0/(NN_genus[genus+1]*SQR(cc_genus[genus+1]));
  return W;
 };
 return 0.0;
}

DECLARE_ACTION_DO(split)
{
 RETURN_IF_FALSE(X.nel<X.max_nel-2, ERR_STACK_OVERFLOW);
 int A, a;
 if(max_A_choices<X.top)
 {
  SAFE_REALLOC(A_choice_probs, double, X.top);
  logs_Write(1, "\n Changing max_A_choices from %i to %i at step %i...", max_A_choices, X.top, step_number);
  max_A_choices = X.top;
 };
 double sA = 0.0;
 for(A=0; A<X.top; A++)
 {
  A_choice_probs[A] = 0.5*(double)((X.len[A] + 1)*(X.len[A] + 2));
  sA += A_choice_probs[A];
 };
 for(A=0; A<X.top; A++)
  A_choice_probs[A] /= sA;
 A = rand_choice(A_choice_probs, X.top);
 (*data_in) = A;
 
 if(max_a_choices<(X.len[A]+1))
 {
  SAFE_REALLOC(a_choice_probs, double, (X.len[A]+1));
  logs_Write(1, "\n Changing max_a_choices from %i to %i at step %i...", max_a_choices, (X.len[A]+1), step_number);
  max_a_choices = (X.len[A]+1);
 };
 
 double sa = 0.0;
 for(a=1; a<=X.len[A]+1; a++)
 {
  a_choice_probs[a-1] = (double)(X.len[A] + 2 - a);
  sa += a_choice_probs[a-1];
 };
 for(a=1; a<=X.len[A]+1; a++)
  a_choice_probs[a-1] /= sa;
 
 a = rand_choice(a_choice_probs, X.len[A]+1)+1;
 
 X.len[A]       += 2 - a;
 for(int i=A+1; i<X.top; i++) //Now we have to shift the rest of stack elements...
  X.start[i] += 2 - a; //X.len[i] does not change!!! 
 
 X.start[ X.top] = (X.top>0? X.start[X.top-1] + X.len[X.top-1] : 0);
 X.top++; 
 X.len[X.top-1]  = a;
 
 X.nel          += 2;
 
 genus ++;
 
 return ACTION_SUCCESS;
}

DECLARE_ACTION_UNDO(split)
{
 RETURN_IF_FALSE(X.top>1, ERR_WRONG_STATE);

 int a    =  X.len[X.top-1];
 int A    =  (*data_in);
 X.top    --;
 X.len[A] -= (2 - a);
 for(int i=A+1; i<X.top; i++)
  X.start[i] -= (2 - a);
 X.nel    -= 2;
 genus    --;

 return ACTION_SUCCESS;
}


/************** Action fetcher *****************/
int my_action_fetcher(t_action_data** action_list, double** amplitude_list, int list_length)
{
 int nact = 0, adata = 0, res;
 double ampl;
 
 logs_Write((step_number%mc_reporting_interval==0? 1 : 2), "Step %08i:\t X.top = %i, X.nel = %i, H.top = %i, H.nel = %i, genus = %i", step_number, X.top, X.nel, H.top, H.nel, genus);
 
 if(check_stack)
 {
  res = check_stack_consistency(&X, "X");
  if(res!=0) system("PAUSE");
  res = check_stack_consistency(&H, "H");
  if(res!=0) system("PAUSE");
 }; 
 
 FETCH_ACTION(         create, 0, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(    evolve_line, 1, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(  evolve_vertex, 2, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(           join, 3, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 FETCH_ACTION(          split, 4, (*action_list), (*amplitude_list), list_length, nact, adata, ampl);
 
 return nact;
}
