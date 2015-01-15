#include "metropolis.h"

int overflow_count = 0;

static inline double A_create_sc()
{
 if(stop[ns]<=ns)
 {
  return 1.0/(fabs(lambda)*NN*cc);
 }
 else
 {
  logs_WriteError("stop[%i] = %i, cannot create new element");
  return 0.0;
 }; 
}

static inline double A_increase_n_sc()
{
 return 1.0/(fabs(lambda)*cc);
}

static inline double A_decrease_n_sc()
{
 if(stack[ns][stop[ns]]>1)
  return cc/fabs(lambda);
 return 0.0;
}

static inline double A_join_sc()
{
 if(stop[ns]>0)
  return NN;
 return 0.0; 
}

int mc_step_sc()
{
 int accepted = 0, i, todo;
 double a, alpha, A[4], p[4];
 //First calculate, where we can go from the current state and the corresponding transition weights - they are necessary for calculating <nA>
 A[0]   = A_create_sc();
 A[1]   = A_increase_n_sc();
 A[2]   = A_decrease_n_sc();
 A[3]   = A_join_sc();
 nA[ns] = A[0] + A[1] + A[2] + A[3];
 //Updating the mean nA, its dispersion and the mean sign
 anA   += nA[ns];
 dnA   += SQR(nA[ns]);
 msign += asign[ns];
 //Decide whether we move forward or backward in the order of expansion
 if(ns<maxn-1)
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
 if(a <= pplus)
 {//Attaching new element to the sequence
  //Acceptance probability
  alpha  = MIN(1.0, (1.0 - pplus)/pplus*nA[ns]);
  a      = rand_double(0.0, 1.0);      
  if(a<alpha)
  {//Accept and perform the forward step
   accepted = 1;
   //Copy the stack to the next step
   memcpy(stack[ns+1], stack[ns], sizeof(int)*(stop[ns]+1));
   stop[ns+1]  = stop[ns];
   asign[ns+1] = asign[ns];
   //Now calculate the normalized probabilities of different steps
   for(i=0; i<4; i++)
    p[i] = A[i]/nA[ns];
   //And increase the order counter 
   ns++; 
   //Finally, nontrivial transformations of the stack
   todo = rand_choice(p, 4);
   //Modify the stack according to the selection
   if(todo==0) //Push new element into the stack
   {
    stop[ns]++;
    stack[ns][stop[ns]] = 1; 
    asign[ns] *= SIGN(lambda);     
   };
   if(todo==1) //Increase n by 1
   {
    stack[ns][stop[ns]] ++;
    asign[ns] *= SIGN(lambda);
   };
   if(todo==2 && stack[ns][stop[ns]]>1) //Decrease n by 1 and flip sign
   {
    stack[ns][stop[ns]] --;
    asign[ns] *= -1*SIGN(lambda);
   };
   if(todo==3 && stop[ns]>0) //Join the two topmost and flip sign 
   {
    stack[ns][stop[ns]-1] += stack[ns][stop[ns]];
    asign[ns] *= -1;
    stop[ns]--;          
   };     
  };
 }
 else
 {//Step backward
  if(ns>0)
  {
   alpha  = MIN( 1.0, pplus/((1.0 - pplus)*nA[ns-1]) ); //Acceptance probability
   a      = rand_double(0.0, 1.0);
   if(a<alpha)
   {
    accepted = 1;
    ns--;
   }; 
  }
  else
  {
   accepted = 1;
   stack[ns][stop[ns]] = 1;
   asign[ns] = 1;
  };
 };
 return accepted; //Return 1 if some new state is accepted, otherwise 0
}

int mc_step_wc()
{
 int accepted = 0, i, k, todo;
 double a, alpha;
 double *A = NULL;
 double *p = NULL;
 //First calculate, where we can go from the current state and the corresponding transition weights - they are necessary for calculating <nA>
 //Determining the number of choices
 int ndecr    = (stack[ns][stop[ns]] - 1); //We can decrease the current element by any number down to 1
 int nchoices = 3 + ndecr;
 SAFE_MALLOC(A, double, nchoices);
 SAFE_MALLOC(p, double, nchoices);
 //Non-normalized transition probabilities
 A[0]   = 1.0/(NN*cc);         //Create "1"
 A[1]   = (stop[ns]>0? NN/cc : 0.0);               //Join two sequences
 A[2]   = 2.0/cc;              //Increase by one
 nA[ns] = A[0] + A[1] + A[2];
 double fct = cc*SQR(alpha_wc);
 for(k=1; k<=ndecr; k++)
 {
  A[2 + k] = 2.0*fct*((double)(4*k + 4) + lambda)/lambda;
  nA[ns] += A[2 + k];
  fct *= cc*alpha_wc;
 };
 //Updating the mean nA, its dispersion and the mean sign
 anA   += nA[ns];
 dnA   += SQR(nA[ns]);
 maxnA  = MAX(maxnA, nA[ns]);
 msign += asign[ns];
 //Decide whether we move forward or backward in the order of expansion
 if(ns<maxn-1)
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
 if(a <= pplus)
 {
  //Acceptance probability
  alpha  = MIN(1.0, (1.0 - pplus)/pplus*nA[ns]);
  a      = rand_double(0.0, 1.0);      
  if(a<alpha)
  {//Accept and perform the forward step
   accepted = 1;
   //Copy the stack to the next step
   memcpy(stack[ns+1], stack[ns], sizeof(int)*(stop[ns]+1));
   stop[ns+1]  = stop[ns];
   asign[ns+1] = asign[ns];
   //Now calculate the probabilities of different steps
   for(i=0; i<nchoices; i++)
    p[i] = A[i]/nA[ns];
   //And increase the order counter 
   ns++; 
   //Finally, nontrivial transformations of the stack
   todo = rand_choice(p, nchoices);
   //Modify the stack according to the selection
   if(todo==0) //Push new element into the stack
   {
    stop[ns]++;
    stack[ns][stop[ns]] = 1;      
   };
   if(todo==1 && stop[ns]>0) //Join the two topmost and flip sign 
   {
    stack[ns][stop[ns]-1] += (stack[ns][stop[ns]] + 1);
    stop[ns]--;          
   };     
   if(todo==2) //Increase n by 1
   {
    stack[ns][stop[ns]] ++;
   };
   if(todo>2 && stack[ns][stop[ns]]>(todo-2)) //Decrease n by (todo-2) and flip sign if (todo-2) is odd
   {
    stack[ns][stop[ns]] -= (todo-2);
    asign[ns] *= ((todo-2)%2==0? -1 : 1);
   };
  };
 }
 else
 {//Step backward
  if(ns>0)
  {
   alpha  = MIN( 1.0, pplus/((1.0 - pplus)*nA[ns-1]) ); //Acceptance probability
   a      = rand_double(0.0, 1.0);
   if(a<alpha)
   {
    accepted = 1;
    ns--;
   }; 
  }
  else
  {
   accepted = 1;
   stack[ns][stop[ns]] = 1;
   asign[ns] = 1;
  };
 };
 SAFE_FREE(A);
 SAFE_FREE(p);
 return accepted; //Return 1 if some new state is accepted, otherwise 0
}
