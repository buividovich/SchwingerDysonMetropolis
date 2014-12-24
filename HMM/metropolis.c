#include "metropolis.h"

int overflow_count = 0;

double A_create()
{
 return 1.0/NN;
}

double A_evolve_line()
{
 return 2.0/cc;
}

double A_evolve_vertex()
{
 if(stack[ns][stop[ns]]>1)
  return fabs(lambda)*cc;
 return 0.0;       
}

double A_join()
{
 if(stop[ns]>0)
  return NN/SQR(cc);
 return 0.0; 
}

int mc_step()
{
 int accepted = 0, i, todo;
 double a, alpha, A[4], p[4];
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
  A[0]   = A_create();
  A[1]   = A_evolve_line();
  A[2]   = A_evolve_vertex();
  A[3]   = A_join();
  nA[ns] = A[0] + A[1] + A[2] + A[3];
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
   };
   if(todo==1) //Add new line to the diagram
   {
    stack[ns][stop[ns]] ++;
   };
   if(todo==2) //Create the new vertex
   {
    stack[ns][stop[ns]] --;
    asign[ns] *= SIGN(lambda);
   };
   if(todo==3)
   {
    stack[ns][stop[ns]-1] += stack[ns][stop[ns]] + 1;
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
