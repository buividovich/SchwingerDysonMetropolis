#include "my_stack.h"

int **stack = NULL;
int   *stop = NULL;
double  *nA = NULL; //Here sums of probabilities will be saved - in order to avoid recalculations
int  *asign = NULL;
int      ns = 0;

void init_stack()
{
 int i;
 SAFE_MALLOC(stack,   int*, maxn);
 SAFE_MALLOC( stop,    int, maxn);
 SAFE_MALLOC(   nA, double, maxn);
 SAFE_MALLOC(asign,    int, maxn);
 for(i=0; i<maxn; i++)
 {
  SAFE_MALLOC(stack[i], int, (i+1));
  stop[i] = 0;
  nA[i]   = 0.0;
 };
 ns = 0;
 //Init the first element in the sta
 stack[ns][stop[ns]] = 1; 
 asign[ns]           = 1;
}

void free_stack()
{
 int i;
 for(i=0; i<maxn; i++)
  SAFE_FREE(stack[i]);
 SAFE_FREE(stack); 
 SAFE_FREE(stop);
 SAFE_FREE(nA);
 SAFE_FREE(asign);
 ns = 0;
}

//TODO: not so "mem-greedy" algorithm? Is it possible?
