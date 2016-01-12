#include "recursion_commons.h"

static double*  G[MAX_M][MAX_M];

static double   D0[MAXLS];          //LS-sized bare kinetic term
static double   G0[MAXLS];          //LS-sized bare propagator Values of the free propagator
static double   G0pq[MAXLS2];       //LS*LS-sized array, for the contact terms in the action

void init_recursion_commons();
void free_recursion_commons();

#ifdef FIX2
#define ils    (0.5)
#define LS2    (4)
#define twoLS  (4)
#define mLS    (2)
static const uint     LSn[MAX_MX2] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 0};
static const uint     LS2n[MAX_M]  = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824};
static const uint     LS2n1[MAX_M] = {0, 2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152, 8388608, 33554432, 134217728, 536870912};
static const double   iLSn[MAX_M]  = {1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625, 0.001953125, 0.0009765625, 0.00048828125, 0.000244140625, 0.0001220703125, 0.00006103515625, 0.000030517578125};
#else
static double   ils          = 0.0;       //Inverse LS
static uint     LS2          = 0;         //LS squared
static uint     twoLS        = 0;         //two times LS
static uint     mLS          = 0;         //Local static copy of LS
static uint     LSn[MAX_MX2];
static uint     LS2n[MAX_M];
static uint     LS2n1[MAX_M];
static double   iLSn[MAX_M];
#endif

void init_recursion_commons()
{
 if(mmax>=MAX_M)
 { 
  mmax = MAX_M;
  logs_WriteWarning("The input value of mmax is too large, reset to %i", MAX_M);
 };       
 
#ifndef FIX2
 mLS   = LS;
 twoLS = 2*LS;   
 ils   = 1.0/(double)LS;
 LS2   = LS*LS;
 
 LSn[0] = 1;
 for(int i=1; i<MAX_MX2; i++)
  LSn[i] = LSn[i-1]*LS;
    
 LS2n[0] = 1;
 iLSn[0] = 1.0;
 for(int i=1; i<MAX_M; i++)
 {
   LS2n[i] = LS2n[i-1]*LS2; 
  iLSn[i]  = iLSn[i-1]/(double)LS;
 };
 LS2n1[0] = 0;
 LS2n1[1] = LS;
 for(int i=2; i<MAX_M; i++)
  LS2n1[i] = LS2n1[i-1]*LS2;
#endif

 //For cleanness, setting all address elements to NULL
 for(int n=0; n<MAX_M; n++)
  for(int m=0; m<MAX_M; m++)
   G[n][m] = NULL;

 //Allocating dynamic memory to contain the pre-calculated propagators    
 uint memory_allocated = 0;
 for(int m=mmin_prc; m<=mmax; m++)
 {
  uint psp_dim = LS2;
  for(int n=1; n<=mmax-m+1; n++) 
  {
   G[m][n] = (double *)malloc(psp_dim*sizeof(double));
   if(G[m][n]!=NULL)
   {
    memory_allocated += sizeof(double)*psp_dim;
    psp_dim *= LS2; 
   }
   else
    logs_WriteErrorAndTerminate("\nUnable to allocate %u bytes of memory at n = %i, m = %i. Exiting.", psp_dim*sizeof(double), n, m); 
  };
 };
 
 logs_WriteParameter(0, "Allocated buffer memory", "%2.4lf Mb", (double)memory_allocated/(double)(1024*1024));
}

void free_recursion_commons()
{
 for(int m=0; m<=mmax; m++)
  for(int n=1; n<=mmax-m+1; n++) 
   SAFE_FREE(G[n][m]);
}
