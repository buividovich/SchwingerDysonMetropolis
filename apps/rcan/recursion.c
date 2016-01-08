#include "recursion.h" 

static double*     prc_Buffer    = NULL;
static double**    prc_Address[MAX_M][MAX_M]; //Address for pre-calculated quantities
static uint        prc_Counter   = 0;

static double   D0[MAXLS];                //LS-sized bare kinetic term
static double   G0[MAXLS];                //LS-sized bare propagator Values of the free propagator
static double   G0pq[MAXLS2];             //LS*LS-sized array, for the contact terms in the action
static double   ms1;                //sigma parameters controlling the thimble structure
static double   ms2s3;              //ratio of s2 to s3

#ifdef FIX2
#define ils    (0.5)
#define LS2    (4)
#define twoLS  (4)
#define mLS    (2)
static const uint     LSn[MAX_MX2] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, 0};
static const uint     LS2n[MAX_M]  = {1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824};
static const uint     LS2n1[MAX_M] = {0, 2, 8, 32, 128, 512, 2048, 8192, 32768, 131072, 524288, 2097152, 8388608, 33554432, 134217728, 536870912};
static const double   iLSn[MAX_M]  = {1.0000000000000000000, 0.50000000000000000000, 0.25000000000000000000, \
0.12500000000000000000, 0.062500000000000000000, 0.031250000000000000000, \
0.015625000000000000000, 0.0078125000000000000000, 0.0039062500000000000000, \
0.0019531250000000000000, 0.00097656250000000000000, \
0.00048828125000000000000, 0.00024414062500000000000, \
0.00012207031250000000000, 0.000061035156250000000000, \
0.000030517578125000000000};
#else
static double   ils          = 0.0;       //Inverse LS
static uint     LS2          = 0;
static uint     twoLS        = 0;
static uint     mLS          = 0;         //Local static copy of LS
static uint     LSn[MAX_MX2];
static uint     LS2n[MAX_M];
static uint     LS2n1[MAX_M];
static double   iLSn[MAX_M];
#endif

void printf_momentum_str(uint P, int n, char* prefix, char* separator, char* suffix)
{
 uint mP = P;
 fprintf(stdout, "%s", prefix);
 for(int i=0; i<2*n-1; i++)
 {
  fprintf(stdout, "%1u%s", mP%LS, separator);
  mP /= LS;
 };
 fprintf(stdout, "%1u%s\n", mP%LS, suffix);
 fflush(stdout);
}

void print_buffer_usage()
{
 logs_WriteParameter(0, "\t Buffer usage (doubles)", "%i of %i\t ( %2.2lf%% ) \n", prc_Counter, MAX_BUFFER, 100.0*(double)prc_Counter/(double)MAX_BUFFER);
 for(int m=mmin_prc; m<=mmax; m++)
 {
  int cnt_used  = 0;
  int cnt_alloc = 0;
  int psp_dim   = LS2;
  for(int n=1; n<=mmax-m+1; n++) 
  {
   for(int P=0; P<psp_dim; P++)
    if((prc_Address[n-1][m])[P]!=NULL)
     cnt_used ++;
   cnt_alloc += psp_dim; 
   psp_dim *= LS2;
  };
  logs_WriteParameter(0, "\t Buffer usage (addresses)", "% 8i of % 8i ( %2.2lf%% ) at m = % 2i", cnt_used, cnt_alloc, 100.0*(double)cnt_used/(double)cnt_alloc, m);
 };
}

static double mG(uint P, int n, int m)
{
 if(m<mmin_prc || (prc_Address[n-1][m])[P]==NULL)
 { 
  double res0 = ils*G0pq[P%LS2];
  if(m<=0)
  {
   uint mP = P;
   for(int i=1; i<n; i++)
   {
    mP   /= LS2;
    res0 *= G0pq[mP%LS2];
   };
   res0 *= iLSn[n-1];
  }
  else
  {
   res0 *= (n>1? mG(P/LS2, n-1, m) : 0.0);
   double res1 = 0.0, res2 = 0.0, res3 = 0.0;
    
   //Extract p1 - necessary everywhere
   uint p1 = P%mLS;
  
   //Now go the terms which join two sequences in one
   for(int A=2; A<=n; A++)
   {
    int A1 = A - 1;
    int A2 = n - A + 1;
    //First split P into two parts
    uint mP2  =   P/LS2n[A1];                  //This is the sequence p_A q_A ... p_n q_n 
    uint mP1a = ((P%LS2n[A1])/mLS)*mLS;        //This is the sequence 0 q_1 ... p_{A-1} q_{A-1}
    uint mP1b = ((P%LS2n[A ])/mLS)*mLS;        //This is the sequence 0 q_1 ... p_A q_A
    //Extract pA
    uint pA   = mP2%mLS;
    //Strip off pA and prepend 0 in place of pAt - now mP2 is the sequence 0 q_A p_{A+1} q_{A+1} ... p_n q_n
    mP2 = (mP2/mLS)*mLS;
    //Now goes the loop over internal momenta - for the second and the third terms in SDs
    for(uint pAt=0; pAt<mLS; pAt++)
    {
     uint p1t0  = (p1 + (mLS - pAt))%mLS; //p1 - pAt
     //Second term in SDs - flip momenta
     uint p1t   = (p1t0 + pA)%mLS;
     uint mP    = mP1a + p1t; //mP is the sequence p1t q_1 ... p_{A-1} q_{A-1}
     for(int m1=0; m1<m; m1++)
      res1 += mG(mP, A1, m1)*mG(mP2, A2, m - m1 - 1);
     //Third term in SDs - join sequences with a vertex
     for(uint qAt=0; qAt<mLS; qAt++)
     {
      p1t = (p1t0 + (mLS - qAt))%mLS; //p1t = p1 - pAt - qAt
      mP  = mP1b + p1t; //mP is the sequence p1t q_1 ... p_A q_A
      mP  = (mP%LS2n1[A]) + qAt*LS2n1[A];
      for(int m1=0; m1<m; m1++)
       res2 += D0[qAt]*mG(mP, A, m1)*mG(mP2, A2, m - m1 - 1);
     };
     mP2 ++; //mP2 is now the sequence pAt q_A p_{A+1} q_{A+1} ... p_n q_n 
    };
   };// <- End of loop over A
    
   //Last term in SDs - a scalar vertex 
   uint mP = (P/mLS)*mLS; //Now mP is the sequence of the form 0 q_1 p_2 q_2 ... p_n q_n
   for(uint p2t=0; p2t<mLS; p2t++) 
   {
    for(uint q1t=0; q1t<mLS; q1t++) 
    { 
     uint p1t   = (p1 + (twoLS - p2t - q1t))%mLS;
     res3 += D0[q1t]*mG(LS2*mP + mLS*q1t + p1t, n+1, m-1);
    };
    mP ++; //mP is the sequence p2t q_1 p_2 q_2 ... p_n q_n         
   };
    
   //Finally, multiply by G0(P0) 
   res1 *= ils*ms2s3; 
   res2 *= -ms2s3; 
   res3 *= -ms2s3*ms1;
   
   res0 += (res1 + res2 + res3)*G0[p1];
  };
  //At this point res0 contains the result - add it to the collection of pre-calculated data
  if(m>=mmin_prc && prc_Counter<MAX_BUFFER)
  {
   if(n>mmax-m+1)
   {
    logs_WriteError("Wrong n=%i for m=%i", n, m); //TODO: remove after debugging
    exit(EXIT_FAILURE);
   }; 
   prc_Buffer[prc_Counter]  = res0;
   (prc_Address[n-1][m])[P] = &(prc_Buffer[prc_Counter]);
   prc_Counter ++;
  };
  
  if(prc_Counter>=MAX_BUFFER)
  {
   logs_WriteError("Buffer overflow");
   exit(EXIT_FAILURE);
  };
  
  return res0; 
 }
 else
 {
  return *((prc_Address[n-1][m])[P]);
 };
}

void getG2(double* G2, int m)
{
 for(int p=0; p<LS; p++)
  G2[p]  = mG(p*LS + (LS-p)%LS, 1, m);
}

void testG4(int m)
{
 for(uint P=0; P<LSn[4]; P++)
 {
  char* mS = NULL;
  double res = mG(P, 2, m);
  get_momentum_str(P, 2, &mS, ", ");
  if(total_momentum(P, 2))
  {
   if(fabs(res)>0.0)
    logs_WriteError("\t [%s]:\t%+2.8lf", mS, res);  
  } 
  else 
  {
   logs_Write(0, "\t [%s]:\t%+2.8lf", mS, res);
  }; 
  SAFE_FREE(mS);
 };    
}

uint total_momentum(uint P, int n)
{
 uint mP = P, res = 0;
 for(int i=0; i<2*n; i++)
 {
  res += mP%LS;
  mP /= LS;
 };
 return (res%LS);
}

void get_momentum_str(uint P, int n, char** s, char* separator)
{
 uint mP = P;
 for(int i=0; i<2*n-1; i++)
 {
  sprintf_append(s, "%1u%s", mP%LS, separator);
  mP /= LS;
 };
 sprintf_append(s, "%1u", mP%LS);
}

void init_recursion()
{
 if(mmax>=MAX_M)
 { 
  mmax = MAX_M;
  logs_WriteWarning("The input value of mmax is too large, reset to %i", MAX_M);
 };       
 
 ms1    = (double)s1;
 ms2s3  = (double)s2s3;

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
 
 //Initializing the bare propagator and the bare momentum term
 for(int p=0; p<LS; p++)
 {
  D0[p] = 4.0*sin(M_PI*(double)p/(double)LS);
  G0[p] = 1.0/(m2 + D0[p]);
  for(int q=0; q<LS; q++)
   G0pq[LS*q + p] = (p==q? G0[p] : 0.0);
 };     

 //For cleanness, setting all address elements to NULL
 for(int n=0; n<MAX_M; n++)
  for(int m=0; m<MAX_M; m++)
   prc_Address[n][m] = NULL;

 //Allocating dynamic memory to contain the pre-calculated propagators    
 int memory_allocated = 0;
 for(int m=mmin_prc; m<=mmax; m++)
 {
  int psp_dim = LS2;
  for(int n=1; n<=mmax-m+1; n++) 
  {
   SAFE_MALLOC((prc_Address[n-1][m]), double*, psp_dim);
   memory_allocated += sizeof(double*)*psp_dim;
   for(int P=0; P<psp_dim; P++)
    (prc_Address[n-1][m])[P] = NULL;
   psp_dim *= LS2; 
  };
  //Increase the size
 };
 
 logs_WriteParameter(0, "Allocated memory for addresses", "%2.4lf Mb", (double)memory_allocated/(double)(1024*1024));
 
 SAFE_MALLOC(prc_Buffer, double, MAX_BUFFER);
 memory_allocated += sizeof(double)*MAX_BUFFER;
 prc_Counter = 0;
 
 logs_WriteParameter(0, "Allocated memory", "%2.2lf Mb", (double)memory_allocated/(double)(1024*1024));
}

void free_recursion()
{
 for(int m=mmin_prc; m<=mmax; m++)
  for(int n=0; n<=mmax-m; n++)
   SAFE_FREE((prc_Address[n][m]));
 SAFE_FREE(prc_Buffer);  
}
