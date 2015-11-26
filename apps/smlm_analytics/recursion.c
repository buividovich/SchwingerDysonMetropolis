#include "recursion.h" 

static double* prcG[MAX_M][MAX_M];
static int*    prcF[MAX_M][MAX_M]; //Flag for precalculated quantities, set to one if the corresponding element of prcG is already calculated

static double   D0[MAXLS2];         //LS-sized bare kinetic term
static double   G0[MAXLS2];         //LS-sized bare propagator Values of the free propagator
static double   ils          = 0.0;          //Inverse LS
static int      LS2          = 0;
static uint     LSn[MAX_MX2];
static double   iLSn[MAX_MX2];
static double   ms1;                //sigma parameters controlling the thimble structure
static double   ms2s3;              //ratio of s2 to s3

//Pack all momenta into a single unsigned integer
//Quick ways to pack and unpack, get all the components?

//But - we never have to really unpack!!!
//Unpack only the lowest few components

//How to code for LS=2?
double G(uint P, int n, int m)
{
 static uint mP;
 //static uint mP1;
 //static uint mP2;
 static double res0; //The very first contact term which is there for all the equations
 static double res1;
 static int i;
 //static int A;

 mP   = P;
 res0 = ils*G0[mP%LS2];
 res1 = 0.0;
 
 if(m<=0)
 {
  for(i=1; i<n; i++)
  {
   mP   /= LS2;
   res0 *= G0[mP%LS2];
  };
  res0 *= iLSn[n-1]; 
 }
 else
 {
  mP  /= LS2;
  res0 *= G(mP, n-1, m);
    
  //Second term in SD equations - flip momenta
  //for(int A=2; A<=n; A++)
  // for(int p1t=0; p1t<2; p1t++)
    
  
  //Finally, multiply by G0(P0) 
  res1 *= G0[P%LS2];
 };
 
 return (res0 + res1);
}

void init_recursion()
{
 if(mmax>=MAX_M)
 { 
  mmax = MAX_M;
  logs_WriteWarning("The input value of mmax is too large, reset to %i", MAX_M);
 };       
     
 ils = 1.0/(double)LS;
 LS2 = LS*LS;
 
 ms1    = (double)s1;
 ms2s3  = (double)s2s3;
 
  LSn[0] = 1;
 iLSn[0] = 1.0;
 for(int i=1; i<2*(mmax+1); i++)
 {
   LSn[i] *= LS;
  iLSn[i] /= (double)LS;
 }; 
 
 //Initializing the bare propagator and the bare momentum term
 for(int p=0; p<LS; p++)
  for(int q=0; q<LS; q++)
  {
   D0[LS*q + p] = (p==q? 4.0*sin(M_PI*(double)p/(double)LS) : 0.0);
   G0[LS*q + p] = (p==q?            1.0/(m2 + D0[LS*q + p]) : 0.0);
  };    
 
 //Allocating dynamic memory to contain the pre-calculated propagators    
 int memory_allocated = 0;
 int psp_dim = LS2;
 for(int n=0; n<=mmax; n++)
 {
  for(int m=0; m<=mmax; m++)
  {
   SAFE_MALLOC(prcG[n][m], double, psp_dim);
   SAFE_MALLOC(prcF[n][m],    int, psp_dim);
   memory_allocated += sizeof(double)*psp_dim + sizeof(int)*psp_dim;
   for(int P=0; P<psp_dim; P++)
    prcF[n][m] = 0;
  };
  //Increase the size
  psp_dim *= LS2;
 };
 logs_WriteParameter(0, "Allocated memory", "%2.2lf Mb", (double)memory_allocated/(double)(1024*1024));
}

void free_recursion()
{
 for(int n=0; n<MAX_M; n++)
  for(int m=0; m<MAX_M; m++)
  {
   SAFE_FREE(prcG[n][m]);
   SAFE_FREE(prcF[n][m]);
  };    
}
