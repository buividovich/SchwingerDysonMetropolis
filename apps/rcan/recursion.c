#include "recursion.h"

static t_real*  G[MAX_M][MAX_M];    //Storage for the expansion coefficients

static t_real   cosines[MAX_LS];
static t_real     sines[MAX_LS];
static t_real   D0[MAX_LS];          //LS-sized bare kinetic term
static t_real   G0[MAX_LS];          //LS-sized bare propagator Values of the free propagator
static t_real   G0pq[MAX_LS2];       //LS*LS-sized array, for the contact terms in the action

static uint     LSn[MAX_MX2];
static uint     LS2n[MAX_M];
static uint     LS2n1[MAX_M];
static t_real   iLSn[MAX_M];

#ifdef FIX2
 #define ils    (0.5)
 #define LS2    (4)
 #define twoLS  (4)
 #define mLS    (2)
#else
static t_real   ils          = 0.5;       //Inverse LS
static uint     LS2          = 4;         //LS squared
static uint     twoLS        = 4;         //two times LS
static uint     mLS          = 2;         //Local static copy of LS
#endif

void init_lattice_constants()
{
 if(mmax>=MAX_M)
 { 
  mmax = MAX_M;
  logs_WriteWarning("The input value of mmax is too large, reset to %i", MAX_M);
 };       
 
#ifndef FIX2
 mLS   = LS;
 twoLS = 2*LS;   
 ils   = 1.0/(t_real)LS;
 LS2   = LS*LS;
#endif
 
 LSn[0] = 1;
 for(int i=1; i<MAX_MX2; i++)
  LSn[i] = LSn[i-1]*LS;
    
 LS2n[0] = 1;
 iLSn[0] = 1.0;
 for(int i=1; i<MAX_M; i++)
 {
   LS2n[i] = LS2n[i-1]*LS2; 
  iLSn[i]  = iLSn[i-1]/(t_real)LS;
 };
 LS2n1[0] = 0;
 LS2n1[1] = LS;
 for(int i=2; i<MAX_M; i++)
  LS2n1[i] = LS2n1[i-1]*LS2;
  
 for(int m=0; m<LS; m++)
 {
    sines[m] = ACALL(sin)(2.0*M_PI*(t_real)m/(t_real)LS);
  cosines[m] = ACALL(cos)(2.0*M_PI*(t_real)m/(t_real)LS);
 }; 
}

static inline uint total_momentum(uint P, int n) //TODO: again, n here is UNDIVIDED BY TWO!!!
{
 uint mP = P, res = 0; //TODO: use here quick total!!!
 for(int i=0; i<n; i++)
 {
  res += mP%mLS;
  mP /= mLS;
 };
 return (res%mLS);
}

void unpack_momenta(uint P, int n, uint* ps)
{
 uint mP = P;
 for(int i=0; i<n; i++)
 {
  ps[i] = mP%mLS;
  mP /= mLS;
 };
}

uint pack_momenta(int n, uint* ps)
{
 uint P = 0, F = 1;
 for(int i=0; i<2*n; i++)
 {
  P += ps[i]*F;
  F *= LS;
 };
 return P;
}
//TODO: fix the ranges of MAX_M for n = (1 ... maxm+1) and m!!!

void init_kinematics(t_real bare_mass)
{
 //Initializing the bare propagator and the bare momentum term
 for(int p=0; p<LS; p++)
 {
  D0[p] = (t_real)meff_sq + 4.0*ACALL(sin)(M_PI*(t_real)p/(t_real)LS);
  G0[p] = 1.0/(bare_mass + D0[p]);
  for(int q=0; q<LS; q++)
   G0pq[LS*q + p] = ((p+q)%LS==0? G0[p]/(t_real)LS : 0.0);
 };
}

void alloc_recursion_unpacked()
{
 //For cleanness, setting all address elements to NULL
 for(int n=1; n<=MAX_MP1; n++)
  for(int m=0; m<=MAX_M; m++)
   G[n][m] = NULL;

 //Allocating dynamic memory to contain the pre-calculated propagators    
 uint memory_allocated = 0;
 for(int m=mmin_prc; m<=mmax; m++)
 {
  uint psp_dim = LS2;
  for(int n=1; n<=mmax-m+1; n++) 
  {
   G[m][n] = (t_real *)malloc(psp_dim*sizeof(t_real));
   if(G[m][n]!=NULL)
   {
    memory_allocated += sizeof(t_real)*psp_dim;
    psp_dim *= LS2; 
   }
   else
    logs_WriteErrorAndTerminate("\nUnable to allocate %u bytes of memory at n = %i, m = %i. Exiting.", psp_dim*sizeof(t_real), n, m); 
  };
 };
 
 logs_WriteParameter(0, "Allocated buffer memory", "%2.4lf Mb", (double)memory_allocated/(double)(1024*1024));
 logs_Write(0, "");
}

//TODO: save the results up to m'th order at every iteration step!!!
/***************** PACKED RECURSION START ***********************/
uint* abook[MAX_MP2];

void alloc_recursion_packed()
{
 //For cleanness, setting all address elements to NULL
 for(int n=1; n<=MAX_M+1; n++)
 { 
  for(int m=0; m<=MAX_M; m++)
   G[n][m] = NULL;
  abook[n] = NULL;
 };  
 //Count physically distinct sequences of momenta...  
 for(int n=1; n<=mmax+1; n++)
 {
  //Init address book counters...
  uint abk_cnt = 0; //Number of elements in the n'th address book
  for(uint tP=0; tP<LS2n1[n]; tP++) //Loop over conserved momenta only
  {
   uint p0 = (mLS - total_momentum(tP, 2*n-1));
   uint P   = p0 + mLS*tP;
   //Now perform the loop of cyclic shifts to find the minimal number + count multiplicity
   uint sP  = P, mP = P; int cP   = 0;
   for(int s=0; s<2*n; s++)
   {
    if(sP!=P) cP++;
    mP = MIN(mP, sP);
    sP = sP/mLS + LS2n1[n]*(sP%mLS); //Shift by 1 bit
   };
   //Now add mP to the address book, if it has not been added yet...
   //We assume that the address book is already ordered in ascending order ad try to maintain it...
   uint apos = abk_cnt-1;
   while(apos>=0 && abook[apos]>mP)
   {
    apos --;             
   };
   abook[apos] = mP; 
  };
 };  
}

static inline t_real get_G(int m, int n, uint P)
{
 return G[m][n][P];
}

static inline void set_G(int m, int n, uint P, t_real value)
{
 G[m][n][P] = value;
}

void test_recursion()
{
 t_real max_err_momentum = 0.0;
 t_real max_err_cyclic   = 0.0;
 
 for(int m=0; m<=mmax; m++)
  for(int n=1; n<=mmax-m+1; n++)
  {
   logs_Write(0, "Checking the m = %02i, n = %02i terms...", m, n);
   for(uint P=0; P<LS2n[n]; P++)
    if(total_momentum(P, 2*n)!=0)
    {
     t_real err_momentum = ACALL(fabs)(get_G(m,n,P));
     max_err_momentum = MAX(max_err_momentum, err_momentum);
     if(err_momentum>1.0E-8 && logs_noise_level>1)
     {
      char* mstr = NULL;
      get_momentum_str(P, 2*n, &mstr, ", ");
      logs_WriteWarning("Momentum non-conservation at m = %02i, n = %02i, P = %u > [%s]", m, n, P, mstr);
     };
    }
    else
    {
     //Testing the full cyclic symmetry
     uint sP = P;
     for(int i=0; i<2*n; i++)
     {
      sP = sP/mLS + LSn[2*n-1]*(sP%mLS);
      t_real err_cyclic = ACALL(fabs)(get_G(m,n,P)-get_G(m,n,sP));
      max_err_cyclic = MAX(max_err_cyclic, err_cyclic);
      if(err_cyclic>1.0E-8 && logs_noise_level>1)
      {
       char* mstr1 = NULL; char* mstr2 = NULL;
       get_momentum_str( P, 2*n, &mstr1, ", ");
       get_momentum_str(sP, 2*n, &mstr2, ", ");
       logs_WriteWarning("Cyclic symmetry violated at m = %02i, n = %02i, P=[%s](=%+2.4E), sP=[%s](=%+2.4E)", m, n,  mstr1, (double)(get_G(m,n,P)), mstr2, (double)(get_G(m,n,sP)));
      };
     };
    };  
  };
 logs_Write(0, "");
 logs_WriteParameter(0, "Max. err for non-conserved momenta",     "%2.4E", (double)max_err_momentum);
 logs_WriteParameter(0, "Max. err for the full cyclic symmetry",  "%2.4E", (double)max_err_cyclic);
}

void free_recursion()
{
 for(int m=0; m<=mmax; m++)
  for(int n=1; n<=mmax-m+1; n++) 
   SAFE_FREE(G[n][m]);
}

void printf_momentum_str(uint P, int n, char* prefix, char* separator, char* suffix)
{
 uint mP = P;
 fprintf(stdout, "%s", prefix);
 for(int i=0; i<n-1; i++)
 {
  fprintf(stdout, "%1u%s", mP%LS, separator);
  mP /= LS;
 };
 fprintf(stdout, "%1u%s", mP%LS, suffix);
 fflush(stdout);
}

void get_momentum_str(uint P, int n, char** s, char* separator) //TODO: control the use of get_momentum_str, now n is the number of momenta, UNDIVIDED BY TWO
{
 uint mP = P;
 for(int i=0; i<n-1; i++)
 {
  sprintf_append(s, "%1u%s", mP%LS, separator);
  mP /= LS;
 };
 sprintf_append(s, "%1u", mP%LS);
}

//TODO: packed/unpacked as a flag, probably at compilation time?

/***************************************************************************/
/************ GENERALIZED STRONG-COUPLING EXPANSION  ***********************/
/***************************************************************************/

void get_Gxy_generalized_sc(double* res)
{
 DECLARE_AND_MALLOC(sG, t_real, mLS);
 t_real factor = (t_real)lambda;
 for(int m=0; m<=mmax; m++)
 {
  for(int x=0; x<mLS; x++)
   sG[x] = 0;
  for(uint p=0; p<mLS; p++)
  {
   uint P = LS*((LS-p)%LS) + p;
   double k = 2.0*M_PI*(t_real)p/(t_real)LS;
   for(int x=0; x<mLS; x++)
    sG[x] += get_G(m,1,P)*ACALL(cos)(k*(double)x);
  };
  for(int x=0; x<mLS; x++) //TODO: here - a potential loss of precision 
   res[m*mLS + x] = (double)(factor*sG[x] + (m>0? (t_real)(res[(m-1)*mLS + x]) : 0.0)); 
  factor *= (t_real)lambda;
 };
 SAFE_FREE(sG);
}

void init_generalized_sc()
{
 init_lattice_constants();
 init_kinematics(lambda);
 alloc_recursion_unpacked();     
}

void run_generalized_sc()
{
 for(int m=0; m<=mmax; m++)
 {
  for(int n=1; n<=mmax-m+1; n++)
  {
   logs_Write(0, "Generating the m = %02i, n = %02i terms...", m, n);
   for(uint P=0; P<LS2n[n]; P++)
   {
    if(P%65536==0) logs_WriteStatus(1, P, LS2n[n], 1024, "");
    
    uint p1 = P%mLS;
    //First contribution - contact terms
    t_real c1 = 0.0;
    if(n>1)
    {
     for(int A=2; A<n; A++)
     {
      uint qA = (P/LS2n1[A])%mLS;
      if((p1+qA)%mLS==0)
      {
       uint P1 = mLS*((P/mLS)%LS2n1[A-1]) + (P/LS2n[A-1])%mLS;
       uint P2 = P/LS2n[A];
       for(int m1=0; m1<=m; m1++)
        c1 += get_G(m1,A-1,P1)*get_G(m-m1,n-A,P2);
      }; 
     };
     c1 *= ils*G0[p1];
     c1 += G0pq[P%LS2]*get_G(m,n-1,P/LS2);
     uint P1 = mLS*((P/mLS)%LS2n1[n-1]) + (P/LS2n[n-1])%mLS;
     c1 += G0pq[(P/LS2n1[n])*mLS + p1]*get_G(m,n-1,P1);
    }
    else
    {
     if(m==0)
      c1 += G0pq[P];   
    }; //End of if(n>1)
    t_real c2 = 0.0; t_real c3 = 0.0;
    if(m>0)
    {
     //Second contribution - double-trace vertex
     for(int A=2; A<=n; A++)
     {
      uint P1 = mLS*((P/mLS)%LS2n1[A-1]); //This is now the sequence 0, q_1, ..., p_{A-1}, q_{A-1} 
      uint P2 = mLS*(P/LS2n1[A]);       //This is now the sequence 0, q_A, ..., p_n, q_n
      uint pA = (P/LS2n[A-1])%mLS;
      for(uint p1t=0; p1t<mLS; p1t++)
      {
       uint pAt = (p1 + pA + (mLS - p1t))%mLS;
       for(int m1=0; m1<m; m1++)
        c2 += get_G(m1,A-1,p1t + P1)*get_G(m-m1-1,n-A+1,pAt + P2);
      }; 
     };
     c2 *= ils*G0[p1];
     //Third contribution - scalar vertex
     uint P1 = (P/mLS)*LS2n1[2];
     for(uint p1t=0; p1t<mLS; p1t++)
      for(uint q1t=0; q1t<mLS; q1t++)
      {
       uint p2t = (p1 + (mLS - p1t) + (mLS - q1t))%mLS;
       c3 += D0[q1t]*get_G(m-1,n+1,p1t + mLS*q1t + LS2*p2t + P1);
      };
     c3 *= G0[p1]; 
    }; //End of if(m>0) 
    //Summing up all the contributions
    set_G(m, n, P, c1 - c2 + c3); //-c2 is important!!! It is the origin of the sign problem!!!
   }; //End of loop over P      
  }; //End of loop over n    
 }; //End of loop over m    
}

/***************************************************************************/
/************ STEREOGRAPHIC PROJECTION *************************************/
/***************************************************************************/
static t_real alpha_stereographic = 0.0;
static t_real    m0_stereographic = 0.0;

void init_stereographic()
{
 meff_sq = 0.0;
 alpha_stereographic = 0.125*(t_real)lambda;
    m0_stereographic = 0.25*(t_real)lambda;
 init_lattice_constants();
 init_kinematics(m0_stereographic);
 alloc_recursion_unpacked();     
}

static uint qsbuf[MAX_MX2];
//Vertex function V(q_0, q_1, ..., q_{2 n_q}), nq = 1, 2, ..., Q = q_1 + LS*q_2 + LS^2*q_3 + ... 
static inline t_real get_vertex(uint Q, int nq)
{
 t_real res = m0_stereographic;
 unpack_momenta(Q, 2*nq+1, qsbuf);
 for(int m=0; m<=2*nq; m++)
 {
  uint Qt = qsbuf[m]; res += D0[Qt];
  int s = -1;
  for(int l=1; l<=2*nq-m; l++)
  {
   Qt += qsbuf[m+l];
   res += s*D0[(Qt%mLS)];
   s *= -1;
  };
 };
 return res;
}

static inline void GammaXQ(uint Q, int nq, t_real* Gamma)
{
 for(uint x=0; x<mLS; x++)
  Gamma[x] = 0.0;
 uint q = 0, mQ = Q;
 for(int k=0; k<2*nq-1; k++)
 {
  q  += mQ%mLS;
  mQ /= mLS;
  int s = (k%2? +1 : -1);
  for(uint x=0; x<mLS; x++)
   Gamma[x] += s*cosines[(q*x)%mLS]; //s*ACALL(cos)(qphys*(t_real)x);
 };
}

void get_Gxy_stereographic(double* Gxy, double* Gx)
{
 DECLARE_AND_MALLOC( gamma, t_real, mLS);
 DECLARE_AND_MALLOC(  tGxy, t_real, (mmax+1)*mLS);
 DECLARE_AND_MALLOC(   tGx, t_real, (mmax+1)    );
 
 for(int i=0; i<=mmax*mLS; i++)
  tGxy[i] = 0;
 for(int i=0; i<=mmax; i++)
  tGx[i] = 0; 
 
 for(int n=1; n<=mmax+1; n++)
  for(uint P=0; P<LS2n[n]; P++) //TODO: here - sum only over conserved momenta...
  {
   GammaXQ(P, n, gamma);
   for(int m=0; m<=(mmax-n+1); m++)
   {
    t_real fc = (n%2==0? 1 : -1)*ACALL(pow)(alpha_stereographic, (t_real)(n+m))*get_G(m, n, P);
    for(int mmax1=m+n-1; mmax1<=mmax; mmax1++)
     tGx[mmax1] += fc;
    for(int x=0; x<mLS; x++)
    {
     t_real fc1 = fc*gamma[x];
     for(int mmax1=m+n-1; mmax1<=mmax; mmax1++)
      tGxy[mmax1*mLS + x] += fc1;
    };
   };
  };
 
 for(int mmax1=0; mmax1<=mmax; mmax1++)
 {
  tGx[mmax1] = 1.0 + 2.0*tGx[mmax1];
  Gx[mmax1] = (double)tGx[mmax1];
  for(int x=0; x<mLS; x++)
   Gxy[mmax1*mLS + x] = (double)(2.0*tGx[mmax1] - 1.0 + 4.0*tGxy[mmax1*mLS + x]);
 };

 SAFE_FREE(gamma);
 SAFE_FREE(tGxy);
 SAFE_FREE(tGx);
}

//TODO: automatic initialization of LSn arrays!
//TODO: for vertex and Gamma - pre-computed array of cosines?
//TODO: pre-computed vertices?

void run_stereographic()
{
 for(int m=0; m<=mmax; m++)
 {
  for(int n=1; n<=mmax-m+1; n++)
  {
   logs_Write(0, "Generating the m = %02i, n = %02i terms...", m, n);
   for(uint P=0; P<LS2n[n]; P++)
   {
    if(P%65536==0) logs_WriteStatus(1, P, LS2n[n], 1024, "");
    uint p0 = P%mLS;
    //First contribution - contact terms
    t_real c1 = 0.0;
    if(n>1)
    {
     for(int A=1; A<=n-2; A++)
     {
      uint pA = (P/LS2n1[A+1])%mLS;
      if((p0+pA)%mLS==0)
      {
       uint P1 = (P/mLS)%LS2n[A];
       uint P2 =  P/LS2n[A+1]; 
       for(int k=0; k<=m; k++)
        c1 += get_G(k,A,P1)*get_G(m-k,n-A-1,P2);
      }; 
     };
     uint p1   = (P%LS2)/mLS;
     if((p0+p1)%mLS==0)
      c1 += get_G(m,n-1,P/LS2);
     p1 = P/LS2n1[n];
     if((p0+p1)%mLS==0)
      c1 += get_G(m,n-1,(P/mLS)%LS2n[n-1]);
     c1 *= ils*G0[p0];
    }
    else
    {
     if(m==0)
      c1 += G0pq[P];   
    }; //End of if(n>1)
    //Second contribution - higher-order vertices
    t_real c2 = 0.0;
    for(int k=1; k<=m; k++)
    {
     uint P1 = (P/mLS)*LS2n1[k+1];
     int   s = (k%2==0? -1 : 1);
     for(uint sQ=0; sQ<LS2n[k]; sQ++) //sequence q1 ... q_{2 k}
     {
      //TODO: both total_momentum and get_vertex require a loop of unpacking... - can one optimize something here?        
      uint q0   = (p0 + (mLS - total_momentum(sQ, 2*k)))%mLS; //This ensures momentum conservation
      uint Q    = q0 + mLS*sQ;
      c2 += s*get_vertex(Q, k)*get_G(m-k,n+k,Q + P1);
     };
    };
    c2 *= G0[p0];
    //Summing up all the contributions
    set_G(m, n, P, c1 + c2); 
   }; //End of loop over P      
  }; //End of loop over n    
 }; //End of loop over m    
}

/*
//For testing purposes, the content of the relevant array entries
 printf("\n");
 for(int n=1; n<=mmax+1; n++)
  for(uint P=0; P<LS2n[n]; P++)
  {
   printf("G%i", 2*n);
   printf_momentum_str(P, 2*n, "[", ",", "] = ");
   for(int m=0; m<=(mmax-n+1); m++)
    printf("%2.4E ", get_G(m, n, P));
   printf("\n"); 
  };
 printf("\n");  
*/
