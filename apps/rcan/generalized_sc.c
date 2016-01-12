#include "generalized_sc.h"
#include "recursion_core.h"

void init_generalized_sc()
{
 init_recursion_commons();
 //Initializing the bare propagator and the bare momentum term
 for(int p=0; p<LS; p++)
 {
  D0[p] = meff_sq + 4.0*sin(M_PI*(double)p/(double)LS);
  G0[p] = 1.0/(lambda + D0[p]);
  for(int q=0; q<LS; q++)
   G0pq[LS*q + p] = ((p+q)%LS==0? G0[p]/(double)LS : 0.0);
 };    
}

void  run_generalized_sc()
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
    double c1 = 0.0;
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
        c1 += G[m1][A-1][P1]*G[m-m1][n-A][P2];
      }; 
     };
     c1 *= ils*G0[p1];
     c1 += G0pq[P%LS2]*G[m][n-1][P/LS2];
     uint P1 = mLS*((P/mLS)%LS2n1[n-1]) + (P/LS2n[n-1])%mLS;
     c1 += G0pq[(P/LS2n1[n])*mLS + p1]*G[m][n-1][P1];
    }
    else
    {
     if(m==0)
      c1 += G0pq[P];   
    }; //End of if(n>1)
    double c2 = 0.0; double c3 = 0.0;
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
        c2 += G[m1][A-1][p1t + P1]*G[m-m1-1][n-A+1][pAt + P2];
      }; 
     };
     c2 *= ils*G0[p1];
     //Third contribution - scalar vertex
     uint P1 = (P/mLS)*LS2n1[2];
     for(uint p1t=0; p1t<mLS; p1t++)
      for(uint q1t=0; q1t<mLS; q1t++)
      {
       uint p2t = (p1 + (mLS - p1t) + (mLS - q1t))%mLS;
       c3 += D0[q1t]*G[m-1][n+1][p1t + mLS*q1t + LS2*p2t + P1];
      };
     c3 *= G0[p1]; 
    }; //End of if(m>0) 
    //Summing up all the contributions
    G[m][n][P] = c1 - c2 + c3; //-c2 is important!!! It is the origin of the sign problem!!!
   }; //End of loop over P      
  }; //End of loop over n    
 }; //End of loop over m    
}

void test_generalized_sc()
{
 double max_err_momentum = 0.0;
 double max_err_cyclic   = 0.0;
 double max_err_pq       = 0.0;
 //uint ps[MAX_MX2];
 //First testing the momentum conservation
 for(int m=0; m<=mmax; m++)
  for(int n=1; n<=mmax-m+1; n++)
  {
   logs_Write(0, "Checking the m = %02i, n = %02i terms...", m, n);
   for(uint P=0; P<LS2n[n]; P++)
    if(total_momentum(P, n)!=0)
    {
     max_err_momentum = MAX(max_err_momentum, fabs(G[m][n][P]));
     /*if(fabs(G[m][n][P])>1.0E-8)
     {
      char* mstr = NULL;
      get_momentum_str(P, n, &mstr, ", ");
      logs_WriteError("Momentum non-conservation at m = %02i, n = %02i, P = %u > [%s]", m, n, P, mstr);
     };*/
    }
    else
    {
     //Testing the full cyclic symmetry
     uint sP = P;
     for(int i=0; i<2*n; i++)
     {
      sP = sP/LS + LSn[2*n-1]*(sP%LS);
      max_err_cyclic = MAX(max_err_cyclic, fabs(G[m][n][P]-G[m][n][sP]));
      if(fabs(G[m][n][P]-G[m][n][sP])>1.0E-8)
      {
       char* mstr1 = NULL; char* mstr2 = NULL;
       get_momentum_str( P, n, &mstr1, ", ");
       get_momentum_str(sP, n, &mstr2, ", ");
       logs_WriteError("Cyclic symmetry violated at m = %02i, n = %02i, P=[%s], sP=[%s]", m, n, mstr1, mstr2);
       system("PAUSE");
      };
     };
     //Testing symmetry for the replacement of all P by Q   
    /* unpack_momenta(P, n, ps);
     for(int i=0; i<2*n; i+=2)
     {
      uint tmp = ps[i];
      ps[i]    = ps[i+1];
      ps[i+1]  = tmp;
     };
     sP = pack_momenta(n, ps); 
     double err_pq = fabs(G[m][n][P]-G[m][n][sP]);
     max_err_pq = MAX(max_err_pq, err_pq);*/
     /*if(err_pq>1.0E-8)
     {
      char* mstr1 = NULL;
      char* mstr2 = NULL;
      get_momentum_str( P, n, &mstr1, ", ");
      get_momentum_str(sP, n, &mstr2, ", ");
      logs_WriteError("P<->Q symmetry violated at m = %02i, n = %02i (deviation %2.4E), P = [%s], sP = [%s]", m, n, err_pq, mstr1, mstr2);
      system("PAUSE");
     };*/
    };  
  };
 logs_Write(0, "");
 logs_WriteParameter(0, "Max. err for non-conserved momenta",     "%2.4E", max_err_momentum);
 logs_WriteParameter(0, "Max. err for the full cyclic symmetry",  "%2.4E", max_err_cyclic);
 logs_WriteParameter(0, "Max. err for P<->Q symmetry",            "%2.4E", max_err_pq);
}

double get_generalized_sc(int m, int n, uint P)
{
 return G[m][n][P];      
}

void free_generalized_sc()
{
 free_recursion_commons();
}
