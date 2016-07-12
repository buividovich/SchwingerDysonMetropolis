#include "test_rand_momentum.h"

void test_rand_momentum()
{
 logs_Write(0, "Testing the fast_random_choice for the case N=2...");
 
 for(int ipr=0; ipr<20; ipr++)
 {
  t_frc_data FD; 
  int nprobs = rand_int(5) + 1;
  DECLARE_AND_MALLOC(probs, double, nprobs);
  double sprobs = 0.0;
  for(int i=0; i<nprobs; i++)
  {
   probs[i] = rand_double(0.0, 1.0);
   sprobs += probs[i];
  };
  for(int i=0; i<nprobs; i++)
   probs[i] /= sprobs;
 
  init_fast_rand_choice(probs, nprobs, &FD);
  DECLARE_AND_MALLOC(mprobs, int, nprobs);
 
  int ntrials = 100000000;
  for(int itrial=0; itrial<ntrials; itrial++)
  {
   int m = fast_rand_choice(&FD);
   mprobs[m] ++;
   if(m<0 || m>=nprobs) logs_WriteError("Wrong value m=%i generated", m);
  };
  double err = 0.0; 
  for(int i=0; i<nprobs; i++)
   err += SQR(probs[i] - (double)mprobs[i]/(double)ntrials);
  logs_Write(0, "nprobs = %i, Error %2.4E\n", nprobs, sqrt(err));
 };
} 
 
 /*logs_Write(0, "Testing the rand_momentum...");
 int m, mc0 = 0, mc1 = 0, m0 = 0, m1 = 1;
 for(int itrial=0; itrial<10000; itrial++)
 {
  rand_momentum(&P, &m);
  if(m==0) mc0++;
  if(m==1) mc1++;
  if(m<0 || m>1) logs_WriteError("Wrong value m=%i generated", m);
 }; 
 logs_Write(0, "mc0 = %i, mc1 = %i, mc0/mc1 = %2.2lf, expected %2.2lf\n\n", mc0, mc1, (double)mc0/(double)mc1, lat_propagator(&m0, P.mass_sq)/lat_propagator(&m1, P.mass_sq) );
 return EXIT_SUCCESS;*/
