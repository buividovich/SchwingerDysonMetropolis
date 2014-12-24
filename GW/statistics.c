#include "statistics.h"

double                   maxnA    = 0; //Max. value of nA
double                     anA    = 0; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
double                     dnA    = 0; //Expectation value of its square - necessary to estimate the error!
double                   msign    = 0.0; //Expectation value of the sign - also necessary for recovering the correct normalization
int                        aac    = 0;
double                     ans    = 0.0;
double                     astop  = 0.0;
int                   *G_hist[2]  = {NULL, NULL};
double                  *ns_hist  = NULL;

void init_mc_stat()
{
 aac   = 0;
 ans   = 0.0;
 SAFE_MALLOC(ns_hist, double, maxn);
 set_zero_double(ns_hist, maxn);
}

void gather_mc_stat()
{
 ans         += (double)(ns+1);
 astop       += (double)(stop[ns]+1);
 if(ns>=0 && ns<maxn)
  ns_hist[ns] += 1.0;
}

void process_mc_stat()
{
 //Summarizing the post-run properties of the MC process
 double acceptance_rate      = (double)aac/(double)nmc;
 double mean_recursion_depth =         ans/(double)nmc;
 double mean_stack_size      =       astop/(double)nmc; 
 double mean_nA              =         anA/(double)nmc;
 double err_nA               = sqrt((dnA/(double)nmc - SQR(mean_nA))/(double)(nmc-1));
 double mean_sign            =       msign/(double)nmc;
 logs_Write(0, "\nStatistics on the MC process (over %i steps): ", nmc);
 logs_Write(0, "\t Acceptance rate         : \t %2.4lf", acceptance_rate);
 logs_Write(0, "\t Mean recursion depth    : \t %2.4lf", mean_recursion_depth);
 logs_Write(0, "\t Mean stack size         : \t %2.4lf", mean_stack_size);
 logs_Write(0, "\t Mean A rescaling factor : \t %2.4lf +/- %2.4lf", mean_nA, err_nA);
 logs_Write(0, "\t Max A rescaling factor  : \t %2.4lf", maxnA);
 logs_Write(0, "\t Mean config. sign       : \t %2.4lf", mean_sign);
 logs_Write(0, "\n");
 //Saving the statistical characteristics of the MC process
 if(mc_stat_file!=NULL)
 {
  FILE* fmcstat = fopen(mc_stat_file, "a");
  fprintf(fmcstat, "%2.4E %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E\n", lambda, NN, cc, acceptance_rate, mean_recursion_depth, mean_stack_size, mean_nA, err_nA);
  fclose(fmcstat);
 };
 if(ns_hist_file!=NULL)
 {
  int ins;
  FILE* nshfile = fopen(ns_hist_file, "w");
  for(ins=0; ins<maxn; ins++)
   if(ns_hist[ins]>0.0)
    fprintf(nshfile, "%i %2.4E %2.4E\n", ins, ns_hist[ins]/(double)nmc, sqrt(ns_hist[ins])/(double)nmc);
  fclose(nshfile);
 }; 
}

void free_mc_stat()
{
 SAFE_FREE(ns_hist);     
}

void init_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G_hist[s], int, maxg);
  set_zero_int(G_hist[s], maxg);
 };
 maxnA = 0.0;
 anA   = 0.0; 
 dnA   = 0.0;
 msign = 0.0;
}

void gather_observable_stat()
{
 //Histograms for observables    
 int gn = stack[ns][stop[ns]];
 ASSERT(gn<=0);
 if(gn<=maxg)
  G_hist[(asign[ns]>0? 0 : 1)][gn-1]++;
}

void process_observable_stat_sc()
{
 int ign;
 if(observables_file!=NULL)
 {
  FILE *ofile = fopen(observables_file, "a");
  fprintf(ofile, "%2.4lf ", lambda);
  for(ign=0; ign<maxg; ign++)
  {
   double rescaling_factor = NN*pow(cc,(double)(ign+1));
   //Here we restore the normalization of our solution. TODO: take into account the error of anA
   double source_norm          = 1.0/(lambda*NN*cc);
   double normalization_factor = source_norm/(1.0 - anA/(double)nmc);
   normalization_factor = normalization_factor/(1.0 + normalization_factor*msign/(double)nmc);
   //Calculate the observables
   double G  = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)nmc;
   double dG = sqrt((double)(G_hist[0][ign] + G_hist[1][ign]))/(double)nmc;
   G  *= rescaling_factor*normalization_factor;
   dG *= rescaling_factor*normalization_factor;
   //SP characterizes the strength of the sign problem
   double SP = 0.0;
   if(G_hist[0][ign]>0 || G_hist[1][ign]>0)
    SP = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)(G_hist[0][ign] + G_hist[1][ign]);
   fprintf(ofile, "%2.4E %2.4E %2.4E ", G, dG, SP);
  };
  fprintf(ofile, "\n");
  fclose(ofile);
 };
}

void process_observable_stat_wc()
{
 int ign;
 //Here we restore the normalization of our solution. TODO: take into account the error of anA
 double source_norm          = 1.0/(NN*cc);
 double normalization_factor = source_norm/(1.0 - anA/(double)nmc);
 normalization_factor = normalization_factor/(1.0 + normalization_factor*msign/(double)nmc);
 //Resummed observables
 double G1  = 1.0;
 double dG1 = 0.0;
 double fct = -2.0*alpha_wc;
 for(ign=0; ign<maxg; ign++)
 {
  double rescaling_factor = NN*pow(cc,(double)(ign+1));
  //Calculate the observables
  double G  = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)nmc;
  double dG = sqrt((double)(G_hist[0][ign] + G_hist[1][ign]))/(double)nmc;
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
  G1 += fct*G;
  dG1 += SQR(fct*dG);
  fct *= -1.0*alpha_wc;
 };
 dG1 = sqrt(dG1);
 logs_Write(0, ">>>>>>>> G1 = %2.4E +/- %2.4E, should be %2.4E\n", G1, dG1, 1.0 - 0.25*lambda);
 if(observables_file!=NULL)
 {
  FILE *ofile = fopen(observables_file, "a");
  fprintf(ofile, "%2.4lf ", lambda);
  fprintf(ofile, "%2.4E %2.4E", G1, dG1);
  fprintf(ofile, "\n");
  fclose(ofile);
 };
}


void free_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
  SAFE_FREE(G_hist[s]);
}
