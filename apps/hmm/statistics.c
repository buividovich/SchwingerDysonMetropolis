#include "statistics.h"

double                     astop = 0.0;
int**                     G_hist = NULL;

int*                    gmn_hist = NULL;
int             actual_max_genus = 0;
int             actual_max_top   = 0;
int             actual_max_nel   = 0;

double IS(double a, double k)
{
 return 0.5*sqrt(M_PI)*pow(a, 2.0 + 2.0*k)*tgamma(k + 0.5)/tgamma(k + 2.0);
}

double a2(double l)
{
 return (1.0 - sqrt(1.0 - 12.0*l))/(6.0*l);  
}

double xm(double l)
{
 return 2.0*sqrt(a2(lambda));      
}

double G_analytic(double l, int n)
{
 return (1.0/M_PI)*IS(xm(l), (double)n)*(0.5 - l*a2(l)) - (l/(2.0*M_PI))*IS(xm(l), (double)(n+1));
}

void init_observable_stat()
{
 SAFE_MALLOC(G_hist, int*, 2*(max_genus+1));
 for(int g=0; g<(max_genus+1); g++)
  for(int s=0; s<2; s++)
  {
   SAFE_MALLOC(G_hist[2*g + s], int, max_correlator_order);
   for(int i=0; i<max_correlator_order; i++)
    G_hist[2*g + s][i] = 0;
  };
  
 SAFE_MALLOC(gmn_hist, int, gmn_hist_vol);
 for(int i=0; i<gmn_hist_vol; i++)
  gmn_hist[i] = 0;
  
 actual_max_genus = 0;
 actual_max_top   = 0;
 actual_max_nel   = 0;
 
 init_stack_statistics(max_stack_nel);
}

void free_observable_stat()
{
 for(int g=0; g<(max_genus+1); g++)
  for(int s=0; s<2; s++)
   SAFE_FREE(G_hist[2*g + s]);
 SAFE_FREE(G_hist);
 SAFE_FREE(gmn_hist);
 free_stack_statistics(); 
}

void gather_observable_stat()
{
 if(genus<gmn_hist_max_genus && X.top<gmn_hist_max_traces && X.nel<gmn_hist_max_nel)
  gmn_hist[GMN_HIST_INDEX(genus, X.top, X.nel)] ++;
 
 if(X.top==1)
 {
  int gn = X.len[X.top-1]/2; 
  ASSERT(gn<=0);
  if(gn<=max_correlator_order && genus<=max_genus)
   G_hist[2*genus + (asign[ns]>0? 0 : 1)][gn-1]++;
 };  
 actual_max_genus = MAX(actual_max_genus, genus);
 actual_max_top   = MAX(actual_max_top  , X.top);
 actual_max_nel   = MAX(actual_max_nel  , X.nel);
 
 gather_stack_statistics(&X); 
}

void process_observable_stat() //It is assumed that process_mc_stat was already called!!!
{
 FILE* ofile = NULL;
 
 if(stack_stat_file!=NULL)
 {
  FILE* f = fopen(stack_stat_file, "w");
  if(f==NULL)
   logs_WriteError("Cannot open the file %s for writing", stack_stat_file);
  else
   print_stack_statistics(f);
  fclose(f);
 };
 
 if(logs_noise_level>=1)
  print_stack_statistics(stdout);
 
 double source_norm          = my_source_norm();
  
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);

 if(observables_file!=NULL)
 { 
  ofile = fopen(observables_file, "a");
  if(ofile!=NULL)
   fprintf(ofile, "%+2.4E ", lambda);
  else
   logs_WriteError("Could not open the file %s for writing", observables_file); 
 };
 
 logs_Write(0, "\nSTATISTICS ON OBSERVABLES:");
 logs_WriteParameter(0, "Max. genus   reached", "%i", actual_max_genus);
 logs_WriteParameter(0, "Max. traces  reached", "%i",   actual_max_top);
 logs_WriteParameter(0, "Max. no. el. reached", "%i",   actual_max_nel);
 
 if(gmn_hist_file!=NULL)
 {
  FILE* f = fopen(gmn_hist_file, "w");
  for(int ig=0; ig<gmn_hist_max_genus; ig++)
   for(int it=0; it<gmn_hist_max_traces; it++)
    for(int in=0; in<gmn_hist_max_nel;    in++)
    {
     int i = GMN_HIST_INDEX(ig, it, in);
     if(gmn_hist[i]>0)
      fprintf(f, "%i %i %i %2.4E\n", ig, it, in, (double)gmn_hist[i]/(double)nmc);      
    };
  fclose(f);
 };
 
 for(int g=0; g<=MIN(max_genus, actual_max_genus); g++)
 {
  logs_Write(0, "Genus %i", g); 
  
  for(int ign=0; ign<max_correlator_order; ign++)
  {
   double rescaling_factor     = f_genus[g]*NN_genus[g]*pow(cc_genus[g],(double)(2*ign+2));
    
   double G  =      (double)(G_hist[2*g + 0][ign] - G_hist[2*g + 1][ign] )/(double)nmc;
   double dG = sqrt((double)(G_hist[2*g + 0][ign] + G_hist[2*g + 1][ign]))/(double)nmc;
   
   double ac_time_estimate = MAX(mean_rt, MAX(mean_rt2, mean_rt4));
    
   G  *= rescaling_factor*normalization_factor;
   dG *= rescaling_factor*normalization_factor*sqrt(ac_time_estimate);
  
   //double G0 = G_analytic(lambda, ign+1);
   //double G0diff  = 100.0*(G0 - G)/G0; 
   //double staterr = 100.0*dG/G;
   
   //Commented output is for simulations with nontrivial lambda 
   //SP characterizes the strength of the sign problem
   double SP = (double)(G_hist[2*g + 0][ign] - G_hist[2*g + 1][ign])/(double)(G_hist[2*g + 0][ign] + G_hist[2*g + 1][ign]);
   if(ofile!=NULL && ign==1)
    fprintf(ofile, "%i %+2.4E %+2.4E ", g, G, dG); 
   logs_Write(0, "G_%02i:\t %+2.4E +/- %+2.4E, \t sp = %+2.4E (%i occurences)", ign+1, G, dG, SP, G_hist[2*g + 0][ign] +  G_hist[2*g + 1][ign]); 
  };
 };
  
 if(ofile!=NULL)
 {
  fprintf(ofile, "\n");  
  fclose(ofile);   
 };
}
