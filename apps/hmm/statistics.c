#include "statistics.h"

double                     astop  = 0.0;
int                   *G_hist[2]  = {NULL, NULL};

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
 int s, i;
 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G_hist[s], int, max_correlator_order);
  for(i=0; i<max_correlator_order; i++)
   G_hist[s][i] = 0;
 };
}

void gather_observable_stat()
{
 int gn = X.len[X.top-1]/2; 
 ASSERT(gn<=0);
 if(gn<=max_correlator_order)
  G_hist[(asign[ns]>0? 0 : 1)][gn-1]++;
}

void process_observable_stat() //It is assumed that process_mc_stat was already called!!!
{
 int ign;
 FILE* ofile = NULL;
 
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 normalization_factor = normalization_factor/(1.0 + mean_sign*normalization_factor);
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);
 
 if(observables_file!=NULL)
 { 
  ofile = fopen(observables_file, "a");
  if(ofile!=NULL)
   fprintf(ofile, "%2.4E ", lambda);
  else
   logs_WriteError("Could not open the file %s for writing", observables_file); 
 };
  
 for(ign=0; ign<max_correlator_order; ign++)
 {
  double rescaling_factor     = NN*pow(cc,(double)(ign));
    
  double G  = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)nmc;
  double dG = sqrt((double)(G_hist[0][ign] + G_hist[1][ign]))/(double)nmc;
    
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
  
  double G0 = G_analytic(lambda, ign+1);
    
  //SP characterizes the strength of the sign problem
  double SP = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)(G_hist[0][ign] + G_hist[1][ign]);
  if(ofile!=NULL)
   fprintf(ofile, "%2.4E %2.4E %2.4E %2.4E ", G, dG, SP, G0);
  logs_Write(0, "G_%i:\t %2.4E +/- %2.4E,\t sp = %2.4E", ign+1, G, dG, SP);
 };
 if(ofile!=NULL)
 {
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

