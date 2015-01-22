#include "statistics.h"

double                     astop  = 0.0;
int                   *G_hist[2]  = {NULL, NULL};

void init_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G_hist[s], int, max_correlator_order);
  set_zero_int(G_hist[s], max_correlator_order);
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
 
 double G1  = 1.0;
 double dG1 = 0.0;
 double fct = -2.0*alpha_wc;
 for(ign=0; ign<max_correlator_order; ign++)
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
  ofile = fopen(observables_file, "a");
  if(ofile!=NULL)
  {
   fprintf(ofile, "%2.4lf ", lambda);
   fprintf(ofile, "%2.4E %2.4E", G1, dG1);
   fprintf(ofile, "\n");
   fclose(ofile);
  }
  else
   logs_WriteError("Cannot open the file %s for writing", observables_file); 
 };
}

void free_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
  SAFE_FREE(G_hist[s]);
}

